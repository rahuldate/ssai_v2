"""
Report / dossier export module.

Pulls ONLY values that were actually computed in this session (structure
descriptors, WoE result, QRA result). If a module hasn't been run, its
section says so explicitly rather than inventing a placeholder number.
This is a decision-support summary export, not a certified regulatory
submission — the disclaimer is baked into the document itself.
"""

import streamlit as st
import datetime
import io

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors


def _section(story, styles, title, lines):
    story.append(Paragraph(title, styles["heading"]))
    if not lines:
        story.append(Paragraph("Not yet run in this session.", styles["body"]))
    else:
        for line in lines:
            story.append(Paragraph(line, styles["body"]))
    story.append(Spacer(1, 8))


def build_pdf(compound_name, assessor, framework_note) -> bytes:
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)

    base = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle("t", parent=base["Heading1"], fontSize=17, textColor=colors.HexColor("#0d6efd")),
        "heading": ParagraphStyle("h", parent=base["Heading2"], fontSize=12, spaceBefore=12, spaceAfter=4),
        "body": ParagraphStyle("b", parent=base["Normal"], fontSize=9.5, leading=13, spaceAfter=3),
        "warn": ParagraphStyle("w", parent=base["Normal"], fontSize=9, leading=12, textColor=colors.HexColor("#a00000")),
    }

    story = []
    story.append(Paragraph("Skin Sensitization Decision-Support Summary", styles["title"]))
    story.append(Paragraph(f"<b>Compound:</b> {compound_name}", styles["body"]))
    story.append(Paragraph(f"<b>Prepared by:</b> {assessor}", styles["body"]))
    story.append(Paragraph(f"<b>Generated:</b> {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles["body"]))
    story.append(Paragraph(f"<b>Framework note:</b> {framework_note}", styles["body"]))
    story.append(Paragraph(
        "This document reflects calculations performed on the inputs entered by the "
        "preparer during this session. It is a decision-support summary only, is NOT "
        "a validated or certified regulatory submission, and must be reviewed by a "
        "qualified toxicologist before any regulatory or safety use.",
        styles["warn"]
    ))
    story.append(Spacer(1, 10))

    # Structure section
    desc = st.session_state.get("structure_descriptors")
    lines = []
    if desc:
        for k, v in desc.items():
            lines.append(f"{k}: {v}")
        alerts = st.session_state.get("structure_alerts", {})
        flagged = [k for k, v in alerts.items() if v]
        lines.append("Structural alerts flagged: " + (", ".join(flagged) if flagged else "None"))
    _section(story, styles, "1. Structural & Physicochemical Analysis", lines)

    # WoE section
    woe = st.session_state.get("woe_result")
    lines = []
    if woe:
        lines.append(f"Positive assays: {woe['n_positive']} / {woe['n_total']}")
        lines.append(f"Classification: {woe['classification']}")
        for name, detail, call in woe["details"]:
            lines.append(f"{name}: {detail} → {call}")
    _section(story, styles, "2. Assay-Based Weight-of-Evidence (2-out-of-3)", lines)

    # Read-across section
    ra = st.session_state.get("read_across_result")
    lines = []
    if ra:
        lines.append(f"Call: {ra['call']}")
        lines.append(f"Applicability domain cutoff (Tanimoto): {ra['domain_cutoff']}")
        for n in ra["top_k"]:
            lines.append(
                f"Neighbor: {n['name']} (similarity {n['similarity']:.3f}, "
                f"sensitizer={n['sensitizer']}, potency={n['potency_category']})"
            )
    _section(story, styles, "3. Read-Across (Analog-Based) Prediction", lines)

    # Potency section
    potency = st.session_state.get("potency_result")
    lines = []
    if potency:
        if potency["path"] == "A":
            lines.append(f"Path: LLNA EC3-based")
            lines.append(f"EC3: {potency['ec3_pct']}% (cutoff {potency['cutoff']}%)")
            lines.append(f"Category: {potency['category']}")
        else:
            lines.append("Path: In-vitro battery proxy score (screening-level, not GHS-official)")
            lines.append(f"DPRA points: {potency['dpra_points']}, h-CLAT points: {potency['hclat_points']}, "
                          f"KeratinoSens points: {potency['keratinosens_points']}")
            lines.append(f"Total score: {potency['total_score']}")
            lines.append(f"Result: {potency['category']}")
    _section(story, styles, "4. Potency Categorization (GHS 1A/1B)", lines)

    # QRA section
    qra = st.session_state.get("qra_result")
    lines = []
    if qra:
        lines.append(f"NESIL: {qra['nesil']} µg/cm² (source: {qra['pod_source']})")
        lines.append(f"SAF total: {qra['SAF_total']:.1f}x")
        lines.append(f"AEL: {qra['AEL_ug_cm2_day']:.3f} µg/cm²/day")
        lines.append(f"Consumer Exposure (CEL): {qra['CEL_per_day_ug_cm2']:.3f} µg/cm²/day")
        lines.append(f"Margin of Safety: {qra['MoS']:.2f}")
    _section(story, styles, "5. Dermal Sensitization QRA", lines)

    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()


def render_dossier_module():
    st.markdown("#### 📄 Decision-Support Summary Export")
    st.caption(
        "Compiles whatever you've actually computed in the Structure, WoE, and QRA "
        "tabs this session. Sections you haven't run will say so — nothing is filled in."
    )

    compound_name = st.text_input("Compound name", value=st.session_state.get("compound_name", ""))
    assessor = st.text_input("Prepared by", value="")
    framework_note = st.text_input(
        "Framework note (free text — e.g. intended regulatory context)",
        value=""
    )

    status_rows = {
        "Section": ["Structural analysis", "Assay-based WoE", "Read-across", "Potency (1A/1B)", "Dermal QRA"],
        "Status": [
            "✅ Computed" if st.session_state.get("structure_descriptors") else "— Not run",
            "✅ Computed" if st.session_state.get("woe_result") else "— Not run",
            "✅ Computed" if st.session_state.get("read_across_result") else "— Not run",
            "✅ Computed" if st.session_state.get("potency_result") else "— Not run",
            "✅ Computed" if st.session_state.get("qra_result") else "— Not run",
        ],
    }
    st.dataframe(status_rows, use_container_width=True, hide_index=True)

    pdf_bytes = build_pdf(compound_name or "Unnamed compound", assessor or "Unspecified", framework_note or "Not specified")
    st.download_button(
        "📥 Download Summary PDF",
        data=pdf_bytes,
        file_name=f"skin_sens_summary_{(compound_name or 'compound').replace(' ', '_')}.pdf",
        mime="application/pdf",
        type="primary",
        use_container_width=True,
    )
