"""
Assay-based weight-of-evidence (WoE) module.

Implements a transparent "2-out-of-3" Integrated Testing Strategy (ITS) style
defined approach, combining the three in-vitro/in-chemico OECD-guideline
assays commonly used for skin sensitization hazard identification:

  - DPRA   (OECD TG 442C) — peptide reactivity (Molecular Initiating Event)
  - KeratinoSens (OECD TG 442D) — keratinocyte Nrf2/ARE activation
  - h-CLAT (OECD TG 442E) — dendritic cell activation (CD86/CD54)

IMPORTANT: the cutoff values below are pre-filled EXAMPLES only, not verified
regulatory text. OECD test guidelines are periodically updated. Before using
this for any real decision, replace the defaults with the current cutoffs
from your lab's SOP / the current OECD TG442C/D/E text. Every cutoff here is
a plain, editable number in the UI — nothing is hidden.

The final call is a simple, auditable vote count (>=2 of 3 positive =>
sensitizer; >=2 of 3 negative => non-sensitizer; otherwise inconclusive /
needs expert review), not a black-box probability.
"""

import streamlit as st


def call_dpra(depletion_pct: float, cutoff: float) -> bool:
    """True = positive (reactive)."""
    return depletion_pct >= cutoff


def call_keratinosens(imax_fold: float, imax_cutoff: float, ec_conc_uM: float, ec_cutoff_uM: float) -> bool:
    """Positive if induction exceeds fold cutoff AND the effective concentration
    triggering that induction is at/below the potency cutoff."""
    return imax_fold >= imax_cutoff and ec_conc_uM <= ec_cutoff_uM


def call_hclat(cd86_rfi: float, cd86_cutoff: float, cd54_rfi: float, cd54_cutoff: float) -> bool:
    """Positive if either marker crosses its cutoff (either-marker rule, as
    used in the standard h-CLAT prediction model)."""
    return cd86_rfi >= cd86_cutoff or cd54_rfi >= cd54_cutoff


def render_woe_module():
    st.markdown("#### 🧪 Assay-Based Weight-of-Evidence (2-out-of-3 Defined Approach)")
    st.caption(
        "Combines DPRA / KeratinoSens / h-CLAT results with a transparent vote-count "
        "rule (OECD GL 497 style defined approach). Cutoffs are editable examples — "
        "confirm against your current SOP before relying on the output."
    )

    with st.expander("⚠️ Cutoff values used (edit before real use)", expanded=False):
        st.write(
            "These are placeholder example thresholds so the calculator is usable "
            "out of the box. They are **not** guaranteed to match the current OECD "
            "TG442C/D/E text — check the guideline documents directly."
        )

    run_dpra = st.checkbox("Include DPRA (TG 442C)", value=True)
    run_ks = st.checkbox("Include KeratinoSens (TG 442D)", value=True)
    run_hclat = st.checkbox("Include h-CLAT (TG 442E)", value=True)

    votes = []
    detail_rows = []

    if run_dpra:
        st.markdown("##### DPRA")
        c1, c2 = st.columns(2)
        depletion = c1.number_input("Mean Cys+Lys peptide depletion (%)", 0.0, 100.0, 20.0, 0.1)
        cutoff = c2.number_input("Positive cutoff (% depletion)", 0.0, 100.0, 6.38, 0.01)
        result = call_dpra(depletion, cutoff)
        votes.append(result)
        detail_rows.append(("DPRA", f"{depletion}% depletion (cutoff {cutoff}%)", "Positive" if result else "Negative"))

    if run_ks:
        st.markdown("##### KeratinoSens")
        c1, c2, c3, c4 = st.columns(4)
        imax = c1.number_input("Max fold induction (Imax)", 0.0, 50.0, 2.0, 0.1)
        imax_cut = c2.number_input("Imax cutoff", 0.0, 50.0, 1.5, 0.1)
        ec = c3.number_input("EC(Imax) conc. (µM)", 0.0, 5000.0, 150.0, 1.0)
        ec_cut = c4.number_input("Potency cutoff (µM, ≤)", 0.0, 5000.0, 1000.0, 1.0)
        result = call_keratinosens(imax, imax_cut, ec, ec_cut)
        votes.append(result)
        detail_rows.append(("KeratinoSens", f"Imax {imax}x @ {ec} µM", "Positive" if result else "Negative"))

    if run_hclat:
        st.markdown("##### h-CLAT")
        c1, c2, c3, c4 = st.columns(4)
        cd86 = c1.number_input("CD86 RFI (%)", 0.0, 1000.0, 160.0, 1.0)
        cd86_cut = c2.number_input("CD86 cutoff (%, ≥)", 0.0, 1000.0, 150.0, 1.0)
        cd54 = c3.number_input("CD54 RFI (%)", 0.0, 1000.0, 180.0, 1.0)
        cd54_cut = c4.number_input("CD54 cutoff (%, ≥)", 0.0, 1000.0, 200.0, 1.0)
        result = call_hclat(cd86, cd86_cut, cd54, cd54_cut)
        votes.append(result)
        detail_rows.append(("h-CLAT", f"CD86 {cd86}% / CD54 {cd54}%", "Positive" if result else "Negative"))

    st.markdown("---")

    if not votes:
        st.info("Select at least one assay to compute a call.")
        return

    n_pos = sum(votes)
    n_total = len(votes)

    st.markdown("##### Per-assay calls")
    st.dataframe(
        {
            "Assay": [r[0] for r in detail_rows],
            "Input": [r[1] for r in detail_rows],
            "Call": [r[2] for r in detail_rows],
        },
        use_container_width=True,
        hide_index=True,
    )

    if n_total < 2:
        classification = "Inconclusive — need at least 2 assays for the 2-out-of-3 rule"
    elif n_pos >= 2:
        classification = "Sensitizer (≥2/{} assays positive)".format(n_total)
    elif (n_total - n_pos) >= 2:
        classification = "Non-sensitizer (≥2/{} assays negative)".format(n_total)
    else:
        classification = "Inconclusive — assays split; needs expert review"

    st.metric("Positive assays", f"{n_pos} / {n_total}")
    if "Sensitizer" in classification and "Non" not in classification:
        st.error(f"**Classification:** {classification}")
    elif "Non-sensitizer" in classification:
        st.success(f"**Classification:** {classification}")
    else:
        st.warning(f"**Classification:** {classification}")

    st.session_state["woe_result"] = {
        "n_positive": n_pos,
        "n_total": n_total,
        "classification": classification,
        "details": detail_rows,
    }

    st.caption(
        "This is a rule-based defined approach for hazard identification only. "
        "It does not by itself establish potency category (1A/1B) — that requires "
        "additional potency-weighted approaches (e.g. IDESI, RhE-based) per OECD GL 497."
    )
