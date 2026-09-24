"""
GHS/CLP skin sensitization potency categorization (1A vs 1B).

Two independent, transparent paths are offered — use whichever data you
actually have:

  PATH A — Point-of-departure based (the actual GHS/CLP legal criterion):
      LLNA EC3 <= cutoff  -> Category 1A (strong sensitizer)
      LLNA EC3 >  cutoff  -> Category 1B
    (Human data or GPMT thresholds can substitute for LLNA per GHS Annex I —
    this tool only implements the LLNA EC3 numeric path since that is the
    one most consistently cited.)

    CUTOFF PROVENANCE: the default 2% EC3 cutoff is corroborated by three
    independent authoritative sources (not just one secondary citation):
      - NTP (2011), evaluating the LLNA for regulatory use: "the LLNA can be
        used to categorize substances as strong sensitizers when the
        estimated concentration that produces a positive LLNA result (EC3)
        is <= 2%."
      - An EU Scientific Committee on Consumer Safety working-group table
        (citing Basketter et al. 2005, SCCP 2005, CLP guidance v4 2013, and
        GHS 5th revision 2013) gives the CLP/GHS potency bands as EC3 >0.2-<=2%
        -> "strong, 1A" and EC3 >2% -> "moderate, 1B".
      - A SenzaGen scientific poster, citing Annex 2 of the Supporting
        Document to OECD Guideline 497, states the GHS classification
        threshold as 500 ug/cm2, "corresponding to an LLNA EC3 of 2%".
    This is corroborated secondary citation of the CLP Regulation (EC) No
    1272/2008 Annex I criterion, not a direct excerpt of that legal text
    itself -- still worth being precise about, since it's a meaningfully
    different evidentiary bar than pulling the primary document verbatim.
    Confirm against the current Annex I text (or your local equivalent)
    before a real classification & labelling decision, since legal texts
    are periodically revised.

  PATH B — In-vitro-battery proxy score (NO animal/human PoD available):
    A simplified points-based scheme loosely based on the published
    Integrated Testing Strategy (ITS) concept (Urbisch et al. 2015): each
    assay contributes 0-3 points from potency-related bands, points are
    summed, and the total maps to a screening-level category. This is a
    SCREENING PROXY, not a GHS-official classification. UNLIKE the DPRA
    band edges in the WoE module (which are verified against the primary
    DPRA validation report), the specific point-band numbers here have NOT
    been independently verified against Urbisch et al. 2015's own published
    table -- they are placeholder, editable defaults. Recalibrate against
    the primary literature / your SOP before relying on this path.

IMPORTANT: GHS/CLP criteria and OECD guideline cutoffs are periodically
revised. Verify all thresholds against the current legal text (CLP
Regulation (EC) No 1272/2008, Annex I, or your local equivalent) before
using this for classification & labelling decisions.
"""

import streamlit as st


def classify_from_llna(ec3_pct: float, cutoff_1a: float = 2.0) -> str:
    if ec3_pct <= 0:
        return "Invalid EC3"
    return "Category 1A (strong sensitizer)" if ec3_pct <= cutoff_1a else "Category 1B"


def its_band_score(value: float, bands: list) -> int:
    """bands = [(threshold, points), ...] ascending threshold; returns points for first band value < threshold,
    else max points."""
    for threshold, points in bands:
        if value < threshold:
            return points
    return bands[-1][1]


def render_potency_module():
    st.markdown("#### 🏷️ Potency Categorization (GHS 1A / 1B)")
    st.caption(
        "Two paths: use real animal/human point-of-departure data if you have it (Path A), "
        "or a screening-level in-vitro proxy score if you don't (Path B). Path B is NOT "
        "a GHS-official classification."
    )

    path = st.radio("Which data do you have?", [
        "Path A — LLNA EC3 (or equivalent PoD) available",
        "Path B — Only in-vitro battery data (DPRA/KeratinoSens/h-CLAT), no PoD",
    ])

    if path.startswith("Path A"):
        with st.expander("📖 Cutoff provenance — what's verified vs. not", expanded=False):
            st.markdown(
                "The 2% EC3 cutoff is corroborated by three independent authoritative "
                "sources: NTP (2011)'s LLNA regulatory evaluation, an EU Scientific "
                "Committee on Consumer Safety working-group table (citing Basketter et al. "
                "2005, SCCP 2005, CLP guidance v4, and GHS 5th revision), and a SenzaGen "
                "poster citing Annex 2 of the Supporting Document to OECD GL 497 — all "
                "three independently state EC3 ≤ 2% corresponds to GHS/CLP Category 1A. "
                "This is corroborated secondary citation of the CLP Regulation (EC) No "
                "1272/2008 Annex I criterion, not a verbatim excerpt of that legal text — "
                "confirm against the current Annex I text before a real classification & "
                "labelling decision, since legal texts are periodically revised."
            )
        c1, c2 = st.columns(2)
        ec3 = c1.number_input("LLNA EC3 (%)", 0.0, 100.0, 5.0, 0.1)
        cutoff = c2.number_input("1A/1B cutoff (% EC3, ≤ = 1A)", 0.0, 100.0, 2.0, 0.1,
                                  help="Corroborated by 3 independent sources — see provenance note above. "
                                       "Still confirm against your current Annex I text before real use.")
        category = classify_from_llna(ec3, cutoff)
        st.session_state["potency_result"] = {"path": "A", "ec3_pct": ec3, "cutoff": cutoff, "category": category}

        if category == "Invalid EC3":
            st.error("Enter an EC3 greater than 0% to classify this compound.")
        elif "1A" in category:
            st.error(f"**{category}** (EC3 {ec3}% ≤ cutoff {cutoff}%)")
        else:
            st.warning(f"**{category}** (EC3 {ec3}% > cutoff {cutoff}%)")

    else:
        st.markdown("##### Screening-level ITS-style proxy score (editable bands)")
        st.caption(
            "Default bands loosely follow the published DPRA/h-CLAT/KeratinoSens "
            "potency-banding concept, but — unlike the DPRA band edges in the WoE tab, "
            "which are verified against the primary DPRA validation report — these "
            "specific point-band numbers have NOT been independently verified against "
            "Urbisch et al. 2015's own published table. Edit them to match your source "
            "literature before trusting the output."
        )

        c1, c2 = st.columns(2)
        with c1:
            dpra_depletion = st.number_input("DPRA mean Cys+Lys depletion (%)", 0.0, 100.0, 20.0, 0.1)
            dpra_bands_txt = st.text_input("DPRA bands (upper%:points, ascending)", "6.38:0,22.62:1,42.47:2,999:3")
        with c2:
            hclat_ec = st.number_input("h-CLAT CD86 EC150 (µM, lower = more potent)", 0.0, 5000.0, 500.0, 1.0)
            hclat_bands_txt = st.text_input("h-CLAT bands (upper µM:points, descending potency)", "10:3,100:2,1000:1,999999:0")

        c3, c4 = st.columns(2)
        with c3:
            ks_ec = st.number_input("KeratinoSens EC1.5 (µM, lower = more potent)", 0.0, 5000.0, 150.0, 1.0)
            ks_bands_txt = st.text_input("KeratinoSens bands (upper µM:points, descending potency)", "10:3,100:2,1000:1,999999:0")

        def parse_bands(txt):
            bands = []
            for part in txt.split(","):
                thr, pts = part.split(":")
                bands.append((float(thr), int(pts)))
            return sorted(bands, key=lambda x: x[0])

        try:
            dpra_pts = its_band_score(dpra_depletion, parse_bands(dpra_bands_txt))
            hclat_pts = its_band_score(hclat_ec, parse_bands(hclat_bands_txt))
            ks_pts = its_band_score(ks_ec, parse_bands(ks_bands_txt))
        except Exception:
            st.error("Could not parse band definitions — use format 'threshold:points,threshold:points,...'")
            return

        total = dpra_pts + hclat_pts + ks_pts
        st.markdown("##### Component scores")
        st.dataframe(
            {"Assay": ["DPRA", "h-CLAT", "KeratinoSens"], "Points": [dpra_pts, hclat_pts, ks_pts]},
            use_container_width=True, hide_index=True,
        )

        c1, c2 = st.columns(2)
        cutoff_1a = c1.number_input("Total score ≥ this → proxy 1A (strong)", 0, 9, 5)
        cutoff_1b = c2.number_input("Total score ≥ this (and < 1A cutoff) → proxy 1B", 0, 9, 2)

        if total >= cutoff_1a:
            proxy = "Proxy: STRONG sensitizer signal (≈1A-range)"
            st.error(f"**{proxy}** — total score {total}")
        elif total >= cutoff_1b:
            proxy = "Proxy: WEAK-to-moderate sensitizer signal (≈1B-range)"
            st.warning(f"**{proxy}** — total score {total}")
        else:
            proxy = "Proxy: LOW/no sensitizer signal"
            st.success(f"**{proxy}** — total score {total}")

        st.session_state["potency_result"] = {
            "path": "B", "dpra_points": dpra_pts, "hclat_points": hclat_pts,
            "keratinosens_points": ks_pts, "total_score": total, "category": proxy,
        }
        st.caption(
            "This proxy score has NOT been validated against an external dataset in this "
            "tool. Treat it as a triage signal to prioritize which compounds need real "
            "PoD data (LLNA/human), not as a classification decision on its own."
        )
