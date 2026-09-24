"""
Assay-based weight-of-evidence (WoE) module.

Implements a transparent "2-out-of-3" Integrated Testing Strategy (ITS) style
defined approach, combining the three in-vitro/in-chemico OECD-guideline
assays commonly used for skin sensitization hazard identification:

  - DPRA   (OECD TG 442C) — peptide reactivity (Molecular Initiating Event)
  - KeratinoSens (OECD TG 442D) — keratinocyte Nrf2/ARE activation
  - h-CLAT (OECD TG 442E) — dendritic cell activation (CD86/CD54)

CUTOFF PROVENANCE (checked against the current adopted OECD test guideline
text directly, not memory or secondary sources -- see README for the exact
source documents):

  - DPRA: the 6.38 / 22.62 / 42.47 %-depletion band edges are the actual
    published DPRA prediction model thresholds (minimal/low/moderate/high
    reactivity), confirmed against the ECVAM DPRA Validation Study Report.

  - KeratinoSens: the 1.5-fold positive threshold is correct, but the
    current OECD TG 442D text also defines a statistically-derived
    BORDERLINE range of 1.35-1.67-fold around it, and requires the EC1.5
    concentration to fall below the IC30 (not an arbitrary µM ceiling) for
    a positive call to be valid. Both are modeled here now -- a prior
    version of this tool used an arbitrary "<=1000 uM" proxy, which did not
    reflect the actual guideline and has been corrected.

  - h-CLAT: THIS TOOL PREVIOUSLY USED CD86>=150% / CD54>=200%, WHICH ARE THE
    ORIGINAL 2006-2007 LABORATORY VALUES, NOT THE CURRENT GUIDELINE. The
    adopted OECD TG 442E text (2023) uses CD86 >184% and CD54 >255% (at
    cell viability >=50%) for a positive call. This has been corrected.
    The current guideline also defines a statistically-derived borderline
    range for h-CLAT, which is NOT modeled here because its exact numeric
    bounds were not directly verified against primary text -- this is
    flagged rather than guessed at. Note also that OECD TG 442E's own
    internal "2 out of 3" rule is about replicate RUNS of h-CLAT itself
    (2 of 3 independent runs agreeing), which is a different "2 out of 3"
    from this tool's cross-ASSAY vote (DPRA vs KeratinoSens vs h-CLAT) --
    the single CD86/CD54 input below is assumed to already be the lab's
    finalized per-chemical h-CLAT call, not a single replicate.

The final cross-assay call is a simple, auditable vote count (>=2 of 3
positive => sensitizer; >=2 of 3 negative => non-sensitizer; otherwise
inconclusive / needs expert review), not a black-box probability. Borderline
or inconclusive per-assay results are excluded from the vote, same as a
missing/unlabeled result, rather than being forced into positive or negative.
"""

import streamlit as st


def call_dpra(depletion_pct: float, cutoff: float) -> bool:
    """True = positive (reactive). Default cutoff (6.38%) is the real minimal/low
    reactivity boundary from the published DPRA prediction model."""
    return depletion_pct >= cutoff


def call_keratinosens(imax_fold: float, ec1_5_uM: float, ic30_uM: float, max_tested_uM: float,
                       borderline_low: float = 1.35, borderline_high: float = 1.67) -> str:
    """Returns 'positive', 'negative', 'borderline', or 'inconclusive'.

    Modeled on the actual OECD TG 442D prediction model: the statistically-derived
    borderline range (1.35-1.67-fold) around the 1.5-fold threshold, and the
    requirement that a positive call's EC1.5 concentration be below the IC30
    (i.e. induction happens at a non-cytotoxic concentration) rather than only
    at a cytotoxic one. A negative result is only conclusive if testing reached
    either cytotoxicity or the 1000 uM (or equivalent) ceiling described in the
    guideline; otherwise it is reported as inconclusive rather than negative.
    """
    if imax_fold >= borderline_high:
        if ec1_5_uM is not None and ic30_uM is not None and ec1_5_uM >= ic30_uM:
            return "inconclusive"  # induction only reached at a cytotoxic concentration
        return "positive"
    elif imax_fold < borderline_low:
        reached_ceiling = (max_tested_uM is not None and max_tested_uM >= 1000) or \
                           (ic30_uM is not None and max_tested_uM is not None and max_tested_uM >= ic30_uM)
        if not reached_ceiling:
            return "inconclusive"
        return "negative"
    else:
        return "borderline"


def call_hclat(cd86_rfi: float, cd54_rfi: float, viability_pct: float,
                cd86_cutoff: float = 184.0, cd54_cutoff: float = 255.0,
                min_viability_pct: float = 50.0) -> str:
    """Returns 'positive', 'negative', or 'invalid'.

    Default cutoffs (CD86 >184%, CD54 >255%) are the current OECD TG 442E
    positive-call thresholds, not the older 150%/200% values sometimes cited
    in secondary literature. A result below the minimum viability is reported
    as 'invalid' rather than forced into a call, per the guideline's own
    validity requirement. This function does not model the guideline's
    separate statistically-derived borderline range (not verified here --
    see module docstring)."""
    if viability_pct < min_viability_pct:
        return "invalid"
    if cd86_rfi > cd86_cutoff or cd54_rfi > cd54_cutoff:
        return "positive"
    return "negative"


def render_woe_module():
    st.markdown("#### 🧪 Assay-Based Weight-of-Evidence (2-out-of-3 Defined Approach)")
    st.caption(
        "Combines DPRA / KeratinoSens / h-CLAT results with a transparent vote-count "
        "rule (OECD GL 497 style defined approach). Default cutoffs below are the "
        "current OECD prediction-model thresholds — see the 'cutoff provenance' note "
        "for exactly what was verified against primary text and what wasn't."
    )

    with st.expander("📖 Cutoff provenance — what's verified vs. not", expanded=False):
        st.markdown(
            "- **DPRA** 6.38 / 22.62 / 42.47% depletion bands: verified against the "
            "ECVAM DPRA Validation Study Report's own prediction-model figure.\n"
            "- **KeratinoSens**: 1.5-fold threshold plus the 1.35–1.67-fold borderline "
            "range, and the EC1.5-must-be-below-IC30 validity rule, verified against "
            "the current OECD TG 442D text.\n"
            "- **h-CLAT**: CD86 >184% / CD54 >255% (viability ≥50%) verified against "
            "the current OECD TG 442E text — **this replaces the older 150%/200% "
            "values** that appear in a lot of secondary literature and were used in "
            "an earlier version of this tool. The guideline's own h-CLAT borderline "
            "range is **not** modeled here since its exact numeric bounds weren't "
            "directly verified — treat a result close to the cutoff with extra caution."
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
        is_positive = call_dpra(depletion, cutoff)
        votes.append("yes" if is_positive else "no")
        detail_rows.append(("DPRA", f"{depletion}% depletion (cutoff {cutoff}%)",
                             "Positive" if is_positive else "Negative"))

    if run_ks:
        st.markdown("##### KeratinoSens")
        c1, c2 = st.columns(2)
        imax = c1.number_input("Max fold induction (Imax)", 0.0, 50.0, 2.0, 0.1)
        ec1_5 = c2.number_input("EC1.5 concentration (µM)", 0.0, 5000.0, 150.0, 1.0)
        c3, c4 = st.columns(2)
        ic30 = c3.number_input("IC30 — conc. causing 30% viability loss (µM)", 0.0, 5000.0, 800.0, 1.0)
        max_tested = c4.number_input("Maximum concentration actually tested (µM)", 0.0, 5000.0, 1000.0, 1.0)
        ks_call = call_keratinosens(imax, ec1_5, ic30, max_tested)
        votes.append(ks_call if ks_call in ("positive", "negative") else "excluded")
        label = {"positive": "Positive", "negative": "Negative",
                 "borderline": "Borderline (excluded from vote)",
                 "inconclusive": "Inconclusive (excluded from vote)"}[ks_call]
        detail_rows.append(("KeratinoSens", f"Imax {imax}x, EC1.5={ec1_5}µM, IC30={ic30}µM", label))

    if run_hclat:
        st.markdown("##### h-CLAT")
        c1, c2, c3 = st.columns(3)
        cd86 = c1.number_input("CD86 RFI (%)", 0.0, 1000.0, 190.0, 1.0)
        cd54 = c2.number_input("CD54 RFI (%)", 0.0, 1000.0, 210.0, 1.0)
        viability = c3.number_input("Cell viability (%)", 0.0, 100.0, 80.0, 1.0)
        hclat_call = call_hclat(cd86, cd54, viability)
        votes.append(hclat_call if hclat_call in ("positive", "negative") else "excluded")
        label = {"positive": "Positive", "negative": "Negative",
                 "invalid": "Invalid (viability <50% — excluded from vote)"}[hclat_call]
        detail_rows.append(("h-CLAT", f"CD86 {cd86}% / CD54 {cd54}% (viability {viability}%)", label))

    st.markdown("---")

    if not votes:
        st.info("Select at least one assay to compute a call.")
        return

    usable_votes = [v for v in votes if v in ("yes", "positive", "no", "negative")]
    n_pos = sum(1 for v in usable_votes if v in ("yes", "positive"))
    n_total = len(usable_votes)
    n_excluded = len(votes) - n_total

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
    if n_excluded:
        st.caption(f"{n_excluded} assay result(s) excluded from the vote (borderline/inconclusive/invalid).")

    if n_total < 2:
        classification = "Inconclusive — need at least 2 usable assay results for the 2-out-of-3 rule"
    elif n_pos >= 2:
        classification = "Sensitizer (≥2/{} usable assays positive)".format(n_total)
    elif (n_total - n_pos) >= 2:
        classification = "Non-sensitizer (≥2/{} usable assays negative)".format(n_total)
    else:
        classification = "Inconclusive — assays split; needs expert review"

    st.metric("Positive assays (of usable)", f"{n_pos} / {n_total}")
    if "Sensitizer" in classification and "Non" not in classification:
        st.error(f"**Classification:** {classification}")
    elif "Non-sensitizer" in classification:
        st.success(f"**Classification:** {classification}")
    else:
        st.warning(f"**Classification:** {classification}")

    st.session_state["woe_result"] = {
        "n_positive": n_pos,
        "n_total": n_total,
        "n_excluded": n_excluded,
        "classification": classification,
        "details": detail_rows,
    }

    st.caption(
        "This is a rule-based defined approach for hazard identification only. "
        "It does not by itself establish potency category (1A/1B) — that requires "
        "additional potency-weighted approaches (e.g. IDESI, RhE-based) per OECD GL 497."
    )
