"""
Assay-based weight-of-evidence (WoE) module.

Implements the "2 out of 3" (2o3) Defined Approach -- one of exactly three
DAs formally adopted in OECD Guideline No. 497, "Defined Approaches for Skin
Sensitisation" (DASS GL, published June 2021). GL 497 defines precisely
three DAs, not a loose family of "497-style" methods:
  - 2o3 (this module): hazard-only (sensitizer/non-sensitizer), using
    DPRA + KeratinoSens + h-CLAT. Traces to Bauch et al. (2012) and
    Urbisch et al. (2015).
  - ITSv1: GHS potency sub-categorization (1A/1B/NC) using DPRA + h-CLAT +
    an in-silico prediction from Derek Nexus (proprietary). Takenouchi
    et al. (2015).
  - ITSv2: same as ITSv1 but with the OECD QSAR Toolbox instead of Derek
    Nexus as the in-silico component.
This tool implements 2o3 here (verified against the current adopted OECD
test guideline text for each assay's own cutoff, see below). It does NOT
implement ITSv1 or ITSv2 -- both require either Derek Nexus or the OECD QSAR
Toolbox as their in-silico component, neither of which is available in this
environment. The Potency tab's "Path B" is explicitly labeled as an
independent, unverified proxy inspired by the general scoring concept, not
an implementation of GL 497's actual ITSv1/ITSv2 -- see potency.py's
docstring for that distinction.

The three assays combined here:
  - DPRA   (OECD TG 442C) — peptide reactivity (Molecular Initiating Event)
  - KeratinoSens (OECD TG 442D) — keratinocyte Nrf2/ARE activation
  - h-CLAT (OECD TG 442E) — dendritic cell activation (CD86/CD54)

INPUT CURATION NOTE (a practice worth adopting from other tools in this
space, e.g. LabMol's PredSkin): this module does not standardize salts,
counter-ions, or tautomers in any SMILES entered elsewhere in this app.
Predictions and structural-alert checks are sensitive to the exact form
submitted -- curate inputs (remove counter-ions, select a canonical
tautomer) before relying on results for a salt or mixture, the same way
this app already flags per-compound where relevant (e.g. chlorpromazine
tested as its HCl salt, with the free base SMILES used for structure
matching).

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

  - h-CLAT LogP>3.5 domain limitation (verified against primary TG 442E
    text): "Test chemicals with a Log Kow greater than 3.5 tend to produce
    false negative results... negative results with test chemicals with a
    Log Kow greater than 3.5 should not be considered." A NEGATIVE h-CLAT
    result is therefore excluded from the vote (not silently trusted) when
    the compound's LogP (pulled from the Structure tab, if run) exceeds 3.5
    -- a quantified, real effect: sensitivity drops from 94% to 52% for
    these compounds per OECD's own supporting analysis
    (ENV/CBC/MONO(2025)2, Annex 5). A POSITIVE result at high LogP is
    unaffected, per the same guideline text.

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
                min_viability_pct: float = 50.0, logp: float = None,
                logp_domain_limit: float = 3.5) -> str:
    """Returns 'positive', 'negative', 'unreliable_negative', or 'invalid'.

    Default cutoffs (CD86 >184%, CD54 >255%) are the current OECD TG 442E
    positive-call thresholds, not the older 150%/200% values sometimes cited
    in secondary literature. A result below the minimum viability is reported
    as 'invalid' rather than forced into a call, per the guideline's own
    validity requirement. This function does not model the guideline's
    separate statistically-derived borderline range (not verified here --
    see module docstring).

    LogP > 3.5 domain limitation, verified against OECD TG 442E primary text:
    "Test chemicals with a Log Kow greater than 3.5 tend to produce false
    negative results... negative results with test chemicals with a Log Kow
    greater than 3.5 should not be considered." A supporting OECD analysis
    (ENV/CBC/MONO(2025)2, Annex 5, citing Takenouchi et al. 2013) quantifies
    this: h-CLAT sensitivity drops from 94% (LogP<=3.5) to 52% (LogP>3.5) --
    a 42-percentage-point reduction, not a minor caveat. A NEGATIVE result
    with LogP>3.5 is therefore reported as 'unreliable_negative', not
    'negative' -- per the guideline, it should not be treated as evidence of
    non-sensitization at all. A POSITIVE result with high LogP is unaffected
    ("could still be used to support the identification... as a skin
    sensitiser") and is still reported as 'positive'."""
    if viability_pct < min_viability_pct:
        return "invalid"
    if cd86_rfi > cd86_cutoff or cd54_rfi > cd54_cutoff:
        return "positive"
    if logp is not None and logp > logp_domain_limit:
        return "unreliable_negative"
    return "negative"


def render_woe_module():
    st.markdown("#### 🧪 Assay-Based Weight-of-Evidence (2-out-of-3 Defined Approach)")
    st.caption(
        "Implements the \"2 out of 3\" (2o3) Defined Approach — one of exactly three DAs "
        "formally adopted in OECD Guideline No. 497 (June 2021), not a loose \"497-style\" "
        "concept. Default cutoffs below are the current OECD prediction-model thresholds — "
        "see the 'cutoff provenance' note for exactly what was verified against primary text."
    )
    st.info(
        "GL 497 defines three DAs total: 2o3 (hazard only — implemented here), and ITSv1/ITSv2 "
        "(GHS potency sub-categorization, requiring Derek Nexus or the OECD QSAR Toolbox as an "
        "in-silico component — neither available here, so this app does NOT implement the "
        "official ITS). The Potency tab's Path B is an independent, unverified proxy inspired "
        "by the general scoring concept, not GL 497's actual ITS."
    )

    with st.expander("📖 Cutoff provenance — what's verified vs. not", expanded=False):
        st.markdown(
            "- **DPRA** 6.38 / 22.62 / 42.47% mean-depletion bands, AND 13.89 / 23.09 / 98.24% "
            "cysteine-only-depletion bands (used when peptides co-elute with lysine): both verified "
            "against OECD GL 497 (25 June 2025 edition), Table 3.1, and the ECVAM DPRA Validation "
            "Study Report.\n"
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
        depletion_type = st.radio(
            "Which depletion measure do you have?",
            ["Mean Cys+Lys depletion (standard)", "Cysteine-only depletion (used when peptides co-elute)"],
            horizontal=True,
            help="Per OECD GL 497 Table 3.1 footnote: cysteine-only depletion is used specifically "
                 "in the case of co-elution with the lysine peptide — it is a separate, real criterion "
                 "with its own cutoff, not an approximation of the mean-depletion criterion.",
        )
        is_cys_only = depletion_type.startswith("Cysteine-only")
        default_cutoff = 13.89 if is_cys_only else 6.38
        label = "Cysteine-only peptide depletion (%)" if is_cys_only else "Mean Cys+Lys peptide depletion (%)"

        c1, c2 = st.columns(2)
        depletion = c1.number_input(label, 0.0, 100.0, 20.0, 0.1, key=f"dpra_depletion_value_{is_cys_only}")
        cutoff = c2.number_input("Positive cutoff (% depletion)", 0.0, 100.0, default_cutoff, 0.01,
                                  key=f"dpra_cutoff_value_{is_cys_only}",
                                  help="6.38% (mean) or 13.89% (cysteine-only) — both verified against "
                                       "OECD GL 497 Table 3.1 / the ECVAM DPRA Validation Study Report.")
        is_positive = call_dpra(depletion, cutoff)
        votes.append("yes" if is_positive else "no")
        metric_label = "Cys-only" if is_cys_only else "mean"
        detail_rows.append(("DPRA", f"{depletion}% {metric_label} depletion (cutoff {cutoff}%)",
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

        descriptors = st.session_state.get("structure_descriptors")
        logp = descriptors.get("LogP (Crippen)") if descriptors else None
        if logp is not None:
            st.caption(f"Using LogP={logp} from the Structure tab for the TG 442E domain check below.")
        else:
            st.caption("Run the Structure tab first to enable the LogP>3.5 domain check (OECD TG 442E) for a negative result.")

        hclat_call = call_hclat(cd86, cd54, viability, logp=logp)
        votes.append(hclat_call if hclat_call in ("positive", "negative") else "excluded")
        label = {"positive": "Positive", "negative": "Negative",
                 "unreliable_negative": f"Negative, but UNRELIABLE (LogP={logp}>3.5 — excluded from vote per TG 442E)",
                 "invalid": "Invalid (viability <50% — excluded from vote)"}[hclat_call]
        detail_rows.append(("h-CLAT", f"CD86 {cd86}% / CD54 {cd54}% (viability {viability}%)", label))
        if hclat_call == "unreliable_negative":
            st.warning(
                f"h-CLAT negative result excluded: LogP={logp} exceeds the OECD TG 442E domain limit of 3.5. "
                "Per the guideline's own text, negative h-CLAT results above this LogP should not be considered "
                "— sensitivity drops from 94% (LogP≤3.5) to 52% (LogP>3.5) in the supporting OECD analysis "
                "(ENV/CBC/MONO(2025)2, Annex 5, citing Takenouchi et al. 2013)."
            )

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
