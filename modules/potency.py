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

  PATH B — OECD Guideline No. 497 "Integrated Testing Strategy" (ITS)
    Defined Approach for GHS potency sub-categorization.

    VERIFIED directly against the primary text of OECD Guideline No. 497,
    "Defined Approaches on Skin Sensitisation" (25 June 2025 edition),
    Part II, Section 3.1.2, Table 3.1 -- fetched and cross-checked, not a
    secondary paraphrase. (Also independently corroborated by the earlier
    2021 supporting document ENV/CBC/MONO(2021)11, a NICEATM/ICCVAM
    scientist's SCHC 2022 presentation, and Kleinstreuer et al. 2022 in
    Frontiers in Toxicology.) The exact bands from Table 3.1:

      DPRA depletion (%):   <6.38->0, >=6.38 & <22.62->1, >=22.62 & <42.47->2, >=42.47->3
      h-CLAT MIT (ug/mL):   negative->0, >150 & <=5000 (positive)->1, >10 & <=150->2, <=10->3
      In silico prediction: negative->0, positive->1
      Total score (0-7):    0-1 -> Not Classified, 2-5 -> GHS 1B, 6-7 -> GHS 1A

    The DPRA band edges match what this app's WoE tab already independently
    verified against the primary DPRA validation report -- that agreement
    across two separately-verified sources is a good sign, not a
    coincidence: both trace back to the same underlying DPRA prediction
    model. Note h-CLAT here uses MIT (Minimum Induction Threshold, ug/mL --
    the lower of the CD86 EC150 or CD54 EC200 concentrations), which is a
    DIFFERENT metric in DIFFERENT units than the CD86/CD54 %RFI the WoE
    tab's h-CLAT call uses for hazard identification -- they are not
    interchangeable, and this app does not attempt to convert between them.

    IMPORTANT SUBSTITUTION, disclosed rather than hidden: the real ITS (v1
    or v2) uses Derek Nexus or the OECD QSAR Toolbox for the in silico
    component -- this app has neither. In their place, Path B offers using
    this app's own Read-Across tab result as the in silico input. That
    substitution is NOT validated to the same standard as Derek Nexus or
    the OECD QSAR Toolbox -- the SCORING FRAMEWORK above is verified against
    primary text, but plugging in a different (and, per this app's own
    LOO-CV validation, only moderately predictive, MCC~0.3) in silico source
    means this tool's Path B accuracy cannot be assumed to match the real
    ITS's published performance. For reference, real ITSv2 achieves ~80%
    balanced accuracy for hazard and ~67-72% accuracy across GHS categories
    for potency (Strickland 2022) -- that describes the validated method
    with its real in silico component, not this app's substitution.

  PATH C -- kDPRA (kinetic Direct Peptide Reactivity Assay), log(kmax):
    Unlike Path A and Path B, this is a REAL, OECD-ADOPTED, single-assay
    method for 1A determination -- not a proxy or a screening heuristic.
    Adopted into OECD TG 442C in 2021 (WNT-endorsed peer review and
    validation report, April 2021). It measures the rate constant of
    cysteine-peptide depletion across a concentration/time matrix, expressed
    as log(kmax), and is validated specifically to discriminate GHS
    Subcategory 1A from everything else (1B or non-sensitizer) -- it is
    NOT a hazard-identification method on its own and is only meaningful
    for a chemical already identified as a sensitizer by another method
    (DPRA/KeratinoSens/h-CLAT battery, a defined approach, or LLNA/human data).

    THRESHOLD PROVENANCE (verified against the primary validation dataset,
    not a secondary paraphrase): log(kmax) >= -2.0 -> Subcategory 1A. This
    exact threshold and its balanced accuracy of 85% (against LLNA-based
    GHS classification, n=180 chemicals) come directly from Natsch A,
    Haupt T, Wareing B, Landsiedel R, Kolle SN (2020), "Predictivity of the
    kinetic direct peptide reactivity assay (kDPRA) for sensitizer potency
    assessment and GHS subclassification," ALTEX 37(4):652-664 -- the
    validation study OECD's own peer-review report (ENV/CBC/MONO(2021)13)
    is built on. Independently corroborated by DB-ALM Protocol No. 217 and
    the official OECD TG 442C guideline page. This is the single
    best-evidenced threshold in this entire app: a real accuracy figure
    from a real validation study, not a commonly-cited default.

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


def classify_from_kdpra(log_kmax: float, threshold: float = -2.0) -> str:
    """Verified against Natsch et al. (2020), ALTEX 37(4):652-664, the primary
    kDPRA validation study underlying OECD TG 442C's 2021 kDPRA appendix.
    kDPRA only resolves 1A vs. not-1A -- it cannot distinguish 1B from
    non-sensitizer, so "not 1A" is reported as such, not as "1B"."""
    return "Subcategory 1A (strong sensitizer)" if log_kmax >= threshold else "Not Subcategory 1A (1B or non-sensitizer — kDPRA cannot distinguish further)"


def its_band_score(value: float, bands: list) -> int:
    """bands = [(threshold, points), ...] ascending threshold; returns points for first band value < threshold,
    else max points."""
    for threshold, points in bands:
        if value < threshold:
            return points
    return bands[-1][1]


def its_dpra_points(depletion_pct: float) -> int:
    """OECD 497 ITS DPRA scoring band, verified against the primary ITS table
    (same band edges independently verified in this app's WoE tab)."""
    if depletion_pct >= 42.47:
        return 3
    if depletion_pct >= 22.62:
        return 2
    if depletion_pct >= 6.38:
        return 1
    return 0


def its_hclat_points(mit_ug_ml: float, is_positive: bool) -> int:
    """OECD 497 ITS h-CLAT scoring band, using MIT (Minimum Induction
    Threshold, ug/mL) -- NOT the same metric as the WoE tab's h-CLAT
    %RFI-based hazard call. A negative h-CLAT result scores 0 regardless
    of MIT (MIT is only meaningful for a positive result). Per the primary
    table (OECD GL 497, 25 June 2025, Table 3.1), the score=1 band for a
    positive result runs from >150 up to <=5000 ug/mL."""
    if not is_positive:
        return 0
    if mit_ug_ml <= 10:
        return 3
    if mit_ug_ml <= 150:
        return 2
    return 1  # positive, MIT > 150 and <= 5000 ug/mL per Table 3.1


def its_total_to_category(total_score: int) -> str:
    """OECD 497 ITS total-score-to-GHS-category mapping, verified against
    the primary table (6-7 -> 1A, 2-5 -> 1B, 0-1 -> Not Classified)."""
    if total_score >= 6:
        return "GHS 1A (strong sensitizer)"
    if total_score >= 2:
        return "GHS 1B"
    return "Not Classified (non-sensitizer)"


def render_potency_module():
    st.markdown("#### 🏷️ Potency Categorization (GHS 1A / 1B)")
    st.caption(
        "Three paths: real animal/human point-of-departure data (Path A), the OECD 497 ITS "
        "Defined Approach with a disclosed in-silico substitution (Path B), or the validated "
        "kDPRA single-assay method for confirming Subcategory 1A (Path C)."
    )

    path = st.radio("Which data do you have?", [
        "Path A — LLNA EC3 (or equivalent PoD) available",
        "Path B — OECD 497 ITS score (DPRA + h-CLAT + in silico), no PoD",
        "Path C — kDPRA log(kmax) available (validated 1A vs. not-1A only)",
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

    elif path.startswith("Path B"):
        st.markdown("##### OECD 497 ITS score (DPRA + h-CLAT verified; in silico substituted)")
        with st.expander("📖 What's verified vs. substituted here", expanded=True):
            st.markdown(
                "The scoring framework below (DPRA bands, h-CLAT MIT bands, total-score cutoffs "
                "6-7→1A / 2-5→1B / 0-1→NC) is verified against the primary OECD Guideline No. 497 "
                "ITS table (ENV/CBC/MONO(2021)11), corroborated by a NICEATM/ICCVAM scientist's SCHC "
                "2022 presentation and a Frontiers in Toxicology review. **What's substituted**: the "
                "real ITS uses Derek Nexus or the OECD QSAR Toolbox for the in silico component — "
                "this app has neither, so it offers using this app's own Read-Across tab result "
                "instead. That substitution is not validated to the same standard, and per this "
                "app's own LOO-CV validation, its read-across MCC is moderate (~0.3), not the ~80% "
                "balanced accuracy the real ITS achieves with its intended in silico source. Treat "
                "this Path's overall accuracy as unproven, even though its scoring arithmetic is not."
            )

        st.markdown("###### DPRA")
        dpra_depletion = st.number_input("DPRA mean Cys+Lys depletion (%)", 0.0, 100.0, 20.0, 0.1)
        dpra_pts = its_dpra_points(dpra_depletion)

        st.markdown("###### h-CLAT")
        c1, c2 = st.columns(2)
        hclat_positive = c1.checkbox("h-CLAT positive?", value=True)
        hclat_mit = c2.number_input("h-CLAT MIT (µg/mL) — lower of CD86 EC150 / CD54 EC200",
                                     0.0, 20000.0, 500.0, 1.0,
                                     help="Only used if positive above. Note: µg/mL, NOT the %RFI the WoE tab's h-CLAT uses.")
        hclat_pts = its_hclat_points(hclat_mit, hclat_positive)

        st.markdown("###### In silico prediction")
        insilico_source = st.radio(
            "Source", ["Use this app's Read-Across tab result", "Enter manually"], horizontal=True,
        )
        read_across_result = st.session_state.get("read_across_result")
        if insilico_source.startswith("Use"):
            if read_across_result and "SENSITIZER" in read_across_result.get("call", "") and "NON" not in read_across_result.get("call", ""):
                insilico_positive = True
                st.caption("✅ Read-Across tab currently says SENSITIZER → in silico score = 1")
            elif read_across_result:
                insilico_positive = False
                st.caption("ℹ️ Read-Across tab currently does not say SENSITIZER → in silico score = 0")
            else:
                insilico_positive = False
                st.warning("Read-Across tab hasn't been run yet — defaulting to 0 (negative). Run it first for a real input here.")
        else:
            insilico_positive = st.checkbox("In silico prediction positive?", value=False)
        insilico_pts = 1 if insilico_positive else 0

        total = dpra_pts + hclat_pts + insilico_pts
        st.markdown("##### Component scores")
        st.dataframe(
            {"Component": ["DPRA", "h-CLAT", "In silico"], "Points": [dpra_pts, hclat_pts, insilico_pts]},
            use_container_width=True, hide_index=True,
        )

        category = its_total_to_category(total)
        st.metric("Total ITS score", f"{total} / 7")
        if "1A" in category:
            st.error(f"**{category}**")
        elif "1B" in category:
            st.warning(f"**{category}**")
        else:
            st.success(f"**{category}**")

        st.session_state["potency_result"] = {
            "path": "B", "dpra_points": dpra_pts, "hclat_points": hclat_pts,
            "insilico_points": insilico_pts, "total_score": total, "category": category,
        }
        st.caption(
            "Scoring framework verified against OECD GL 497's primary ITS table. Overall accuracy "
            "of this specific implementation is unproven due to the in-silico substitution — see the "
            "provenance note above."
        )

    else:
        st.markdown("##### kDPRA — validated single-assay 1A determination")
        with st.expander("📖 Threshold provenance — the strongest-evidenced number in this app", expanded=True):
            st.markdown(
                "log(kmax) ≥ **-2.0** → GHS Subcategory 1A, with a real, published **85% balanced "
                "accuracy** (n=180 chemicals, vs. LLNA-based GHS classification). Source: Natsch A, "
                "Haupt T, Wareing B, Landsiedel R, Kolle SN (2020), *Predictivity of the kinetic direct "
                "peptide reactivity assay (kDPRA) for sensitizer potency assessment and GHS "
                "subclassification*, ALTEX 37(4):652-664 — the study OECD's own peer-review report "
                "(ENV/CBC/MONO(2021)13, WNT-endorsed April 2021) is built on, and the basis for "
                "kDPRA's 2021 adoption into OECD TG 442C. Independently corroborated by DB-ALM "
                "Protocol No. 217 and the official OECD TG 442C page. Unlike Paths A and B above, "
                "this isn't a commonly-cited default or an unverified proxy — it's a real accuracy "
                "figure from the primary validation study."
            )
        st.warning(
            "kDPRA only resolves 1A vs. not-1A. It does NOT identify whether a chemical is a "
            "sensitizer in the first place — it's meant to be applied to a chemical already "
            "identified as a sensitizer by another method (the WoE tab, read-across, or LLNA/human "
            "data). Applying it to an untested chemical and getting 'not 1A' does not mean "
            "non-sensitizer; it means 'not confirmed as 1A' — could still be 1B, or could be "
            "untested for hazard entirely."
        )

        woe_result = st.session_state.get("woe_result")
        if woe_result and "Sensitizer" in woe_result.get("classification", "") and "Non" not in woe_result.get("classification", ""):
            st.caption("✅ The WoE tab already reports a sensitizer call for this compound — kDPRA is appropriately scoped here.")
        elif woe_result:
            st.caption("ℹ️ The WoE tab's current result for this compound doesn't confirm sensitizer hazard — check that before relying on a kDPRA result.")

        c1, c2 = st.columns(2)
        log_kmax = c1.number_input("log(kmax)", -6.0, 3.0, -2.5, 0.1,
                                    help="Rate constant of cysteine-peptide depletion from the kDPRA assay, per DB-ALM Protocol No. 217 / OECD TG 442C.")
        threshold = c2.number_input("1A threshold (log kmax ≥ this)", -6.0, 3.0, -2.0, 0.1,
                                     help="Verified against Natsch et al. (2020) — see provenance note above.")

        category = classify_from_kdpra(log_kmax, threshold)
        st.session_state["potency_result"] = {"path": "C", "log_kmax": log_kmax, "threshold": threshold, "category": category}

        if "Subcategory 1A" in category and "Not" not in category:
            st.error(f"**{category}** (log kmax {log_kmax} ≥ threshold {threshold})")
        else:
            st.info(f"**{category}** (log kmax {log_kmax} < threshold {threshold})")
