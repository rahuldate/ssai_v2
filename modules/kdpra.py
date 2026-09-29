"""
kDPRA (kinetic Direct Peptide Reactivity Assay) — GHS Subcategory 1A screen.

kDPRA is Appendix III of OECD Test Guideline 442C (adopted 2021). It is
fundamentally different in scope from the DPRA/ADRA appendices of the same
guideline, and that difference is enforced here, not just mentioned:

    Per the OECD TG 442C text itself: "the kDPRA... allows discrimination of
    UN GHS subcategory 1A skin sensitisers from those not categorised as
    subcategory 1A (non-subcategory 1A) i.e., subcategory 1B or no category
    (1) but does not allow to distinguish sensitisers (Category 1) from
    non-sensitisers."

In plain terms: kDPRA answers "is this a STRONG (1A) sensitizer?" -- it does
NOT answer "is this a sensitizer at all?". Using it on a compound with no
prior positive hazard call from somewhere else (the WoE tab, LLNA, human
data) would silently misuse the method. This module requires the user to
confirm hazard has already been established before it will return a real
category; otherwise it explicitly returns "not applicable" rather than a
number that looks like an answer.

THRESHOLD PROVENANCE (verified against primary + peer-reviewed sources):
    - log(kmax) >= -2.0  ->  GHS Subcategory 1A
    - log(kmax) <  -2.0  ->  non-Subcategory 1A (1B or no category)
  Confirmed independently by: the OECD TG 442C guideline text (Appendix III,
  oecd.org); an OECD/ECHA-hosted guidance document stating the same cutoff
  in the same units; and Wareing, Natsch et al. 2020 (ALTEX 37(4), 652-664),
  the peer-reviewed validation study reporting log kmax = -2 as the
  threshold achieving 85% balanced accuracy against LLNA-based GHS 1A calls
  across a 180-chemical database. That 85% figure is the validated method's
  own reported performance -- not something re-derived in this app -- and is
  surfaced in the UI so the person using it sees the real, published
  reliability of the number they're about to rely on.
"""

import streamlit as st


def classify_kdpra(log_kmax: float, already_identified_sensitizer: bool, cutoff: float = -2.0) -> str:
    """Returns a category string. Deliberately returns a non-category message
    (not a guess) when the precondition isn't met."""
    if not already_identified_sensitizer:
        return ("Not applicable — hazard not yet confirmed. kDPRA sub-categorizes potency "
                 "among already-identified sensitizers; it cannot tell you whether a compound "
                 "is a sensitizer at all. Establish hazard first (WoE tab, LLNA, or human data).")
    if log_kmax >= cutoff:
        return "Subcategory 1A (strong sensitizer)"
    return "Non-Subcategory 1A (1B or no category)"


def render_kdpra_module():
    st.markdown("#### ⚡ kDPRA (Kinetic DPRA) — GHS Subcategory 1A Screen")
    st.warning(
        "Scope limit, enforced not just stated: kDPRA sub-categorizes POTENCY among compounds "
        "already identified as sensitizers — it cannot distinguish sensitizer from non-sensitizer "
        "on its own (per OECD TG 442C's own text). Confirm hazard elsewhere first (WoE tab, LLNA, "
        "or human data) before using this result."
    )

    already_sensitizer = st.checkbox(
        "This compound has already been identified as a sensitizer (Category 1) by another "
        "method (WoE tab, LLNA, human data, etc.)",
        value=False,
    )

    c1, c2 = st.columns(2)
    kmax_input_mode = c1.radio("Input as", ["log(kmax) directly", "raw kmax (M⁻¹ min⁻¹)"], horizontal=True)
    if kmax_input_mode.startswith("log"):
        log_kmax = c2.number_input("log(kmax)", -6.0, 3.0, -2.5, 0.1)
    else:
        raw_kmax = c2.number_input("kmax (M⁻¹ min⁻¹)", 0.0000001, 1000.0, 0.01, format="%.7f")
        import math
        log_kmax = math.log10(raw_kmax) if raw_kmax > 0 else float("-inf")
        st.caption(f"log(kmax) = {log_kmax:.2f}")

    cutoff = st.number_input(
        "Subcategory 1A cutoff (log kmax, ≥ = 1A)", -6.0, 3.0, -2.0, 0.1,
        help="Verified against OECD TG 442C Appendix III and Wareing/Natsch et al. 2020 (ALTEX) — "
             "see module docstring for the full citation trail."
    )

    category = classify_kdpra(log_kmax, already_sensitizer, cutoff)
    st.session_state["kdpra_result"] = {
        "path": "kDPRA", "log_kmax": log_kmax, "cutoff": cutoff,
        "already_identified_sensitizer": already_sensitizer, "category": category,
    }
    # Also populate potency_result so this flows into the Integrated Summary / dossier
    # the same way Path A and Path B already do, without duplicating that plumbing.
    st.session_state["potency_result"] = {
        "path": "kDPRA", "log_kmax": log_kmax, "cutoff": cutoff, "category": category,
    }

    if not already_sensitizer:
        st.info(category)
    elif "Subcategory 1A" in category and "Non" not in category:
        st.error(f"**{category}** (log kmax {log_kmax:.2f} ≥ cutoff {cutoff})")
    else:
        st.warning(f"**{category}** (log kmax {log_kmax:.2f} < cutoff {cutoff})")

    st.caption(
        "Published performance of this threshold: 85% balanced accuracy for predicting LLNA-based "
        "GHS 1A calls, across a 180-chemical validation database (Wareing, Natsch et al. 2020, "
        "ALTEX 37(4), 652-664) — a real, reported figure for the validated method itself, not a "
        "claim about this app's implementation of it."
    )
