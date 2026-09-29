"""
Dermal sensitization Quantitative Risk Assessment (QRA) module.

Implements the general structure of the Api et al. (2008) dermal
sensitization QRA framework used by RIFM/IFRA for fragrance ingredients:

    AEL (Acceptable Exposure Level) = NESIL / SAF
    SAF (Sensitization Assessment Factor) = product of user-set sub-factors
    CEL (Consumer Exposure Level)   = (amount per use x concentration) / skin area  [x frequency if needed]
    Margin of Safety (MoS)          = AEL / CEL   (MoS >= 1 => acceptable)

All numbers are either entered directly by the user (NESIL/PoD, exposure
parameters) or are simple, visible multiplications/divisions of those numbers.
The individual SAF sub-factor DEFAULT values below are common textbook
starting points (interspecies=10, inter-individual=5) — they are editable and
should be set per your own SOP / the current RIFM-IFRA QRA methodology
documents, not treated as fixed regulatory constants.
"""

import streamlit as st


def compute_qra(nesil_ug_cm2: float, saf_factors: dict, amount_g: float,
                 concentration_pct: float, skin_area_cm2: float, uses_per_day: float):
    saf_total = 1.0
    for v in saf_factors.values():
        saf_total *= v

    ael = nesil_ug_cm2 / saf_total if saf_total > 0 else float("inf")

    # amount applied (g) * fraction concentration -> mass of ingredient (g) -> µg
    ingredient_mass_ug = amount_g * (concentration_pct / 100.0) * 1_000_000
    cel_per_use = ingredient_mass_ug / skin_area_cm2 if skin_area_cm2 > 0 else float("inf")
    cel_per_day = cel_per_use * uses_per_day

    mos = ael / cel_per_day if cel_per_day > 0 else float("inf")

    return {
        "SAF_total": saf_total,
        "AEL_ug_cm2_day": ael,
        "CEL_per_use_ug_cm2": cel_per_use,
        "CEL_per_day_ug_cm2": cel_per_day,
        "MoS": mos,
    }


def render_qra_module():
    st.markdown("#### 🛡️ Dermal Sensitization QRA (Api et al. 2008 framework)")
    st.caption(
        "Transparent AEL / CEL / Margin-of-Safety calculator. Every output is a "
        "direct arithmetic function of the inputs you provide below — nothing is "
        "pre-filled with compound-specific results."
    )

    st.markdown("##### 1. Point of Departure")
    nesil = st.number_input(
        "NESIL — No Expected Sensitization Induction Level (µg/cm²)",
        min_value=0.0, value=100.0, step=1.0,
        help="Derived from LLNA EC3 (animal) or human NOEL/HRIPT data for this substance."
    )
    pod_source = st.selectbox("Point-of-departure source", ["Human NOEL / HRIPT", "LLNA EC3 (animal)", "Other in-vitro/QSAR estimate"])

    st.markdown("##### 2. Sensitization Assessment Factor (SAF) — editable sub-factors")
    c1, c2 = st.columns(2)
    with c1:
        interspecies = st.number_input("Interspecies factor", 1.0, 100.0,
                                        10.0 if pod_source == "LLNA EC3 (animal)" else 1.0, 0.5)
        interindividual = st.number_input("Inter-individual variability factor", 1.0, 50.0, 5.0, 0.5)
    with c2:
        matrix_factor = st.number_input("Matrix / vehicle factor", 1.0, 20.0, 1.0, 0.5,
                                         help="Accounts for formulation type, e.g. rinse-off vs leave-on, per your SOP.")
        use_factor = st.number_input("Use/exposure condition factor", 1.0, 20.0, 1.0, 0.5,
                                      help="E.g. occlusion, skin condition, or other exposure modifiers per your SOP.")

    saf_factors = {
        "interspecies": interspecies,
        "interindividual": interindividual,
        "matrix": matrix_factor,
        "use": use_factor,
    }

    st.markdown("##### 3. Consumer Exposure Parameters")
    c1, c2, c3, c4 = st.columns(4)
    amount_g = c1.number_input("Amount applied per use (g)", 0.0, 1000.0, 1.0, 0.1)
    conc_pct = c2.number_input("Ingredient concentration in product (%)", 0.0, 100.0, 0.1, 0.01)
    skin_area = c3.number_input("Skin area exposed (cm²)", 1.0, 20000.0, 500.0, 1.0)
    uses_per_day = c4.number_input("Applications per day", 0.0, 20.0, 1.0, 1.0)

    result = compute_qra(nesil, saf_factors, amount_g, conc_pct, skin_area, uses_per_day)
    st.session_state["qra_result"] = {**result, "nesil": nesil, "saf_factors": saf_factors,
                                       "amount_g": amount_g, "conc_pct": conc_pct,
                                       "skin_area_cm2": skin_area, "uses_per_day": uses_per_day,
                                       "pod_source": pod_source}

    st.markdown("---")
    st.markdown("##### Results")
    c1, c2, c3 = st.columns(3)
    c1.metric("SAF (total)", f"{result['SAF_total']:.1f}x")
    c2.metric("AEL", f"{result['AEL_ug_cm2_day']:.3f} µg/cm²/day")
    c3.metric("Consumer Exposure (CEL)", f"{result['CEL_per_day_ug_cm2']:.3f} µg/cm²/day")

    mos = result["MoS"]
    if mos >= 1:
        st.success(f"**Margin of Safety = {mos:.2f}** (≥ 1 → exposure is below the acceptable level, subject to expert review)")
    else:
        st.error(f"**Margin of Safety = {mos:.2f}** (< 1 → exposure exceeds the acceptable level; flag for review/reformulation)")

    st.caption(
        "This calculator reproduces the arithmetic structure of the published dermal "
        "sensitization QRA approach. It does not replace the full RIFM/IFRA QRA2 "
        "methodology documents, which specify additional considerations (e.g. "
        "aggregate exposure across multiple product types) not modeled here."
    )
