"""
Integrated evidence synthesis (IATA-style weight-of-evidence summary).

This module does NOT run any new assay or computation. It reasons ONLY over
what other tabs have actually produced in this session (structure alerts,
WoE result, read-across result, potency result) and:

  1. Reports which evidence streams are actually available -- never implies
     a stream ran when it didn't.
  2. Cross-checks the hazard CALLS (sensitizer / non-sensitizer) from WoE and
     read-across against each other. Agreement is reported as agreement;
     disagreement is reported as disagreement, explicitly and by name -- this
     is the single most useful thing a synthesis step can surface, and it is
     never smoothed over into a false consensus.
  3. Uses structural alerts as CONTEXT, not a vote: a positive hazard call
     with no matching structural alert, or a negative call despite an alert
     being flagged, is noted as worth a closer look -- not silently ignored,
     and not treated as if it changes the call itself.
  4. Reports a QUALITATIVE evidence-strength label (e.g. "single stream
     only", "two streams, concordant") based purely on how many streams ran
     and whether they agree. This is NOT a numeric confidence score or
     probability -- there is no validated way to turn "2 out of 2 agree" into
     a percentage with this data, and presenting one would manufacture false
     precision exactly like the fabricated MCC in the original app.

The output is a reasoning summary for an expert to review, not a
determination on its own.
"""

import streamlit as st


def synthesize(structure_alerts: dict = None, woe_result: dict = None,
                read_across_result: dict = None, potency_result: dict = None) -> dict:
    """Pure function: takes the same result dicts other tabs store in
    st.session_state and returns a structured synthesis. Kept separate from
    Streamlit so it's directly unit-testable."""

    streams_run = []
    if structure_alerts is not None:
        streams_run.append("structure")
    if woe_result is not None:
        streams_run.append("woe")
    if read_across_result is not None:
        streams_run.append("read_across")
    if potency_result is not None:
        streams_run.append("potency")

    # --- Extract a hazard call (yes/no/None) from each stream that gives one ---
    hazard_votes = {}

    if woe_result is not None:
        cls = woe_result.get("classification", "")
        if "Non-sensitizer" in cls:
            hazard_votes["woe"] = "no"
        elif "Sensitizer" in cls:
            hazard_votes["woe"] = "yes"
        else:
            hazard_votes["woe"] = None  # inconclusive

    if read_across_result is not None:
        call = read_across_result.get("call", "")
        if "NON-SENSITIZER" in call:
            hazard_votes["read_across"] = "no"
        elif "SENSITIZER" in call:
            hazard_votes["read_across"] = "yes"
        else:
            hazard_votes["read_across"] = None  # inconclusive / out of domain

    usable_votes = {k: v for k, v in hazard_votes.items() if v is not None}
    flags = []

    if len(usable_votes) == 0:
        if len(streams_run) == 0:
            evidence_strength = "no evidence"
            overall = "No evidence streams have been run yet."
        else:
            evidence_strength = "no usable hazard call"
            overall = ("Evidence streams were run, but none produced a usable hazard call "
                       "(inconclusive or out of applicability domain).")
    elif len(usable_votes) == 1:
        stream, vote = next(iter(usable_votes.items()))
        evidence_strength = "single stream only"
        label = "SENSITIZER" if vote == "yes" else "NON-SENSITIZER"
        overall = f"Only one usable hazard call available ({stream}): {label}. Treat as preliminary."
    else:
        distinct = set(usable_votes.values())
        if len(distinct) == 1:
            label = "SENSITIZER" if "yes" in distinct else "NON-SENSITIZER"
            evidence_strength = f"{len(usable_votes)} streams, concordant"
            overall = f"{label} — {len(usable_votes)} of {len(usable_votes)} usable evidence streams agree ({', '.join(usable_votes.keys())})."
        else:
            evidence_strength = f"{len(usable_votes)} streams, DISCORDANT"
            detail = "; ".join(f"{k}={'SENSITIZER' if v == 'yes' else 'NON-SENSITIZER'}" for k, v in usable_votes.items())
            overall = f"DISCORDANT evidence — expert review required. {detail}."
            flags.append("Hazard-call streams disagree with each other; do not average or silently pick one.")

    # --- Structural alert context (never a vote, always contextual) ---
    structural_note = None
    if structure_alerts is not None:
        any_alert = any(structure_alerts.values())
        flagged = [k for k, v in structure_alerts.items() if v]
        if any_alert and usable_votes and all(v == "no" for v in usable_votes.values()):
            structural_note = (
                f"Structural alert(s) present ({', '.join(flagged)}) despite a negative hazard call — "
                "worth a closer mechanistic look; a structural alert alone doesn't override an assay-based "
                "or read-across call, but a real reactive group with no positive evidence is a combination "
                "worth double-checking rather than dismissing."
            )
        elif not any_alert and usable_votes and all(v == "yes" for v in usable_votes.values()):
            structural_note = (
                "Positive hazard call(s) with none of the 8 checked structural alerts triggered — the "
                "reactive mechanism (if any) may fall outside this tool's simplified alert set, or the "
                "positive call may be driven by structural similarity rather than an identified mechanism. "
                "Not a reason to doubt the call, just a reason not to expect a mechanistic explanation from "
                "the Structure tab alone."
            )
        elif any_alert:
            structural_note = f"Structural alert(s) present: {', '.join(flagged)}."
        else:
            structural_note = "No structural alerts triggered from the 8 checked patterns."

    # --- Potency context ---
    potency_note = None
    if potency_result is not None:
        cat = potency_result.get("category")
        path = potency_result.get("path")
        if path == "A":
            potency_note = f"Potency (Path A, LLNA EC3): {cat}."
        elif path == "B":
            potency_note = f"Potency (Path B, OECD 497 ITS score — in silico component substituted): {cat}."
        elif path == "C":
            potency_note = f"Potency (Path C, kDPRA — validated 1A vs. not-1A only): {cat}."
        # A "1A" call means the category string names 1A AND doesn't negate it
        # (kDPRA's negative result is literally "Not Subcategory 1A...", which
        # would false-positive on a naive substring check for "1A").
        is_1a_call = bool(cat) and "1A" in str(cat) and "Not" not in str(cat)
        if usable_votes and all(v == "no" for v in usable_votes.values()) and is_1a_call:
            flags.append("Potency tab reports a strong-sensitizer category (1A) while hazard-call "
                         "streams say non-sensitizer — this is a direct contradiction worth resolving "
                         "before relying on either.")

    return {
        "streams_run": streams_run,
        "hazard_votes": hazard_votes,
        "usable_votes": usable_votes,
        "evidence_strength": evidence_strength,
        "overall": overall,
        "structural_note": structural_note,
        "potency_note": potency_note,
        "flags": flags,
    }


def render_synthesis_module():
    st.markdown("#### 🧩 Integrated Evidence Synthesis")
    st.caption(
        "Reasons over whatever the other tabs have actually computed this session — runs nothing new. "
        "Reports agreement and disagreement between evidence streams explicitly; never averages or "
        "silently resolves a disagreement into a false consensus, and never assigns a numeric confidence "
        "score, since there's no validated way to turn 'streams agree' into a percentage with this data."
    )

    structure_alerts = st.session_state.get("structure_alerts")
    woe_result = st.session_state.get("woe_result")
    read_across_result = st.session_state.get("read_across_result")
    potency_result = st.session_state.get("potency_result")

    result = synthesize(structure_alerts, woe_result, read_across_result, potency_result)
    st.session_state["synthesis_result"] = result

    st.markdown("##### Which evidence streams have been run")
    all_streams = ["structure", "woe", "read_across", "potency"]
    st.dataframe(
        {
            "Stream": all_streams,
            "Status": ["✅ Run" if s in result["streams_run"] else "— Not run" for s in all_streams],
        },
        use_container_width=True, hide_index=True,
    )

    st.markdown("##### Overall")
    if "DISCORDANT" in result["evidence_strength"]:
        st.error(result["overall"])
    elif "concordant" in result["evidence_strength"]:
        st.success(result["overall"])
    elif result["evidence_strength"] in ("no evidence", "no usable hazard call"):
        st.info(result["overall"])
    else:
        st.warning(result["overall"])

    st.caption(f"Evidence strength: {result['evidence_strength']}")

    if result["flags"]:
        st.markdown("##### Flags")
        for f in result["flags"]:
            st.warning(f)

    if result["structural_note"]:
        st.markdown("##### Structural context")
        st.caption(result["structural_note"])

    if result["potency_note"]:
        st.markdown("##### Potency context")
        st.caption(result["potency_note"])
