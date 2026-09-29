"""
Batch / multi-compound screening.

Scope, deliberately limited: this runs Structure analysis (RDKit descriptors
+ structural alerts) and Read-Across against a chosen reference set for
every compound in an uploaded list, plus GHS Path-A potency categorization
if an ec3_pct column is supplied. It does NOT batch-run the assay-based WoE
tab -- that needs several numeric columns per assay (depletion %, EC1.5,
IC30, CD86/CD54 RFI, viability) that are realistically only available for a
handful of compounds actively being tested, not a triage-sized batch list.

This is a triage tool: run it BEFORE deciding which compounds are worth
actual assay time, not a replacement for the single-compound tabs. Each
compound gets the same per-compound disclaimers as the rest of this app
(structural alerts are heuristic flags; read-across depends entirely on
reference-set relevance).

A bad SMILES in one row must not crash the whole batch -- each row is
processed independently and errors are reported per-row, not raised.
"""

import io
import streamlit as st
import pandas as pd

from modules.structure import analyze_smiles
from modules.potency import classify_from_llna
from modules.read_across import compute_read_across, _load_combined_reference

MAX_BATCH_ROWS = 500


def run_batch_screening(df: pd.DataFrame, ref_df: pd.DataFrame, k: int = 3,
                         domain_cutoff: float = 0.4) -> pd.DataFrame:
    """Pure function: takes a compound list and a reference set, returns a
    results DataFrame. One row in, one row out, always -- a failure on one
    compound becomes an error message in that row, never a raised exception
    that kills the rest of the batch."""
    cols = {c.strip().lower(): c for c in df.columns}
    if "name" not in cols or "smiles" not in cols:
        raise ValueError("Input must have at least 'name' and 'smiles' columns.")
    has_ec3 = "ec3_pct" in cols

    results = []
    for _, row in df.iterrows():
        name = row[cols["name"]]
        smiles = row[cols["smiles"]]
        out = {"name": name, "smiles": smiles, "error": None}

        if pd.isna(smiles) or not str(smiles).strip():
            out["error"] = "Empty SMILES"
            results.append(out)
            continue

        mol, descriptors, alerts, err, was_standardized = analyze_smiles(str(smiles))
        if err:
            out["error"] = err
            results.append(out)
            continue

        out["salt_stripped"] = was_standardized
        out["mw"] = descriptors.get("Molecular Weight (g/mol)")
        out["logp"] = descriptors.get("LogP (Crippen)")
        flagged = [a for a, present in alerts.items() if present]
        out["n_alerts"] = len(flagged)
        out["alerts_flagged"] = "; ".join(flagged) if flagged else "None"

        ra_result, ra_err = compute_read_across(str(smiles), ref_df, k=k, domain_cutoff=domain_cutoff)
        if ra_err:
            out["read_across_call"] = f"Error: {ra_err}"
            out["nearest_neighbor"] = None
            out["nearest_neighbor_similarity"] = None
        else:
            out["read_across_call"] = ra_result["call"]
            top = ra_result["ranked"][0] if ra_result["ranked"] else None
            out["nearest_neighbor"] = top["name"] if top else None
            out["nearest_neighbor_similarity"] = round(top["similarity"], 3) if top else None

        if has_ec3:
            ec3_val = row[cols["ec3_pct"]]
            if pd.notna(ec3_val):
                try:
                    out["potency_category"] = classify_from_llna(float(ec3_val))
                except (TypeError, ValueError):
                    out["potency_category"] = "Invalid EC3"
            else:
                out["potency_category"] = None

        results.append(out)

    return pd.DataFrame(results)


def render_batch_screening_module():
    st.markdown("#### 📋 Batch / Multi-Compound Screening")
    st.warning(
        "Triage only: this runs Structure + Read-Across (and GHS Path-A potency if you "
        "supply an ec3_pct column) across a whole list at once. It does NOT batch-run the "
        "assay-based WoE tab, which needs per-assay numeric data unlikely to exist for a "
        "batch-sized list. Use this to prioritize which compounds are worth actual assay "
        "time — not as a replacement for the single-compound tabs."
    )

    uploaded = st.file_uploader(
        "Upload compound list (CSV with columns: name, smiles, optional ec3_pct)",
        type=["csv"],
    )
    if uploaded is None:
        st.info("Upload a CSV to run a batch screen.")
        return

    try:
        df = pd.read_csv(uploaded)
    except Exception as e:
        st.error(f"Could not read that file as a CSV: {e}")
        return

    if len(df) > MAX_BATCH_ROWS:
        st.error(f"This file has {len(df)} rows; batch screening is capped at {MAX_BATCH_ROWS} "
                 "per run to keep this responsive. Split it into smaller files.")
        return

    st.caption(f"{len(df)} compounds loaded.")
    st.dataframe(df.head(10), use_container_width=True, hide_index=True)
    if len(df) > 10:
        st.caption(f"(showing first 10 of {len(df)} rows)")

    st.markdown("##### Reference set for read-across")
    ref_choice = st.radio(
        "Which reference set?",
        ["Combined real reference set (54 compounds)", "Upload a different reference CSV"],
        horizontal=True,
    )
    if ref_choice.startswith("Combined"):
        ref_df = _load_combined_reference()
    else:
        ref_upload = st.file_uploader("Reference CSV (columns: name, smiles, sensitizer)", type=["csv"], key="batch_ref")
        if ref_upload is None:
            st.info("Upload a reference CSV to continue.")
            return
        try:
            ref_df = pd.read_csv(ref_upload)
        except Exception as e:
            st.error(f"Could not read reference file: {e}")
            return

    c1, c2 = st.columns(2)
    k = c1.slider("Nearest neighbors (k)", 1, 10, 3, key="batch_k")
    domain_cutoff = c2.slider("Applicability domain cutoff (Tanimoto)", 0.0, 1.0, 0.4, 0.05, key="batch_cutoff")

    if st.button("▶ Run batch screen", type="primary"):
        with st.spinner(f"Screening {len(df)} compounds..."):
            try:
                results = run_batch_screening(df, ref_df, k=k, domain_cutoff=domain_cutoff)
            except ValueError as e:
                st.error(str(e))
                return

        st.session_state["batch_screening_result"] = results

        n_errors = results["error"].notna().sum()
        st.success(f"Screened {len(results)} compounds ({n_errors} with errors).")
        st.dataframe(results, use_container_width=True, hide_index=True)

        csv_bytes = results.to_csv(index=False).encode("utf-8")
        st.download_button(
            "📥 Download results CSV",
            data=csv_bytes,
            file_name="batch_screening_results.csv",
            mime="text/csv",
            use_container_width=True,
        )
