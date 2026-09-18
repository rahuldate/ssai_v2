"""
Read-across (analog-based) prediction module.

This is NOT a claim of "validated QSAR toolbox equivalence." It is the same
underlying technique the OECD QSAR Toolbox itself uses for read-across:

    1. Compute a molecular fingerprint for the target compound (RDKit Morgan
       fingerprint, radius 2, 2048 bits — the standard ECFP4-equivalent).
    2. Compute Tanimoto similarity against every compound in a REFERENCE SET
       that the user supplies (CSV: name, smiles, sensitizer, potency_category,
       optional ec3_pct / notes).
    3. Rank nearest neighbors. If the closest neighbor(s) are below a
       similarity threshold, the target is flagged OUT OF the applicability
       domain of that reference set and no read-across call is made.
    4. The "prediction" is a simple, auditable majority vote among the
       k nearest in-domain neighbors — never a black-box probability.

Validity of the output depends entirely on the quality and relevance of the
reference set you provide. The bundled demo set is a handful of
well-documented, unambiguous reference chemicals for illustration only —
replace it with your own curated, source-cited dataset before relying on
this for anything real.
"""

import io
import os
import streamlit as st
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NICEATM_REFERENCE_PATH = os.path.join(_THIS_DIR, "..", "data", "niceatm_reference_18.csv")

# Tiny illustrative demo set (6 compounds) — purely for a zero-setup smoke test.
DEMO_REFERENCE_CSV = """name,smiles,sensitizer,potency_category,ec3_pct,notes
DNCB (positive control),O=[N+]([O-])c1ccc(Cl)c(c1)[N+]([O-])=O,yes,1A,0.05,Classic strong LLNA positive control
Cinnamaldehyde,O=CC=CC1=CC=CC=C1,yes,1B,7.9,Moderate sensitizer; alpha-beta-unsaturated aldehyde
Eugenol,COc1cc(CC=C)ccc1O,yes,1B,10.0,Weak-to-moderate fragrance sensitizer
Isopropanol,CC(O)C,no,NC,,Common negative control solvent
Lactic acid,CC(O)C(=O)O,no,NC,,Common negative control (non-sensitizer)
Glycerol,C(C(CO)O)O,no,NC,,Common negative control
"""


def _fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def compute_read_across(target_smiles: str, ref_df: pd.DataFrame, k: int = 3, domain_cutoff: float = 0.4):
    if not RDKIT_AVAILABLE:
        return None, "RDKit not available."

    target_mol = Chem.MolFromSmiles(target_smiles)
    if target_mol is None:
        return None, "Could not parse target SMILES."

    target_fp = _fingerprint(target_mol)

    rows = []
    for _, r in ref_df.iterrows():
        ref_mol = Chem.MolFromSmiles(str(r["smiles"]))
        if ref_mol is None:
            continue
        ref_fp = _fingerprint(ref_mol)
        sim = DataStructs.TanimotoSimilarity(target_fp, ref_fp)
        rows.append({
            "name": r["name"],
            "similarity": sim,
            "sensitizer": r.get("sensitizer", ""),
            "potency_category": r.get("potency_category", ""),
            "ec3_pct": r.get("ec3_pct", ""),
        })

    if not rows:
        return None, "No valid reference structures could be parsed."

    ranked = sorted(rows, key=lambda x: x["similarity"], reverse=True)
    top_k = ranked[:k]
    in_domain = [n for n in top_k if n["similarity"] >= domain_cutoff]

    result = {
        "ranked": ranked,
        "top_k": top_k,
        "in_domain": in_domain,
        "domain_cutoff": domain_cutoff,
    }

    if not in_domain:
        result["call"] = "OUT OF APPLICABILITY DOMAIN — no reference neighbor met the similarity cutoff"
        return result, None

    votes = [n["sensitizer"] for n in in_domain]
    pos_votes = sum(1 for v in votes if str(v).lower() == "yes")
    neg_votes = len(votes) - pos_votes

    if pos_votes > neg_votes:
        call = f"Read-across suggests SENSITIZER (majority of {len(in_domain)} in-domain neighbors positive)"
    elif neg_votes > pos_votes:
        call = f"Read-across suggests NON-SENSITIZER (majority of {len(in_domain)} in-domain neighbors negative)"
    else:
        call = "Read-across INCONCLUSIVE (neighbor vote tied)"

    result["call"] = call
    return result, None


def render_read_across_module():
    st.markdown("#### 🔬 Read-Across (Analog-Based) Prediction")
    st.warning(
        "This is a similarity/read-across engine, not a validated QSAR model. Two bundled "
        "reference sets are offered: a small set of REAL, published, source-cited compounds "
        "(NICEATM/ICCVAM) that is heavily weighted toward hard/edge cases, and a tiny "
        "illustrative demo set for smoke-testing. Neither is large or representative enough "
        "to be a real applicability domain on its own — upload your own curated set for "
        "production use."
    )

    if not RDKIT_AVAILABLE:
        st.error("RDKit is not installed.")
        return

    st.markdown("##### Reference set")
    ref_choice = st.radio(
        "Which reference set do you want to use?",
        [
            "NICEATM/ICCVAM 18-compound reference set (real, published, cited)",
            "Tiny illustrative demo set (6 compounds, for a quick smoke test)",
            "Upload my own reference CSV",
        ],
    )

    if ref_choice.startswith("NICEATM"):
        ref_df = pd.read_csv(NICEATM_REFERENCE_PATH)
        with st.expander("📖 Source & important caveats for this reference set", expanded=True):
            st.markdown(
                "- **Source:** Strickland J, Choksi N, Allen D, Casey W. *In Silico Predictions "
                "of Skin Sensitization Using OECD QSAR Toolbox.* NICEATM SOT 2015 Poster "
                "(U.S. NIEHS/NTP, public domain). Compounds drawn from that poster's Tables 3 & 4, "
                "with CAS numbers, LLNA results, and (where positive) LLNA EC3 values as reported there.\n"
                "- **Every SMILES here has been cross-checked**: the RDKit-computed molecular weight "
                "for each compound matches the source table's reported molecular weight.\n"
                "- **Important limitation:** these 18 substances were specifically the ones a structural "
                "read-across method *misclassified* in that study — i.e. this is a curated set of "
                "**hard/discordant cases**, not a representative random sample. It is genuinely useful "
                "for stress-testing your applicability domain logic, but a small, edge-case-heavy set "
                "like this is not a substitute for a large, representative validated database.\n"
                "- Potency categories (1A/1B) shown here were derived by applying this tool's own "
                "EC3 ≤ 2% → 1A rule to the source's reported EC3 values — they are not independently "
                "confirmed GHS classifications from a regulatory body."
            )
    elif ref_choice.startswith("Tiny"):
        ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
        st.caption(f"Using tiny illustrative demo set — {len(ref_df)} compounds, for smoke-testing only.")
    else:
        uploaded = st.file_uploader(
            "Upload reference CSV (columns: name, smiles, sensitizer [yes/no], potency_category [1A/1B/NC], ec3_pct optional)",
            type=["csv"],
        )
        if uploaded is None:
            st.info("Upload a CSV to continue, or switch to one of the bundled sets above.")
            return
        ref_df = pd.read_csv(uploaded)

    required_cols = {"name", "smiles", "sensitizer"}
    if not required_cols.issubset(set(c.lower() for c in ref_df.columns)):
        st.error(f"Reference CSV must include at least these columns: {required_cols}")
        return

    st.dataframe(ref_df, use_container_width=True, hide_index=True)

    st.markdown("##### Target compound")
    target_smiles = st.text_input("Target SMILES", value=st.session_state.get("compound_smiles", "O=CC=CC1=CC=CC=C1"))
    k = st.slider("Number of nearest neighbors (k)", 1, 10, 3)
    domain_cutoff = st.slider(
        "Applicability domain similarity cutoff (Tanimoto)", 0.0, 1.0, 0.4, 0.05,
        help="Below this similarity, no reference compound is considered close enough to read across from."
    )

    if not target_smiles.strip():
        return

    result, error = compute_read_across(target_smiles, ref_df, k=k, domain_cutoff=domain_cutoff)
    if error:
        st.error(error)
        return

    st.markdown("##### Nearest reference neighbors")
    st.dataframe(
        pd.DataFrame(result["top_k"])[["name", "similarity", "sensitizer", "potency_category", "ec3_pct"]]
        .round({"similarity": 3}),
        use_container_width=True,
        hide_index=True,
    )

    if "OUT OF" in result["call"]:
        st.error(result["call"])
    elif "SENSITIZER" in result["call"] and "NON" not in result["call"]:
        st.warning(result["call"])
    else:
        st.info(result["call"])

    st.session_state["read_across_result"] = result
    st.caption(
        "Read-across confidence depends on reference-set relevance and coverage. "
        "Always confirm mechanistic plausibility (shared structural alerts) before accepting a call."
    )
