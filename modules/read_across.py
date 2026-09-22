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
    from rdkit.Chem import rdFingerprintGenerator
    _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NICEATM_REFERENCE_PATH = os.path.join(_THIS_DIR, "..", "data", "niceatm_reference_18.csv")
ECVAM_REFERENCE_PATH = os.path.join(_THIS_DIR, "..", "data", "ecvam_dpra_hclat_24.csv")
OECD_PROFICIENCY_PATH = os.path.join(_THIS_DIR, "..", "data", "oecd_tg442d_proficiency_6.csv")

# Tiny illustrative demo set (6 compounds) — purely for a zero-setup smoke test.
DEMO_REFERENCE_CSV = """name,smiles,sensitizer,potency_category,ec3_pct,notes
DNCB (positive control),O=[N+]([O-])c1ccc(Cl)c(c1)[N+]([O-])=O,yes,1A,0.05,Classic strong LLNA positive control
Cinnamaldehyde,O=CC=CC1=CC=CC=C1,yes,1B,7.9,Moderate sensitizer; alpha-beta-unsaturated aldehyde
Eugenol,COc1cc(CC=C)ccc1O,yes,1B,10.0,Weak-to-moderate fragrance sensitizer
Isopropanol,CC(O)C,no,NC,,Common negative control solvent
Lactic acid,CC(O)C(=O)O,no,NC,,Common negative control (non-sensitizer)
Glycerol,C(C(CO)O)O,no,NC,,Common negative control
"""


@st.cache_data
def _load_csv(path: str) -> pd.DataFrame:
    """Cached so repeated Streamlit reruns (which happen on every widget
    interaction) don't re-read and re-parse the same bundled CSV every time."""
    return pd.read_csv(path)


@st.cache_data
def _load_combined_reference() -> pd.DataFrame:
    return pd.concat([
        _load_csv(NICEATM_REFERENCE_PATH),
        _load_csv(ECVAM_REFERENCE_PATH),
        _load_csv(OECD_PROFICIENCY_PATH),
    ], ignore_index=True)


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Lowercase and strip column names so a CSV with 'SMILES' or ' Sensitizer '
    still matches the columns this module indexes by. Without this, the
    required-column check (which lowercases for comparison) can pass while
    the actual lookups below still use the literal lowercase names and crash."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df


def _fingerprint(mol):
    # Uses the current, non-deprecated MorganGenerator API (GetMorganFingerprintAsBitVect
    # is being phased out by RDKit and floods logs with deprecation warnings).
    return _MORGAN_GEN.GetFingerprint(mol)


def compute_read_across(target_smiles: str, ref_df: pd.DataFrame, k: int = 3, domain_cutoff: float = 0.4):
    if not RDKIT_AVAILABLE:
        return None, "RDKit not available."

    target_mol = Chem.MolFromSmiles(target_smiles)
    if target_mol is None:
        return None, "Could not parse target SMILES."

    ref_df = _normalize_columns(ref_df)
    if not {"name", "smiles", "sensitizer"}.issubset(ref_df.columns):
        return None, "Reference data is missing one of the required columns: name, smiles, sensitizer."

    target_fp = _fingerprint(target_mol)

    rows = []
    for _, r in ref_df.iterrows():
        smiles_val = r.get("smiles")
        if pd.isna(smiles_val):
            continue
        ref_mol = Chem.MolFromSmiles(str(smiles_val))
        if ref_mol is None:
            continue
        ref_fp = _fingerprint(ref_mol)
        sim = DataStructs.TanimotoSimilarity(target_fp, ref_fp)
        sensitizer_val = r.get("sensitizer")
        rows.append({
            "name": r.get("name", "(unnamed)"),
            "similarity": sim,
            # Missing/blank sensitizer values are kept as "unknown" rather than
            # silently counted as a negative vote further down.
            "sensitizer": str(sensitizer_val).strip().lower() if pd.notna(sensitizer_val) else "unknown",
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
        "potency_note": None,
    }

    if not in_domain:
        result["call"] = "OUT OF APPLICABILITY DOMAIN — no reference neighbor met the similarity cutoff"
        return result, None

    votable = [n for n in in_domain if n["sensitizer"] in ("yes", "no")]
    if not votable:
        result["call"] = "INCONCLUSIVE — in-domain neighbors have no usable sensitizer label"
        return result, None

    pos_votes = sum(1 for n in votable if n["sensitizer"] == "yes")
    neg_votes = len(votable) - pos_votes

    if pos_votes > neg_votes:
        call = f"Read-across suggests SENSITIZER (majority of {len(votable)} in-domain, labeled neighbors positive)"
    elif neg_votes > pos_votes:
        call = f"Read-across suggests NON-SENSITIZER (majority of {len(votable)} in-domain, labeled neighbors negative)"
    else:
        call = "Read-across INCONCLUSIVE (neighbor vote tied)"

    result["call"] = call

    # Potency read-across: report what potency category the in-domain positive
    # neighbors carry, so the person sees more than a binary yes/no when the
    # data supports it. This does NOT vote/average categories -- it just
    # surfaces what's actually in the data, since collapsing 1A/1B into a
    # single "majority" would imply a precision the small reference set can't support.
    pos_categories = [
        n["potency_category"] for n in votable
        if n["sensitizer"] == "yes" and str(n.get("potency_category", "")).strip() not in ("", "nan")
    ]
    if pos_categories:
        from collections import Counter
        counts = Counter(pos_categories)
        result["potency_note"] = (
            "Potency categories among positive in-domain neighbors: " +
            ", ".join(f"{cat} (x{n})" for cat, n in counts.most_common())
        )
    else:
        result["potency_note"] = None

    return result, None


def render_read_across_module():
    st.markdown("#### 🔬 Read-Across (Analog-Based) Prediction")
    st.warning(
        "This is a similarity/read-across engine, not a validated QSAR model. Bundled reference "
        "sets are offered at increasing size: a real 24-compound official validation set spanning "
        "the full potency range, a real 18-compound hard-case set, both combined (42), and a tiny "
        "demo set. Even combined, 42 compounds is not a large or representative applicability "
        "domain — upload your own curated set for production use."
    )

    if not RDKIT_AVAILABLE:
        st.error("RDKit is not installed.")
        return

    st.markdown("##### Reference set")
    ref_choice = st.radio(
        "Which reference set do you want to use?",
        [
            "Combined real reference set (48 compounds: NICEATM + ECVAM/JRC + OECD TG442D)",
            "NICEATM/ICCVAM 18-compound set (hard/edge cases, real & cited)",
            "ECVAM/JRC 24-compound set (official DPRA/h-CLAT validation set, real & cited)",
            "OECD TG 442D 6-compound set (official KeratinoSens proficiency substances, real & cited)",
            "Tiny illustrative demo set (6 compounds, for a quick smoke test)",
            "Upload my own reference CSV",
        ],
    )

    if ref_choice.startswith("Combined"):
        ref_df = _load_combined_reference()
        with st.expander("📖 Source & important caveats for this reference set", expanded=True):
            st.markdown(
                "Combines three independently sourced, real reference sets (see the individual "
                "set descriptions below for full citations). Together they span the full LLNA "
                "potency range — extreme, strong, moderate, weak, and non-sensitizer — rather "
                "than being weighted toward edge cases alone. Still only 48 compounds: a real "
                "but modest-sized reference set, not a substitute for a large validated database.\n\n"
                "**Known cross-source discordance, kept rather than resolved:** salicylic acid appears "
                "in the NICEATM set as an LLNA-positive weak sensitizer, while the official OECD "
                "TG 442D proficiency-substance table lists it as a non-sensitizer. Both are real, "
                "cited sources that genuinely disagree — this is left visible rather than silently "
                "picking one, since that disagreement is itself useful information for anyone doing "
                "a careful weight-of-evidence assessment."
            )
    elif ref_choice.startswith("NICEATM"):
        ref_df = _load_csv(NICEATM_REFERENCE_PATH)
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
                "**hard/discordant cases**, not a representative random sample."
            )
    elif ref_choice.startswith("ECVAM"):
        ref_df = _load_csv(ECVAM_REFERENCE_PATH)
        with st.expander("📖 Source & important caveats for this reference set", expanded=True):
            st.markdown(
                "- **Source:** European Commission Joint Research Centre / ECVAM. *Direct Peptide "
                "Reactivity Assay, human Cell Line Activation Test, Myeloid U937 Skin Sensitisation "
                "Test — Phase III Pre-validation Study, Chemical Selection Report* (2012, public JRC "
                "TSAR archive). These are the 24 chemicals officially selected to validate DPRA/h-CLAT, "
                "spanning the full LLNA potency range: extreme, strong, moderate, weak, and "
                "non-sensitizer — not an edge-case set.\n"
                "- **Every SMILES here has been cross-checked**: the RDKit-computed molecular weight "
                "matches known values for each compound.\n"
                "- **Known limitations flagged in the data itself:** nickel chloride and xylene are "
                "LLNA/human discordant cases included deliberately by the original study as "
                "performance-standard reference chemicals — see each row's notes column. Kathon CG "
                "was tested as a commercial mixture; only its active ingredient (CMIT) is represented "
                "here as a single structure. Chlorpromazine was tested as its HCl salt; the SMILES "
                "shown is the free base."
            )
    elif ref_choice.startswith("OECD"):
        ref_df = _load_csv(OECD_PROFICIENCY_PATH)
        with st.expander("📖 Source & important caveats for this reference set", expanded=True):
            st.markdown(
                "- **Source:** OECD Test Guideline No. 442D *In Vitro Skin Sensitisation: ARE-Nrf2 "
                "Luciferase Test Method* (adopted 2015), Annex 2, Table 1 — the official 'Proficiency "
                "Substances' a laboratory must correctly classify before it's considered competent to "
                "run the KeratinoSens assay. This is a regulatory test-guideline document itself, not "
                "a secondary study.\n"
                "- Only the 6 of the original 10 proficiency substances not already present in the "
                "other bundled sets are included here, to avoid duplicate entries in the combined set "
                "(isopropanol, glycerol, 2-mercaptobenzothiazole, and salicylic acid are the overlaps — "
                "see the Combined set's notes for the salicylic acid discordance in particular).\n"
                "- 4-Methylaminophenol was tested as its sulfate salt; the SMILES shown is the free base."
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
        try:
            ref_df = pd.read_csv(uploaded)
        except Exception as e:
            st.error(f"Could not read that file as a CSV: {e}")
            return

    ref_df = _normalize_columns(ref_df)
    required_cols = {"name", "smiles", "sensitizer"}
    if not required_cols.issubset(set(ref_df.columns)):
        st.error(f"Reference CSV must include at least these columns: {sorted(required_cols)}")
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

    if result.get("potency_note"):
        st.caption(result["potency_note"])

    st.session_state["read_across_result"] = result
    st.caption(
        "Read-across confidence depends on reference-set relevance and coverage. "
        "Always confirm mechanistic plausibility (shared structural alerts) before accepting a call."
    )
