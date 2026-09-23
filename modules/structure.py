"""
Structural & physicochemical analysis module.

Everything in this module is computed directly from the RDKit molecule
object for whatever SMILES the user enters. Nothing is hardcoded per-compound.
Structural alerts are simplified SMARTS-based heuristics commonly used as a
first-pass screen (comparable in spirit to OECD QSAR Toolbox alerts) — they
are NOT a validated skin-sensitization prediction on their own and should be
read as flags for expert review, not verdicts.
"""

import streamlit as st

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw, Lipinski, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Simplified structural alert SMARTS patterns.
# References (for the underlying chemistry, not for exact regulatory thresholds):
#   - Michael acceptors / Schiff base formers / acyl transfer agents / SN2 alkylators
#     are the four reactive-mechanism domains used throughout the skin-sensitization
#     structural-alert literature (e.g. Aptula & Roberts 2006; OECD QSAR Toolbox).
STRUCTURAL_ALERTS = {
    "Michael acceptor (α,β-unsaturated carbonyl)": "[CX3]=[CX3][CX3]=[OX1]",
    "Aldehyde (Schiff base former)": "[CX3H1](=O)[#6]",
    "Ketone (Schiff base former)": "[#6][CX3](=O)[#6]",
    "Ester (acyl transfer agent)": "[CX3](=O)[OX2H0][#6]",
    "Anhydride (acyl transfer agent)": "[CX3](=O)[OX2][CX3](=O)",
    "Alkyl halide (SN2 alkylating agent)": "[CX4][Cl,Br,I]",
    "Epoxide (SN2 alkylating agent)": "[OX2r3]1[#6r3][#6r3]1",
    "Isocyanate": "[NX2]=[CX2]=[OX1]",
}


def analyze_smiles(smiles: str):
    """Return (mol, descriptors dict, alerts dict, error) for a SMILES string."""
    if not RDKIT_AVAILABLE:
        return None, {}, {}, "RDKit is not installed in this environment."

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, {}, {}, "Could not parse this as a valid SMILES string."

    descriptors = {
        "Molecular Weight (g/mol)": round(Descriptors.MolWt(mol), 2),
        "LogP (Crippen)": round(Crippen.MolLogP(mol), 2),
        "TPSA (Ų)": round(Descriptors.TPSA(mol), 2),
        "H-Bond Donors": Lipinski.NumHDonors(mol),
        "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Aromatic Rings": Lipinski.NumAromaticRings(mol),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 2),
    }

    lipinski_violations = sum([
        descriptors["Molecular Weight (g/mol)"] > 500,
        descriptors["LogP (Crippen)"] > 5,
        descriptors["H-Bond Donors"] > 5,
        descriptors["H-Bond Acceptors"] > 10,
    ])
    descriptors["Lipinski Violations"] = lipinski_violations

    alerts = {}
    for name, smarts in STRUCTURAL_ALERTS.items():
        patt = Chem.MolFromSmarts(smarts)
        alerts[name] = bool(patt is not None and mol.HasSubstructMatch(patt))

    return mol, descriptors, alerts, None


def render_structure_module():
    st.markdown("#### 🧬 Structural & Physicochemical Analysis")
    st.caption(
        "Real RDKit computation on the SMILES you enter. Structural alerts are "
        "simplified SMARTS-based flags for expert review, not a standalone "
        "sensitization prediction."
    )

    if not RDKIT_AVAILABLE:
        st.error(
            "RDKit is not installed. Run `pip install rdkit` — this module cannot "
            "compute anything meaningful without it."
        )
        return

    default_smiles = st.session_state.get("compound_smiles", "O=CC=CC1=CC=CC=C1")
    compound_name = st.text_input(
        "Compound name (for your own reference — not used in the calculation)",
        value=st.session_state.get("compound_name", "Cinnamaldehyde"),
    )
    smiles = st.text_input("SMILES", value=default_smiles)

    if not smiles.strip():
        st.info("Enter a SMILES string to analyze.")
        return

    mol, descriptors, alerts, error = analyze_smiles(smiles)

    if error:
        st.error(error)
        return

    st.session_state["compound_smiles"] = smiles
    st.session_state["compound_name"] = compound_name
    st.session_state["structure_descriptors"] = descriptors
    st.session_state["structure_alerts"] = alerts

    col1, col2 = st.columns([1, 1.2], gap="medium")

    with col1:
        st.markdown("##### 2D Structure")
        img = Draw.MolToImage(mol, size=(380, 300))
        st.image(img, use_container_width=True)

    with col2:
        st.markdown("##### Computed Physicochemical Properties")
        st.dataframe(
            {"Property": list(descriptors.keys()), "Value": list(descriptors.values())},
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("##### Structural Alerts (simplified SMARTS screen)")
    flagged = [name for name, present in alerts.items() if present]
    alert_rows = {
        "Alert": list(alerts.keys()),
        "Present": ["✅ Yes" if v else "— No" for v in alerts.values()],
    }
    st.dataframe(alert_rows, use_container_width=True, hide_index=True)

    if flagged:
        st.warning(
            "Flagged reactive-mechanism domain(s): " + ", ".join(flagged) +
            ". These are heuristic flags only — confirm with mechanistic/assay data."
        )
    else:
        st.info("No reactive-mechanism structural alerts matched for this substructure set.")
