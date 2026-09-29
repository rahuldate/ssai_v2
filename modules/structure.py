
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import streamlit as st

# Expanded mechanistic domain SMARTS including SNAr (Nucleophilic Aromatic Substitution for DNCB-like compounds)
SKIN_SENS_ALERTS = {
    "Nucleophilic Aromatic Substitution (SNAr / DNCB-like)": "c1cc([Cl,Br,I,F])c(cc1[N+](=O)[O-])[N+](=O)[O-]",
    "Michael Acceptor (alpha,beta-unsaturated carbonyl)": "[C,c]=[C,c]-[C,c]=O",
    "Aldehyde / Schiff Base Former": "[CX3H1](=O)[#6]",
    "SN2 / Aliphatic Halide": "[CX4][Cl,Br,I]",
    "Epoxide / Ring-Opening": "C1OC1",
    "Acyl Transfer / Anhydride": "[CX3](=O)[O][CX3](=O)"
}

def get_compound_identifiers(smiles: str):
    """Provides basic metadata placeholders or lookups for common names / CAS."""
    # Dictionary mapping for common reference test substances
    known_db = {
        "c1cc(c(cc1[N+](=O)[O-])[N+](=O)[O-])Cl": {"name": "2,4-Dinitrochlorobenzene (DNCB)", "cas": "97-00-7"},
        "O=CC=CC1=CC=CC=C1": {"name": "Cinnamaldehyde", "cas": "104-55-2"},
        "CC(=O)OC1=CC=CC=C1C(=O)O": {"name": "Aspirin", "cas": "50-78-2"}
    }
    cleaned_smiles = smiles.strip()
    return known_db.get(cleaned_smiles, {"name": "Custom / User-Defined Molecule", "cas": "N/A (Check inventory database)"})

def render_mol_preview(smiles: str, size: tuple = (300, 300)):
    """Generates a 2D image preview and displays identifiers like CAS and Common Name."""
    if not smiles:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        # Display metadata
        info = get_compound_identifiers(smiles)
        st.markdown(f"**Common Name:** {info["name"]}")
        st.markdown(f"**CAS Number:** {info["cas"]}")
        
        # Draw image
        img = Draw.MolToImage(mol, size=size)
        buf = BytesIO()
        img.save(buf, format="PNG")
        st.image(buf.getvalue(), caption=f"Structure preview: {smiles}", width=size[0])
        return mol
    else:
        st.error("Invalid SMILES string: RDKit could not parse the chemical structure.")
        return None

def check_mechanistic_alerts(smiles: str):
    """Scans a target SMILES string against known protein-binding structural alerts including SNAr."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ["Invalid SMILES structure"]
    
    matched_alerts = []
    for alert_name, smarts in SKIN_SENS_ALERTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_alerts.append(alert_name)
            
    return matched_alerts if matched_alerts else ["No classic protein-binding alerts flagged (may act via non-covalent or precursor mechanisms)"]
