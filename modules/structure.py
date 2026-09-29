
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import streamlit as st

SKIN_SENS_ALERTS = {
    "Nucleophilic Aromatic Substitution (SNAr / DNCB-like)": "c1cc([Cl,Br,I,F])c(cc1[N+](=O)[O-])[N+](=O)[O-]",
    "Michael Acceptor (alpha,beta-unsaturated carbonyl)": "[C,c]=[C,c]-[C,c]=O",
    "Aldehyde / Schiff Base Former": "[CX3H1](=O)[#6]",
    "SN2 / Aliphatic Halide": "[CX4][Cl,Br,I]",
    "Epoxide / Ring-Opening": "C1OC1",
    "Acyl Transfer / Anhydride": "[CX3](=O)[O][CX3](=O)"
}

def get_compound_identifiers(smiles: str):
    known_db = {
        "c1cc(c(cc1[N+](=O)[O-])[N+](=O)[O-])Cl": {"name": "2,4-Dinitrochlorobenzene (DNCB)", "cas": "97-00-7"},
        "O=CC=CC1=CC=CC=C1": {"name": "Cinnamaldehyde", "cas": "104-55-2"},
        "CC(=O)OC1=CC=CC=C1C(=O)O": {"name": "Aspirin", "cas": "50-78-2"}
    }
    cleaned_smiles = smiles.strip()
    return known_db.get(cleaned_smiles, {"name": "Custom / User-Defined Molecule", "cas": "N/A"})

def render_mol_preview(smiles: str, size: tuple = (300, 300)):
    if not smiles:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        info = get_compound_identifiers(smiles)
        c_name = info["name"]
        c_cas = info["cas"]
        st.markdown(f"**Common Name:** {c_name}")
        st.markdown(f"**CAS Number:** {c_cas}")
        
        img = Draw.MolToImage(mol, size=size)
        buf = BytesIO()
        img.save(buf, format="PNG")
        st.image(buf.getvalue(), caption=f"Structure preview: {smiles}", width=size[0])
        return mol
    else:
        st.error("Invalid SMILES string: RDKit could not parse the chemical structure.")
        return None

def check_mechanistic_alerts(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ["Invalid SMILES structure"]
    
    matched_alerts = []
    for alert_name, smarts in SKIN_SENS_ALERTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_alerts.append(alert_name)
            
    return matched_alerts if matched_alerts else ["No classic protein-binding alerts flagged"]

def render_structure_module():
    st.subheader("Chemical Structure & Mechanistic Alert Screening")
    user_smiles = st.text_input("Enter SMILES string:", value="c1cc(c(cc1[N+](=O)[O-])[N+](=O)[O-])Cl")
    if user_smiles:
        mol = render_mol_preview(user_smiles)
        if mol:
            st.markdown("##### 🧬 Mechanistic Domain & Structural Alerts")
            alerts = check_mechanistic_alerts(user_smiles)
            for alert in alerts:
                if "Invalid" in alert:
                    st.error(alert)
                elif "No classic" in alert:
                    st.info(alert)
                else:
                    st.warning(f"⚠️ Flagged Alert: {alert}")
