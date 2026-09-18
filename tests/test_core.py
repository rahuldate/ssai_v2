import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import io
import pandas as pd
from modules.structure import analyze_smiles
from modules.woe import call_dpra, call_keratinosens, call_hclat
from modules.qra import compute_qra
from modules.read_across import compute_read_across, DEMO_REFERENCE_CSV
from modules.potency import classify_from_llna, its_band_score


def test_structure_computes_real_descriptors_for_different_molecules():
    # Two different molecules must yield different molecular weights —
    # this fails immediately if the module falls back to hardcoded values.
    mol1, desc1, alerts1, err1 = analyze_smiles("O=CC=CC1=CC=CC=C1")  # cinnamaldehyde
    mol2, desc2, alerts2, err2 = analyze_smiles("CCO")  # ethanol
    assert err1 is None and err2 is None
    assert desc1["Molecular Weight (g/mol)"] != desc2["Molecular Weight (g/mol)"]
    assert desc1["Molecular Weight (g/mol)"] > 100  # cinnamaldehyde ~132
    assert desc2["Molecular Weight (g/mol)"] < 60   # ethanol ~46


def test_structure_flags_aldehyde_alert_for_cinnamaldehyde():
    mol, desc, alerts, err = analyze_smiles("O=CC=CC1=CC=CC=C1")
    assert err is None
    assert alerts["Aldehyde (Schiff base former)"] is True
    assert alerts["Michael acceptor (α,β-unsaturated carbonyl)"] is True


def test_structure_rejects_invalid_smiles():
    mol, desc, alerts, err = analyze_smiles("not_a_smiles!!!")
    assert mol is None
    assert err is not None


def test_dpra_call_respects_cutoff():
    assert call_dpra(20.0, cutoff=6.38) is True
    assert call_dpra(2.0, cutoff=6.38) is False


def test_keratinosens_requires_both_conditions():
    assert call_keratinosens(imax_fold=2.0, imax_cutoff=1.5, ec_conc_uM=100, ec_cutoff_uM=1000) is True
    # High induction but only at a concentration above the potency cutoff -> negative
    assert call_keratinosens(imax_fold=2.0, imax_cutoff=1.5, ec_conc_uM=5000, ec_cutoff_uM=1000) is False


def test_hclat_either_marker_rule():
    assert call_hclat(cd86_rfi=160, cd86_cutoff=150, cd54_rfi=50, cd54_cutoff=200) is True
    assert call_hclat(cd86_rfi=50, cd86_cutoff=150, cd54_rfi=50, cd54_cutoff=200) is False


def test_qra_margin_of_safety_arithmetic():
    result = compute_qra(
        nesil_ug_cm2=100.0,
        saf_factors={"a": 10.0, "b": 5.0},
        amount_g=1.0,
        concentration_pct=0.1,
        skin_area_cm2=500.0,
        uses_per_day=1.0,
    )
    assert result["SAF_total"] == 50.0
    assert result["AEL_ug_cm2_day"] == 2.0
    # ingredient mass = 1g * 0.001 = 0.001g = 1000 µg; /500cm2 = 2 µg/cm2
    assert abs(result["CEL_per_day_ug_cm2"] - 2.0) < 1e-9
    assert abs(result["MoS"] - 1.0) < 1e-9


def test_qra_higher_exposure_lowers_margin_of_safety():
    low = compute_qra(100.0, {"a": 10.0}, amount_g=1.0, concentration_pct=0.1, skin_area_cm2=500.0, uses_per_day=1.0)
    high = compute_qra(100.0, {"a": 10.0}, amount_g=5.0, concentration_pct=0.1, skin_area_cm2=500.0, uses_per_day=1.0)
    assert high["MoS"] < low["MoS"]


def test_read_across_identical_molecule_gets_similarity_1():
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    # Cinnamaldehyde is in the demo set — querying with its own SMILES must
    # produce a top neighbor similarity of (approximately) 1.0.
    result, err = compute_read_across("O=CC=CC1=CC=CC=C1", ref_df, k=3, domain_cutoff=0.4)
    assert err is None
    assert result["ranked"][0]["name"] == "Cinnamaldehyde"
    assert result["ranked"][0]["similarity"] > 0.99


def test_read_across_flags_out_of_domain_for_unrelated_structure():
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    # A totally unrelated, unusual structure should not have any close neighbor
    # in a 6-compound demo set at a strict cutoff.
    result, err = compute_read_across(
        "C1CC2CCC1CC2N3C(=O)c4ccccc4C3=O", ref_df, k=3, domain_cutoff=0.9
    )
    assert err is None
    assert "OUT OF" in result["call"]


def test_llna_potency_classification_boundary():
    assert classify_from_llna(2.0, cutoff_1a=2.0) == "Category 1A (strong sensitizer)"
    assert classify_from_llna(2.01, cutoff_1a=2.0) == "Category 1B"


def test_its_band_score_picks_correct_band():
    bands = [(6.38, 0), (22.62, 1), (42.47, 2), (999, 3)]
    assert its_band_score(3.0, bands) == 0
    assert its_band_score(10.0, bands) == 1
    assert its_band_score(30.0, bands) == 2
    assert its_band_score(90.0, bands) == 3


def test_niceatm_reference_set_loads_and_all_smiles_valid():
    import os
    from rdkit import Chem
    path = os.path.join(os.path.dirname(__file__), "..", "data", "niceatm_reference_18.csv")
    df = pd.read_csv(path)
    assert len(df) == 18
    invalid = [r["name"] for _, r in df.iterrows() if Chem.MolFromSmiles(r["smiles"]) is None]
    assert invalid == []


def test_niceatm_reference_set_mw_matches_source_table():
    # Cross-check: RDKit-computed MW should match the molecular weights
    # reported in the source poster table for a sample of entries.
    import os
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    path = os.path.join(os.path.dirname(__file__), "..", "data", "niceatm_reference_18.csv")
    df = pd.read_csv(path).set_index("name")
    expected_mw = {
        "Pyridine": 79.10,
        "Propyl gallate": 212.20,
        "Benzoyl peroxide": 242.23,
        "Aniline": 93.13,
    }
    for name, mw in expected_mw.items():
        mol = Chem.MolFromSmiles(df.loc[name, "smiles"])
        assert abs(Descriptors.MolWt(mol) - mw) < 0.1


def test_read_across_works_against_niceatm_reference_set():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "niceatm_reference_18.csv")
    ref_df = pd.read_csv(path)
    # Query with propyl gallate's own SMILES -> should match itself with similarity ~1.0
    result, err = compute_read_across("CCCOC(=O)c1cc(O)c(O)c(O)c1", ref_df, k=3, domain_cutoff=0.4)
    assert err is None
    assert result["ranked"][0]["name"] == "Propyl gallate"
    assert result["ranked"][0]["similarity"] > 0.99
