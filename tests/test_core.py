import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import io
import pandas as pd
from modules.structure import analyze_smiles
from modules.woe import call_dpra, call_keratinosens, call_hclat
from modules.qra import compute_qra
from modules.read_across import compute_read_across, DEMO_REFERENCE_CSV
from modules.potency import classify_from_llna, its_band_score
from modules.dossier import build_pdf, _safe_filename


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


def test_keratinosens_prediction_model_matches_verified_oecd_thresholds():
    # Verified against primary OECD TG 442D text: 1.5-fold threshold, with a
    # statistically-derived borderline range of 1.35-1.67-fold, and validity
    # gated on EC1.5 falling below (not above) the IC30.
    # Clearly positive: above the borderline range, non-cytotoxic.
    assert call_keratinosens(imax_fold=2.0, ec1_5_uM=100, ic30_uM=800, max_tested_uM=1000) == "positive"
    # Clearly negative: below the borderline range, ceiling reached.
    assert call_keratinosens(imax_fold=1.0, ec1_5_uM=None, ic30_uM=800, max_tested_uM=1000) == "negative"
    # Inside the statistically-derived borderline range -> borderline, not forced either way.
    assert call_keratinosens(imax_fold=1.5, ec1_5_uM=100, ic30_uM=800, max_tested_uM=1000) == "borderline"
    # Above threshold, but induction only reached at a cytotoxic concentration -> inconclusive.
    assert call_keratinosens(imax_fold=2.0, ec1_5_uM=900, ic30_uM=800, max_tested_uM=1000) == "inconclusive"
    # Below the borderline range, but testing never reached a real ceiling -> inconclusive, not negative.
    assert call_keratinosens(imax_fold=1.0, ec1_5_uM=None, ic30_uM=None, max_tested_uM=200) == "inconclusive"


def test_hclat_prediction_model_matches_verified_oecd_thresholds():
    # Verified against primary OECD TG 442E text (2023): CD86 >184% or CD54 >255%
    # at viability >=50% -- not the older 150%/200% values from secondary literature.
    assert call_hclat(cd86_rfi=190, cd54_rfi=50, viability_pct=80) == "positive"
    assert call_hclat(cd86_rfi=50, cd54_rfi=260, viability_pct=80) == "positive"
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80) == "negative"
    # Below minimum viability -> invalid, not forced into positive/negative.
    assert call_hclat(cd86_rfi=190, cd54_rfi=260, viability_pct=30) == "invalid"


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


def test_read_across_handles_mixed_case_columns():
    # Bug regression: a reference CSV with 'SMILES'/'Sensitizer' (different case)
    # used to pass the column check but then KeyError on lookup.
    ref_df = pd.DataFrame({
        "Name": ["Test Compound"],
        "SMILES": ["O=CC=CC1=CC=CC=C1"],
        "Sensitizer": ["yes"],
    })
    result, err = compute_read_across("O=CC=CC1=CC=CC=C1", ref_df, k=1, domain_cutoff=0.4)
    assert err is None
    assert result["ranked"][0]["name"] == "Test Compound"


def test_read_across_missing_sensitizer_label_does_not_silently_count_as_negative():
    # Bug regression: a blank/NaN sensitizer value used to be treated as "not yes",
    # i.e. silently counted as a negative vote.
    ref_df = pd.DataFrame({
        "name": ["Unlabeled analog"],
        "smiles": ["O=CC=CC1=CC=CC=C1"],
        "sensitizer": [None],
    })
    result, err = compute_read_across("O=CC=CC1=CC=CC=C1", ref_df, k=1, domain_cutoff=0.4)
    assert err is None
    assert "INCONCLUSIVE" in result["call"]
    assert "NON-SENSITIZER" not in result["call"]


def test_potency_zero_ec3_is_flagged_invalid_not_silently_categorized():
    assert classify_from_llna(0.0) == "Invalid EC3"
    assert classify_from_llna(-1.0) == "Invalid EC3"


def test_dossier_pdf_survives_special_characters_in_user_text():
    # Bug regression: '&', '<', '>' in a compound name or note used to crash
    # ReportLab's Paragraph XML parser.
    pdf_bytes = build_pdf(
        compound_name="R&D Compound <Test> 5>3",
        assessor="Rahul & Team",
        framework_note="Uses < and > and & in the same note",
    )
    assert isinstance(pdf_bytes, bytes)
    assert pdf_bytes[:4] == b"%PDF"


def test_safe_filename_strips_unsafe_characters():
    assert _safe_filename("Compound / X & Y?.csv") == "Compound_X_Y_csv"
    assert _safe_filename("") == "compound"
    assert _safe_filename("Cinnamaldehyde") == "Cinnamaldehyde"


def test_ecvam_reference_set_loads_and_all_smiles_valid():
    import os
    from rdkit import Chem
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_hclat_24.csv")
    df = pd.read_csv(path)
    assert len(df) == 24
    invalid = [r["name"] for _, r in df.iterrows() if Chem.MolFromSmiles(r["smiles"]) is None]
    assert invalid == []


def test_ecvam_reference_set_mw_matches_known_values():
    import os
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_hclat_24.csv")
    df = pd.read_csv(path).set_index("name")
    expected_mw = {
        "Chlorpromazine": 318.86,
        "Imidazolidinyl urea": 388.29,
        "2-Mercaptobenzothiazole (MBT)": 167.25,
        "1,4-Benzoquinone": 108.09,
    }
    for name, mw in expected_mw.items():
        mol = Chem.MolFromSmiles(df.loc[name, "smiles"])
        assert abs(Descriptors.MolWt(mol) - mw) < 0.2


def test_ecvam_reference_set_spans_full_potency_range_not_just_edge_cases():
    # Unlike the NICEATM hard-case set, this set should include every potency band.
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_hclat_24.csv")
    df = pd.read_csv(path)
    ec3 = pd.to_numeric(df["ec3_pct"], errors="coerce").dropna()
    assert (ec3 < 1).any()      # extreme/strong
    assert (ec3 >= 10).any()    # weak
    assert (df["sensitizer"] == "no").any()  # non-sensitizers present


def test_combined_reference_sets_have_no_duplicate_names():
    import os
    p1 = os.path.join(os.path.dirname(__file__), "..", "data", "niceatm_reference_18.csv")
    p2 = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_hclat_24.csv")
    combined = pd.concat([pd.read_csv(p1), pd.read_csv(p2)], ignore_index=True)
    assert len(combined) == 42
    assert combined["name"].duplicated().sum() == 0


def test_read_across_works_against_ecvam_reference_set():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_hclat_24.csv")
    ref_df = pd.read_csv(path)
    result, err = compute_read_across("C=O", ref_df, k=1, domain_cutoff=0.1)
    assert err is None
    assert result["ranked"][0]["name"] == "Formaldehyde"


def test_oecd_proficiency_set_loads_and_all_smiles_valid():
    import os
    from rdkit import Chem
    path = os.path.join(os.path.dirname(__file__), "..", "data", "oecd_tg442d_proficiency_6.csv")
    df = pd.read_csv(path)
    assert len(df) == 6
    invalid = [r["name"] for _, r in df.iterrows() if Chem.MolFromSmiles(r["smiles"]) is None]
    assert invalid == []


def test_oecd_proficiency_set_mw_matches_known_values():
    import os
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    path = os.path.join(os.path.dirname(__file__), "..", "data", "oecd_tg442d_proficiency_6.csv")
    df = pd.read_csv(path).set_index("name")
    expected_mw = {
        "Methyldibromo glutaronitrile (MDBGN)": 265.94,
        "2,4-Dinitrochlorobenzene (DNCB)": 202.55,
        "Lactic acid": 90.08,
    }
    for name, mw in expected_mw.items():
        mol = Chem.MolFromSmiles(df.loc[name, "smiles"])
        assert abs(Descriptors.MolWt(mol) - mw) < 0.1


def test_all_three_real_reference_sets_combined_have_no_duplicate_names():
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    combined = pd.concat([
        pd.read_csv(os.path.join(base, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(base, "oecd_tg442d_proficiency_6.csv")),
    ], ignore_index=True)
    assert len(combined) == 48
    assert combined["name"].duplicated().sum() == 0


def test_read_across_works_against_oecd_proficiency_set():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "oecd_tg442d_proficiency_6.csv")
    ref_df = pd.read_csv(path)
    result, err = compute_read_across("N#CC(Br)(CBr)CCC#N", ref_df, k=1, domain_cutoff=0.1)
    assert err is None
    assert result["ranked"][0]["name"] == "Methyldibromo glutaronitrile (MDBGN)"


def test_read_across_reports_potency_note_for_positive_neighbors():
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    # Querying with cinnamaldehyde's own SMILES should surface its 1B potency
    # category among the positive-neighbor potency note.
    result, err = compute_read_across("O=CC=CC1=CC=CC=C1", ref_df, k=3, domain_cutoff=0.2)
    assert err is None
    assert result["potency_note"] is not None
    assert "1B" in result["potency_note"]


def test_read_across_potency_note_is_none_when_out_of_domain():
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    result, err = compute_read_across(
        "C1CC2CCC1CC2N3C(=O)c4ccccc4C3=O", ref_df, k=3, domain_cutoff=0.9
    )
    assert err is None
    assert result["potency_note"] is None


def test_fingerprint_uses_non_deprecated_api_and_still_matches_self():
    # Regression guard: after switching from GetMorganFingerprintAsBitVect to
    # rdFingerprintGenerator, self-similarity must still be ~1.0.
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    result, err = compute_read_across("CCO", ref_df, k=1, domain_cutoff=0.9)
    assert err is None
    assert result["ranked"][0]["name"] == "Ethanol"
    assert result["ranked"][0]["similarity"] > 0.99


def test_demo_set_has_no_compounds_overlapping_the_real_reference_sets():
    # Regression guard: the demo set used to reuse isopropanol, glycerol,
    # lactic acid, and DNCB from the real sets, which was harmless to any
    # calculation (the app never merges demo with real sets) but confusing
    # to browse -- the same compound appearing to be "duplicated" across
    # different bundled lists. Compare by canonical structure, not name,
    # since the same molecule can have different display names.
    from rdkit import Chem
    import os

    demo_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    real_df = pd.concat([
        pd.read_csv(os.path.join(base, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(base, "oecd_tg442d_proficiency_6.csv")),
    ], ignore_index=True)

    def canon(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else None

    demo_structures = {canon(s) for s in demo_df["smiles"]}
    real_structures = {canon(s) for s in real_df["smiles"]}
    assert demo_structures.isdisjoint(real_structures)


def test_ecvam_transfer_set_loads_and_all_smiles_valid():
    import os
    from rdkit import Chem
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_transfer_qualification_6.csv")
    df = pd.read_csv(path)
    assert len(df) == 6
    invalid = [r["name"] for _, r in df.iterrows() if Chem.MolFromSmiles(r["smiles"]) is None]
    assert invalid == []


def test_ecvam_transfer_set_mw_matches_known_values():
    import os
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_transfer_qualification_6.csv")
    df = pd.read_csv(path).set_index("name")
    expected_mw = {
        "6-Methylcoumarin": 160.17,
        "Diethyl maleate": 172.18,
        "1-Butanol": 74.12,
    }
    for name, mw in expected_mw.items():
        mol = Chem.MolFromSmiles(df.loc[name, "smiles"])
        assert abs(Descriptors.MolWt(mol) - mw) < 0.1


def test_all_four_real_reference_sets_combined_have_no_duplicate_names():
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    combined = pd.concat([
        pd.read_csv(os.path.join(base, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(base, "oecd_tg442d_proficiency_6.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_transfer_qualification_6.csv")),
    ], ignore_index=True)
    assert len(combined) == 54
    assert combined["name"].duplicated().sum() == 0


def test_ecvam_transfer_set_adds_non_sensitizers_not_just_sensitizers():
    # This set was added specifically to fix a diagnosed sensitizer/non-sensitizer
    # imbalance (33 vs 15 before). Guard that it actually contributes negatives.
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_transfer_qualification_6.csv")
    df = pd.read_csv(path)
    assert (df["sensitizer"] == "no").sum() >= 2


def test_read_across_works_against_ecvam_transfer_set():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "ecvam_dpra_transfer_qualification_6.csv")
    ref_df = pd.read_csv(path)
    result, err = compute_read_across("CCCCO", ref_df, k=1, domain_cutoff=0.4)
    assert err is None
    assert result["ranked"][0]["name"] == "1-Butanol"
    assert result["ranked"][0]["similarity"] > 0.99


def test_alert_vector_and_similarity_helpers():
    from rdkit import Chem
    from modules.read_across import _alert_vector, _alert_similarity, _combined_score

    # Cinnamaldehyde triggers Michael-acceptor + aldehyde alerts; ethanol triggers none.
    cinnamaldehyde = Chem.MolFromSmiles("O=CC=CC1=CC=CC=C1")
    ethanol = Chem.MolFromSmiles("CCO")
    vec_cinnamaldehyde = _alert_vector(cinnamaldehyde)
    vec_ethanol = _alert_vector(ethanol)

    assert any(vec_cinnamaldehyde)
    assert not any(vec_ethanol)
    # Both-absent pairs must NOT score as similar (this is the documented design
    # choice; the combined-mode validation shows keeping this at 0.0 -- not 1.0 --
    # is itself not enough to make the combined score competitive, but flipping it
    # to 1.0 would be strictly worse, rewarding any two unrelated inert molecules).
    assert _alert_similarity(vec_ethanol, vec_ethanol) == 0.0
    # A molecule compared to itself, when it DOES trigger alerts, must be 1.0.
    assert _alert_similarity(vec_cinnamaldehyde, vec_cinnamaldehyde) == 1.0

    assert _combined_score(tanimoto=0.8, alert_sim=0.0, alpha=1.0) == 0.8
    assert _combined_score(tanimoto=0.8, alert_sim=1.0, alpha=0.0) == 1.0
    assert abs(_combined_score(tanimoto=0.8, alert_sim=0.4, alpha=0.5) - 0.6) < 1e-9


def test_combined_similarity_mode_is_opt_in_and_does_not_change_default_behavior():
    # Regression guard: adding similarity_mode/alpha parameters must not change
    # a single call's behavior when the caller doesn't pass them.
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    default_result, err1 = compute_read_across("O=CC=CC1=CC=CC=C1", ref_df, k=3, domain_cutoff=0.2)
    explicit_tanimoto_result, err2 = compute_read_across(
        "O=CC=CC1=CC=CC=C1", ref_df, k=3, domain_cutoff=0.2, similarity_mode="tanimoto"
    )
    assert err1 is None and err2 is None
    assert default_result["call"] == explicit_tanimoto_result["call"]
    assert default_result["ranked"][0]["similarity"] == explicit_tanimoto_result["ranked"][0]["similarity"]


def test_combined_mode_runs_without_error_even_though_validated_as_worse():
    # combined mode is a real, callable code path (documented as not-recommended
    # pending a better design) -- it must not crash even though it isn't wired
    # into the main UI.
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    result, err = compute_read_across(
        "O=CC=CC1=CC=CC=C1", ref_df, k=3, domain_cutoff=0.2, similarity_mode="combined", alpha=0.5
    )
    assert err is None
    assert "tanimoto" in result["ranked"][0]
    assert "alert_similarity" in result["ranked"][0]
