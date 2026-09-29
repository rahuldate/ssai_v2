import pytest
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import io
import pandas as pd
from modules.structure import analyze_smiles, standardize_mol
from modules.woe import call_dpra, call_keratinosens, call_hclat
from modules.qra import compute_qra
from modules.read_across import compute_read_across, DEMO_REFERENCE_CSV
from modules.potency import classify_from_llna, its_band_score, classify_from_kdpra, its_dpra_points, its_hclat_points, its_total_to_category
from modules.dossier import build_pdf, _safe_filename


def test_structure_computes_real_descriptors_for_different_molecules():
    # Two different molecules must yield different molecular weights —
    # this fails immediately if the module falls back to hardcoded values.
    mol1, desc1, alerts1, err1, std1 = analyze_smiles("O=CC=CC1=CC=CC=C1")  # cinnamaldehyde
    mol2, desc2, alerts2, err2, std2 = analyze_smiles("CCO")  # ethanol
    assert err1 is None and err2 is None
    assert desc1["Molecular Weight (g/mol)"] != desc2["Molecular Weight (g/mol)"]
    assert desc1["Molecular Weight (g/mol)"] > 100  # cinnamaldehyde ~132
    assert desc2["Molecular Weight (g/mol)"] < 60   # ethanol ~46


def test_structure_flags_aldehyde_alert_for_cinnamaldehyde():
    mol, desc, alerts, err, std = analyze_smiles("O=CC=CC1=CC=CC=C1")
    assert err is None
    assert alerts["Aldehyde (Schiff base former)"] is True
    assert alerts["Michael acceptor (α,β-unsaturated carbonyl)"] is True


def test_structure_rejects_invalid_smiles():
    mol, desc, alerts, err, std = analyze_smiles("not_a_smiles!!!")
    assert mol is None
    assert err is not None


def test_dpra_call_respects_cutoff():
    assert call_dpra(20.0, cutoff=6.38) is True
    assert call_dpra(2.0, cutoff=6.38) is False


def test_dpra_cysteine_only_depletion_uses_its_own_verified_cutoff():
    # Verified against OECD GL 497 Table 3.1: cysteine-only depletion cutoff is
    # 13.89%, distinct from the mean Cys+Lys cutoff of 6.38%. A value between
    # the two cutoffs must be classified differently depending on which
    # depletion measure it represents.
    value_between_cutoffs = 10.0
    assert call_dpra(value_between_cutoffs, cutoff=6.38) is True   # positive under the mean-depletion cutoff
    assert call_dpra(value_between_cutoffs, cutoff=13.89) is False  # negative under the cysteine-only cutoff
    # Boundary check for the cysteine-only cutoff itself.
    assert call_dpra(13.89, cutoff=13.89) is True
    assert call_dpra(13.88, cutoff=13.89) is False


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


def test_hclat_logp_domain_limitation_verified_against_primary_tg442e_text():
    # A negative result at LogP<=3.5 is a normal, usable negative.
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=2.0) == "negative"
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=3.5) == "negative"
    # A negative result at LogP>3.5 is unreliable per OECD TG 442E and must NOT
    # be reported as a normal, trustworthy negative.
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=3.6) == "unreliable_negative"
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=6.0) == "unreliable_negative"
    # A POSITIVE result at high LogP is unaffected by this limitation (per the
    # same guideline text: "positive results... could still be used").
    assert call_hclat(cd86_rfi=190, cd54_rfi=50, viability_pct=80, logp=6.0) == "positive"
    # No LogP supplied at all (Structure tab not run) -> falls back to a normal
    # negative rather than silently assuming high LogP.
    assert call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=None) == "negative"


def test_hclat_unreliable_negative_is_excluded_from_2o3_vote():
    # An unreliable negative must be excluded from the vote, same as an
    # invalid/borderline/inconclusive result -- never silently counted as a
    # trustworthy negative vote in the 2-of-3 rule.
    result = call_hclat(cd86_rfi=100, cd54_rfi=100, viability_pct=80, logp=5.0)
    assert result not in ("positive", "negative")


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


from modules.synthesis import synthesize


def test_synthesis_no_evidence():
    result = synthesize(None, None, None, None)
    assert result["streams_run"] == []
    assert result["evidence_strength"] == "no evidence"
    assert "No evidence streams" in result["overall"]


def test_synthesis_single_stream_is_preliminary():
    woe = {"classification": "Sensitizer (≥2/3 assays positive)"}
    result = synthesize(None, woe, None, None)
    assert result["evidence_strength"] == "single stream only"
    assert "preliminary" in result["overall"]
    assert "SENSITIZER" in result["overall"]


def test_synthesis_concordant_streams():
    woe = {"classification": "Sensitizer (≥2/3 assays positive)"}
    read_across = {"call": "Read-across suggests SENSITIZER (majority of 3 in-domain, labeled neighbors positive)"}
    result = synthesize(None, woe, read_across, None)
    assert "concordant" in result["evidence_strength"]
    assert "SENSITIZER" in result["overall"]
    assert result["flags"] == []


def test_synthesis_discordant_streams_are_flagged_not_averaged():
    woe = {"classification": "Non-sensitizer (≥2/3 assays negative)"}
    read_across = {"call": "Read-across suggests SENSITIZER (majority of 3 in-domain, labeled neighbors positive)"}
    result = synthesize(None, woe, read_across, None)
    assert "DISCORDANT" in result["evidence_strength"]
    assert "DISCORDANT" in result["overall"]
    assert len(result["flags"]) == 1
    assert "disagree" in result["flags"][0]


def test_synthesis_structural_alert_context_flags_concern_despite_negative_call():
    woe = {"classification": "Non-sensitizer (≥2/3 assays negative)"}
    alerts = {"Michael acceptor (α,β-unsaturated carbonyl)": True, "Aldehyde (Schiff base former)": False}
    result = synthesize(alerts, woe, None, None)
    assert "worth a closer" in result["structural_note"]
    assert "Michael acceptor" in result["structural_note"]


def test_synthesis_no_alert_context_for_positive_call_without_mechanism():
    woe = {"classification": "Sensitizer (≥2/3 assays positive)"}
    alerts = {k: False for k in ["Michael acceptor (α,β-unsaturated carbonyl)", "Aldehyde (Schiff base former)"]}
    result = synthesize(alerts, woe, None, None)
    assert "outside this tool's simplified alert set" in result["structural_note"]


def test_synthesis_potency_contradiction_flagged():
    woe = {"classification": "Non-sensitizer (≥2/3 assays negative)"}
    read_across = {"call": "Read-across suggests NON-SENSITIZER (majority of 3 in-domain, labeled neighbors negative)"}
    potency = {"path": "A", "category": "Category 1A (strong sensitizer)"}
    result = synthesize(None, woe, read_across, potency)
    assert any("Potency tab reports a strong-sensitizer" in f for f in result["flags"])


def test_synthesis_inconclusive_streams_produce_no_usable_call():
    woe = {"classification": "Inconclusive — assays split; needs expert review"}
    result = synthesize(None, woe, None, None)
    assert result["evidence_strength"] == "no usable hazard call"


from modules.batch_screening import run_batch_screening


def test_batch_screening_runs_multiple_compounds():
    batch_df = pd.DataFrame({
        "name": ["Cinnamaldehyde", "Ethanol"],
        "smiles": ["O=CC=CC1=CC=CC=C1", "CCO"],
    })
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.4)
    assert len(results) == 2
    assert pd.isna(results.loc[0, "error"])
    assert results.loc[0, "n_alerts"] >= 1  # cinnamaldehyde has real alerts
    assert results.loc[1, "n_alerts"] == 0  # ethanol has none


def test_batch_screening_one_bad_smiles_does_not_crash_the_batch():
    batch_df = pd.DataFrame({
        "name": ["Good compound", "Bad compound", "Another good compound"],
        "smiles": ["CCO", "not_a_valid_smiles!!!", "O=CC=CC1=CC=CC=C1"],
    })
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.4)
    assert len(results) == 3  # all three rows present, none dropped
    assert pd.isna(results.loc[0, "error"])
    assert not pd.isna(results.loc[1, "error"])  # bad row reports an error
    assert pd.isna(results.loc[2, "error"])  # third row unaffected by second row's failure


def test_batch_screening_empty_smiles_reported_not_crashed():
    batch_df = pd.DataFrame({"name": ["Missing structure"], "smiles": [None]})
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.4)
    assert len(results) == 1
    assert results.loc[0, "error"] == "Empty SMILES"


def test_batch_screening_requires_name_and_smiles_columns():
    bad_df = pd.DataFrame({"compound": ["X"], "structure": ["CCO"]})
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    with pytest.raises(ValueError):
        run_batch_screening(bad_df, ref_df)


def test_batch_screening_applies_potency_when_ec3_column_present():
    batch_df = pd.DataFrame({
        "name": ["Strong one", "Weak one", "No data"],
        "smiles": ["O=CC=CC1=CC=CC=C1", "CCO", "CC(O)C"],
        "ec3_pct": [1.0, 5.0, None],
    })
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.4)
    assert results.loc[0, "potency_category"] == "Category 1A (strong sensitizer)"
    assert results.loc[1, "potency_category"] == "Category 1B"
    assert pd.isna(results.loc[2, "potency_category"])


def test_batch_screening_works_against_real_combined_reference_set():
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    ref_df = pd.concat([
        pd.read_csv(os.path.join(base, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(base, "oecd_tg442d_proficiency_6.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_transfer_qualification_6.csv")),
    ], ignore_index=True)
    batch_df = pd.DataFrame({"name": ["Formaldehyde query"], "smiles": ["C=O"]})
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.1)
    assert results.loc[0, "nearest_neighbor"] == "Formaldehyde"


def test_kdpra_classification_boundary():
    # Verified against Natsch et al. (2020): log(kmax) >= -2.0 -> 1A.
    assert classify_from_kdpra(-2.0) == "Subcategory 1A (strong sensitizer)"
    assert classify_from_kdpra(-1.5) == "Subcategory 1A (strong sensitizer)"
    assert "Not Subcategory 1A" in classify_from_kdpra(-2.1)
    assert "Not Subcategory 1A" in classify_from_kdpra(-5.0)


def test_kdpra_not_1a_result_does_not_falsely_trigger_synthesis_contradiction_flag():
    # Regression guard: "Not Subcategory 1A" contains the substring "1A", which
    # a naive check would misread as a positive 1A call.
    woe = {"classification": "Non-sensitizer (≥2/3 assays negative)"}
    potency = {"path": "C", "log_kmax": -3.0, "threshold": -2.0,
               "category": classify_from_kdpra(-3.0)}
    result = synthesize(None, woe, None, potency)
    assert not any("strong-sensitizer category (1A)" in f for f in result["flags"])


def test_kdpra_positive_1a_result_does_trigger_synthesis_contradiction_flag():
    woe = {"classification": "Non-sensitizer (≥2/3 assays negative)"}
    potency = {"path": "C", "log_kmax": -1.0, "threshold": -2.0,
               "category": classify_from_kdpra(-1.0)}
    result = synthesize(None, woe, None, potency)
    assert any("strong-sensitizer category (1A)" in f for f in result["flags"])


def test_kdpra_potency_note_labels_path_c_correctly():
    potency = {"path": "C", "log_kmax": -1.0, "threshold": -2.0, "category": classify_from_kdpra(-1.0)}
    result = synthesize(None, None, None, potency)
    assert "Path C, kDPRA" in result["potency_note"]


def test_dossier_pdf_builds_with_kdpra_path_c_result():
    potency_result = {"path": "C", "log_kmax": -1.5, "threshold": -2.0,
                       "category": classify_from_kdpra(-1.5)}
    import streamlit as st

    class FakeSessionState(dict):
        def get(self, k, default=None):
            return super().get(k, default)
    st.session_state = FakeSessionState({"potency_result": potency_result})
    pdf_bytes = build_pdf("Test Compound", "Rahul", "kDPRA test")
    assert pdf_bytes[:4] == b"%PDF"


def test_atompair_is_the_validated_default_fingerprint():
    import inspect
    sig = inspect.signature(compute_read_across)
    assert sig.parameters["fingerprint_type"].default == "atompair"
    assert sig.parameters["k"].default == 7
    assert sig.parameters["domain_cutoff"].default == 0.2


def test_atompair_and_morgan_fingerprints_both_self_match_but_can_rank_differently():
    # Both fingerprint types must still correctly identify a molecule as most
    # similar to itself -- this must hold regardless of which fingerprint is used.
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    for fp_type in ["morgan", "atompair", "maccs"]:
        result, err = compute_read_across(
            "O=CC=CC1=CC=CC=C1", ref_df, k=1, domain_cutoff=0.2, fingerprint_type=fp_type
        )
        assert err is None
        assert result["ranked"][0]["name"] == "Cinnamaldehyde"
        assert result["ranked"][0]["similarity"] > 0.99


def test_combined_mode_alert_failure_confirmed_with_both_fingerprints():
    # Regression guard for the corrected scope note: the alert-blend failure
    # (specificity collapsing) must be empirically confirmed with AtomPair
    # (the current default), not just assumed to carry over from Morgan.
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "data")
    combined = pd.concat([
        pd.read_csv(os.path.join(base, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(base, "oecd_tg442d_proficiency_6.csv")),
        pd.read_csv(os.path.join(base, "ecvam_dpra_transfer_qualification_6.csv")),
    ], ignore_index=True)

    for fp_type in ["morgan", "atompair"]:
        tn = 0
        for i, row in combined.iterrows():
            true_label = str(row["sensitizer"]).strip().lower()
            if true_label != "no":
                continue
            rest = combined.drop(index=i).reset_index(drop=True)
            result, err = compute_read_across(
                row["smiles"], rest, k=5, domain_cutoff=0.3,
                similarity_mode="combined", alpha=0.7, fingerprint_type=fp_type
            )
            if err is None and "NON-SENSITIZER" in result["call"]:
                tn += 1
        assert tn == 0  # confirms specificity=0.000 finding holds for both fingerprints


def test_standardize_mol_strips_salt_and_flags_the_change():
    from rdkit import Chem
    # Chlorpromazine HCl represented with an explicit chloride counter-ion.
    mol_with_salt = Chem.MolFromSmiles("CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21.Cl")
    stripped, changed = standardize_mol(mol_with_salt)
    assert changed is True
    # The stripped structure should be the free base only -- no more disconnected fragments.
    assert len(Chem.GetMolFrags(stripped)) == 1


def test_standardize_mol_does_not_flag_a_change_for_a_plain_molecule():
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")  # ethanol, no salt present
    stripped, changed = standardize_mol(mol)
    assert changed is False


def test_standardize_mol_falls_back_gracefully_for_simple_inorganic_salts():
    from rdkit import Chem
    # Sodium chloride: stripping "everything" would leave nothing -- must not
    # silently return an empty molecule.
    mol = Chem.MolFromSmiles("[Na+].[Cl-]")
    stripped, changed = standardize_mol(mol)
    assert stripped is not None
    assert stripped.GetNumAtoms() > 0


def test_analyze_smiles_reports_was_standardized_flag():
    # Chlorpromazine HCl -- salt should be stripped and flagged.
    mol, desc, alerts, err, was_standardized = analyze_smiles("CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21.Cl")
    assert err is None
    assert was_standardized is True

    # Plain compound -- nothing to strip, flag should be False.
    mol2, desc2, alerts2, err2, was_standardized2 = analyze_smiles("CCO")
    assert err2 is None
    assert was_standardized2 is False


def test_batch_screening_reports_salt_stripped_column():
    batch_df = pd.DataFrame({
        "name": ["Chlorpromazine HCl", "Ethanol"],
        "smiles": ["CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21.Cl", "CCO"],
    })
    ref_df = pd.read_csv(io.StringIO(DEMO_REFERENCE_CSV))
    results = run_batch_screening(batch_df, ref_df, k=1, domain_cutoff=0.4)
    assert bool(results.loc[0, "salt_stripped"]) is True
    assert bool(results.loc[1, "salt_stripped"]) is False


def test_its_dpra_points_match_primary_oecd_gl497_table_3_1():
    # Verified directly against OECD GL 497 (25 June 2025), Table 3.1.
    assert its_dpra_points(50.0) == 3   # >= 42.47
    assert its_dpra_points(42.47) == 3  # boundary, inclusive
    assert its_dpra_points(30.0) == 2   # >= 22.62, < 42.47
    assert its_dpra_points(22.62) == 2  # boundary, inclusive
    assert its_dpra_points(10.0) == 1   # >= 6.38, < 22.62
    assert its_dpra_points(6.38) == 1   # boundary, inclusive
    assert its_dpra_points(3.0) == 0    # < 6.38


def test_its_hclat_points_match_primary_oecd_gl497_table_3_1():
    assert its_hclat_points(mit_ug_ml=500, is_positive=False) == 0  # negative always scores 0
    assert its_hclat_points(mit_ug_ml=5, is_positive=True) == 3     # <= 10
    assert its_hclat_points(mit_ug_ml=10, is_positive=True) == 3    # boundary, inclusive
    assert its_hclat_points(mit_ug_ml=100, is_positive=True) == 2   # > 10, <= 150
    assert its_hclat_points(mit_ug_ml=150, is_positive=True) == 2   # boundary, inclusive
    assert its_hclat_points(mit_ug_ml=4000, is_positive=True) == 1  # > 150, <= 5000


def test_its_total_to_category_matches_primary_oecd_gl497_table_3_1():
    assert its_total_to_category(7) == "GHS 1A (strong sensitizer)"
    assert its_total_to_category(6) == "GHS 1A (strong sensitizer)"  # boundary
    assert its_total_to_category(5) == "GHS 1B"
    assert its_total_to_category(2) == "GHS 1B"  # boundary
    assert its_total_to_category(1) == "Not Classified (non-sensitizer)"
    assert its_total_to_category(0) == "Not Classified (non-sensitizer)"


def test_its_path_b_total_score_reproduces_a_known_official_example():
    # A strong sensitizer scoring maximum on both KE1 and KE3, with a positive
    # in silico call, should hit the top of the 1A band (score 7).
    dpra = its_dpra_points(50.0)          # 3
    hclat = its_hclat_points(5, True)     # 3
    insilico = 1                          # positive
    total = dpra + hclat + insilico
    assert total == 7
    assert its_total_to_category(total) == "GHS 1A (strong sensitizer)"
    # The guideline text notes a chemical can still reach 1A (score >= 6) even
    # without an in silico prediction, since DPRA(3) + h-CLAT(3) = 6 alone.
    assert its_total_to_category(dpra + hclat) == "GHS 1A (strong sensitizer)"
