# Skin Sensitization Decision-Support Tool

A Streamlit calculator for skin sensitization screening: structural analysis,
assay-based hazard classification, read-across prediction, GHS potency
categorization, dermal exposure/risk assessment, and a decision-support
summary export.

**This is a decision-support calculator, not a validated regulatory tool.**
Every number it shows is computed from what you enter — nothing is
pre-filled or simulated. See "What this tool is (and isn't)" below before
using it for anything real.

## Quick start

```bash
pip install -r requirements.txt
streamlit run app.py
```

Run the tests:

```bash
python3 -m pytest tests/ -v
```

## Project layout

```
app.py                          Streamlit entry point / tab navigation
modules/
  structure.py                  RDKit-based structural & physicochemical analysis
  woe.py                        2-out-of-3 assay-based weight-of-evidence (DPRA/KeratinoSens/h-CLAT)
  read_across.py                Similarity-based (Tanimoto) read-across prediction
  potency.py                    GHS 1A/1B potency categorization (LLNA path + in-vitro proxy path)
  synthesis.py                  Integrated evidence synthesis: cross-checks WoE/read-across/structure/potency for agreement or disagreement
  qra.py                        Dermal sensitization QRA (Api et al. 2008 framework)
  dossier.py                    Compiles a PDF summary from whatever was actually computed
data/
  niceatm_reference_18.csv      Real reference compounds, see Data provenance below
  ecvam_dpra_hclat_24.csv       Real reference compounds, see Data provenance below
  oecd_tg442d_proficiency_6.csv Real reference compounds, see Data provenance below
tests/
  test_core.py                  Output-checking tests (not just "does it run")
requirements.txt                Python dependencies
packages.txt                    System-level (apt) dependencies — see Deployment notes
```

## What this tool is (and isn't)

- **Structure tab**: real RDKit computation on whatever SMILES you enter.
  Structural alerts are simplified SMARTS-based flags for expert review, not
  a standalone prediction.
- **WoE tab**: a transparent vote-count rule combining DPRA/KeratinoSens/
  h-CLAT results (OECD GL 497 "defined approach" style). Unlike an earlier
  version of this tool, the DPRA, KeratinoSens, and h-CLAT cutoffs are now
  verified against primary OECD test-guideline text, not assumed defaults —
  see "Cutoff & criterion provenance" below for exactly what was checked.
- **Read-across tab**: real Tanimoto similarity (RDKit Morgan fingerprints)
  against a reference set. It is explicitly **not** presented as equivalent
  to a validated QSAR toolbox — see Data provenance below for exactly what
  the bundled reference data is and isn't.
- **Potency tab**: GHS 1A/1B categorization via three paths — a real LLNA EC3
  (Path A, cutoff corroborated by three independent authoritative sources),
  the real OECD GL 497 ITS scoring framework (Path B, verified directly
  against the primary Guideline text's Table 3.1 — though the in-silico
  component is disclosedly substituted with this app's own Read-Across
  result, not Derek Nexus/OECD QSAR Toolbox as the real ITS specifies), or
  kDPRA (Path C) — a real, OECD-adopted
  single-assay method for confirming Subcategory 1A, with the
  best-evidenced threshold in this entire app (85% balanced accuracy from
  its primary validation study, not a commonly-cited default). kDPRA only
  resolves 1A vs. not-1A and assumes hazard has already been established by
  another method.
- **QRA tab**: the real arithmetic of the Api et al. (2008) dermal
  sensitization QRA framework (AEL = NESIL / SAF; MoS = AEL / CEL).
- **Integrated Summary tab**: reasons over whatever the other tabs actually
  produced this session (never runs anything new). Cross-checks WoE and
  read-across hazard calls against each other and reports agreement or
  disagreement explicitly -- a disagreement is never averaged or silently
  resolved into a false consensus. Uses structural alerts as context, not a
  vote. Deliberately has NO numeric confidence score: there is no validated
  way to turn "2 of 2 streams agree" into a percentage with this data, and
  manufacturing one would repeat the fabricated-MCC mistake this whole
  rebuild started by fixing.
- **Export tab**: a PDF that pulls only what was actually computed in-session.
  Steps you skip are labeled "Not run," never filled with a placeholder.

## Cutoff & criterion provenance

A running record of which numeric thresholds in this tool have been checked
against a primary source, versus which are still placeholder defaults. This
distinction matters more than it might seem: "verified" here always means a
specific document was fetched and the number cross-checked against it in
this project's history, not that the number merely sounds right.

| Threshold | Status | Source |
|---|---|---|
| DPRA 6.38 / 22.62 / 42.47% mean-depletion bands, AND 13.89 / 23.09 / 98.24% cysteine-only bands (co-elution case) | Verified (primary) | OECD GL 497 (25 June 2025), Table 3.1; ECVAM DPRA Validation Study Report |
| KeratinoSens 1.5-fold threshold + 1.35-1.67-fold borderline range + EC1.5-vs-IC30 validity rule | Verified (primary) | Current OECD TG 442D text |
| h-CLAT CD86 >184% / CD54 >255% (viability >=50%) | Verified (primary) | Current OECD TG 442E text (2023) — supersedes the older 150%/200% values still common in secondary literature |
| GHS/CLP 1A/1B cutoff: LLNA EC3 <= 2% -> 1A | Corroborated (3 independent secondary sources) | NTP (2011); EU SCCS working-group table; SenzaGen poster citing OECD GL 497 Annex 2 — not a verbatim excerpt of CLP Annex I itself |
| kDPRA: log(kmax) >= -2.0 -> Subcategory 1A | **Verified against the primary validation study**, with a real published accuracy figure (85% balanced accuracy, n=180) | Natsch et al. (2020), ALTEX 37(4):652-664 — the study OECD's own 2021 peer-review report and TG 442C kDPRA adoption are built on; corroborated by DB-ALM Protocol No. 217 and the official OECD TG 442C page. The strongest-evidenced threshold in this app. |
| QRA: AEL = NESIL / SAF, MoS = AEL / CEL | Structural citation | Api et al. (2008) dermal sensitization QRA framework — the arithmetic structure is real; SAF sub-factor values are user-supplied, not hardcoded |
| Path B (OECD 497 ITS scoring framework) | **Verified directly against the primary Guideline text** (25 June 2025 edition, Table 3.1) for the scoring arithmetic; in-silico component substituted (unverified) | OECD Guideline No. 497, Part II Section 3.1.2, Table 3.1 — fetched and cross-checked directly, not a secondary paraphrase; also corroborated by the 2021 supporting document ENV/CBC/MONO(2021)11, Strickland (NICEATM/ICCVAM, SCHC 2022), and Kleinstreuer et al. (2022, Frontiers in Toxicology). Real ITS uses Derek Nexus/OECD QSAR Toolbox for the in-silico input; this app substitutes its own (moderately predictive, MCC~0.3) Read-Across result instead — disclosed explicitly in-app, not hidden |
| OECD GL 497 scope | 2o3 (WoE tab) implemented and cutoff-verified per assay; ITS (Potency Path B) implemented with verified scoring arithmetic but a disclosed in-silico substitution; SARA-ICE (Part III, added in the 2025 guideline update) NOT implemented — it's a separate Bayesian point-of-departure model, distinct from this app's kDPRA Path C | OECD Guideline No. 497 (25 June 2025 edition) now defines three DA families: 2o3 and ITS (Parts I/II, both present in this app in some form) and SARA-ICE (Part III, not implemented) — see "What this tool is (and isn't)" above |
| Read-across reference data | See "Data provenance" below | Four bundled real, cited datasets (54 compounds); LOO-CV validated, AtomPair fingerprint gives MCC=+0.306 (moderate, real) vs. Morgan's near-chance +0.164 |

## Data provenance

Read-across needs a reference set to compare against. Three real,
source-cited sets are bundled (48 compounds total, no name overlaps),
plus a 6-compound smoke-test demo set. None of these are a substitute for
a large validated database — see the in-app caveats for each set.

| File | Compounds | Source | Character |
|---|---|---|---|
| `niceatm_reference_18.csv` | 18 | Strickland et al. 2015 NICEATM SOT poster (NIEHS/NTP, public domain) | Hard/discordant cases a structural read-across method misclassified — not a representative sample |
| `ecvam_dpra_hclat_24.csv` | 24 | EU JRC/ECVAM DPRA/h-CLAT/MUSST Phase III chemical-selection report (public JRC TSAR archive) | Officially selected to span the full LLNA potency range: extreme → non-sensitizer |
| `oecd_tg442d_proficiency_6.csv` | 6 | OECD Test Guideline 442D, Annex 2, Table 1 (official regulatory text) | KeratinoSens proficiency substances; 4 of the original 10 overlap with the sets above and were excluded here to avoid duplicates |
| `ecvam_dpra_transfer_qualification_6.csv` | 6 | Same ECVAM DPRA Validation Study Report as above, Table 7 (transfer-lab qualification chemicals) | Added specifically to fix a diagnosed sensitizer/non-sensitizer imbalance (33 vs 15 before); contributes 2 more non-sensitizers and 4 more sensitizers |

Every SMILES in every bundled file has been validated (RDKit can parse it)
and, wherever the source gave a molecular weight, cross-checked against
RDKit's computed MW as an independent correctness check.

**A known, deliberately-preserved discordance:** the NICEATM set lists
salicylic acid as a weak LLNA-positive sensitizer; the official OECD TG 442D
proficiency table lists it as a non-sensitizer. Both are real, cited
sources. Rather than silently picking one, the combined reference set keeps
both and flags the disagreement — that's useful information for a careful
assessment, not a bug to be hidden.

To use your own reference data instead, upload a CSV with columns `name,
smiles, sensitizer` (yes/no) and optionally `potency_category, ec3_pct`.

## Predictive validation status (read-across)

Data validation (SMILES parse, MW cross-checks) is not the same thing as
model validation (does it predict correctly?) -- conflating those two was
the core problem with the pre-rebuild version of this app. Run
`python3 validate_read_across.py` for a reproducible leave-one-out
cross-validation of the read-across engine.

### Fingerprint choice matters more than k or the similarity cutoff

A head-to-head comparison across three fingerprint types (Morgan/ECFP4,
MACCS keys, RDKit AtomPair), each swept across k in {3,5,7,9} and cutoff in
{0.2,0.25,0.3,0.35,0.4} on the 54-compound combined set, found a real,
reproducible difference:

| fingerprint | best (k, cutoff) | n scored | accuracy | sensitivity | specificity | MCC |
|---|---|---|---|---|---|---|
| Morgan/ECFP4 (the original default) | k=5, cutoff=0.25 | 33 | 0.606 | 0.700 | 0.462 | +0.164 |
| MACCS keys | k=5, cutoff=0.2 | 52 | 0.577 | 0.714 | 0.294 | +0.009 |
| **AtomPair** | **k=7, cutoff=0.2** | **42** | **0.667** | **0.679** | **0.643** | **+0.306** |

**AtomPair is now the shipped default** (`fingerprint_type="atompair"`,
`k=7`, `domain_cutoff=0.2` -- see `compute_read_across`'s defaults). This is
a real, moderate positive correlation -- not chance, and a meaningful jump
from the ~0 (or negative) MCC the original Morgan-fingerprint default gave.
By conventional correlation-strength interpretation MCC=0.306 is "moderate,"
not "strong": treat this as a genuinely better triage signal, not a solved
prediction problem. It also scores fewer compounds (42 of 54, ~78%) than a
looser cutoff would -- the improvement comes partly from being pickier about
what counts as "in domain," which is an honest trade, not free lunch.

Similarity-weighted voting (weighting each neighbor's vote by its similarity
rather than simple majority count) was also tested against this same
AtomPair/k=7/cutoff=0.2 configuration and **underperformed** plain majority
voting (MCC +0.202 vs +0.306) -- so majority voting remains the default.

### Experiment: combining similarity with structural-alert overlap

`compute_read_across` supports an opt-in `similarity_mode="combined"` that
blends fingerprint similarity with Jaccard similarity over the Structure
tab's 8 SMARTS alerts. Tested with both the original Morgan fingerprint and
the new AtomPair default -- **harmful in both cases** (MCC -0.16 to -0.47,
specificity collapses to 0.000 in every tested configuration). Diagnosed
cause: 68% of sensitizers and 65% of non-sensitizers in this dataset trigger
zero of the 8 SMARTS alerts, so the alert-similarity term mostly contributes
a uniform drag on non-sensitizer pairs while occasionally boosting
sensitizer pairs that do share a triggered alert -- nothing ever outranks a
sensitizer. **Not wired into the app's UI** as a result. The code path is
implemented and tested so a future attempt starts from a known failure mode;
a revised design would need a richer alert set and/or different handling of
the both-absent case before it's worth re-testing.

Re-run `validate_read_across.py` and update this section whenever the
reference data or similarity logic changes materially -- every number above
was independently reproduced by re-running the script, not carried forward
from an earlier claim.

## Deployment notes (Streamlit Community Cloud)

Two gotchas we've hit in practice, in case the app fails to start or a
tab shows "RDKit is not installed" despite `requirements.txt` looking correct:

1. **Python version.** Streamlit Cloud has, at times, defaulted new
   deployments to a Python version newer than RDKit has wheels for. If the
   build log shows a very recent Python version and RDKit still won't import,
   delete and redeploy the app, explicitly picking Python 3.11 or 3.12 in
   "Advanced settings" at deploy time (this can't be changed on an existing
   deployment — only at creation).
2. **Missing system library for `rdkit.Chem.Draw`.** RDKit's PyPI wheel needs
   `libXrender.so.1` at the OS level to render 2D structures, and Streamlit
   Cloud's base image doesn't include it. This repo's `packages.txt` installs
   it (and a few related X11/font libraries) automatically — if you fork this
   repo, make sure `packages.txt` stays in the root alongside `requirements.txt`.

## Running the tests

`tests/test_core.py` checks actual computed outputs (not just "the function
didn't crash") — e.g., two different molecules must have different molecular
weights, higher exposure must lower the margin of safety, a reference set's
SMILES must all be RDKit-parseable, and combining reference sets must not
produce name collisions. Run with:

```bash
python3 -m pytest tests/ -v
```
