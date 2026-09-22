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
- **WoE tab**: a transparent, editable-cutoff vote-count rule combining
  DPRA/KeratinoSens/h-CLAT results (OECD GL 497 "defined approach" style).
  Cutoffs are pre-filled examples — verify against your current SOP.
- **Read-across tab**: real Tanimoto similarity (RDKit Morgan fingerprints)
  against a reference set. It is explicitly **not** presented as equivalent
  to a validated QSAR toolbox — see Data provenance below for exactly what
  the bundled reference data is and isn't.
- **Potency tab**: GHS 1A/1B categorization from a real LLNA EC3 (Path A), or
  a screening-level in-vitro proxy score (Path B) that is *not* a GHS-official
  classification — it's a triage signal only.
- **QRA tab**: the real arithmetic of the Api et al. (2008) dermal
  sensitization QRA framework (AEL = NESIL / SAF; MoS = AEL / CEL).
- **Export tab**: a PDF that pulls only what was actually computed in-session.
  Steps you skip are labeled "Not run," never filled with a placeholder.

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
cross-validation of the read-across engine. As of the last run:

| k | similarity cutoff | n scored | accuracy | sensitivity | specificity | MCC |
|---|---|---|---|---|---|---|
| 1 | 0.4 | 25 | 0.440 | 0.429 | 0.455 | -0.116 |
| 3 | 0.4 | 18 | 0.500 | 0.545 | 0.429 | -0.025 |
| 3 | 0.3 | 26 | 0.538 | 0.625 | 0.400 | 0.025 |
| 5 | 0.3 | 26 | 0.577 | 0.688 | 0.400 | 0.089 |

**MCC near zero means close to chance performance.** This is not a
statement that read-across as a method is weak -- it's a measurement of
*this specific 48-compound bundled set*, which is small and, for 18 of the
48 compounds, deliberately adversarial (selected because a structural
method misclassified them in the source study). Re-run this script and
update this table whenever the reference data changes materially.

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
