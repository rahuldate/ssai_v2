"""
Leave-one-out cross-validation of the read-across engine against the bundled
real reference sets.

This is the FIRST predictive-performance validation run on this tool. Prior
work only validated the DATA (SMILES parse correctly, MW matches known
values) -- this validates the MODEL (does read-across actually predict
correctly on held-out compounds?). Those are different things; conflating
them was the exact mistake the original (pre-rebuild) app made when it
displayed a fabricated MCC of 0.81.

Run with: python3 validate_read_across.py

Interpretation of results as of the last run (see README.md for full
discussion): MCC ranges from about -0.12 to +0.09 across parameter choices
-- i.e. close to chance performance. This is expected given (a) n=48 is a
small sample for LOO-CV, and (b) 18 of the 48 compounds were deliberately
selected in their source study as cases a structural method misclassified,
so a meaningful fraction of the data is adversarial to structure-based
read-across by construction. This script should be re-run (and its output
in README.md updated) whenever the reference data changes.
"""

import os
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")  # silence RDKit's noisy deprecation/info logging for a clean report
from modules.read_across import compute_read_across

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def load_combined():
    return pd.concat([
        pd.read_csv(os.path.join(BASE, "niceatm_reference_18.csv")),
        pd.read_csv(os.path.join(BASE, "ecvam_dpra_hclat_24.csv")),
        pd.read_csv(os.path.join(BASE, "oecd_tg442d_proficiency_6.csv")),
        pd.read_csv(os.path.join(BASE, "ecvam_dpra_transfer_qualification_6.csv")),
    ], ignore_index=True)


def loo_validate(combined: pd.DataFrame, k: int, domain_cutoff: float,
                  similarity_mode: str = "tanimoto", alpha: float = 0.7,
                  fingerprint_type: str = "morgan", vote_mode: str = "majority",
                  standardize: bool = False):
    y_true, y_pred = [], []
    out_of_domain = inconclusive = 0
    for i, row in combined.iterrows():
        true_label = str(row["sensitizer"]).strip().lower()
        if true_label not in ("yes", "no"):
            continue
        rest = combined.drop(index=i).reset_index(drop=True)
        result, err = compute_read_across(row["smiles"], rest, k=k, domain_cutoff=domain_cutoff,
                                           similarity_mode=similarity_mode, alpha=alpha,
                                           fingerprint_type=fingerprint_type, vote_mode=vote_mode,
                                           standardize=standardize)
        if err is not None:
            continue
        call = result["call"]
        if "OUT OF" in call:
            out_of_domain += 1
            continue
        if "INCONCLUSIVE" in call:
            inconclusive += 1
            continue
        pred = "yes" if ("SENSITIZER" in call and "NON" not in call) else "no"
        y_true.append(true_label)
        y_pred.append(pred)
    return y_true, y_pred, out_of_domain, inconclusive


def compute_metrics(y_true, y_pred):
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == "yes" and p == "yes")
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == "no" and p == "no")
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == "no" and p == "yes")
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == "yes" and p == "no")
    n = tp + tn + fp + fn
    accuracy = (tp + tn) / n if n else float("nan")
    sensitivity = tp / (tp + fn) if (tp + fn) else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) else float("nan")
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = ((tp * tn - fp * fn) / denom) if denom else float("nan")
    return {"n": n, "tp": tp, "tn": tn, "fp": fp, "fn": fn,
            "accuracy": accuracy, "sensitivity": sensitivity,
            "specificity": specificity, "mcc": mcc}


def main():
    combined = load_combined()
    print(f"Combined reference set: {len(combined)} compounds\n")

    print("=== Fingerprint comparison (majority vote, wide (k, cutoff) grid, real shipped code) ===\n")
    grid = [(k, c) for k in [3, 5, 7, 9] for c in [0.2, 0.25, 0.3, 0.35, 0.4]]
    best = {}
    for fp_type in ["morgan", "maccs", "atompair"]:
        print(f"--- fingerprint_type={fp_type} ---")
        for k, cutoff in grid:
            y_true, y_pred, ood, inc = loo_validate(combined, k, cutoff, fingerprint_type=fp_type)
            m = compute_metrics(y_true, y_pred)
            if m["n"] < 20:  # skip configs that score too few compounds to be meaningful
                continue
            print(f"  k={k} cutoff={cutoff}: n={m['n']:3d}  acc={m['accuracy']:.3f}  "
                  f"sens={m['sensitivity']:.3f}  spec={m['specificity']:.3f}  MCC={m['mcc']:+.3f}")
            if fp_type not in best or m["mcc"] > best[fp_type][0]:
                best[fp_type] = (m["mcc"], k, cutoff, m)
        print()

    print("=== Best config per fingerprint type ===")
    for fp_type, (mcc, k, cutoff, m) in best.items():
        print(f"  {fp_type}: k={k} cutoff={cutoff} -> MCC={mcc:+.3f} (n={m['n']}, "
              f"acc={m['accuracy']:.3f}, sens={m['sensitivity']:.3f}, spec={m['specificity']:.3f})")

    print("\n=== Vote mode comparison: majority vs. similarity-weighted (atompair, best k/cutoff) ===\n")
    best_fp, (_, best_k, best_cutoff, _) = max(best.items(), key=lambda kv: kv[1][0])
    for vote_mode in ["majority", "weighted"]:
        y_true, y_pred, ood, inc = loo_validate(combined, best_k, best_cutoff,
                                                 fingerprint_type=best_fp, vote_mode=vote_mode)
        m = compute_metrics(y_true, y_pred)
        print(f"  vote_mode={vote_mode}: n={m['n']:3d}  acc={m['accuracy']:.3f}  "
              f"sens={m['sensitivity']:.3f}  spec={m['specificity']:.3f}  MCC={m['mcc']:+.3f}")

    print("\n=== Mode: combined (tanimoto + structural-alert overlap) -- kept for reference, known worse ===\n")
    print("--- with fingerprint_type=morgan (the original combination this was first validated against) ---\n")
    for alpha in [0.7, 0.5, 0.3]:
        for k, cutoff in [(3, 0.3), (5, 0.3)]:
            y_true, y_pred, ood, inc = loo_validate(combined, k, cutoff, similarity_mode="combined",
                                                     alpha=alpha, fingerprint_type="morgan")
            m = compute_metrics(y_true, y_pred)
            print(f"alpha={alpha} k={k} cutoff={cutoff}: n_scored={m['n']}, out_of_domain={ood}, inconclusive={inc}")
            print(f"  accuracy={m['accuracy']:.3f}  sensitivity={m['sensitivity']:.3f}  "
                  f"specificity={m['specificity']:.3f}  MCC={m['mcc']:.3f}\n")

    print("--- with fingerprint_type=atompair (the current default -- does the same failure persist?) ---\n")
    for alpha in [0.7, 0.5, 0.3]:
        for k, cutoff in [(3, 0.3), (5, 0.3)]:
            y_true, y_pred, ood, inc = loo_validate(combined, k, cutoff, similarity_mode="combined",
                                                     alpha=alpha, fingerprint_type="atompair")
            m = compute_metrics(y_true, y_pred)
            print(f"alpha={alpha} k={k} cutoff={cutoff}: n_scored={m['n']}, out_of_domain={ood}, inconclusive={inc}")
            print(f"  accuracy={m['accuracy']:.3f}  sensitivity={m['sensitivity']:.3f}  "
                  f"specificity={m['specificity']:.3f}  MCC={m['mcc']:.3f}\n")

    print("=== Experiment: salt/counter-ion standardization before fingerprinting (PRED-SKIN-style curation) ===\n")
    print("At the current validated-best config (atompair, majority vote):\n")
    for standardize in [False, True]:
        for k, cutoff in [(7, 0.2), (5, 0.2), (9, 0.2), (7, 0.25)]:
            y_true, y_pred, ood, inc = loo_validate(combined, k, cutoff, fingerprint_type="atompair",
                                                     standardize=standardize)
            m = compute_metrics(y_true, y_pred)
            print(f"standardize={standardize} k={k} cutoff={cutoff}: n={m['n']:3d}  acc={m['accuracy']:.3f}  "
                  f"sens={m['sensitivity']:.3f}  spec={m['specificity']:.3f}  MCC={m['mcc']:+.3f}")
        print()


if __name__ == "__main__":
    main()
