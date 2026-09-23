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
<<<<<<< HEAD
        pd.read_csv(os.path.join(BASE, "ecvam_dpra_transfer_qualification_6.csv")),
=======
>>>>>>> a1c85051a37a249f5dc90a2a3f2a639590693e3f
    ], ignore_index=True)


def loo_validate(combined: pd.DataFrame, k: int, domain_cutoff: float):
    y_true, y_pred = [], []
    out_of_domain = inconclusive = 0
    for i, row in combined.iterrows():
        true_label = str(row["sensitizer"]).strip().lower()
        if true_label not in ("yes", "no"):
            continue
        rest = combined.drop(index=i).reset_index(drop=True)
        result, err = compute_read_across(row["smiles"], rest, k=k, domain_cutoff=domain_cutoff)
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
    for k, cutoff in [(1, 0.4), (3, 0.4), (3, 0.3), (5, 0.3)]:
        y_true, y_pred, ood, inc = loo_validate(combined, k, cutoff)
        m = compute_metrics(y_true, y_pred)
        print(f"k={k} cutoff={cutoff}: n_scored={m['n']}, out_of_domain={ood}, inconclusive={inc}")
        print(f"  TP={m['tp']} TN={m['tn']} FP={m['fp']} FN={m['fn']}")
        print(f"  accuracy={m['accuracy']:.3f}  sensitivity={m['sensitivity']:.3f}  "
              f"specificity={m['specificity']:.3f}  MCC={m['mcc']:.3f}\n")


if __name__ == "__main__":
    main()
