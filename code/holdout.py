#!/usr/bin/env python3
"""Table 7 — temporal holdout (pre-reg s.5.5).

Trained on structure-resolved PRE edges; tested on structure-resolved POST
edges. Honest threshold = the frozen sweep on TRAINING scores; oracle
threshold = the same sweep on TEST scores.
Output: results/table7.json; deposited scores in
data/scores_holdout_{struct,struct_tgt}.csv.gz (target_id, mol_id, era, score).
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from metrics import auc_rank, filter_metrics, sweep  # noqa: E402

DESCRIPTORS = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
               "NumRotatableBonds", "RingCount", "NumAromaticRings",
               "FractionCSP3", "HeavyAtomCount", "QED", "LipinskiViolations"]

desc = pd.read_csv(os.path.join(REBUILD, "data", "descriptors.csv"))
pre = pd.read_csv(os.path.join(REBUILD, "data", "edges_pre.csv"))
post = pd.read_csv(os.path.join(REBUILD, "data", "edges_post.csv"))
tr = pre.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
te = post.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
print(f"train (struct pre): {len(tr)}; test (struct post): {len(te)}")
targets = sorted(tr.target_id.unique())


def design(df, with_target):
    X = df[DESCRIPTORS].values.astype(float)
    if with_target:
        oh = np.zeros((len(df), len(targets)))
        idx = df.target_id.map({t: i for i, t in enumerate(targets)}).values
        oh[np.arange(len(df)), idx] = 1.0
        X = np.hstack([X, oh])
    return X


def run_arm(with_target, tag):
    Xtr_raw = design(tr, with_target)
    Xte_raw = design(te, with_target)
    nd = len(DESCRIPTORS)
    mu = Xtr_raw[:, :nd].mean(axis=0)
    sd = Xtr_raw[:, :nd].std(axis=0)
    sd[sd == 0] = 1.0
    Xtr = Xtr_raw.copy(); Xte = Xte_raw.copy()
    Xtr[:, :nd] = (Xtr_raw[:, :nd] - mu) / sd
    Xte[:, :nd] = (Xte_raw[:, :nd] - mu) / sd
    m = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)
    m.fit(Xtr, tr.label.values)
    s_tr = m.predict_proba(Xtr)[:, 1]
    s_te = m.predict_proba(Xte)[:, 1]
    honest = sweep(s_tr, tr.label.values)
    keep_te = (s_te >= honest["cut"]) if honest["direction"] == ">=" \
        else (s_te <= honest["cut"])
    fm_honest = filter_metrics(te.label.values, keep_te)
    oracle = sweep(s_te, te.label.values)
    res = {"auc_train": auc_rank(s_tr, tr.label.values),
           "auc_test": auc_rank(s_te, te.label.values),
           "honest": {"train_cut": honest["cut"],
                      "direction": honest["direction"],
                      "train_rho": honest["rho"], "train_sigma": honest["sigma"],
                      "train_bits": honest["bits"], **fm_honest},
           "oracle": oracle,
           "n_train": int(len(tr)), "n_test": int(len(te))}
    dep = pd.concat([
        tr[["target_id", "mol_id"]].assign(era="pre", score=s_tr),
        te[["target_id", "mol_id"]].assign(era="post", score=s_te)],
        ignore_index=True)
    with gzip.open(os.path.join(REBUILD, "data",
                                f"scores_holdout_{tag}.csv.gz"), "wt") as f:
        dep.to_csv(f, index=False)
    return res


out = {"structure_only": run_arm(False, "struct"),
       "structure_target": run_arm(True, "struct_tgt"),
       "train_basis": "structure-resolved pre edges",
       "test_basis": "structure-resolved post edges"}
with open(os.path.join(REBUILD, "results", "table7.json"), "w") as f:
    json.dump(out, f, indent=2)
for arm, r in out.items():
    if isinstance(r, dict) and "auc_train" in r:
        h, o = r["honest"], r["oracle"]
        print(f"{arm}: AUCtr={r['auc_train']:.4f} AUCte={r['auc_test']:.4f} "
              f"honest rho={h['rho']:.4f} sigma={h['sigma']:.4f} "
              f"bits={h['bits']:+.3f} | oracle bits={o['bits']:+.3f}")
