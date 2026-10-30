#!/usr/bin/env python3
"""Tables 5 (descriptors + in-sample oracle) and 6 (canonical rules),
pre-reg ss.5.3-5.4. Basis: structure-resolved post edges.

Input : data/edges_post.csv, data/descriptors.csv
Output: results/table5.json, results/table6.json,
        data/scores_oracle.csv.gz (deposited oracle scores per edge).
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

post = pd.read_csv(os.path.join(REBUILD, "data", "edges_post.csv"))
desc = pd.read_csv(os.path.join(REBUILD, "data", "descriptors.csv"))
sp = post.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
y = sp.label.values
print(f"structure-resolved post edges: {len(sp)} of {len(post)}")

t5 = {"basis": {"struct_post_edges": int(len(sp)),
                "struct_post_true": int((y == 1).sum()),
                "struct_post_wrong": int((y == 0).sum())},
      "descriptors": {}}
for name in DESCRIPTORS:
    v = sp[name].values.astype(float)
    sw = sweep(v, y)
    t5["descriptors"][name] = {"auc": auc_rank(v, y), **sw}

X = sp[DESCRIPTORS].values.astype(float)
mu, sd = X.mean(axis=0), X.std(axis=0)
sd[sd == 0] = 1.0
Xs = (X - mu) / sd
model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)
model.fit(Xs, y)
score = model.predict_proba(Xs)[:, 1]
sw = sweep(score, y)
t5["oracle"] = {"auc": auc_rank(score, y), **sw,
                "model": "LogisticRegression(C=1.0, lbfgs, max_iter=1000) "
                         "on 12 standardized descriptors, in-sample"}
with open(os.path.join(REBUILD, "results", "table5.json"), "w") as f:
    json.dump(t5, f, indent=2)
sp_out = sp[["target_id", "mol_id"]].copy()
sp_out["score"] = score
with gzip.open(os.path.join(REBUILD, "data", "scores_oracle.csv.gz"),
               "wt") as f:
    sp_out.to_csv(f, index=False)

rules = {
    "lipinski_ro5": (sp.LipinskiViolations.values == 0),
    "mw_le500": (sp.MolWt.values <= 500),
    "qed_ge05": (sp.QED.values >= 0.5),
    "veber": (sp.NumRotatableBonds.values <= 10) & (sp.TPSA.values <= 140),
    "leadlike": (sp.MolWt.values <= 350) & (sp.MolLogP.values <= 3),
}
t6 = {k: filter_metrics(y, keep) for k, keep in rules.items()}
with open(os.path.join(REBUILD, "results", "table6.json"), "w") as f:
    json.dump(t6, f, indent=2)
for k, v in t5["descriptors"].items():
    print(f"{k:20s} AUC={v['auc']:.4f} rho={v['rho']:.4f} "
          f"sigma={v['sigma']:.4f} bits={v['bits']:+.3f} ({v['direction']}{v['cut']})")
print(f"{'ORACLE':20s} AUC={t5['oracle']['auc']:.4f} rho={t5['oracle']['rho']:.4f} "
      f"sigma={t5['oracle']['sigma']:.4f} bits={t5['oracle']['bits']:+.3f}")
for k, v in t6.items():
    print(f"rule {k:14s} rho={v['rho']:.4f} sigma={v['sigma']:.4f} bits={v['bits']:+.3f}")
