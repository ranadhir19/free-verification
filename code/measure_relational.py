#!/usr/bin/env python3
"""Table 4 (relational constraints) + Section 9.1 feature AUCs (pre-reg s.5.1/5.2).

Input : data/edges_pre.csv, data/edges_post.csv
Output: results/table4.json — per-row filter metrics WITH the underlying
        counts, the four AUCs, per-target pre hit-rates and potent counts.
Conventions adopted verbatim from Cycle 1 rerun_section9.py.
"""
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from metrics import auc_rank, filter_metrics  # noqa: E402

pre = pd.read_csv(os.path.join(REBUILD, "data", "edges_pre.csv"))
post = pd.read_csv(os.path.join(REBUILD, "data", "edges_post.csv"))
y = post.label.values

# pre-cutoff relational features
mol_potent_counts = (pre[pre.label == 1].groupby("mol_id")
                     .target_id.nunique())
mol_seen = set(pre.mol_id)
tgt_hitrate = pre.groupby("target_id").label.mean()
tgt_potent = pre[pre.label == 1].groupby("target_id").size()

cnt = post.mol_id.map(mol_potent_counts).fillna(0).astype(int).values
seen = post.mol_id.isin(mol_seen).values
hr = post.target_id.map(tgt_hitrate).values
tp = post.target_id.map(tgt_potent).fillna(0).astype(int).values

rows = {}
rows["potent_ge1"] = filter_metrics(y, cnt >= 1)
rows["potent_ge2"] = filter_metrics(y, cnt >= 2)
rows["potent_ge3"] = filter_metrics(y, cnt >= 3)
rows["mol_seen"] = filter_metrics(y, seen)
rows["tgt_hr_ge50"] = filter_metrics(y, hr >= 0.5)
rows["tgt_hr_ge70"] = filter_metrics(y, hr >= 0.7)
rows["both"] = filter_metrics(y, (cnt >= 1) & (hr >= 0.5))

aucs = {
    "mol_potent_pre": auc_rank(cnt, y),
    "mol_seen_pre": auc_rank(seen.astype(float), y),
    "tgt_potent_pre": auc_rank(tp, y),
    "tgt_hitrate_pre": auc_rank(hr, y),
}
out = {"rows": rows, "aucs": aucs,
       "tgt_hitrate_pre": {t: float(v) for t, v in tgt_hitrate.items()},
       "tgt_potent_pre": {t: int(v) for t, v in tgt_potent.items()},
       "basis": {"post_edges": int(len(post)),
                 "post_true": int((y == 1).sum()),
                 "post_wrong": int((y == 0).sum())}}
with open(os.path.join(REBUILD, "results", "table4.json"), "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps({k: (round(v["rho"], 4), round(v["sigma"], 4),
                      round(v["p_prime"], 4), round(v["bits"], 3))
                  for k, v in rows.items()}, indent=1))
print({k: round(v, 4) for k, v in aucs.items()})
