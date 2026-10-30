#!/usr/bin/env python3
"""Rebuild corpus construction — PREREGISTRATION_REBUILD.md ss.1-3.

Input : data/chembl36_activities.csv.gz, data/chembl36_targets.csv
        (the pinned ChEMBL 36 extract deposited in this package)
Output: data/edges_pre.csv, data/edges_post.csv,
        results/targets8.json, results/corpus_stats.json
Edge = (target, molecule) pair; era by first-appearance year;
pre potency = max pChEMBL over activities with doc_year < 2018;
post potency = max pChEMBL over all activities of the pair;
label = 1 iff potency >= 6.0. Universe: SINGLE PROTEIN, Homo sapiens.
Targets: top 8 by pre-edge count desc, tie-break post-edge count desc,
then target_chembl_id ascending.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)
# The deposited copies of the ChEMBL 36 extract (self-contained deposit;
# SHA-256 in README.md). build_corpus reads these copies, so this package
# regenerates the corpus offline.
ACT = os.path.join(REBUILD, "data", "chembl36_activities.csv.gz")
TGT = os.path.join(REBUILD, "data", "chembl36_targets.csv")
CUTOFF = 2018
THETA = 6.0

tg = pd.read_csv(TGT)
universe = set(tg[(tg.target_type == "SINGLE PROTEIN") &
                  (tg.organism == "Homo sapiens")].target_id)
pref = dict(zip(tg.target_id, tg.pref_name))
print(f"universe targets: {len(universe)}", flush=True)

d = pd.read_csv(ACT, dtype={"target_id": str, "mol_id": str})
pool_rows = len(d)
d = d[d.pchembl.notna() & d.doc_year.notna()].copy()
d = d[d.target_id.isin(universe)]
d["doc_year"] = d.doc_year.astype(np.int32)
d["pchembl"] = d.pchembl.astype(np.float64)
print(f"pool rows: {pool_rows}; rows in universe: {len(d)}", flush=True)

# pair aggregates
g = d.groupby(["target_id", "mol_id"], sort=True)
pairs = g.agg(min_year=("doc_year", "min"),
              max_all=("pchembl", "max")).reset_index()
pre_act = d[d.doc_year < CUTOFF]
pre_max = pre_act.groupby(["target_id", "mol_id"], sort=True)["pchembl"].max()
pairs = pairs.merge(pre_max.rename("max_pre"),
                    on=["target_id", "mol_id"], how="left")
pairs["era"] = np.where(pairs.min_year < CUTOFF, "pre", "post")
pairs["potency"] = np.where(pairs.era == "pre", pairs.max_pre, pairs.max_all)
assert pairs.potency.notna().all()
pairs["label"] = (pairs.potency >= THETA).astype(np.int8)
print(f"pairs total: {len(pairs)}", flush=True)

# target ranking over the universe (frozen mechanical rule)
cnt = pairs.groupby("target_id").agg(
    pre_n=("era", lambda s: int((s == "pre").sum())),
    post_n=("era", lambda s: int((s == "post").sum()))).reset_index()
cnt = cnt.sort_values(["pre_n", "post_n", "target_id"],
                      ascending=[False, False, True]).reset_index(drop=True)
top8 = cnt.head(8)
targets8 = top8.target_id.tolist()
print("top 8 targets:")
for _, r in top8.iterrows():
    print(f"  {r.target_id}  pre={r.pre_n} post={r.post_n}  {pref.get(r.target_id)}")

corpus = pairs[pairs.target_id.isin(targets8)].copy()
pre = corpus[corpus.era == "pre"][["target_id", "mol_id", "potency", "label"]]
post = corpus[corpus.era == "post"][["target_id", "mol_id", "potency", "label"]]
pre.to_csv(os.path.join(REBUILD, "data", "edges_pre.csv"), index=False)
post.to_csv(os.path.join(REBUILD, "data", "edges_post.csv"), index=False)

stats = {
    "release": "ChEMBL 36",
    "pool_rows_extract": int(pool_rows),
    "universe_targets": int(len(universe)),
    "pairs_total_universe": int(len(pairs)),
    "cutoff": CUTOFF, "theta": THETA,
    "pre_edges": int(len(pre)),
    "post_edges": int(len(post)),
    "post_potent": int(post.label.sum()),
    "post_wrong": int((post.label == 0).sum()),
    "p_error": float((post.label == 0).mean()),
    "pre_potent": int(pre.label.sum()),
    "pre_wrong": int((pre.label == 0).sum()),
    "unique_mols_pre": int(pre.mol_id.nunique()),
    "unique_mols_post": int(post.mol_id.nunique()),
}
with open(os.path.join(REBUILD, "results", "corpus_stats.json"), "w") as f:
    json.dump(stats, f, indent=2)
t8 = {"rule": "top 8 by pre-edge count desc; tie-break post-edge count desc, "
              "then target_chembl_id asc; universe = SINGLE PROTEIN, Homo sapiens",
      "targets": [{"target_id": t, "pref_name": pref.get(t),
                   "pre_edges": int(cnt.loc[cnt.target_id == t, "pre_n"].iloc[0]),
                   "post_edges": int(cnt.loc[cnt.target_id == t, "post_n"].iloc[0])}
                  for t in targets8],
      "ranking_top12": [{"target_id": r.target_id,
                         "pref_name": pref.get(r.target_id),
                         "pre_edges": int(r.pre_n), "post_edges": int(r.post_n)}
                        for _, r in cnt.head(12).iterrows()]}
with open(os.path.join(REBUILD, "results", "targets8.json"), "w") as f:
    json.dump(t8, f, indent=2)
print(json.dumps(stats, indent=2))
