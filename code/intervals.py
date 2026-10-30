#!/usr/bin/env python3
"""Section 9.5 intervals (Table 8), the bar, and the cluster bootstrap
(pre-reg ss.5.6-5.7).

Output: results/table8.json, results/bar.json, results/bootstrap.json.
All intervals from the rebuild's ACTUAL 2x2 counts on each quantity's
basis set. Cluster bootstrap: B=2000, seed 20261030, clusters = molecule
and (separately) target.
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from metrics import (auc_rank, cluster_bootstrap_bits, filter_metrics,  # noqa: E402
                     hanley_mcneil_se, wilson)
from bar import model_bar, static_bar  # noqa: E402

R = os.path.join(REBUILD, "results")
D = os.path.join(REBUILD, "data")
cs = json.load(open(os.path.join(R, "corpus_stats.json")))
t4 = json.load(open(os.path.join(R, "table4.json")))
t5 = json.load(open(os.path.join(R, "table5.json")))
t6 = json.load(open(os.path.join(R, "table6.json")))
t7 = json.load(open(os.path.join(R, "table7.json")))

pre = pd.read_csv(os.path.join(D, "edges_pre.csv"))
post = pd.read_csv(os.path.join(D, "edges_post.csv"))
desc = pd.read_csv(os.path.join(D, "descriptors.csv"))


def bits_ci(m):
    """Wilson CIs for rho (x = rejected wrong) and sigma (x = retained
    true), endpoint transform (printed method, pre-reg s.5.7)."""
    rlo, rhi = wilson(m["tot_wrong"] - m["ret_wrong"], m["tot_wrong"])
    slo, shi = wilson(m["ret_true"], m["tot_true"])
    lo = math.log2(slo / (1 - rlo)) if rlo < 1 and slo > 0 else float("nan")
    hi = math.log2(shi / (1 - rhi)) if rhi < 1 and shi > 0 else float("nan")
    return {"bits_lo": lo, "bits_hi": hi,
            "rho_ci": [rlo, rhi], "sigma_ci": [slo, shi]}


# ---- the bar ----
p = cs["p_error"]
bar = {"p": p, "static": static_bar(p)}
bar["kappa_sweep"] = {str(k): model_bar(p, k)
                      for k in (0.0, 0.019, 0.038, 0.057)}
bar["headline"] = bar["kappa_sweep"]["0.038"]
with open(os.path.join(R, "bar.json"), "w") as f:
    json.dump(bar, f, indent=2)

# ---- p interval ----
plo, phi = wilson(cs["post_wrong"], cs["post_edges"])

# ---- Table 8 rows ----
t8 = {"p": {"value": p, "ci": [plo, phi], "n": cs["post_edges"]}}
for key, row in (("potent_ge1", "potent_ge1"), ("potent_ge3", "potent_ge3"),
                 ("tgt_hr_ge70", "tgt_hr_ge70")):
    m = t4["rows"][row]
    t8[key] = {**{k: m[k] for k in ("rho", "sigma", "bits")}, **bits_ci(m)}
for key, row in (("rule_ro5", "lipinski_ro5"), ("rule_mw", "mw_le500"),
                 ("rule_qed", "qed_ge05"), ("rule_veber", "veber"),
                 ("rule_leadlike", "leadlike")):
    m = t6[row]
    t8[key] = {**{k: m[k] for k in ("rho", "sigma", "bits")}, **bits_ci(m)}
m = t5["oracle"]
t8["oracle"] = {**{k: m[k] for k in ("rho", "sigma", "bits", "auc")},
                **bits_ci(m)}
for key, arm in (("holdout_struct", "structure_only"),
                 ("holdout_struct_tgt", "structure_target")):
    h = t7[arm]["honest"]
    t8[key] = {**{k: h[k] for k in ("rho", "sigma", "bits")}, **bits_ci(h),
               "auc_test": t7[arm]["auc_test"]}
y = post.label.values
n1, n0 = int((y == 1).sum()), int((y == 0).sum())
mol_seen_auc = t4["aucs"]["mol_seen_pre"]
se = hanley_mcneil_se(mol_seen_auc, n1, n0)
t8["auc_mol_seen"] = {"value": mol_seen_auc,
                      "ci": [mol_seen_auc - 1.96 * se, mol_seen_auc + 1.96 * se]}
with open(os.path.join(R, "table8.json"), "w") as f:
    json.dump(t8, f, indent=2)

# ---- cluster bootstrap ----
# Table 4 keep masks on the full post set (features as in measure_relational)
mol_potent_counts = (pre[pre.label == 1].groupby("mol_id")
                     .target_id.nunique())
cnt = post.mol_id.map(mol_potent_counts).fillna(0).astype(int).values
tgt_hitrate = pre.groupby("target_id").label.mean()
hr = post.target_id.map(tgt_hitrate).values
masks = {
    "potent_ge1": cnt >= 1,
    "potent_ge3": cnt >= 3,
    "tgt_hr_ge70": hr >= 0.7,
}
sp = post.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
ys = sp.label.values
masks_s = {
    "rule_ro5": sp.LipinskiViolations.values == 0,
    "rule_mw": sp.MolWt.values <= 500,
    "rule_qed": sp.QED.values >= 0.5,
    "rule_veber": (sp.NumRotatableBonds.values <= 10) & (sp.TPSA.values <= 140),
    "rule_leadlike": (sp.MolWt.values <= 350) & (sp.MolLogP.values <= 3),
}
orc = pd.read_csv(os.path.join(D, "scores_oracle.csv.gz"))
sp = sp.merge(orc, on=["target_id", "mol_id"], how="left")
oc = t5["oracle"]
masks_s["oracle"] = (sp.score.values >= oc["cut"]) if oc["direction"] == ">=" \
    else (sp.score.values <= oc["cut"])
for key, arm, tag in (("holdout_struct", "structure_only", "struct"),
                      ("holdout_struct_tgt", "structure_target", "struct_tgt")):
    sc = pd.read_csv(os.path.join(D, f"scores_holdout_{tag}.csv.gz"))
    scp = sc[sc.era == "post"].merge(
        sp[["target_id", "mol_id"]].assign(_i=np.arange(len(sp))),
        on=["target_id", "mol_id"], how="inner").sort_values("_i")
    h = t7[arm]["honest"]
    masks_s[key] = (scp.score.values >= h["train_cut"]) \
        if h["direction"] == ">=" else (scp.score.values <= h["train_cut"])

boot = {}
for key, keep in masks.items():
    boot[key] = {
        "by_molecule": cluster_bootstrap_bits(y, keep, post.mol_id.values),
        "by_target": cluster_bootstrap_bits(y, keep, post.target_id.values)}
for key, keep in masks_s.items():
    boot[key] = {
        "by_molecule": cluster_bootstrap_bits(ys, keep, sp.mol_id.values),
        "by_target": cluster_bootstrap_bits(ys, keep, sp.target_id.values)}
with open(os.path.join(R, "bootstrap.json"), "w") as f:
    json.dump(boot, f, indent=2)

print(f"p = {p:.4f} [{plo:.3f}, {phi:.3f}]")
print(f"static bar: rho*={bar['static']['rho_star']:.4f} bits={bar['static']['bits']:.3f}")
h = bar["headline"]
print(f"model bar (k=0.038): rho*={h['rho_star']:.4f} r={h['r']:.4f} "
      f"pc'={h['pc_prime']:.4f} bits={h['bits']:.3f}")
for k, v in bar["kappa_sweep"].items():
    print(f"  kappa={k}: bits={v['bits']:.3f}")
def _fmt(b):
    if b["lo"] is None:
        return f"[degenerate, valid {b['valid']}]"
    return f"[{b['lo']:+.3f}, {b['hi']:+.3f}] (valid {b['valid']})"


for k, v in boot.items():
    print(f"boot {k:18s} mol {_fmt(v['by_molecule'])} "
          f"tgt {_fmt(v['by_target'])}")
