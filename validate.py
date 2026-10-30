#!/usr/bin/env python3
"""Rebuild validation (PREREGISTRATION_REBUILD.md s.7).

Recomputes EVERY number printed in paper/Free_Verification_rebuilt.tex
from the deposited artifacts in this directory and diffs them:

  stage A  corpus: rebuild edges from data/chembl36_*.csv via the frozen
           design (code/build_corpus.py logic, re-derived here) and
           compare against the deposited data/edges_*.csv exactly;
           verify the top-8 target rule and the corpus counts.
  stage B  Tables 4-8, AUCs, bar: recompute from edges/descriptors and
           the deposited model scores; refit the three logistic models
           and require AUC agreement within 2e-3 (solver reproducibility);
           recompute the cluster bootstrap with the frozen seed and
           require exact agreement with results/bootstrap.json.
  stage C  tex diff: extract the printed values from
           paper/Free_Verification_rebuilt.tex at anchored locations and
           compare with the recomputation at printed precision.
  stage D  near-miss negatives (silent agreement = FAIL):
           (i) theta 6.0 -> 6.5 must change the post potent count;
           (ii) dropping the first selected target must change counts;
           (iii) shuffled labels (seed 7) must destroy the Table 4
           potent>=1 signal (|bits| < 0.2).

Exit 0 and 'OVERALL PASS' only if every stage passes.
Run: python validate.py   (from this directory; see README.md)
"""
import gzip
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, "paper", "Free_Verification_rebuilt.tex")
sys.path.insert(0, os.path.join(HERE, "code"))
from metrics import (auc_rank, cluster_bootstrap_bits, filter_metrics,  # noqa: E402
                     hanley_mcneil_se, sweep, wilson)
from bar import model_bar, static_bar  # noqa: E402

DESCRIPTORS = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
               "NumRotatableBonds", "RingCount", "NumAromaticRings",
               "FractionCSP3", "HeavyAtomCount", "QED", "LipinskiViolations"]
CUTOFF, THETA = 2018, 6.0
FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL.append(name)


def close(a, b, tol):
    return a is not None and b is not None and abs(a - b) <= tol


# ---------------------------------------------------------------- stage A
print("== stage A: corpus re-derivation from the deposited extract ==")
tg = pd.read_csv(os.path.join(HERE, "data", "chembl36_targets.csv"))
universe = set(tg[(tg.target_type == "SINGLE PROTEIN") &
                  (tg.organism == "Homo sapiens")].target_id)
d = pd.read_csv(os.path.join(HERE, "data", "chembl36_activities.csv.gz"),
                dtype={"target_id": str, "mol_id": str})
d = d[d.pchembl.notna() & d.doc_year.notna()]
d = d[d.target_id.isin(universe)].copy()
d["doc_year"] = d.doc_year.astype(np.int32)
pairs = d.groupby(["target_id", "mol_id"], sort=True).agg(
    min_year=("doc_year", "min"), max_all=("pchembl", "max")).reset_index()
pm = d[d.doc_year < CUTOFF].groupby(["target_id", "mol_id"], sort=True)[
    "pchembl"].max()
pairs = pairs.merge(pm.rename("max_pre"), on=["target_id", "mol_id"],
                    how="left")
pairs["era"] = np.where(pairs.min_year < CUTOFF, "pre", "post")
pairs["potency"] = np.where(pairs.era == "pre", pairs.max_pre, pairs.max_all)
pairs["label"] = (pairs.potency >= THETA).astype(np.int8)
cnt = pairs.groupby("target_id").agg(
    pre_n=("era", lambda s: int((s == "pre").sum())),
    post_n=("era", lambda s: int((s == "post").sum()))).reset_index()
cnt = cnt.sort_values(["pre_n", "post_n", "target_id"],
                      ascending=[False, False, True]).reset_index(drop=True)
top8 = cnt.head(8).target_id.tolist()
t8 = json.load(open(os.path.join(HERE, "results", "targets8.json")))
check("top-8 target rule", [t["target_id"] for t in t8["targets"]] == top8,
      f"-> {top8}")
corpus = pairs[pairs.target_id.isin(top8)]
for era in ("pre", "post"):
    dep = pd.read_csv(os.path.join(HERE, "data", f"edges_{era}.csv"))
    mine = (corpus[corpus.era == era][["target_id", "mol_id", "potency",
                                       "label"]]
            .sort_values(["target_id", "mol_id"]).reset_index(drop=True))
    dep = dep.sort_values(["target_id", "mol_id"]).reset_index(drop=True)
    same = len(dep) == len(mine) and \
        (dep.target_id.values == mine.target_id.values).all() and \
        (dep.mol_id.values == mine.mol_id.values).all() and \
        np.allclose(dep.potency.values, mine.potency.values) and \
        (dep.label.values == mine.label.values).all()
    check(f"edges_{era} exact re-derivation", same, f"({len(dep)} edges)")

pre = pd.read_csv(os.path.join(HERE, "data", "edges_pre.csv"))
post = pd.read_csv(os.path.join(HERE, "data", "edges_post.csv"))
cs = json.load(open(os.path.join(HERE, "results", "corpus_stats.json")))
check("corpus counts", cs["pre_edges"] == len(pre)
      and cs["post_edges"] == len(post)
      and cs["post_potent"] == int(post.label.sum())
      and cs["post_wrong"] == int((post.label == 0).sum())
      and close(cs["p_error"], float((post.label == 0).mean()), 1e-12),
      f"pre={len(pre)} post={len(post)} p={cs['p_error']:.6f}")

# ---------------------------------------------------------------- stage B
print("== stage B: tables, AUCs, bar, bootstrap recomputation ==")
y = post.label.values
mol_pc = pre[pre.label == 1].groupby("mol_id").target_id.nunique()
cntv = post.mol_id.map(mol_pc).fillna(0).astype(int).values
seenv = post.mol_id.isin(set(pre.mol_id)).values
tgt_hr = pre.groupby("target_id").label.mean()
hrv = post.target_id.map(tgt_hr).values
tgt_pot = pre[pre.label == 1].groupby("target_id").size()
tpv = post.target_id.map(tgt_pot).fillna(0).astype(int).values
masks4 = {"potent_ge1": cntv >= 1, "potent_ge2": cntv >= 2,
          "potent_ge3": cntv >= 3, "mol_seen": seenv,
          "tgt_hr_ge50": hrv >= 0.5, "tgt_hr_ge70": hrv >= 0.7,
          "both": (cntv >= 1) & (hrv >= 0.5)}
t4 = json.load(open(os.path.join(HERE, "results", "table4.json")))
for k, keep in masks4.items():
    m = filter_metrics(y, keep)
    r = t4["rows"][k]
    ok = all(m[f] == r[f] for f in ("ret_true", "ret_wrong", "tot_true",
                                    "tot_wrong")) and \
        close(m["bits"], r["bits"], 1e-9) if math.isfinite(r["bits"]) \
        else (m["ret_true"] == r["ret_true"] == 0)
    check(f"table4 {k}", bool(ok),
          f"bits mine={m['bits'] if not math.isfinite(m['bits']) else round(m['bits'],4)} "
          f"stored={r['bits']}")
aucs = {"mol_potent_pre": auc_rank(cntv, y),
        "mol_seen_pre": auc_rank(seenv.astype(float), y),
        "tgt_potent_pre": auc_rank(tpv, y),
        "tgt_hitrate_pre": auc_rank(hrv, y)}
for k, v in aucs.items():
    check(f"auc {k}", close(v, t4["aucs"][k], 1e-9), f"{v:.6f}")

desc = pd.read_csv(os.path.join(HERE, "data", "descriptors.csv"))
sp = post.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
ys = sp.label.values
spr = pre.merge(desc, on="mol_id", how="inner").reset_index(drop=True)
t5 = json.load(open(os.path.join(HERE, "results", "table5.json")))
check("structure basis", t5["basis"]["struct_post_edges"] == len(sp),
      f"{len(sp)} edges")
for name in DESCRIPTORS:
    sw = sweep(sp[name].values.astype(float), ys)
    st = t5["descriptors"][name]
    check(f"table5 {name}", close(sw["bits"], st["bits"], 1e-9)
          and close(auc_rank(sp[name].values.astype(float), ys),
                    st["auc"], 1e-9),
          f"bits={sw['bits']:+.4f}")
t6 = json.load(open(os.path.join(HERE, "results", "table6.json")))
rules = {"lipinski_ro5": sp.LipinskiViolations.values == 0,
         "mw_le500": sp.MolWt.values <= 500,
         "qed_ge05": sp.QED.values >= 0.5,
         "veber": (sp.NumRotatableBonds.values <= 10) & (sp.TPSA.values <= 140),
         "leadlike": (sp.MolWt.values <= 350) & (sp.MolLogP.values <= 3)}
for k, keep in rules.items():
    m = filter_metrics(ys, keep)
    check(f"table6 {k}", close(m["bits"], t6[k]["bits"], 1e-9)
          and m["ret_true"] == t6[k]["ret_true"], f"bits={m['bits']:+.4f}")

# model refits (solver reproducibility) + score-based recomputation
from sklearn.linear_model import LogisticRegression  # noqa: E402


def std_fit(df):
    X = df[DESCRIPTORS].values.astype(float)
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd, mu, sd


orc_dep = pd.read_csv(os.path.join(HERE, "data", "scores_oracle.csv.gz"))
sp2 = sp.merge(orc_dep, on=["target_id", "mol_id"], how="left")
oc = t5["oracle"]
keep = (sp2.score.values >= oc["cut"]) if oc["direction"] == ">=" \
    else (sp2.score.values <= oc["cut"])
m = filter_metrics(ys, keep)
check("table5 oracle (deposited scores)",
      close(m["bits"], oc["bits"], 1e-9)
      and close(auc_rank(sp2.score.values, ys), oc["auc"], 1e-9),
      f"bits={m['bits']:+.4f}")
oracle_keep = keep  # mask over the structure-resolved post edges (sp order)
Xs, _, _ = std_fit(sp)
refit = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000).fit(Xs, ys)
check("table5 oracle (refit AUC)",
      close(auc_rank(refit.predict_proba(Xs)[:, 1], ys), oc["auc"], 2e-3),
      f"refit auc={auc_rank(refit.predict_proba(Xs)[:, 1], ys):.4f} "
      f"stored={oc['auc']:.4f}")

t7 = json.load(open(os.path.join(HERE, "results", "table7.json")))
targets = sorted(spr.target_id.unique())


def design(df, with_target, mu, sd):
    X = df[DESCRIPTORS].values.astype(float)
    Xn = (X - mu) / sd
    if with_target:
        oh = np.zeros((len(df), len(targets)))
        idx = df.target_id.map({t: i for i, t in enumerate(targets)}).values
        oh[np.arange(len(df)), idx] = 1.0
        Xn = np.hstack([Xn, oh])
    return Xn


Xtr_raw = spr[DESCRIPTORS].values.astype(float)
mu7, sd7 = Xtr_raw.mean(axis=0), Xtr_raw.std(axis=0)
sd7[sd7 == 0] = 1.0
for arm, tag, wt in (("structure_only", "struct", False),
                     ("structure_target", "struct_tgt", True)):
    dep = pd.read_csv(os.path.join(HERE, f"data/scores_holdout_{tag}.csv.gz"))
    dtr = dep[dep.era == "pre"].merge(
        spr[["target_id", "mol_id", "label"]], on=["target_id", "mol_id"])
    dte = dep[dep.era == "post"].merge(
        sp[["target_id", "mol_id", "label"]], on=["target_id", "mol_id"])
    h = t7[arm]["honest"]
    keep_te = (dte.score.values >= h["train_cut"]) \
        if h["direction"] == ">=" else (dte.score.values <= h["train_cut"])
    fm = filter_metrics(dte.label.values, keep_te)
    ok = close(fm["bits"], h["bits"], 1e-9) \
        and close(auc_rank(dtr.score.values, dtr.label.values),
                  t7[arm]["auc_train"], 1e-9) \
        and close(auc_rank(dte.score.values, dte.label.values),
                  t7[arm]["auc_test"], 1e-9)
    check(f"table7 {arm} (deposited scores)", bool(ok),
          f"honest bits={fm['bits']:+.4f}")
    Xtr = design(spr, wt, mu7, sd7)
    Xte = design(sp, wt, mu7, sd7)
    rf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)
    rf.fit(Xtr, spr.label.values)
    check(f"table7 {arm} (refit AUCs)",
          close(auc_rank(rf.predict_proba(Xtr)[:, 1], spr.label.values),
                t7[arm]["auc_train"], 2e-3)
          and close(auc_rank(rf.predict_proba(Xte)[:, 1], sp.label.values),
                    t7[arm]["auc_test"], 2e-3), "")

barj = json.load(open(os.path.join(HERE, "results", "bar.json")))
p = cs["p_error"]
check("bar static", close(static_bar(p)["bits"], barj["static"]["bits"], 1e-12)
      and close(static_bar(p)["rho_star"], barj["static"]["rho_star"], 1e-12),
      f"bits={barj['static']['bits']:.4f}")
for k in ("0.0", "0.019", "0.038", "0.057"):
    mb = model_bar(p, float(k))
    check(f"bar kappa={k}", close(mb["bits"], barj["kappa_sweep"][k]["bits"],
                                  1e-9), f"bits={mb['bits']:.4f}")

boot = json.load(open(os.path.join(HERE, "results", "bootstrap.json")))


def lohi_ok(b, s):
    if s["lo"] is None:
        return b["lo"] is None and b["hi"] is None
    return close(b["lo"], s["lo"], 1e-12) and close(b["hi"], s["hi"], 1e-12)
boot_masks = dict(masks4)
boot_masks = {k: boot_masks[k] for k in ("potent_ge1", "potent_ge3",
                                         "tgt_hr_ge70")}
for k, keep in boot_masks.items():
    for cl_name, cl in (("by_molecule", post.mol_id.values),
                        ("by_target", post.target_id.values)):
        b = cluster_bootstrap_bits(y, keep, cl)
        s = boot[k][cl_name]
        check(f"bootstrap {k} {cl_name}",
              lohi_ok(b, s) and b["valid"] == s["valid"],
              f"[{b['lo']:+.3f}, {b['hi']:+.3f}] valid={b['valid']}"
              if b["lo"] is not None else
              f"degenerate valid={b['valid']}")
spk = dict(rules)
spk["oracle"] = oracle_keep  # sp2 is sp merged with scores in sp order
for arm, tag in (("holdout_struct", "struct"), ("holdout_struct_tgt", "struct_tgt")):
    dep = pd.read_csv(os.path.join(HERE, f"data/scores_holdout_{tag}.csv.gz"))
    dte = dep[dep.era == "post"].merge(
        sp[["target_id", "mol_id"]].assign(_i=np.arange(len(sp))),
        on=["target_id", "mol_id"], how="inner").sort_values("_i")
    h = t7[{"struct": "structure_only", "struct_tgt": "structure_target"}[tag]]["honest"]
    spk[arm] = (dte.score.values >= h["train_cut"]) \
        if h["direction"] == ">=" else (dte.score.values <= h["train_cut"])
for k, keepm in spk.items():
    key = {"lipinski_ro5": "rule_ro5", "mw_le500": "rule_mw",
           "qed_ge05": "rule_qed", "veber": "rule_veber",
           "leadlike": "rule_leadlike"}.get(k, k)
    for cl_name, cl in (("by_molecule", sp.mol_id.values),
                        ("by_target", sp.target_id.values)):
        b = cluster_bootstrap_bits(ys, keepm, cl)
        s = boot[key][cl_name]
        check(f"bootstrap {key} {cl_name}",
              lohi_ok(b, s) and b["valid"] == s["valid"],
              f"[{b['lo']:+.3f}, {b['hi']:+.3f}] valid={b['valid']}"
              if b["lo"] is not None else
              f"degenerate valid={b['valid']}")

# ---------------------------------------------------------------- stage C
print("== stage C: printed values in Free_Verification_rebuilt.tex ==")
tex = open(TEX).read()
t8j = json.load(open(os.path.join(HERE, "results", "table8.json")))
bar_full = json.load(open(os.path.join(HERE, "results", "bar.json")))


def num(s):
    return float(s.replace("{,}", "").replace(",", ""))


def row_cells(label, occurrence=1):
    """All numeric cells on the tex line(s) containing label."""
    lines = [ln for ln in tex.splitlines() if label in ln]
    assert len(lines) >= occurrence, f"label not found: {label}"
    ln = lines[occurrence - 1]
    cells = [c.strip() for c in ln.split("&")]
    return cells


def floats(s):
    return [float(x) for x in re.findall(r"-?\d+\.\d+", s)]


def cmp_cells(name, got, want, tol):
    tols = tol if isinstance(tol, (list, tuple)) else [tol] * len(want)
    ok = len(got) == len(want) and all(abs(g - w) <= t
                                       for g, w, t in zip(got, want, tols))
    check(f"tex {name}", ok, f"got={got} want={[round(w, 4) for w in want]}")


# Table 4 rows: cells after label are rho, sigma, p', bits
for key, label in (("potent_ge1", "potent binder of $\\geq$1 target"),
                   ("potent_ge2", "potent binder of $\\geq$2 targets"),
                   ("mol_seen", "molecule seen at all"),
                   ("tgt_hr_ge50", "target hit-rate pre-cutoff $\\geq$ 0.5"),
                   ("tgt_hr_ge70", "target hit-rate pre-cutoff $\\geq$ 0.7"),
                   ("both", "both: potent molecule")):
    m = t4["rows"][key]
    cells = row_cells(label)
    cmp_cells(f"table4 {key}", floats(" & ".join(cells[1:])),
              [m["rho"], m["sigma"], m["p_prime"], m["bits"]],
              [6e-5, 6e-5, 6e-5, 6e-4])
cells = row_cells("potent binder of $\\geq$3 targets")
m = t4["rows"]["potent_ge3"]
check("tex table4 potent_ge3 (-inf)",
      "$-\\infty$" in cells[-1] and abs(floats(cells[1])[0] - m["rho"]) < 6e-5
      and abs(floats(cells[2])[0] - m["sigma"]) < 6e-5, f"cells={cells[1:]}")
# Table 5 rows
for name in DESCRIPTORS:
    mm = t5["descriptors"][name]
    cells = row_cells(f"{name} &")
    cmp_cells(f"table5 {name}", floats(" & ".join(cells[1:])),
              [mm["auc"], mm["rho"], mm["sigma"], mm["bits"]],
              [6e-5, 6e-5, 6e-5, 6e-4])
cells = row_cells("all descriptors (oracle, in-sample) &")
mm = t5["oracle"]
cmp_cells("table5 oracle", floats(" & ".join(cells[1:])),
          [mm["auc"], mm["rho"], mm["sigma"], mm["bits"]],
          [6e-5, 6e-5, 6e-5, 6e-4])
# Table 6 rows
for key, label in (("lipinski_ro5", "Lipinski Ro5 pass"),
                   ("mw_le500", "MW $\\leq$ 500 &"),
                   ("qed_ge05", "QED $\\geq$ 0.5"),
                   ("veber", "Veber (RotB"),
                   ("leadlike", "lead-like (MW")):
    mm = t6[key]
    cells = row_cells(label)
    cmp_cells(f"table6 {key}", floats(" & ".join(cells[1:])),
              [mm["rho"], mm["sigma"], mm["bits"]], [6e-5, 6e-5, 6e-4])
# Table 7 rows
for arm, label in (("structure_only", "structure only &"),
                   ("structure_target", "structure + target &")):
    r7 = t7[arm]
    cells = row_cells(label)
    cmp_cells(f"table7 {arm}", floats(" & ".join(cells[1:])),
              [r7["auc_train"], r7["auc_test"], r7["honest"]["rho"],
               r7["honest"]["sigma"], r7["honest"]["bits"],
               r7["oracle"]["bits"]],
              [6e-5, 6e-5, 6e-5, 6e-5, 6e-4, 6e-4])
# Table 8 rows
for key, label in (("rule_ro5", "Ro5 &"), ("rule_mw", "MW $\\leq$ 500 &"),
                   ("rule_qed", "QED $\\geq$ 0.5 &"), ("rule_veber", "Veber &"),
                   ("rule_leadlike", "lead-like &"),
                   ("potent_ge1", "potent binder $\\geq$1 target (Table 4)"),
                   ("tgt_hr_ge70", "target hit-rate $\\geq$ 0.7 (Table 4)"),
                   ("oracle", "all descriptors, oracle (Table 5)"),
                   ("holdout_struct_tgt", "holdout structure + target, honest"),
                   ("holdout_struct", "holdout structure only, honest")):
    v = t8j[key]
    cells = row_cells(label, 2 if key == "rule_mw" else 1)
    got = floats(" & ".join(cells[1:]))
    cmp_cells(f"table8 {key}", got, [v["bits"], v["bits_lo"], v["bits_hi"]],
              6e-4)
cells = row_cells("potent binder $\\geq$3 targets (Table 4)")
check("tex table8 potent_ge3 (-inf)",
      "$-\\infty$" in cells[1] and "$[-\\infty," in cells[2]
      and abs(floats(cells[2])[0] - t8j["potent_ge3"]["bits_hi"]) < 6e-4,
      f"cells={cells[1:]}")
cells = row_cells("error rate $p$ &")
cmp_cells("table8 p", floats(" & ".join(cells[1:])),
          [t8j["p"]["value"], *t8j["p"]["ci"]], 6e-4)
cells = row_cells("AUC, mol seen pre")
cmp_cells("table8 auc", floats(" & ".join(cells[1:])),
          [t8j["auc_mol_seen"]["value"], *t8j["auc_mol_seen"]["ci"]], 6e-4)

# in-text anchors
def grab(name, pattern, want, tol):
    mm = re.search(pattern, tex, re.S)
    if not mm:
        check(f"tex {name}", False, "pattern not found")
        return
    got = [num(g) for g in mm.groups()]
    tols = tol if isinstance(tol, (list, tuple)) else [tol] * len(want)
    ok = len(got) == len(want) and all(abs(g - w) <= t
                                       for g, w, t in zip(got, want, tols))
    check(f"tex {name}", ok, f"got={got} want={[round(w, 4) for w in want]}")


grab("abstract counts",
     r"hypotheses \(([\d{},]+) genuine,\s*([\d{},]+) spurious; measured error rate \$p = ([\d.]+)\$\)",
     [cs["post_potent"], cs["post_wrong"], cs["p_error"]], 6e-5)
grab("bar static",
     r"i\.e\.\\ \$\\rho\^\\star = ([\d.]+)\$,\s*or \$([\d.]+)\$ bits",
     [bar_full["static"]["rho_star"], bar_full["static"]["bits"]],
     [6e-5, 6e-4])
grab("bar model",
     r"to \$\\rho\^\\star = ([\d.]+)\$, i\.e\.\\ \$([\d.]+)\$ bits",
     [bar_full["headline"]["rho_star"], bar_full["headline"]["bits"]],
     [6e-5, 6e-4])
grab("s5 instance",
     r"instance: \$\\rho = ([\d.]+)\$ raises \$p\$ from \$([\d.]+)\$ to \$([\d.]+)\$",
     [t4["rows"]["potent_ge1"]["rho"], cs["p_error"],
      t4["rows"]["potent_ge1"]["p_prime"]], 6e-5)
grab("s4 aside", r"retaining ([\d.]+)\\% of its truths",
     [round(t4["rows"]["potent_ge1"]["sigma"] * 100, 2)], 0.006)
grab("s9.4 best honest",
     r"The best honest value is ([\d.]+) bits, against the oracle bound of ([\d.]+)",
     [t7["structure_target"]["honest"]["bits"],
      t7["structure_target"]["oracle"]["bits"]], 6e-4)
grab("s9.1 aucs",
     r"mol\s+potent pre ([\d.]+), mol seen pre ([\d.]+), tgt potent pre\s+([\d.]+), tgt hitrate pre ([\d.]+)",
     [t4["aucs"]["mol_potent_pre"], t4["aucs"]["mol_seen_pre"],
      t4["aucs"]["tgt_potent_pre"], t4["aucs"]["tgt_hitrate_pre"]], 6e-5)
for tid in top8:
    check(f"tex reproducibility lists {tid}", tid in tex, "")
check("tex provenance (release/cutoff/threshold)",
      "ChEMBL release 36" in tex and "cutoff 2018" in tex
      and "pChEMBL $\\geq 6.0$" in tex, "")

# ---------------------------------------------------------------- stage D
print("== stage D: near-miss negative checks ==")
alt = int((post.potency >= 6.5).sum())
check("near-miss (i) theta 6.5 changes post potent count",
      alt != cs["post_potent"], f"6.5 -> {alt} vs {cs['post_potent']}")
drop = post[post.target_id != top8[0]]
check("near-miss (ii) dropping first target changes counts",
      len(drop) != len(post) and int(drop.label.sum()) != cs["post_potent"],
      f"without {top8[0]}: {len(drop)} edges")
# Amended per the PROGRESS.md deviation note: the frozen single-seed
# bound |bits| < 0.2 is inside the shuffle null's own SD (0.320) for this
# 105-edge filter. Multi-seed version of the same check: the shuffle
# distribution must centre on 0 and the true value must lie >= 5 SD away.
shuf_bits = []
for s in range(1, 51):
    rng = np.random.default_rng(s)
    shuf_bits.append(filter_metrics(rng.permutation(y),
                                    masks4["potent_ge1"])["bits"])
shuf_bits = np.array(shuf_bits)
obs = t4["rows"]["potent_ge1"]["bits"]
z = (obs - shuf_bits.mean()) / shuf_bits.std(ddof=1)
check("near-miss (iii) shuffled labels destroy potent>=1 signal",
      abs(shuf_bits.mean()) < 0.15 and z <= -5,
      f"shuffle mean={shuf_bits.mean():+.4f} sd={shuf_bits.std(ddof=1):.4f} "
      f"observed z={z:.1f}")

print()
if FAIL:
    print(f"OVERALL FAIL ({len(FAIL)} checks): {FAIL}")
    sys.exit(1)
print("OVERALL PASS")
