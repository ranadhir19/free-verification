# Free Verification — replication package

Replication package for the paper:

> **How Much Is Domain Knowledge Worth? Pricing Hypothesis-Space Constraints
> in Bits and Auditing Them on a Dated Corpus**
> Ranadhir Ghosh, FIS Global, Jacksonville, Florida, USA
> (ranadhir.ghosh@fisglobal.com; ORCID 0009-0001-5439-4258)
> Submitted to *Machine Learning* (Springer), 30 October 2026.

The manuscript's Data Availability statement says the package accompanies the
manuscript; this repository is that package. It contains everything needed to
regenerate, offline, every number printed in Section 9 of the paper: the
measurement code, the deposited data artifacts, the frozen pre-registration,
and a validator that recomputes each printed value from the artifacts and
diffs it against the manuscript source.

**Headline result.** On a dated ChEMBL 36 corpus (55,559 pre-cutoff and
23,710 post-cutoff (target, molecule) hypotheses; measured error rate
*p* = 0.2770), a structural constraint must supply 1.85–1.99 bits of evidence
to rescue the discovery loop at its operating point. None does:
molecule-level relational priors are at chance, a target's pre-cutoff
hit-rate earns +0.857 bits, all five canonical drug-likeness rules are
harmful in point estimate (decisive for four of the five), and the best
honest learned constraint (structure + target identity, leakage-free
temporal holdout) earns +1.206 bits — 0 of 2 honest arms clear the bar.

## Layout

| path | contents |
|---|---|
| `paper/Free_Verification_rebuilt.tex` | the submitted manuscript source (validator stage C diffs its printed values) |
| `PREREGISTRATION_REBUILD.md` | the design, frozen before execution (SHA-256 below; `.sha256` file alongside) |
| `code/build_corpus.py` | corpus construction from the pinned extract → edge lists (pre-reg §§1–3) |
| `code/extract_smiles.py` | provenance: canonical SMILES extraction from the ChEMBL chemreps release file → `data/smiles.csv` |
| `code/descriptors.py` | 12 RDKit descriptors per molecule → `data/descriptors.csv` |
| `code/measure_relational.py` | Table 4 (relational constraints) and the §9.1 feature AUCs |
| `code/measure_structure.py` | Tables 5–6 and the deposited oracle model scores |
| `code/holdout.py` | Table 7 (leakage-free temporal holdout) and the deposited holdout scores |
| `code/intervals.py` | Table 8 intervals, the bar, and the cluster bootstrap (§9.5) |
| `code/metrics.py`, `code/bar.py` | shared metric definitions and the bar solver — used by every measurement script *and* by `validate.py`, so printed numbers and their recomputation share one definition source |
| `data/` | the deposited artifacts (below) |
| `results/` | the deposited result JSONs the manuscript's tables were set from (`table4.json` … `table8.json`, `bar.json`, `bootstrap.json`, `corpus_stats.json`, `targets8.json`, `descriptors_meta.json`) |
| `validate.py` | full recomputation and diff (see below) |

## Regenerating Section 9

**One-time join (GitHub copy only).** `data/chembl36_activities.csv.gz`
(28,256,850 bytes) exceeds GitHub's 25 MB web-upload limit, so in this
repository it is stored as two byte parts,
`data/chembl36_activities.csv.gz.part00` and `.part01` (14,128,425 bytes
each). Reassemble it before running anything:

```sh
bash join_data.sh
```

The script concatenates the parts and verifies the SHA-256
(`1583f55cbb72676ff2a8f7a7083fcc27c3f2f5825162170d237445f81be45af7`)
before continuing; the reassembled file is byte-identical to the
deposited artifact.

From this directory, with the requirements installed:

```sh
python validate.py
```

`validate.py` exits 0 and prints `OVERALL PASS` only if all four stages
pass (125 checks):

- **Stage A** re-derives the corpus from the deposited ChEMBL 36 extract
  under the frozen design and requires the deposited edge files to match
  exactly; it also re-verifies the mechanical top-8 target rule.
- **Stage B** recomputes Tables 4–8, the AUCs, the bar and the cluster
  bootstrap from the edge lists, descriptors and deposited model scores;
  refits the three logistic models (AUC agreement within 2e-3); and
  requires the bootstrap at the frozen seed to reproduce
  `results/bootstrap.json` exactly.
- **Stage C** extracts every printed Section 9 value from
  `paper/Free_Verification_rebuilt.tex` at anchored locations and matches
  it against the recomputation at printed precision.
- **Stage D** near-miss negatives: a potency threshold of 6.5 instead of
  6.0 must change the counts; dropping the first selected target must
  change the counts; shuffled labels must destroy the Table 4 signal.

To regenerate the artifacts themselves rather than verify them, run the
pipeline in order (each script overwrites its deposited outputs):

```sh
python code/build_corpus.py        # extract -> edges_pre/post.csv
python code/descriptors.py         # smiles.csv -> descriptors.csv (needs RDKit)
python code/measure_relational.py  # -> results/table4.json
python code/measure_structure.py   # -> results/table5.json, table6.json, scores_oracle
python code/holdout.py             # -> results/table7.json, holdout scores
python code/intervals.py           # -> results/table8.json, bar.json, bootstrap.json
```

`code/extract_smiles.py` is provenance code: its input was the ChEMBL 36
`chembl_36_chemreps.txt.gz` release file (SHA-256 verified against the
release checksums before use), which is not redistributed here; its
deposited output, `data/smiles.csv`, is what the pipeline consumes.

## Requirements

Python 3.12 with the pinned versions in `requirements.txt`
(numpy 1.26.4, pandas 2.1.4, scipy 1.11.4, scikit-learn 1.9.1,
RDKit 2026.9.1). RDKit is needed only for `code/descriptors.py`;
`validate.py` and the measurement scripts run without it.

## Data provenance

All study data derive from **ChEMBL release 36** (Mendez et al., *Nucleic
Acids Research* 2019; EBI), a public resource. The full PostgreSQL dump was
verified against the release checksums (dump SHA-256
`370f3160…e7343`) before extraction; the structure file
`chembl_36_chemreps.txt.gz` likewise (SHA-256 `f9c2fa73…7cc4`). Only the
pinned extract the study was built from is deposited here — not the raw
dumps.

Deposited artifacts (SHA-256):

| file | SHA-256 |
|---|---|
| `data/chembl36_activities.csv.gz` | `1583f55cbb72676ff2a8f7a7083fcc27c3f2f5825162170d237445f81be45af7` |
| `data/chembl36_targets.csv` | `046b3e0fb6782e65834810ffafe3a342453b8e9d9e06f707e160d7abcf07357b` |
| `data/edges_pre.csv` | `30c1b50a3fde5cf25331acca5223f6ef2b9d5a7ef26114c366c51b2d02b82a84` |
| `data/edges_post.csv` | `c0c9fb96583f66823769ee1db35c3414b700a12add7883e3599c100a1e2db23d` |
| `data/smiles.csv` | `6e5f51f415813f6ce33b14fccb4dbc876382a256e4a1b03cbd2b176c665d108f` |
| `data/descriptors.csv` | `f7be4ada670f29f715b9e405d37af838e0c674009586f6834a04bf46289aa9e6` |
| `data/scores_oracle.csv.gz` | `64008cd26c7bbff40bad7057bdaafb13cb2fedf3bb8b5c294b83276caaae8eb5` |
| `data/scores_holdout_struct.csv.gz` | `2ac16dd9084acbbda5f6c98c60aa41fef3a55a0b903527c495ca347515685d99` |
| `data/scores_holdout_struct_tgt.csv.gz` | `2f29ca0377efbab45318da8d48ba6f62eb25e635766155f860a11a5359204de4` |
| `PREREGISTRATION_REBUILD.md` | `ec2688e5da773fb627d9b656a095a4d3c9d39651c2434ccbd408b82e64302e4e` |

The design in brief (full text in the pre-registration): activities with a
pChEMBL value and a document year; potency threshold θ = 6.0; document-year
cutoff 2018; an edge is a (target, molecule) pair whose era is its
first-appearance year, with pre-cutoff potency the maximum over pre-cutoff
activities only; universe = human single-protein targets; targets = the
eight with the most pre-cutoff edges (ties: post-cutoff count, then ChEMBL
ID) — CHEMBL4822 (BACE1), CHEMBL279 (VEGFR2), CHEMBL240 (hERG),
CHEMBL4296 (Nav1.7), CHEMBL253 (CB2), CHEMBL217 (D2R), CHEMBL203 (EGFR),
CHEMBL2971 (JAK2).

## Provenance note

The study reported in the manuscript was executed once, under the
pre-registration deposited here, after the design facts of an earlier run
of the same study proved unrecoverable — two documented reconstruction
attempts (roughly 8,000 candidate designs) failed to reproduce the earlier
run's printed counts, so the study was re-executed under a design frozen
in advance and the manuscript rewritten to the re-execution's verified
numbers. This package is the deposit of that re-execution; it makes no
claim of continuity with the earlier run. One disclosed deviation: the
pre-registration's literal single-seed shuffle bound for validator stage D
proved to lie inside the shuffle null's own noise for the filter in
question, and the check was amended to the same intent in a stronger
multi-seed form (the multi-seed form is what `validate.py` stage D
implements); the pre-registration text itself was not altered.

## License

Code and documentation: MIT (see `LICENSE`). The deposited data artifacts
derive from ChEMBL, which is distributed under a Creative Commons
Attribution-ShareAlike 3.0 Unported license; the derived artifacts here
are shared under the same terms.
