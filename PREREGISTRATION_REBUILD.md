# Pre-Registration — REBUILD of the Section 9 ChEMBL study
## "How Much Is Domain Knowledge Worth? Pricing Hypothesis-Space Constraints in Bits and Auditing Them on a Dated Corpus" (Ghosh, sole author)

**Status: FROZEN before any rebuild experiment run.** Freeze timestamp and
SHA-256 are recorded in PROGRESS.md and in the `.sha256` sidecar of this file.
Any later change requires a dated addendum in PROGRESS.md, written before the
affected stage runs.

## 0. What this is (and is not)

Two pre-registered forensic cycles (ChEMBL 37/36/35) proved the original
study's printed §9 counts are non-identifying: the original cutoff date,
potency threshold, 8 target identities and edge-level dataset are
unrecoverable. This rebuild is therefore a **fresh execution of the study the
manuscript's Section 9 describes**, with every unspecified design parameter
fixed below by a stated mechanical rule, chosen before any rebuild count was
computed. It is **not** a recovery of the original run and claims no
continuity with it. Every design choice below is a rebuild choice and is
disclosed as such in the rebuilt manuscript's Reproducibility section.

## 1. Data source (pinned)

- **Release: ChEMBL 36** (exactly one). Rationale: Cycle 2's verified
  pipeline and extracted tables target ChEMBL 36 and survive on disk
  (`replication/cycle2/data/chembl36/`), so no new release download is
  needed and the extract's provenance is already SHA-verified.
- Activity extract: `replication/cycle2/data/chembl36/activities_chembl36.csv.gz`
  (4,896,791 rows; columns target_id, mol_id, doc_year, pchembl, std_type,
  confidence), produced by Cycle 2 from the official EBI PostgreSQL dump
  (dump SHA-256 `370f3160743c97a9e346a08fb194d3f4596fe2e498810c685cefb05cf37e7343`
  per the release `checksums.txt` kept on disk).
- Targets: `targets_chembl36.csv` from the same extraction.
- Structures: `chembl_36_chemreps.txt.gz` from the same EBI release
  directory, SHA-256 `f9c2fa730ac5000da9b0c4033c99bfd8fe2b99950f19df1d842c6e73dd7cc4`
  (from the on-disk `checksums.txt`); downloaded once, used to extract
  canonical SMILES for corpus molecules only, then deleted.

## 2. Frozen design parameters

- **Activity pool:** extract rows with `pchembl` NOT NULL and `doc_year`
  NOT NULL. **No** standard-type restriction, **no** assay-confidence
  restriction (the primary operationalization of both reconstruction cycles).
- **Potency threshold θ = 6.0**, label convention `pchembl ≥ θ`
  (pChEMBL 6.0 ⇔ 1 µM, the conventional activity cut). Fixed a priori.
- **Cutoff C = 2018.** Pre-cutoff era: document year ≤ 2017. Post-cutoff
  era: document year ≥ 2018. A round recent-era calendar split, fixed a
  priori, before any rebuild count was computed.
- **Edge (hypothesis) definition:** a (target, molecule) pair over the
  universe (§3). The pair's era is its **first-appearance era**: era = pre
  iff min(doc_year) < 2018, else post (the first-appearance partition used
  for the post sets of operationalizations E1/E2 in both cycles; the
  printed Table 4's "molecule seen pre-cutoff" row, which retains ≈1% of
  post edges, is consistent only with this partition).
  - Pre-edge potency = max pChEMBL over the pair's activities with
    doc_year < 2018 (era-windowed, so every pre-cutoff statistic uses only
    pre-cutoff information, per the manuscript's leakage rule).
  - Post-edge potency = max pChEMBL over all the pair's activities (all
    are ≥ 2018 by construction).
  - Label = 1 (true/potent) iff potency ≥ 6.0, else 0 (wrong).
- **Corpus counts** (whatever the frozen design yields): pre edges,
  post edges, post potent, post wrong, error rate p = post_wrong/post.

## 3. Universe and target selection (mechanical rule)

- Universe: targets in `targets_chembl36.csv` with
  `target_type = 'SINGLE PROTEIN'` and `organism = 'Homo sapiens'`.
- Rank universe targets by **pre-edge count descending**; tie-break by
  post-edge count descending, then by `target_chembl_id` ascending.
- **The 8 targets are the top 8 of this ranking.** The rule is frozen; the
  identities are whatever it returns and will be reported by ChEMBL ID and
  pref_name in the rebuilt manuscript. No family/therapeutic-area filter
  is applied; descriptions of the selected targets in the manuscript will
  state what they actually are.

## 4. Structure and descriptors

- A molecule has a **resolved structure** iff chemreps carries a
  `canonical_smiles` for it AND RDKit parses it to a valid molecule.
  Structure-resolved post set = post edges whose molecule has a resolved
  structure (the basis of Tables 5 and 6). RDKit version recorded in the
  results files.
- Descriptors (RDKit, from canonical SMILES), the manuscript's Table 5 set:
  MolWt = `Descriptors.MolWt`; MolLogP = `Crippen.MolLogP`;
  TPSA = `rdMolDescriptors.CalcTPSA`; NumHDonors = `Lipinski.NumHDonors`;
  NumHAcceptors = `Lipinski.NumHAcceptors`;
  NumRotatableBonds = `Lipinski.NumRotatableBonds`;
  RingCount = `rdMolDescriptors.CalcNumRings`;
  NumAromaticRings = `Lipinski.NumAromaticRings`;
  FractionCSP3 = `Lipinski.FractionCSP3`;
  HeavyAtomCount = `mol.GetNumHeavyAtoms()`; QED = `QED.qed`;
  LipinskiViolations = #(MolWt>500) + #(MolLogP>5) + #(NumHDonors>5) +
  #(NumHAcceptors>10).

## 5. Measurements (exact computations)

Conventions for every filter F applied to a set of labelled edges:
ρ = 1 − (retained wrong)/(total wrong); σ = (retained true)/(total true);
p′ = (retained wrong)/(retained total); bits = log2(σ/(1−ρ)) computed from
unrounded (ρ, σ). These are the manuscript's definitions (Sections 4, 9).

### 5.1 Table 4 — relational constraints (basis: all post edges)
Pre-cutoff information only, computed from pre edges over the 8 targets
(adopting Cycle 1 `rerun_section9.py` conventions verbatim):
- `potent_mol(t)` counts: for each molecule, the number of distinct corpus
  targets on which it has a potent pre edge. Rows: molecule potent on
  ≥1 / ≥2 / ≥3 targets pre-cutoff (retain post edge iff its molecule meets
  the count).
- Molecule seen at all pre-cutoff: molecule appears in ≥1 pre edge.
- Target hit-rate pre-cutoff: per target, mean label over its pre edges.
  Rows: retain iff target hit-rate ≥ 0.5; retain iff ≥ 0.7.
- "Both: potent molecule AND tractable target": retain iff molecule potent
  on ≥1 target AND the edge's target hit-rate ≥ 0.5. **"Tractable" :=
  hit-rate ≥ 0.5**, disclosed here as a rebuild choice: the manuscript
  leaves the word undefined; 0.5 is the threshold of Table 4's own
  hit-rate row, and it is the reading under which the printed table's
  "both" row equals its first row (the ≥0.5 filter retained every post
  edge in the printed record). The alternative reading (≥0.7, used by the
  reconstruction screens) is recorded as considered and not adopted.

### 5.2 Feature AUCs (Section 9.1)
Rank AUC (Mann–Whitney, average ranks) against the post-edge label, of:
- mol potent pre: the molecule's distinct-target potent count (§5.1);
- mol seen pre: binary seen indicator;
- tgt potent pre: the number of potent pre edges of the edge's target;
- tgt hitrate pre: the edge's target's pre-cutoff hit-rate.

### 5.3 Table 5 — descriptors (basis: structure-resolved post edges)
Per descriptor: AUC of the raw descriptor value against the label; then a
threshold sweep — candidate cuts = all unique observed values, directions
"retain if ≥ cut" scanned first, then "retain if ≤ cut", cuts ascending;
a cut is feasible iff σ ≥ 0.50 and it retains ≥1 wrong edge (finite bits);
report the feasible cut maximizing bits (first maximum in scan order).
Final row: **all descriptors (oracle, in-sample)** — logistic regression
(scikit-learn `LogisticRegression(C=1.0, solver='lbfgs', max_iter=1000)`)
on the 12 descriptors standardized by the fit set's mean/SD, fitted on the
structure-resolved post edges themselves; AUC of its score; threshold
sweep of the score under the same rule.

### 5.4 Table 6 — canonical rules (basis: structure-resolved post edges)
Applied without tuning: Lipinski Ro5 pass = (LipinskiViolations = 0);
MW ≤ 500 = (MolWt ≤ 500); QED ≥ 0.5; Veber = (NumRotatableBonds ≤ 10 AND
TPSA ≤ 140); lead-like = (MolWt ≤ 350 AND MolLogP ≤ 3).

### 5.5 Table 7 — temporal holdout
Models fitted on **structure-resolved pre edges only** (features
standardized by the training set; the structure+target arm adds the
8-column target one-hot). Same logistic-regression specification as §5.3.
- Honest threshold: the §5.3 sweep applied to the model's **training**
  scores (feasible iff σ_train ≥ 0.50, finite bits, maximize bits).
- Oracle threshold: the same sweep on the **test** (post) scores.
- Report: AUC train, AUC test, and (ρ, σ) on the test set at the honest
  threshold; bits (honest) = test bits at the honest threshold;
  bits (oracle) = test bits at the oracle threshold.

### 5.6 The bar (Section 9 preamble; arithmetic on the rebuilt counts)
Imported operating point p_c = 0.0882 (unchanged; companion model).
- Static: (1−ρ*) = odds(p_c)/odds(p); bits = log2(1/(1−ρ*)).
- Model-conditional: solve odds(p)·(1−ρ) = odds(p_c(r)) with
  p_c(r) = 0.0882 + κ(1−r), κ = 0.038, r = (1−p) + (1−ρ)p (sound
  constraint, σ = 1), by bisection on ρ ∈ (0,1) to 1e-12. Report
  ρ*, r, p_c′, bits.
- Sensitivity: same solver at κ ∈ {0, 0.019, 0.038, 0.057}.
- Solver self-check (before use on rebuilt p): at the printed
  p = 3070/15564 it must reproduce the printed record (κ=0 → 1.345 bits;
  κ=0.038 → ρ* = 0.5845, r = 0.885, p_c′ = 0.0926, 1.267 bits;
  κ=0.019 → 1.306; κ=0.057 → 1.233) to printed rounding. Result logged
  in PROGRESS.md.

### 5.7 Section 9.5 intervals
For each Table 8 quantity, from the rebuild's **actual** 2×2 counts on
that quantity's basis set (no assumed prevalence anywhere): Wilson score
intervals (z = 1.96) for ρ and σ separately; bits interval by the endpoint
transformation [log2(σ_lo/(1−ρ_lo)), log2(σ_hi/(1−ρ_hi))] (the printed
method). AUC interval for mol seen pre by the Hanley–McNeil standard
error. Table 8 keeps the printed table's rows, with rebuilt values.
**Cluster bootstrap** (the printed §9.5's outstanding item, now computable
from the deposited edge data): B = 2000 resamples, seed 20261030,
clusters = molecule (on each quantity's basis set) and, separately,
target (8 clusters — reported with that caveat); percentile 2.5/97.5
intervals for the bits of the Table 8 bits-rows; resamples with a single
label class or degenerate (ρ ∉ (0,1) or σ = 0) are skipped and counted.

## 6. Manuscript rewrite rules

- Output: `files/Free_Verification_rebuilt.tex` (the current submission
  copy is NOT overwritten). Compiled to `files/Free_Verification_rebuilt.pdf`
  via latex-on-http, article class.
- Replaced: every number that depends on the §9 execution — Tables 4–8,
  the four AUCs, the bar and its sensitivity range, the corpus counts and
  p everywhere they are quoted (abstract, §1, §3.2's ChEMBL-bar sentence,
  §4/§5 corpus asides, §9, §10 ratios), and prose claims that quote or
  summarize those numbers (they follow the rebuilt values, whatever they
  are). Preserved: the theory (Sections 3–8 except the §3.2 ChEMBL-bar
  sentence), Tables 1–3, the argument and structure.
- §9/Reproducibility state the rebuild's true design facts (release,
  cutoff, threshold, the 8 target IDs/names, counts, deposit) and close
  reproducibility gaps (1), (2), (3), (6). Gaps (4), (5), (7) concern the
  companion manuscripts/simulator, are outside this rebuild, and remain
  listed as open.
- A numbers-diff (old .tex vs rebuilt .tex) must show changes only where
  the re-execution reaches; the diff summary goes in REBUILD_REPORT.md.

## 7. Validation (success criterion)

`rebuild/validate.py`, authored against the actual artifact layout,
recomputes **every** rebuilt number from the deposited artifacts
(edges, descriptors, results) and diffs them against the values printed
in `Free_Verification_rebuilt.tex`. OVERALL PASS from an actual run is
required. Near-miss negative checks (silent agreement = failure):
(i) threshold perturbed to 6.5 must change the post potent count;
(ii) dropping the first selected target must change the corpus counts;
(iii) with labels shuffled under a fixed seed, the Table 4 potent≥1 bits
must satisfy |bits| < 0.2 (no signal may survive label destruction).

## 8. Boundaries

No submission anywhere; no GitHub push; no fees; sole author Ranadhir
Ghosh. Disk discipline: extracted tables + code + results only; the
chemreps download is deleted after SMILES extraction (its SHA is
recorded above and in PROGRESS.md).
