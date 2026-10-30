#!/usr/bin/env python3
"""Extract canonical SMILES for corpus molecules from chembl_36_chemreps.

Input : chemreps.tmp.gz in this package's root — a scratch copy of the
        ChEMBL 36 chembl_36_chemreps.txt.gz release file (SHA-256 verified
        against the release checksums before this runs; pass its path as
        argv[1] if it lives elsewhere). corpus = union of mol_ids in
        data/edges_pre.csv and data/edges_post.csv.
Output: data/smiles.csv (mol_id, canonical_smiles) for molecules found.
The chemreps scratch file is not part of the deposit: it was deleted
after extraction (pre-reg s.8), so this step is provenance code — the
deposited output, data/smiles.csv, is what the pipeline consumes.
"""
import csv
import gzip
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REBUILD, "chemreps.tmp.gz")

need = set()
for name in ("edges_pre.csv", "edges_post.csv"):
    need |= set(pd.read_csv(os.path.join(REBUILD, "data", name),
                            usecols=["mol_id"]).mol_id)
print(f"corpus molecules: {len(need)}", flush=True)

found = {}
with gzip.open(SRC, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    i_id = header.index("chembl_id")
    i_sm = header.index("canonical_smiles")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        mid = parts[i_id]
        if mid in need and parts[i_sm]:
            found[mid] = parts[i_sm]
print(f"with canonical_smiles in chemreps: {len(found)}")

out = os.path.join(REBUILD, "data", "smiles.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["mol_id", "canonical_smiles"])
    for mid in sorted(found):
        w.writerow([mid, found[mid]])
print(f"wrote {out}")
