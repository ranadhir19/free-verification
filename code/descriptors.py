#!/usr/bin/env python3
"""RDKit descriptors for corpus molecules (pre-reg s.4).

Input : data/smiles.csv
Output: data/descriptors.csv — one row per molecule that RDKit parses to a
        valid molecule (a molecule with no row has NO resolved structure),
        mol_id + the 12 deposited descriptors. Parse failures are counted
        in results/descriptors_meta.json.
"""
import json
import os

import pandas as pd
from rdkit import RDLogger, __version__ as rdkit_version
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.dirname(HERE)

from rdkit import Chem  # noqa: E402

smiles = pd.read_csv(os.path.join(REBUILD, "data", "smiles.csv"))
rows = []
n_fail = 0
for mid, smi in zip(smiles.mol_id, smiles.canonical_smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        n_fail += 1
        continue
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    rows.append({
        "mol_id": mid, "MolWt": mw, "MolLogP": logp,
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "NumHDonors": hbd, "NumHAcceptors": hba,
        "NumRotatableBonds": Lipinski.NumRotatableBonds(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "NumAromaticRings": Lipinski.NumAromaticRings(mol),
        "FractionCSP3": Lipinski.FractionCSP3(mol),
        "HeavyAtomCount": mol.GetNumHeavyAtoms(),
        "QED": QED.qed(mol), "LipinskiViolations": viol,
    })
out = pd.DataFrame(rows)
out.to_csv(os.path.join(REBUILD, "data", "descriptors.csv"), index=False)
meta = {"rdkit_version": rdkit_version, "smiles_rows": int(len(smiles)),
        "parsed_ok": int(len(out)), "parse_failures": int(n_fail)}
with open(os.path.join(REBUILD, "results", "descriptors_meta.json"), "w") as f:
    json.dump(meta, f, indent=2)
print(meta)
