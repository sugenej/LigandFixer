# MolTriage
**Repair and normalize ligand structure files for docking pipelines.**

Software to repair and normalize broken small molecule files for computational chemistry pipelines.
Handles the broken MOL2 aromaticity, wrong charges, missing bonds, and format conversion failures.
---

## Install
```bash
pip install rdkit
pip install openbabel-wheel   # optional: enables MOL2 + PDBQT output

# from source
pip install -e .
```

---
## Quick start

```bash
# CLI
python ligandfix.py ligand.mol2
python ligandfix.py ligand.sdf --output clean.mol2
python ligandfix.py ligand.pdb --output docking_ready.pdbqt --generate-3d
python ligandfix.py ligand.smi --output fixed.sdf --add-hydrogens --generate-3d
```

```python
# Python API
from ligandfixer import fix_ligand

report = fix_ligand(
    "ligand.mol2",
    output_path="ligand_fixed.sdf",
    add_hydrogens=True,
    generate_3d=True,
    standardize=True,
    validate_charge=True,
)
print(report.summary())
```

---

## What it fixes

| Problem | Fix |
|---|---|
| Broken aromatic bonds | Re-perceives aromaticity (RDKit + MDL models) |
| Valence errors | Relaxed property cache update |
| Kekulization failures | Partial sanitization fallback |
| Duplicate atoms at (0,0,0) | Coordinate-based deduplication |
| Lost stereochemistry | Re-assigns from 3D coordinates |
| Salts and solvents | Largest fragment selection |
| Zwitterions | Charge neutralization |
| Missing hydrogens | Explicit H addition |
| Missing 3D coords | ETKDG + MMFF geometry generation |
| Wrong charges | Gasteiger partial charge validation |

---

## Supported formats

| Format | Read | Write |
|--------|------|-------|
| SDF / MOL | ✓ | ✓ |
| MOL2 | ✓ | ✓ (requires OpenBabel) |
| PDB | ✓ | ✓ |
| PDBQT | ✓ | ✓ (requires OpenBabel) |
| SMILES | ✓ | ✓ |

---

## Charge validation

LigandFixer runs a Gasteiger charge consistency check on every output molecule:

```
Charge check:   ✓ PASS  (formal=0, Gasteiger Σ=+0.002)
```

It warns you if:
- Gasteiger sum diverges from formal charge by > 0.5 e
- Any atom has a suspiciously large partial charge (|q| > 0.8)
- Formal charge is extreme (|q| > 4) — usually a parsing error

---

## Example report

```
============================================================
  LigandFixer Report
============================================================
  Input:   ligand.mol2 (MOL2)
  Output:  ligand_fixed.sdf (SDF)
  Status:  ✓ SUCCESS
  Atoms:   24 → 45
  Formal charge:  0 → 0
  Charge check:   ✓ PASS  (formal=0, Gasteiger Σ=+0.001)

  Fixes applied (4):
    + Re-perceived aromaticity (AROMATICITY_RDKIT)
    + Re-assigned stereochemistry from 3D coordinates
    + Kept largest fragment (removed salts/solvents)
    + Added explicit hydrogens
============================================================
```

---

## CLI reference

```
usage: ligandfix.py [-h] [-o OUTPUT] [--output-format {sdf,mol2,pdb,pdbqt,smiles}]
                    [--add-hydrogens] [--no-hydrogens]
                    [--generate-3d] [--no-standardize]
                    [--no-charge-check] [--quiet]
                    input
```
---


## License

MIT
