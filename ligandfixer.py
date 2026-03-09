"""
LigandFixer — Ligand structure repair and format normalization tool.
Fixes broken SDF, MOL2, PDB, PDBQT, SMILES files for docking pipelines.
"""

import os
import sys
import copy
from dataclasses import dataclass, field
from typing import Optional
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdmolops, SanitizeMol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

# Suppress RDKit noise
RDLogger.DisableLog("rdApp.*")


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class FixReport:
    input_file: str
    input_format: str
    success: bool
    fixes_applied: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    atoms_before: int = 0
    atoms_after: int = 0
    charge_before: Optional[int] = None
    charge_after: Optional[int] = None

    def summary(self) -> str:
        lines = [
            f"\n{'='*55}",
            f"  LigandFixer Report",
            f"{'='*55}",
            f"  Input:   {self.input_file} ({self.input_format.upper()})",
            f"  Status:  {'✓ Success' if self.success else '✗ Failed'}",
            f"  Atoms:   {self.atoms_before} → {self.atoms_after}",
        ]
        if self.charge_before is not None:
            lines.append(f"  Charge:  {self.charge_before:+d} → {self.charge_after:+d}")
        if self.fixes_applied:
            lines.append(f"\n  Fixes applied ({len(self.fixes_applied)}):")
            for f in self.fixes_applied:
                lines.append(f"    + {f}")
        if self.warnings:
            lines.append(f"\n  Warnings:")
            for w in self.warnings:
                lines.append(f"    ! {w}")
        if self.errors:
            lines.append(f"\n  Errors:")
            for e in self.errors:
                lines.append(f"    ✗ {e}")
        lines.append(f"{'='*55}\n")
        return "\n".join(lines)


# ── Format detection ─────────────────────────────────────────────────────────

def detect_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    return {
        ".sdf": "sdf", ".mol": "sdf",
        ".mol2": "mol2",
        ".pdb": "pdb",
        ".pdbqt": "pdbqt",
        ".smi": "smiles", ".smiles": "smiles",
    }.get(ext, "sdf")


# ── Readers ──────────────────────────────────────────────────────────────────

def _read_sdf(path: str) -> tuple[Optional[Chem.Mol], list[str]]:
    warnings = []
    # Try strict first
    supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=False)
    mol = next(iter(supplier), None)
    if mol is not None:
        return mol, warnings
    # Fallback: no sanitize
    warnings.append("Standard SDF parsing failed — attempting lax read")
    supplier = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)
    mol = next(iter(supplier), None)
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol, warnings


def _read_mol2(path: str) -> tuple[Optional[Chem.Mol], list[str]]:
    warnings = []
    mol = Chem.MolFromMol2File(path, sanitize=True, removeHs=False)
    if mol is not None:
        return mol, warnings
    warnings.append("Standard MOL2 parsing failed — attempting lax read")
    mol = Chem.MolFromMol2File(path, sanitize=False, removeHs=False)
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol, warnings


def _read_pdb(path: str) -> tuple[Optional[Chem.Mol], list[str]]:
    warnings = []
    mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=False)
    if mol is not None:
        return mol, warnings
    warnings.append("Standard PDB parsing failed — attempting lax read")
    mol = Chem.MolFromPDBFile(path, sanitize=False, removeHs=False)
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol, warnings


def _read_pdbqt(path: str) -> tuple[Optional[Chem.Mol], list[str]]:
    """Strip PDBQT-specific lines and read as PDB."""
    warnings = ["PDBQT: stripping AutoDock fields to read as PDB"]
    pdb_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM", "CONECT", "END")):
                pdb_lines.append(line[:80])  # truncate PDBQT extras
    pdb_block = "".join(pdb_lines)
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol, warnings


def _read_smiles(path: str) -> tuple[Optional[Chem.Mol], list[str]]:
    with open(path) as fh:
        smi = fh.readline().strip().split()[0]
    mol = Chem.MolFromSmiles(smi)
    return mol, []


def read_molecule(path: str, fmt: Optional[str] = None) -> tuple[Optional[Chem.Mol], str, list[str]]:
    fmt = fmt or detect_format(path)
    readers = {
        "sdf": _read_sdf,
        "mol2": _read_mol2,
        "pdb": _read_pdb,
        "pdbqt": _read_pdbqt,
        "smiles": _read_smiles,
    }
    reader = readers.get(fmt)
    if reader is None:
        return None, fmt, [f"Unsupported format: {fmt}"]
    mol, warnings = reader(path)
    return mol, fmt, warnings


# ── Repair pipeline ──────────────────────────────────────────────────────────

def _fix_valence(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Fix valence errors by adjusting radical/charge."""
    try:
        mol.UpdatePropertyCache(strict=True)
        return mol
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
        report.fixes_applied.append("Relaxed valence update (strict=False)")
    except Exception as e:
        report.warnings.append(f"Valence update failed: {e}")
    return mol


def _fix_aromaticity(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Re-perceive aromaticity using multiple models and pick the best."""
    for model in [
        Chem.AromaticityModel.AROMATICITY_RDKIT,
        Chem.AromaticityModel.AROMATICITY_MDL,
    ]:
        try:
            cp = copy.copy(mol)
            Chem.SetAromaticity(cp, model)
            Chem.SanitizeMol(cp)
            report.fixes_applied.append(
                f"Re-perceived aromaticity ({model.name})"
            )
            return cp
        except Exception:
            continue
    report.warnings.append("Aromaticity re-perception failed — keeping original")
    return mol


def _fix_kekulization(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Attempt to kekulize if molecule has unset aromatic bonds."""
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception:
        report.warnings.append("Kekulization issue — may affect downstream tools")
    return mol


def _remove_duplicate_atoms(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Remove atoms at identical coordinates (common in bad PDB conversions)."""
    conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
    if conf is None:
        return mol
    seen = set()
    to_remove = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        key = (round(pos.x, 3), round(pos.y, 3), round(pos.z, 3))
        if key in seen:
            to_remove.append(atom.GetIdx())
        else:
            seen.add(key)
    if to_remove:
        ed = Chem.RWMol(mol)
        for idx in sorted(to_remove, reverse=True):
            ed.RemoveAtom(idx)
        report.fixes_applied.append(f"Removed {len(to_remove)} duplicate atom(s)")
        return ed.GetMol()
    return mol


def _fix_stereo(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Re-assign stereo from 3D coordinates if available."""
    if mol.GetNumConformers() > 0:
        try:
            Chem.AssignStereochemistryFrom3D(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            report.fixes_applied.append("Re-assigned stereochemistry from 3D coordinates")
        except Exception as e:
            report.warnings.append(f"Stereo assignment failed: {e}")
    else:
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass
    return mol


def _standardize(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Apply RDKit MolStandardize pipeline."""
    try:
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol2 = lfc.choose(mol)
        if mol2.GetNumAtoms() != mol.GetNumAtoms():
            report.fixes_applied.append("Kept largest fragment (removed salts/solvents)")
        mol = mol2
    except Exception:
        pass
    try:
        uncharger = rdMolStandardize.Uncharger()
        mol2 = uncharger.uncharge(mol)
        before = Chem.GetFormalCharge(mol)
        after = Chem.GetFormalCharge(mol2)
        if before != after:
            report.fixes_applied.append(
                f"Neutralized zwitterionic charges ({before:+d} → {after:+d})"
            )
        mol = mol2
    except Exception:
        pass
    return mol


def _add_hydrogens(mol: Chem.Mol, report: FixReport, add_hs: bool) -> Chem.Mol:
    if not add_hs:
        return mol
    try:
        mol = Chem.AddHs(mol, addCoords=True)
        report.fixes_applied.append("Added explicit hydrogens")
    except Exception as e:
        report.warnings.append(f"Could not add hydrogens: {e}")
    return mol


def _embed_3d(mol: Chem.Mol, report: FixReport) -> Chem.Mol:
    """Generate 3D coordinates if missing."""
    if mol.GetNumConformers() > 0:
        return mol
    report.warnings.append("No 3D coordinates found — generating with ETKDG")
    try:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol)
            report.fixes_applied.append("Generated 3D coordinates (ETKDG + MMFF)")
        else:
            report.warnings.append("3D embedding failed — outputting 2D/flat structure")
    except Exception as e:
        report.warnings.append(f"3D generation error: {e}")
    return mol


def _final_sanitize(mol: Chem.Mol, report: FixReport) -> Optional[Chem.Mol]:
    """Final sanitization pass — return None if unrecoverable."""
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        # Try partial sanitize flags
        try:
            flags = (
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                Chem.SanitizeFlags.SANITIZE_SYMMRINGS
            )
            Chem.SanitizeMol(mol, flags)
            report.fixes_applied.append("Partial sanitization (skipped strict valence check)")
            return mol
        except Exception as e2:
            report.errors.append(f"Could not sanitize molecule: {e2}")
            return None


# ── Writers ──────────────────────────────────────────────────────────────────

def write_molecule(mol: Chem.Mol, path: str, fmt: Optional[str] = None) -> bool:
    fmt = fmt or detect_format(path)
    try:
        if fmt == "sdf":
            writer = Chem.SDWriter(path)
            writer.write(mol)
            writer.close()
        elif fmt == "pdb":
            Chem.MolToPDBFile(mol, path)
        elif fmt == "smiles":
            smi = Chem.MolToSmiles(mol)
            with open(path, "w") as fh:
                fh.write(smi + "\n")
        else:
            # mol2 and pdbqt: write as SDF and warn
            writer = Chem.SDWriter(path.replace(f".{fmt}", ".sdf"))
            writer.write(mol)
            writer.close()
            return False
        return True
    except Exception:
        return False


# ── Main fix function ─────────────────────────────────────────────────────────

def fix_ligand(
    input_path: str,
    output_path: Optional[str] = None,
    output_format: Optional[str] = None,
    add_hydrogens: bool = True,
    generate_3d: bool = False,
    standardize: bool = True,
) -> FixReport:
    """
    Main entry point. Reads, repairs, and writes a ligand structure.

    Parameters
    ----------
    input_path     : path to input file
    output_path    : path to write fixed file (default: input_stem_fixed.sdf)
    output_format  : sdf | pdb | smiles (default: sdf)
    add_hydrogens  : add explicit H atoms
    generate_3d    : generate 3D coords if missing
    standardize    : apply charge/fragment standardization
    """
    report = FixReport(
        input_file=input_path,
        input_format=detect_format(input_path),
        success=False,
    )

    # ── Read ────────────────────────────────────────────────────────────────
    mol, fmt, read_warnings = read_molecule(input_path)
    report.warnings.extend(read_warnings)

    if mol is None:
        report.errors.append("Could not parse input file — file may be severely corrupted")
        return report

    report.atoms_before = mol.GetNumAtoms()
    report.charge_before = Chem.GetFormalCharge(mol)

    # ── Repair pipeline ──────────────────────────────────────────────────────
    mol = _fix_valence(mol, report)
    mol = _fix_aromaticity(mol, report)
    mol = _fix_kekulization(mol, report)
    mol = _remove_duplicate_atoms(mol, report)
    mol = _fix_stereo(mol, report)

    if standardize:
        mol = _standardize(mol, report)

    if add_hydrogens:
        mol = _add_hydrogens(mol, report, True)

    if generate_3d:
        mol = _embed_3d(mol, report)

    mol = _final_sanitize(mol, report)
    if mol is None:
        return report

    report.atoms_after = mol.GetNumAtoms()
    report.charge_after = Chem.GetFormalCharge(mol)

    # ── Write ────────────────────────────────────────────────────────────────
    if output_path is None:
        stem = os.path.splitext(input_path)[0]
        out_fmt = output_format or "sdf"
        output_path = f"{stem}_fixed.{out_fmt}"

    out_fmt = output_format or detect_format(output_path)
    ok = write_molecule(mol, output_path, out_fmt)
    if not ok:
        report.warnings.append(
            f"MOL2/PDBQT output not fully supported — saved as SDF instead"
        )

    report.success = True
    return report
