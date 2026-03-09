#!/usr/bin/env python3
"""
LigandFixer CLI

Usage:
  python ligandfix.py input.mol2
  python ligandfix.py input.sdf --output clean.sdf
  python ligandfix.py input.pdb --output fixed.sdf --add-hydrogens --generate-3d
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from ligandfixer import fix_ligand


def main():
    parser = argparse.ArgumentParser(
        prog="ligandfix",
        description="LigandFixer — repair broken ligand files for docking pipelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ligandfix.py ligand.mol2
  python ligandfix.py ligand.sdf --output clean.sdf
  python ligandfix.py ligand.pdb --output fixed.sdf --add-hydrogens --generate-3d
  python ligandfix.py ligand.pdbqt --no-standardize

Supported input formats:  .sdf .mol .mol2 .pdb .pdbqt .smi .smiles
Supported output formats: .sdf .pdb .smi
        """,
    )

    parser.add_argument("input", help="Input ligand file path")
    parser.add_argument("-o", "--output", help="Output file path (default: <input>_fixed.sdf)")
    parser.add_argument(
        "--output-format",
        choices=["sdf", "pdb", "smiles"],
        default="sdf",
        help="Output format (default: sdf)",
    )
    parser.add_argument(
        "--add-hydrogens",
        action="store_true",
        default=True,
        help="Add explicit hydrogens (default: on)",
    )
    parser.add_argument(
        "--no-hydrogens",
        dest="add_hydrogens",
        action="store_false",
        help="Do not add explicit hydrogens",
    )
    parser.add_argument(
        "--generate-3d",
        action="store_true",
        default=False,
        help="Generate 3D coordinates if missing",
    )
    parser.add_argument(
        "--no-standardize",
        action="store_false",
        dest="standardize",
        default=True,
        help="Skip charge/fragment standardization",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress report output",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    report = fix_ligand(
        input_path=args.input,
        output_path=args.output,
        output_format=args.output_format,
        add_hydrogens=args.add_hydrogens,
        generate_3d=args.generate_3d,
        standardize=args.standardize,
    )

    if not args.quiet:
        print(report.summary())

    sys.exit(0 if report.success else 1)


if __name__ == "__main__":
    main()
