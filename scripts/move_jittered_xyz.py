#!/usr/bin/env python3
"""
Move all .xyz files containing '_jittered_' in their basename
into a 'jittered' subfolder within each PDB subdirectory.
"""

import os
from pathlib import Path
import shutil

def move_jittered_files(pdb_dir="PDB"):
    """
    Move jittered XYZ files into organized subfolders.

    Args:
        pdb_dir: Path to PDB directory (default: "PDB")
    """
    pdb_path = Path(pdb_dir)

    if not pdb_path.exists():
        print(f"Error: Directory {pdb_dir} does not exist")
        return

    moved_count = 0

    # Iterate through all subdirectories in PDB/
    for subdir in pdb_path.iterdir():
        if not subdir.is_dir():
            continue

        # Find all .xyz files with '_jittered_' in the name
        jittered_files = [f for f in subdir.glob("*.xyz") if "_jittered_" in f.name]

        if jittered_files:
            # Create jittered subfolder if it doesn't exist
            jittered_folder = subdir / "jittered"
            jittered_folder.mkdir(exist_ok=True)

            # Move each jittered file
            for xyz_file in jittered_files:
                dest = jittered_folder / xyz_file.name
                shutil.move(str(xyz_file), str(dest))
                print(f"Moved: {xyz_file.relative_to(pdb_path)} -> {dest.relative_to(pdb_path)}")
                moved_count += 1

    print(f"\nTotal files moved: {moved_count}")

if __name__ == "__main__":
    move_jittered_files()
