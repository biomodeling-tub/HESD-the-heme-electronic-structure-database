#!/usr/bin/env python3
"""
Test script to verify axial ligand sorting is working correctly.
Compares the sorted and unsorted versions of the processed data.
"""

import pandas as pd
import numpy as np
import os
from repo_paths import resolve_table_input

# Decoding dictionaries for axial ligands (from preprocessor.py)
axial1_decode = {0: 'CYS', 1: 'HIS'}
axial2_decode = {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'}

def decode_axial_ligand(value, position):
    """Decode axial ligand value to string."""
    if pd.isna(value):
        return 'NaN'

    decode_dict = axial1_decode if position == 1 else axial2_decode
    return decode_dict.get(int(value), f'Unknown({value})')

def main():
    print("="*80)
    print("AXIAL LIGAND SORTING TEST")
    print("="*80)

    # Check if files exist
    sorted_file = str(resolve_table_input("processed_output_sorted_axial.csv"))
    unsorted_file = str(resolve_table_input("processed_output_unsorted_axial.csv"))
    main_file = str(resolve_table_input("processed_output.csv"))

    files_exist = []
    for f in [sorted_file, unsorted_file, main_file]:
        exists = os.path.exists(f)
        files_exist.append(exists)
        status = "✓ EXISTS" if exists else "✗ NOT FOUND"
        print(f"{status}: {f}")

    if not any(files_exist):
        print("\nERROR: No CSV files found. Run the preprocessor first.")
        return

    print("\n" + "="*80)

    # Load the data
    df_sorted = None
    df_unsorted = None
    df_main = None

    try:
        if os.path.exists(sorted_file):
            df_sorted = pd.read_csv(sorted_file)
            print(f"Loaded sorted data: {len(df_sorted)} rows, {len(df_sorted.columns)} columns")

        if os.path.exists(unsorted_file):
            df_unsorted = pd.read_csv(unsorted_file)
            print(f"Loaded unsorted data: {len(df_unsorted)} rows, {len(df_unsorted.columns)} columns")

        if os.path.exists(main_file):
            df_main = pd.read_csv(main_file)
            print(f"Loaded main output: {len(df_main)} rows, {len(df_main.columns)} columns")
    except Exception as e:
        print(f"\nERROR loading CSV files: {e}")
        return

    # Check if axial columns exist
    if df_main is not None:
        if 'axial1' not in df_main.columns or 'axial2' not in df_main.columns:
            print("\nERROR: axial1 and/or axial2 columns not found in main output!")
            print(f"Available columns: {df_main.columns.tolist()[:10]}...")
            return

    print("\n" + "="*80)
    print("COMPARING SORTED VS UNSORTED DATA")
    print("="*80)

    if df_unsorted is not None and df_sorted is not None:
        # Compare axial ligands
        unsorted_axial1 = df_unsorted['axial1'].values
        unsorted_axial2 = df_unsorted['axial2'].values
        sorted_axial1 = df_sorted['axial1'].values
        sorted_axial2 = df_sorted['axial2'].values

        # Find rows where swapping occurred
        swapped_mask = (unsorted_axial1 != sorted_axial1) | (unsorted_axial2 != sorted_axial2)
        num_swapped = swapped_mask.sum()

        print(f"\nTotal rows: {len(df_unsorted)}")
        print(f"Rows with swapped axial ligands: {num_swapped}")
        print(f"Percentage swapped: {num_swapped/len(df_unsorted)*100:.2f}%")

        if num_swapped > 0:
            print(f"\n{'='*80}")
            print("FIRST 20 SWAPPED ROWS:")
            print(f"{'='*80}")
            print(f"{'Row':<6} {'Before':<20} {'After':<20} {'File Name':<30}")
            print("-"*80)

            swapped_indices = np.where(swapped_mask)[0][:20]
            for idx in swapped_indices:
                unsorted_combo = f"{decode_axial_ligand(unsorted_axial1[idx], 1)}-{decode_axial_ligand(unsorted_axial2[idx], 2)}"
                sorted_combo = f"{decode_axial_ligand(sorted_axial1[idx], 1)}-{decode_axial_ligand(sorted_axial2[idx], 2)}"

                file_name = df_unsorted.iloc[idx]['file_name'] if 'file_name' in df_unsorted.columns else 'N/A'
                if isinstance(file_name, str) and len(file_name) > 27:
                    file_name = file_name[:27] + "..."

                print(f"{idx:<6} {unsorted_combo:<20} {sorted_combo:<20} {file_name:<30}")

            # Summary of combinations
            print(f"\n{'='*80}")
            print("COMBINATION CHANGES:")
            print(f"{'='*80}")

            changes = {}
            for idx in np.where(swapped_mask)[0]:
                before = f"{decode_axial_ligand(unsorted_axial1[idx], 1)}-{decode_axial_ligand(unsorted_axial2[idx], 2)}"
                after = f"{decode_axial_ligand(sorted_axial1[idx], 1)}-{decode_axial_ligand(sorted_axial2[idx], 2)}"
                key = f"{before} → {after}"
                changes[key] = changes.get(key, 0) + 1

            for change, count in sorted(changes.items(), key=lambda x: x[1], reverse=True):
                print(f"{change:<40} {count:>5} occurrences")
        else:
            print("\n⚠ WARNING: No swaps detected!")
            print("This could mean:")
            print("  1. All axial ligands were already in alphabetical order")
            print("  2. The sorting function didn't execute")
            print("  3. The axial columns don't contain the expected data")

    # Verify main output matches sorted version
    print(f"\n{'='*80}")
    print("VERIFYING MAIN OUTPUT (processed_output.csv)")
    print(f"{'='*80}")

    if df_main is not None and df_sorted is not None:
        if df_main['axial1'].equals(df_sorted['axial1']) and df_main['axial2'].equals(df_sorted['axial2']):
            print("✓ Main output contains SORTED axial ligands")
        else:
            print("✗ Main output does NOT match sorted version!")

            # Check if it matches unsorted
            if df_unsorted is not None:
                if df_main['axial1'].equals(df_unsorted['axial1']) and df_main['axial2'].equals(df_unsorted['axial2']):
                    print("✗ Main output still contains UNSORTED axial ligands!")
                    print("\nPossible issue: preprocessor.py may not be replacing self.df correctly")

    # Show distribution of axial combinations in main output
    if df_main is not None:
        print(f"\n{'='*80}")
        print("AXIAL LIGAND COMBINATIONS IN MAIN OUTPUT")
        print(f"{'='*80}")

        combinations = {}
        for idx, row in df_main.iterrows():
            if pd.notna(row['axial1']) and pd.notna(row['axial2']):
                ligand1 = decode_axial_ligand(row['axial1'], 1)
                ligand2 = decode_axial_ligand(row['axial2'], 2)
                combo = f"{ligand1}-{ligand2}"
                combinations[combo] = combinations.get(combo, 0) + 1

        print(f"\n{'Combination':<20} {'Count':<10} {'Sorted?'}")
        print("-"*45)
        for combo, count in sorted(combinations.items(), key=lambda x: x[1], reverse=True):
            parts = combo.split('-')
            if len(parts) == 2:
                is_sorted = "✓" if parts[0] <= parts[1] else "✗ NOT SORTED"
            else:
                is_sorted = "?"
            print(f"{combo:<20} {count:<10} {is_sorted}")

    print("\n" + "="*80)

if __name__ == "__main__":
    main()
