#!/usr/bin/env python3
"""
Analyze HOH positions in HSD-HOH structures to check if HOH is axial1 or axial2.
This script investigates whether there's inconsistency in axial ligand assignments.
"""

import pandas as pd
import numpy as np

def main():
    # Load the data
    csv_file = "tables/iron_axial_distances.csv"
    df = pd.read_csv(csv_file)

    print(f"Total measurements: {len(df)}")
    print(f"Unique PDB IDs: {df['PDB_ID'].nunique()}")
    print()

    # Find all HSD-HOH structures by checking which PDB IDs have both HSD and HOH
    hsd_hoh_structures = []

    for pdb_id, group in df.groupby('PDB_ID'):
        axial_resnames = set(group['Axial_Resname'].values)

        # Check if this PDB has both HSD and HOH
        if 'HSD' in axial_resnames and 'HOH' in axial_resnames:
            # Should have exactly 2 measurements (one per axial ligand)
            if len(group) == 2:
                hsd_hoh_structures.append(pdb_id)

    print(f"Found {len(hsd_hoh_structures)} structures with HSD-HOH combination")
    print()

    if len(hsd_hoh_structures) == 0:
        print("No HSD-HOH structures found!")
        return

    # Filter to only HSD-HOH structures
    df_hsd_hoh = df[df['PDB_ID'].isin(hsd_hoh_structures)].copy()

    # Analyze positions
    hoh_axial1 = []  # HOH is in position 1
    hoh_axial2 = []  # HOH is in position 2
    hsd_axial1 = []  # HSD is in position 1
    hsd_axial2 = []  # HSD is in position 2

    for pdb_id, group in df_hsd_hoh.groupby('PDB_ID'):
        for _, row in group.iterrows():
            if row['Axial_Resname'] == 'HOH':
                if row['Axial_Ligand'] == 1:
                    hoh_axial1.append({
                        'PDB_ID': pdb_id,
                        'Distance_PDB': row['PDB_Iron_Axial_Distance'],
                        'Distance_XYZ': row['XYZ_Iron_Axial_Distance']
                    })
                elif row['Axial_Ligand'] == 2:
                    hoh_axial2.append({
                        'PDB_ID': pdb_id,
                        'Distance_PDB': row['PDB_Iron_Axial_Distance'],
                        'Distance_XYZ': row['XYZ_Iron_Axial_Distance']
                    })
            elif row['Axial_Resname'] == 'HSD':
                if row['Axial_Ligand'] == 1:
                    hsd_axial1.append({
                        'PDB_ID': pdb_id,
                        'Distance_PDB': row['PDB_Iron_Axial_Distance'],
                        'Distance_XYZ': row['XYZ_Iron_Axial_Distance']
                    })
                elif row['Axial_Ligand'] == 2:
                    hsd_axial2.append({
                        'PDB_ID': pdb_id,
                        'Distance_PDB': row['PDB_Iron_Axial_Distance'],
                        'Distance_XYZ': row['XYZ_Iron_Axial_Distance']
                    })

    # Create summary
    print("="*80)
    print("SUMMARY OF AXIAL LIGAND POSITIONS IN HSD-HOH STRUCTURES")
    print("="*80)
    print()

    print(f"HOH in Axial Position 1: {len(hoh_axial1)} structures")
    print(f"HOH in Axial Position 2: {len(hoh_axial2)} structures")
    print()
    print(f"HSD in Axial Position 1: {len(hsd_axial1)} structures")
    print(f"HSD in Axial Position 2: {len(hsd_axial2)} structures")
    print()

    # This should be the issue: HOH should ALWAYS be in position 2
    if len(hoh_axial1) > 0:
        print("⚠️  WARNING: Found structures where HOH is in axial position 1!")
        print("   HOH should typically be in axial position 2 for HSD-HOH combinations.")
        print()

    # Analyze distances for each case
    print("="*80)
    print("DISTANCE STATISTICS")
    print("="*80)
    print()

    if len(hoh_axial1) > 0:
        hoh_ax1_df = pd.DataFrame(hoh_axial1)
        print("HOH in Axial Position 1:")
        print(f"  PDB Distance - Mean: {hoh_ax1_df['Distance_PDB'].mean():.3f} Å, Std: {hoh_ax1_df['Distance_PDB'].std():.3f} Å")
        print(f"  XYZ Distance - Mean: {hoh_ax1_df['Distance_XYZ'].mean():.3f} Å, Std: {hoh_ax1_df['Distance_XYZ'].std():.3f} Å")
        print(f"  PDB IDs: {sorted(hoh_ax1_df['PDB_ID'].tolist())}")
        print()

    if len(hoh_axial2) > 0:
        hoh_ax2_df = pd.DataFrame(hoh_axial2)
        print("HOH in Axial Position 2:")
        print(f"  PDB Distance - Mean: {hoh_ax2_df['Distance_PDB'].mean():.3f} Å, Std: {hoh_ax2_df['Distance_PDB'].std():.3f} Å")
        print(f"  XYZ Distance - Mean: {hoh_ax2_df['Distance_XYZ'].mean():.3f} Å, Std: {hoh_ax2_df['Distance_XYZ'].std():.3f} Å")
        print(f"  PDB IDs: {sorted(hoh_ax2_df['PDB_ID'].tolist())}")
        print()

    if len(hsd_axial1) > 0:
        hsd_ax1_df = pd.DataFrame(hsd_axial1)
        print("HSD in Axial Position 1:")
        print(f"  PDB Distance - Mean: {hsd_ax1_df['Distance_PDB'].mean():.3f} Å, Std: {hsd_ax1_df['Distance_PDB'].std():.3f} Å")
        print(f"  XYZ Distance - Mean: {hsd_ax1_df['Distance_XYZ'].mean():.3f} Å, Std: {hsd_ax1_df['Distance_XYZ'].std():.3f} Å")
        print(f"  PDB IDs: {sorted(hsd_ax1_df['PDB_ID'].tolist())}")
        print()

    if len(hsd_axial2) > 0:
        hsd_ax2_df = pd.DataFrame(hsd_axial2)
        print("HSD in Axial Position 2:")
        print(f"  PDB Distance - Mean: {hsd_ax2_df['Distance_PDB'].mean():.3f} Å, Std: {hsd_ax2_df['Distance_PDB'].std():.3f} Å")
        print(f"  XYZ Distance - Mean: {hsd_ax2_df['Distance_XYZ'].mean():.3f} Å, Std: {hsd_ax2_df['Distance_XYZ'].std():.3f} Å")
        print(f"  PDB IDs: {sorted(hsd_ax2_df['PDB_ID'].tolist())}")
        print()

    # Analysis of the plotting code behavior
    print("="*80)
    print("IMPACT ON PLOTTING")
    print("="*80)
    print()
    print("In create_actual_distances_histogram:")
    print("  - Subplot 1 (left) shows distances for 'Axial_Ligand == 1'")
    print("  - Subplot 2 (right) shows distances for 'Axial_Ligand == 2'")
    print()
    print(f"Current situation:")
    print(f"  - Subplot 1 will contain {len(hsd_axial1)} HSD distances and {len(hoh_axial1)} HOH distances")
    print(f"  - Subplot 2 will contain {len(hsd_axial2)} HSD distances and {len(hoh_axial2)} HOH distances")
    print()

    if len(hoh_axial1) > 0:
        print("⚠️  ISSUE IDENTIFIED:")
        print(f"  - {len(hoh_axial1)} HOH ligands are in axial position 1")
        print(f"  - These HOH distances will be plotted in the LEFT subplot (Axial 1)")
        print(f"  - But the legend says 'HIS-HOH' (sorted alphabetically)")
        print(f"  - This could make it appear that HIS has longer distances than expected")
        print(f"  - when in fact those are HOH distances!")
        print()
        print("The sorting in line 1227 of calculate_rmsd.py:")
        print("  return '-'.join(sorted([norm_axial1, norm_axial2]))")
        print("  ^^ This hides which ligand is actually in which position!")

    # Check if the low HIS distances mentioned by user are actually from structures
    # where HOH is axial1 and HIS is axial2
    print()
    print("="*80)
    print("EXAMINING THE 'LOW HIS DISTANCE' ISSUE")
    print("="*80)
    print()

    if len(hoh_axial1) > 0 and len(hsd_axial2) > 0:
        print("Structures where HOH is axial1 and HSD is axial2:")
        print("(These HSD distances appear in the RIGHT subplot)")
        print()

        for pdb_id in hoh_axial1:
            # Find corresponding HSD entry
            hsd_entry = [x for x in hsd_axial2 if x['PDB_ID'] == pdb_id['PDB_ID']]
            if hsd_entry:
                print(f"  {pdb_id['PDB_ID']}:")
                print(f"    HOH (axial1): PDB={pdb_id['Distance_PDB']:.3f} Å, XYZ={pdb_id['Distance_XYZ']:.3f} Å")
                print(f"    HSD (axial2): PDB={hsd_entry[0]['Distance_PDB']:.3f} Å, XYZ={hsd_entry[0]['Distance_XYZ']:.3f} Å")

if __name__ == "__main__":
    main()
