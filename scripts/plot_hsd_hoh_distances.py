#!/usr/bin/env python3
"""
Plot iron-axial distances for PDB structures with HSD-HOH axial ligand combinations.
Creates histograms comparing PDB and XYZ distances.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Load the data
    csv_file = "tables/iron_axial_distances.csv"
    df = pd.read_csv(csv_file)

    print(f"Loaded {len(df)} total measurements from {csv_file}")
    print(f"Unique PDB IDs: {df['PDB_ID'].nunique()}")

    # Group by PDB_ID and find structures with HSD-HOH combination
    hsd_hoh_pdb_ids = []

    for pdb_id, group in df.groupby('PDB_ID'):
        axial_resnames = set(group['Axial_Resname'].values)

        # Check if this PDB has both HSD and HOH (and only these two)
        if 'HSD' in axial_resnames and 'HOH' in axial_resnames:
            # Check if it has exactly 2 measurements (one per axial ligand)
            if len(group) == 2:
                hsd_hoh_pdb_ids.append(pdb_id)

    print(f"\nFound {len(hsd_hoh_pdb_ids)} PDB structures with HSD-HOH combination")

    if len(hsd_hoh_pdb_ids) == 0:
        print("No structures found with HSD-HOH combination!")
        return

    # Filter the dataframe for these PDB IDs
    df_filtered = df[df['PDB_ID'].isin(hsd_hoh_pdb_ids)].copy()

    print(f"Total measurements for HSD-HOH structures: {len(df_filtered)}")

    # Separate HSD and HOH measurements
    df_hsd = df_filtered[df_filtered['Axial_Resname'] == 'HSD']
    df_hoh = df_filtered[df_filtered['Axial_Resname'] == 'HOH']

    print(f"\nHSD measurements: {len(df_hsd)}")
    print(f"HOH measurements: {len(df_hoh)}")

    # Get all distance values for histograms
    pdb_distances = df_filtered['PDB_Iron_Axial_Distance'].values
    xyz_distances = df_filtered['XYZ_Iron_Axial_Distance'].values

    # Calculate statistics
    print("\n" + "="*80)
    print("STATISTICS FOR HSD-HOH STRUCTURES")
    print("="*80)

    print("\nPDB Iron-Axial Distances:")
    print(f"  Mean: {pdb_distances.mean():.3f} Å")
    print(f"  Std:  {pdb_distances.std():.3f} Å")
    print(f"  Min:  {pdb_distances.min():.3f} Å")
    print(f"  Max:  {pdb_distances.max():.3f} Å")

    print("\nXYZ Iron-Axial Distances:")
    print(f"  Mean: {xyz_distances.mean():.3f} Å")
    print(f"  Std:  {xyz_distances.std():.3f} Å")
    print(f"  Min:  {xyz_distances.min():.3f} Å")
    print(f"  Max:  {xyz_distances.max():.3f} Å")

    # Create histogram plot
    fig, ax = plt.subplots(figsize=(12, 7))

    # Define bins for histogram
    min_dist = min(pdb_distances.min(), xyz_distances.min())
    max_dist = max(pdb_distances.max(), xyz_distances.max())
    bins = np.linspace(min_dist - 0.1, max_dist + 0.1, 40)

    # Plot histograms with transparency
    ax.hist(pdb_distances, bins=bins, alpha=0.6, label=f'PDB Distances (n={len(pdb_distances)})',
            color='#1f77b4', edgecolor='black')
    ax.hist(xyz_distances, bins=bins, alpha=0.6, label=f'XYZ Distances (n={len(xyz_distances)})',
            color='#ff7f0e', edgecolor='black')

    # Formatting
    ax.set_xlabel('Iron-Axial Distance (Å)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax.set_title('Iron-Axial Ligand Distances for HSD-HOH Combinations\n' +
                 f'({len(hsd_hoh_pdb_ids)} PDB structures)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)

    # Add text box with statistics
    stats_text = (f"PDB: {pdb_distances.mean():.2f} ± {pdb_distances.std():.2f} Å\n"
                  f"XYZ: {xyz_distances.mean():.2f} ± {xyz_distances.std():.2f} Å")
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Save plot
    output_file = "plots/hsd_hoh_iron_axial_distances.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to: {output_file}")

    # Create a second plot: separate histograms for HSD and HOH
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # HSD distances
    pdb_hsd = df_hsd['PDB_Iron_Axial_Distance'].values
    xyz_hsd = df_hsd['XYZ_Iron_Axial_Distance'].values

    bins_hsd = np.linspace(min(pdb_hsd.min(), xyz_hsd.min()) - 0.05,
                           max(pdb_hsd.max(), xyz_hsd.max()) + 0.05, 30)

    ax1.hist(pdb_hsd, bins=bins_hsd, alpha=0.6, label=f'PDB (n={len(pdb_hsd)})',
             color='#1f77b4', edgecolor='black')
    ax1.hist(xyz_hsd, bins=bins_hsd, alpha=0.6, label=f'XYZ (n={len(xyz_hsd)})',
             color='#ff7f0e', edgecolor='black')
    ax1.set_xlabel('Iron-HSD Distance (Å)', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=13, fontweight='bold')
    ax1.set_title('HSD Axial Ligand', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.tick_params(axis='both', labelsize=11)

    stats_hsd = (f"PDB: {pdb_hsd.mean():.2f} ± {pdb_hsd.std():.2f} Å\n"
                 f"XYZ: {xyz_hsd.mean():.2f} ± {xyz_hsd.std():.2f} Å")
    ax1.text(0.02, 0.98, stats_hsd, transform=ax1.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # HOH distances
    pdb_hoh = df_hoh['PDB_Iron_Axial_Distance'].values
    xyz_hoh = df_hoh['XYZ_Iron_Axial_Distance'].values

    bins_hoh = np.linspace(min(pdb_hoh.min(), xyz_hoh.min()) - 0.05,
                           max(pdb_hoh.max(), xyz_hoh.max()) + 0.05, 30)

    ax2.hist(pdb_hoh, bins=bins_hoh, alpha=0.6, label=f'PDB (n={len(pdb_hoh)})',
             color='#1f77b4', edgecolor='black')
    ax2.hist(xyz_hoh, bins=bins_hoh, alpha=0.6, label=f'XYZ (n={len(xyz_hoh)})',
             color='#ff7f0e', edgecolor='black')
    ax2.set_xlabel('Iron-HOH Distance (Å)', fontsize=13, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=13, fontweight='bold')
    ax2.set_title('HOH Axial Ligand', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='both', labelsize=11)

    stats_hoh = (f"PDB: {pdb_hoh.mean():.2f} ± {pdb_hoh.std():.2f} Å\n"
                 f"XYZ: {xyz_hoh.mean():.2f} ± {xyz_hoh.std():.2f} Å")
    ax2.text(0.02, 0.98, stats_hoh, transform=ax2.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig2.suptitle('Iron-Axial Distances by Ligand Type (HSD-HOH Combinations)',
                  fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save second plot
    output_file2 = "plots/hsd_hoh_iron_axial_distances_separated.png"
    plt.savefig(output_file2, dpi=300, bbox_inches='tight')
    print(f"✓ Separated plot saved to: {output_file2}")

    # Show summary of PDB IDs
    print(f"\n{'='*80}")
    print(f"PDB IDs with HSD-HOH combination ({len(hsd_hoh_pdb_ids)} total):")
    print(f"{'='*80}")
    for i, pdb_id in enumerate(sorted(hsd_hoh_pdb_ids), 1):
        if i % 10 == 0:
            print()  # New line every 10 IDs
        print(f"{pdb_id:6s}", end=" ")
    print("\n" + "="*80)

if __name__ == "__main__":
    main()
