#!/usr/bin/env python3
"""
Test script for the new iron displacement barplot function.
This script demonstrates how to create a bar plot showing Fe displacement from the porphyrin plane,
with bars colored by charge-multiplicity states.
"""

import pandas as pd
from calculate_rmsd import RMSDAnalyzer
from repo_paths import resolve_table_input

def main():
    # Initialize analyzer
    analyzer = RMSDAnalyzer(verbose=True)

    # Load iron plane distance data
    print("Loading iron plane distance data...")
    df_plane = pd.read_csv(resolve_table_input("iron_plane_distances.csv"))
    print(f"Loaded {len(df_plane)} structures with iron plane distance data")

    # Load charge-multiplicity data from processed_output.csv
    print("\nLoading charge-multiplicity data...")
    processed_df = pd.read_csv(resolve_table_input("processed_output.csv"))

    # Extract PDB_ID and charge-multiplicity information
    # The file_name format is like "1a2f01.log" where first 4 chars are PDB_ID
    processed_df['PDB_ID'] = processed_df['file_name'].str[:4]

    # Select relevant columns for charge-multiplicity data
    charge_mult_data = processed_df[['PDB_ID', 'charge', 'multiplicity']].drop_duplicates()
    print(f"Loaded charge-multiplicity data for {len(charge_mult_data)} unique PDB IDs")

    # Display available charge-multiplicity combinations
    print("\nAvailable charge-multiplicity combinations:")
    charge_mult_counts = processed_df.groupby(['charge', 'multiplicity']).size().reset_index(name='count')
    for _, row in charge_mult_counts.iterrows():
        print(f"  q={int(row['charge'])}, m={int(row['multiplicity'])}: {row['count']} structures")

    # Create the barplot with default parameters
    print("\n" + "="*80)
    print("Creating iron displacement barplot with charge-multiplicity coloring...")
    print("="*80)
    analyzer.create_iron_displacement_barplot_by_charge_mult(
        df=df_plane,
        charge_mult_data=charge_mult_data,
        output_file="plots/iron_displacement_charge_mult_bars.png",
        bin_width=0.01,  # 0.01 Angstrom bins
        max_distance=0.12  # Maximum distance of 0.12 Angstroms
    )

    # Create a second version with larger bins for comparison
    print("\n" + "="*80)
    print("Creating version with larger bins (0.02 Å)...")
    print("="*80)
    analyzer.create_iron_displacement_barplot_by_charge_mult(
        df=df_plane,
        charge_mult_data=charge_mult_data,
        output_file="plots/iron_displacement_charge_mult_bars_wide.png",
        bin_width=0.02,  # 0.02 Angstrom bins
        max_distance=0.12
    )

    # Create the normalized distribution plot
    print("\n" + "="*80)
    print("Creating NORMALIZED distribution plot...")
    print("This shows what percentage of each charge-multiplicity group falls into each bin")
    print("="*80)
    analyzer.create_iron_displacement_normalized_by_group(
        df=df_plane,
        charge_mult_data=charge_mult_data,
        output_file="plots/iron_displacement_normalized_by_group.png",
        bin_width=0.01,
        max_distance=0.12
    )

    # Create normalized version with wider bins
    print("\n" + "="*80)
    print("Creating NORMALIZED distribution plot with wider bins (0.02 Å)...")
    print("="*80)
    analyzer.create_iron_displacement_normalized_by_group(
        df=df_plane,
        charge_mult_data=charge_mult_data,
        output_file="plots/iron_displacement_normalized_by_group_wide.png",
        bin_width=0.02,
        max_distance=0.12
    )

    print("\n" + "="*80)
    print("Done! Check the plots directory for the generated figures:")
    print("  - iron_displacement_charge_mult_bars.png (absolute counts, stacked)")
    print("  - iron_displacement_charge_mult_bars_wide.png (absolute counts, wider bins)")
    print("  - iron_displacement_normalized_by_group.png (normalized percentages)")
    print("  - iron_displacement_normalized_by_group_wide.png (normalized, wider bins)")
    print("="*80)

if __name__ == "__main__":
    main()
