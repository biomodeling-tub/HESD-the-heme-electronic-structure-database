#!/usr/bin/env python3
"""
Test script to verify the iron-plane analysis fix works correctly.
This script tests the new single histogram with axial ligand coloring.
"""

import pandas as pd
from calculate_rmsd import RMSDAnalyzer

def test_iron_plane_analysis_fix():
    """Test the fixed iron-plane analysis function."""
    
    print("Testing Iron-Plane Analysis Fix")
    print("=" * 50)
    
    # Read the preprocessed data that should only have 5 axial combinations
    try:
        df_plots = pd.read_csv("tables/processed_output.csv")
        print(f"Loaded preprocessed data: {len(df_plots)} structures")
    except FileNotFoundError:
        try:
            df_plots = pd.read_csv("tables/DB.csv")
            print(f"Loaded raw data: {len(df_plots)} structures")
            print("WARNING: Using raw data - you may see more than 5 axial combinations")
        except FileNotFoundError:
            print("ERROR: No data file found. Please ensure processed_output.csv or DB.csv exists.")
            return
    
    # Create PDB_ID column if it doesn't exist
    if 'PDB_ID' not in df_plots.columns:
        df_plots['PDB_ID'] = df_plots['file_name'].str[:4]
    
    # Create RMSDAnalyzer instance
    analyzer = RMSDAnalyzer(verbose=True)
    
    # Get iron-plane distance data
    print("\n1. Processing iron-plane distances...")
    df_plane = analyzer.process_iron_plane_distances()
    
    if df_plane.empty:
        print("ERROR: No iron-plane distance data found")
        return
    
    print(f"   Found iron-plane data for {len(df_plane)} structures")
    
    # Get filtered PDB IDs from the preprocessed dataset
    valid_pdb_ids = df_plots['PDB_ID'].unique()
    plane_pdb_ids = df_plane['PDB_ID'].unique()
    
    # Find intersection of PDB IDs
    common_pdb_ids = set(valid_pdb_ids) & set(plane_pdb_ids)
    print(f"   PDB IDs in preprocessed data: {len(valid_pdb_ids)}")
    print(f"   PDB IDs with iron-plane data: {len(plane_pdb_ids)}")
    print(f"   Common PDB IDs: {len(common_pdb_ids)}")
    
    # Filter iron-plane data to only include preprocessed structures
    df_plane_filtered = df_plane[df_plane['PDB_ID'].isin(valid_pdb_ids)]
    print(f"   Filtered iron-plane data: {len(df_plane_filtered)} structures")
    
    # Get axial data for the filtered PDB IDs
    print("\n2. Processing axial ligand data for filtered PDB IDs...")
    axial_df = analyzer.process_iron_axial_distances(pdb_ids_filter=list(common_pdb_ids))
    
    if axial_df.empty:
        print("ERROR: No axial ligand data found")
        return
    
    print(f"   Found axial data for {len(axial_df)} measurements")
    
    # Count unique axial combinations
    unique_combos = axial_df.groupby('PDB_ID')['Axial_Resname'].apply(
        lambda x: '-'.join(sorted(x.astype(str)))
    ).unique()
    print(f"   Unique axial combinations: {len(unique_combos)}")
    print("   Combinations found:")
    for i, combo in enumerate(sorted(unique_combos)):
        count = sum(axial_df.groupby('PDB_ID')['Axial_Resname'].apply(
            lambda x: '-'.join(sorted(x.astype(str)))
        ) == combo)
        print(f"     {i+1}. {combo}: {count} structures")
    
    # Test the new comprehensive analysis function
    print("\n3. Testing new comprehensive analysis function...")
    try:
        analyzer.create_comprehensive_iron_plane_analysis(
            df_plane_filtered, 
            output_file="plots/test_iron_plane_analysis.png",
            color_by_axial_ligands=True,
            axial_data=axial_df
        )
        print("   ✓ Function executed successfully!")
        print("   ✓ Plot saved as: plots/test_iron_plane_analysis.png")
        
    except Exception as e:
        print(f"   ✗ Function failed with error: {e}")
        import traceback
        traceback.print_exc()
        return
    
    print("\n" + "=" * 50)
    print("TEST SUMMARY:")
    print(f"- Preprocessed structures: {len(df_plots)}")
    print(f"- Iron-plane data: {len(df_plane_filtered)} structures")
    print(f"- Axial combinations found: {len(unique_combos)}")
    print(f"- Expected axial combinations: 5")
    
    if len(unique_combos) == 5:
        print("✓ SUCCESS: Found exactly 5 axial combinations as expected!")
    else:
        print(f"⚠ WARNING: Found {len(unique_combos)} combinations instead of 5")
        print("  This might indicate that the data filtering isn't working as expected")
    
    print("\nThe new function should now:")
    print("- Create a single histogram (not 3 separate plots)")
    print("- Show iron-plane distance differences only")
    print("- Color each bar by axial ligand combination proportions")
    print("- Use only the filtered dataset (not all PDB files)")

if __name__ == "__main__":
    test_iron_plane_analysis_fix()