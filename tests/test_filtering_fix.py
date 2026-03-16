#!/usr/bin/env python3
"""
Test script to demonstrate the filtering fix for axial ligand combinations.
"""

import pandas as pd

def test_filtering_fix():
    """Test how filtering axial data fixes the combination explosion."""
    
    print("Testing Filtering Fix for Axial Ligand Combinations")
    print("=" * 60)
    
    # Simulate the FULL PDB directory data (what process_iron_axial_distances returns)
    # This represents ALL structures in the PDB directory, including filtered ones
    full_axial_data = pd.DataFrame({
        'PDB_ID': ['1abc', '1abc', '2def', '2def', '3ghi', '3ghi', '4jkl', '4jkl', '5mno', '5mno',
                   '6xyz', '6xyz', '7pqr', '7pqr', '8stu', '8stu', '9vwx', '9vwx', '0aaa', '0aaa',
                   '1bbb', '1bbb', '2ccc', '2ccc', '3ddd', '3ddd', '4eee', '4eee', '5fff', '5fff',
                   '6ggg', '6ggg'],
        'Axial_Ligand': [1, 2] * 16,  # 16 PDB_IDs × 2 axial ligands each
        'Axial_Resname': ['HIS', 'HIS', 'HIS', 'MET', 'CYS', 'HOH', 'HIS', 'OXY', 'HIS', 'HOH',
                          'TYR', 'PHE', 'ASP', 'GLU', 'ARG', 'LYS', 'THR', 'SER', 'VAL', 'ILE',
                          'LEU', 'ALA', 'GLY', 'PRO', 'TRP', 'ASN', 'GLN', 'MSE', 'CYM', 'HID',
                          'HIE', 'HIP']
    })
    
    # Simulate the filtered iron-plane data (what your df_plots contains)
    # This represents only the structures that passed your filtering pipeline
    iron_plane_data = pd.DataFrame({
        'PDB_ID': ['1abc', '2def', '3ghi', '4jkl', '5mno'],  # Only 5 structures
        'PDB_Iron_Plane_Distance': [0.1, 0.2, 0.15, 0.25, 0.3],
        'XYZ_Iron_Plane_Distance': [0.12, 0.18, 0.14, 0.28, 0.32],
        'Plane_Distance_Difference': [0.02, -0.02, -0.01, 0.03, 0.02]
    })
    
    print(f"Full PDB directory contains {len(full_axial_data)} axial measurements")
    print(f"Filtered dataset contains {len(iron_plane_data)} structures")
    print()
    
    # Show what happens WITHOUT filtering (the current problem)
    print("WITHOUT FILTERING (current problem):")
    full_combinations = full_axial_data.groupby('PDB_ID')['Axial_Resname'].apply(
        lambda x: '-'.join(sorted(x.astype(str)))
    ).unique()
    print(f"  Unique combinations found: {len(full_combinations)}")
    for i, combo in enumerate(sorted(full_combinations)):
        print(f"    {i+1:2d}. {combo}")
    print()
    
    # Show what happens WITH filtering (the fix)
    print("WITH FILTERING (the fix):")
    valid_pdb_ids = iron_plane_data['PDB_ID'].unique()
    filtered_axial_data = full_axial_data[full_axial_data['PDB_ID'].isin(valid_pdb_ids)]
    
    print(f"  Filtered axial data: {len(full_axial_data)} -> {len(filtered_axial_data)} measurements")
    
    filtered_combinations = filtered_axial_data.groupby('PDB_ID')['Axial_Resname'].apply(
        lambda x: '-'.join(sorted(x.astype(str)))
    ).unique()
    
    print(f"  Unique combinations after filtering: {len(filtered_combinations)}")
    for i, combo in enumerate(sorted(filtered_combinations)):
        print(f"    {i+1}. {combo}")
    print()
    
    print("COMPARISON:")
    print(f"  Before filtering: {len(full_combinations)} combinations")
    print(f"  After filtering:  {len(filtered_combinations)} combinations")
    print(f"  Reduction:        {len(full_combinations) - len(filtered_combinations)} combinations removed")
    print()
    
    if len(filtered_combinations) == 5:
        print("✓ SUCCESS: Filtering reduces to exactly 5 combinations!")
    else:
        print(f"✗ Issue: Still have {len(filtered_combinations)} combinations after filtering")
    
    print()
    print("This explains why you saw 16 combinations:")
    print("- The RMSD analyzer processes ALL PDB files in the directory")
    print("- Your df_plots was filtered to specific structures with specific axial ligands")
    print("- The plotting function was using unfiltered axial data")
    print("- The fix ensures only axial data matching your filtered structures is used")

if __name__ == "__main__":
    test_filtering_fix()