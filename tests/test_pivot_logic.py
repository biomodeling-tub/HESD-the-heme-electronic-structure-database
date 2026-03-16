#!/usr/bin/env python3
"""
Test the pivot logic fix for axial ligand combinations.
"""

import pandas as pd

def test_pivot_logic():
    """Test the fixed pivot logic."""
    
    print("Testing Pivot Logic Fix for Axial Ligand Combinations")
    print("=" * 60)
    
    # Create test axial ligand data (simulating the format from process_iron_axial_distances)
    axial_data = pd.DataFrame({
        'PDB_ID': ['1abc', '1abc', '2def', '2def', '3ghi', '3ghi', '4jkl', '4jkl', '5mno', '5mno'],
        'Axial_Ligand': [1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
        'Axial_Resname': ['HIS', 'HIS', 'HIS', 'MET', 'CYS', 'HOH', 'HIS', 'OXY', 'HIS', 'HOH']
    })
    
    print("Original axial ligand data:")
    print(axial_data)
    print()
    
    print("OLD METHOD (problematic):")
    # The old problematic method
    old_method = axial_data.groupby('PDB_ID').agg({
        'Axial_Resname': lambda x: '-'.join(sorted(x.astype(str)))
    }).reset_index()
    old_method.rename(columns={'Axial_Resname': 'Axial_Combo'}, inplace=True)
    
    print("Old method results:")
    print(old_method)
    print(f"Old method unique combinations: {len(old_method['Axial_Combo'].unique())}")
    for combo in sorted(old_method['Axial_Combo'].unique()):
        print(f"  {combo}")
    print()
    
    print("NEW METHOD (fixed):")
    # The new fixed method
    axial_pivot = axial_data.pivot_table(
        index='PDB_ID', 
        columns='Axial_Ligand', 
        values='Axial_Resname', 
        aggfunc='first'
    ).reset_index()
    
    print("Pivot table result:")
    print(axial_pivot)
    print()
    
    # Create proper axial combinations (axial1-axial2)
    if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
        axial_pivot['Axial_Combo'] = axial_pivot[1].astype(str) + '-' + axial_pivot[2].astype(str)
    elif 1 in axial_pivot.columns:
        axial_pivot['Axial_Combo'] = axial_pivot[1].astype(str) + '-Unknown'
    elif 2 in axial_pivot.columns:
        axial_pivot['Axial_Combo'] = 'Unknown-' + axial_pivot[2].astype(str)
    else:
        axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
    
    new_method = axial_pivot[['PDB_ID', 'Axial_Combo']]
    
    print("New method results:")
    print(new_method)
    print(f"New method unique combinations: {len(new_method['Axial_Combo'].unique())}")
    for combo in sorted(new_method['Axial_Combo'].unique()):
        print(f"  {combo}")
    print()
    
    print("COMPARISON:")
    print(f"Old method: {len(old_method['Axial_Combo'].unique())} combinations")
    print(f"New method: {len(new_method['Axial_Combo'].unique())} combinations")
    print()
    
    expected_combinations = {'HIS-HIS', 'HIS-MET', 'CYS-HOH', 'HIS-OXY', 'HIS-HOH'}
    found_combinations = set(new_method['Axial_Combo'].unique())
    
    print("Expected combinations: ", expected_combinations)
    print("Found combinations:   ", found_combinations)
    print()
    
    if found_combinations == expected_combinations:
        print("✓ SUCCESS: New method produces exactly the expected 5 combinations!")
    else:
        print("✗ ISSUE: New method doesn't match expected combinations")
        print(f"Missing: {expected_combinations - found_combinations}")
        print(f"Extra:   {found_combinations - expected_combinations}")
    
    print()
    print("Summary of the fix:")
    print("- OLD: Grouped by PDB_ID and joined ALL Axial_Resname values")
    print("- NEW: Pivoted by Axial_Ligand to separate axial1 and axial2")
    print("- NEW: Creates proper axial1-axial2 combinations")
    print("- This eliminates malformed combinations and matches the dataset structure")

if __name__ == "__main__":
    test_pivot_logic()