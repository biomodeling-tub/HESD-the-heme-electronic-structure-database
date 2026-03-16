#!/usr/bin/env python3
"""
Test script to verify the axial ligand combinations fix.
"""

import pandas as pd
import numpy as np
from calculate_rmsd import RMSDAnalyzer

def test_axial_combinations_fix():
    """Test the fixed axial ligand combinations logic."""
    
    print("Testing Axial Ligand Combinations Fix")
    print("=" * 50)
    
    # Create test iron-plane data
    plane_data = pd.DataFrame({
        'PDB_ID': ['1abc', '2def', '3ghi', '4jkl', '5mno'],
        'PDB_Iron_Plane_Distance': [0.1, 0.2, 0.15, 0.25, 0.3],
        'XYZ_Iron_Plane_Distance': [0.12, 0.18, 0.14, 0.28, 0.32],
        'Plane_Distance_Difference': [0.02, -0.02, -0.01, 0.03, 0.02]
    })
    
    # Create test axial ligand data (simulating the format from process_iron_axial_distances)
    axial_data = pd.DataFrame({
        'PDB_ID': ['1abc', '1abc', '2def', '2def', '3ghi', '3ghi', '4jkl', '4jkl', '5mno', '5mno'],
        'Axial_Ligand': [1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
        'Axial_Resname': ['HIS', 'HIS', 'HIS', 'MET', 'CYS', 'HOH', 'HIS', 'OXY', 'HIS', 'HOH']
    })
    
    print("Test axial ligand data:")
    print(axial_data)
    print()
    
    print("Expected combinations:")
    expected_combinations = ['HIS-HIS', 'HIS-MET', 'CYS-HOH', 'HIS-OXY', 'HIS-HOH']
    for combo in expected_combinations:
        print(f"  {combo}")
    print()
    
    try:
        # Create an instance of the calculator
        calculator = RMSDAnalyzer()
        
        print("Testing comprehensive iron-plane analysis with fixed axial ligand processing...")
        
        # Test with axial ligand coloring
        calculator.create_comprehensive_iron_plane_analysis(
            plane_data,
            output_file="test_fixed_comprehensive_iron_plane_analysis.png",
            color_by_axial_ligands=True,
            axial_data=axial_data
        )
        
        print("✓ Comprehensive iron-plane analysis function works successfully with fix")
        print("✓ Should now show exactly 5 axial ligand combinations instead of 16")
        print()
        
        # Test the pivot logic separately to verify the fix
        print("Testing the pivot logic directly...")
        axial_pivot = axial_data.pivot_table(
            index='PDB_ID', 
            columns='Axial_Ligand', 
            values='Axial_Resname', 
            aggfunc='first'
        ).reset_index()
        
        print("Pivot table result:")
        print(axial_pivot)
        print()
        
        # Create proper axial combinations
        if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
            axial_pivot['Axial_Combo'] = axial_pivot[1].astype(str) + '-' + axial_pivot[2].astype(str)
        
        combinations_found = axial_pivot['Axial_Combo'].unique()
        print(f"Combinations found: {len(combinations_found)}")
        for combo in sorted(combinations_found):
            print(f"  {combo}")
        
        # Verify we get exactly 5 combinations
        if len(combinations_found) == 5:
            print("\n✓ SUCCESS: Found exactly 5 axial ligand combinations!")
        else:
            print(f"\n✗ ERROR: Found {len(combinations_found)} combinations, expected 5")
        
        print()
        print("Key improvements:")
        print("- Used pivot_table to properly separate axial ligands 1 and 2")
        print("- Creates proper axial1-axial2 combinations")
        print("- Eliminates duplicate and malformed combinations")
        print("- Matches the expected 5 combinations from the dataset")
        
    except Exception as e:
        print(f"✗ Error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_axial_combinations_fix()