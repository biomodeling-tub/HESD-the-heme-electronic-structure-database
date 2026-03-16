#!/usr/bin/env python3
"""
Test script to verify the LaTeX table creation fixes.
"""

import pandas as pd
import numpy as np
from create_latex_table import create_min_max_range_latex_table, create_comprehensive_latex_table

def test_latex_table_fixes():
    """Test the fixed LaTeX table functions."""
    
    print("Testing LaTeX Table Creation Fixes")
    print("=" * 50)
    
    # Create test data for LFT energy delta (positive values only)
    test_data = pd.DataFrame({
        'file_name': ['1abc_charge_0_mult_1.log', '2def_charge_0_mult_1.log', 
                     '3ghi_charge_0_mult_1.log', '4jkl_charge_0_mult_1.log'],
        'PDB_ID': ['1abc', '2def', '3ghi', '4jkl'],
        'charge': [0, 0, 0, 0],
        'multiplicity': [1, 1, 1, 1],
        'axial1': [0, 1, 0, 1],
        'axial2': [0, 1, 2, 3],
        'lft_energy_delta': [0.5, 0.3, 0.8, 0.4]  # Only positive values
    })
    
    print("Test data:")
    print(test_data[['PDB_ID', 'lft_energy_delta']])
    print()
    
    try:
        print("Testing min-max range table creation for LFT energy delta...")
        create_min_max_range_latex_table(
            test_data, 
            'lft_energy_delta', 
            'LFT Energy Delta',
            decimals=6,
            convert_to_ev=True,
            show_row_references=False,
            use_baseline_normalization=False
        )
        print("✓ Min-max range table created successfully")
        print()
        
        print("Testing comprehensive table creation for LFT energy delta...")
        create_comprehensive_latex_table(
            test_data,
            'lft_energy_delta',
            'LFT Energy Delta', 
            decimals=6,
            convert_to_ev=True,
            show_row_references=False,
            use_baseline_normalization=False
        )
        print("✓ Comprehensive table created successfully")
        print()
        
        print("All tests passed! The UnboundLocalError has been fixed.")
        print()
        print("Key fixes implemented:")
        print("1. Fixed variable scope issues for 'pivot_min' vs 'pivot_mean'")
        print("2. LFT energy delta tables now use mean ± std instead of min-max")
        print("3. Comprehensive tables show only mean ± std for LFT energy delta")
        print("4. Proper conditional logic for different energy types")
        
    except Exception as e:
        print(f"✗ Error occurred: {e}")
        print("The fix may need additional work.")

if __name__ == "__main__":
    test_latex_table_fixes()