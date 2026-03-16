#!/usr/bin/env python3
"""
Test script to demonstrate the LFT energy delta changes:
1. Exclusion of negative values in preprocessor.py
2. Positive-only filtering and mean ± std display in create_latex_table.py
"""

import pandas as pd
import numpy as np

def test_lft_changes():
    """Test the LFT energy delta changes."""
    
    print("Testing LFT Energy Delta Changes")
    print("=" * 50)
    
    # Create test data with both positive and negative LFT energy delta values
    test_data = pd.DataFrame({
        'file_name': ['1abc_charge_0_mult_1.log', '2def_charge_0_mult_1.log', 
                     '3ghi_charge_0_mult_1.log', '4jkl_charge_0_mult_1.log'],
        'PDB_ID': ['1abc', '2def', '3ghi', '4jkl'],
        'charge': [0, 0, 0, 0],
        'multiplicity': [1, 1, 1, 1],
        'axial1': [0, 1, 0, 1],
        'axial2': [0, 1, 2, 3],
        'lft_energy_delta': [0.5, -0.2, 0.8, -0.1]  # Mix of positive and negative
    })
    
    print("Original test data:")
    print(test_data[['PDB_ID', 'lft_energy_delta']])
    print()
    
    # Test the filtering behavior
    print("Changes implemented:")
    print()
    
    print("1. PREPROCESSOR.PY:")
    print("   - Added exclude_negative=True flag to calculate_lft_energy_delta()")
    print("   - Negative LFT energy delta values are automatically set to NaN")
    print("   - This effectively excludes them from downstream analysis")
    print()
    
    print("2. CREATE_LATEX_TABLE.PY:")
    print("   - Added positive-only filtering for lft_energy_delta tables")
    print("   - LFT energy delta tables now show only mean ± std (not min/max)")
    print("   - Range tables for LFT energy delta show mean ± std instead of min-max")
    print()
    
    # Simulate the filtering
    filtered_data = test_data.copy()
    filtered_data.loc[filtered_data['lft_energy_delta'] < 0, 'lft_energy_delta'] = np.nan
    
    print("After filtering (negative values → NaN):")
    print(filtered_data[['PDB_ID', 'lft_energy_delta']])
    print()
    
    print("Valid positive values for table creation:")
    valid_data = filtered_data.dropna(subset=['lft_energy_delta'])
    print(valid_data[['PDB_ID', 'lft_energy_delta']])
    print()
    
    if len(valid_data) > 0:
        mean_val = valid_data['lft_energy_delta'].mean()
        std_val = valid_data['lft_energy_delta'].std()
        print(f"Table would show: {mean_val:.3f} ± {std_val:.3f}")
    else:
        print("No valid positive values found")
    
    print()
    print("Key benefits:")
    print("- Ensures only physically meaningful positive LFT values are used")
    print("- Provides cleaner statistical representation (mean ± std)")
    print("- Maintains data integrity by filtering at source (preprocessor)")
    print("- Consistent handling across all table types")

if __name__ == "__main__":
    test_lft_changes()