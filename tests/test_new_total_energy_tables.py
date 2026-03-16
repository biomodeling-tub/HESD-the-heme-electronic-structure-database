#!/usr/bin/env python3
"""
Test script to demonstrate the new total energy table functionality.
"""

import pandas as pd
import numpy as np
from create_latex_table import create_pdb_baseline_latex_table, create_pdb_ground_state_table

def create_test_data():
    """Create test data that mimics real quantum chemistry data."""
    
    # Define test PDB IDs
    pdb_ids = ['1abc', '2def', '3ghi']
    
    # Define charge-multiplicity combinations
    charge_mult_combinations = [
        (0, 1),   # Neutral singlet
        (0, 5),   # Neutral quintet  
        (1, 2),   # +1 doublet
        (1, 6),   # +1 sextet
    ]
    
    # Create test data
    data = []
    for pdb_id in pdb_ids:
        for charge, mult in charge_mult_combinations:
            # Create realistic total energies (in Hartree)
            # The 0-1 state will be our reference
            if charge == 0 and mult == 1:
                base_energy = -2500.0  # Reference energy
            elif charge == 0 and mult == 5:
                base_energy = -2499.8  # Slightly higher
            elif charge == 1 and mult == 2:
                base_energy = -2499.5  # Higher for charged state
            else:  # charge == 1 and mult == 6
                base_energy = -2499.2  # Highest energy
            
            # Add some variation per PDB-ID
            pdb_variation = {'1abc': 0.0, '2def': 0.1, '3ghi': -0.05}[pdb_id]
            total_energy = base_energy + pdb_variation + np.random.normal(0, 0.01)
            
            # Create axial ligand combinations (just for completeness)
            axial1 = 0 if 'a' in pdb_id else 1  # CYS or HIS
            axial2 = hash(pdb_id) % 4  # Random axial2
            
            data.append({
                'file_name': f"{pdb_id}_charge_{charge}_mult_{mult}.log",
                'PDB_ID': pdb_id,
                'charge': charge,
                'multiplicity': mult,
                'total_energies[0]': total_energy,
                'axial1': axial1,
                'axial2': axial2
            })
    
    return pd.DataFrame(data)

def test_new_functionality():
    """Test the new total energy table functionality."""
    
    print("Testing new total energy table functionality")
    print("=" * 60)
    
    # Create test data
    df = create_test_data()
    
    print(f"Created test data with {len(df)} rows")
    print(f"PDB IDs: {sorted(df['PDB_ID'].unique())}")
    print(f"Charge-multiplicity combinations: {sorted(df['charge_mult'].unique())}")
    
    # Display the test data
    print("\nTest data preview:")
    print(df[['PDB_ID', 'charge', 'multiplicity', 'total_energies[0]']].round(6))
    
    # Test the updated baseline table (should now use 0-1 state as reference)
    print("\n" + "=" * 60)
    print("TESTING UPDATED BASELINE TABLE")
    print("(should now use 0-1 state as reference for ALL states)")
    print("=" * 60)
    
    create_pdb_baseline_latex_table(
        df, 
        'total_energies[0]', 
        'Total Energy',
        decimals=6,
        convert_to_ev=True,
        show_row_references=False
    )
    
    # Test the new ground state table
    print("\n" + "=" * 60)
    print("TESTING NEW GROUND STATE TABLE")
    print("(shows gauged energies per PDB-ID with ground states in bold)")
    print("=" * 60)
    
    create_pdb_ground_state_table(
        df,
        'total_energies[0]',
        'Total Energy',
        decimals=6,
        convert_to_ev=True
    )
    
    print("\n" + "=" * 60)
    print("TEST COMPLETED")
    print("=" * 60)
    print("Key changes implemented:")
    print("1. Total energy tables now use the 0-1 state as reference (not average)")
    print("2. All other states are gauged against this single reference energy")
    print("3. New ground state table shows per-PDB-ID energies with lowest in bold")
    print("4. Statistical tables are calculated from these gauged energies")
    print("\nCheck the generated files in the tables/ directory:")
    print("- total_energies_0__pdb_baseline_latex_table.txt")
    print("- total_energies_0__ground_states_table.txt")

if __name__ == "__main__":
    test_new_functionality()