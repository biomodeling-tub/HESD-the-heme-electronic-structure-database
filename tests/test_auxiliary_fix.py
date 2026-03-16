#!/usr/bin/env python3
"""
Test script to verify the auxiliary data loading fix works.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))

from plot_qm_analysis import QMAnalysisPlotter

def test_auxiliary_loading():
    """Test the fixed auxiliary data loading."""
    
    print("Testing auxiliary data loading fix...")
    plotter = QMAnalysisPlotter()
    
    print(f"Auxiliary data loaded: {not plotter.auxiliary_data.empty}")
    if not plotter.auxiliary_data.empty:
        print(f"Auxiliary data shape: {plotter.auxiliary_data.shape}")
        print(f"Auxiliary data columns: {list(plotter.auxiliary_data.columns)}")
        
        # Test if charge and multiplicity columns exist
        charge_cols = [col for col in plotter.auxiliary_data.columns if 'charge' in col.lower()]
        mult_cols = [col for col in plotter.auxiliary_data.columns if 'multiplicity' in col.lower()]
        axial_cols = [col for col in plotter.auxiliary_data.columns if 'axial' in col.lower()]
        
        print(f"Charge columns found: {charge_cols}")
        print(f"Multiplicity columns found: {mult_cols}")
        print(f"Axial columns found: {axial_cols}")
        
        # Test a sample merge
        if 'fe_coordination_analysis' in plotter.tables:
            fe_df = plotter.tables['fe_coordination_analysis']
            print(f"\nTesting merge with fe_coordination_analysis...")
            print(f"Fe coordination shape: {fe_df.shape}")
            
            merged = plotter.merge_with_auxiliary(fe_df.head())
            print(f"Merged shape: {merged.shape}")
            print(f"Merged columns: {list(merged.columns)}")
            
            # Check if charge_multiplicity column was created
            if 'charge_multiplicity' in merged.columns:
                print(f"Success! charge_multiplicity column created")
                print(f"Sample values: {merged['charge_multiplicity'].unique()}")
            else:
                print("Warning: charge_multiplicity column not created")
                
            # Check if axial columns exist
            axial_check = [col for col in merged.columns if col.startswith('axial')]
            if axial_check:
                print(f"Axial columns in merged data: {axial_check}")
            else:
                print("Warning: No axial columns in merged data")

if __name__ == "__main__":
    test_auxiliary_loading()
