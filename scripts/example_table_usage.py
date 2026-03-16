#!/usr/bin/env python3
"""
Example usage of the new flexible table creation system.
This shows how to easily create different types of tables for different energy variables.
"""

import pandas as pd
from create_latex_table import (
    EnergyVariable, TableType, create_energy_table, create_ground_state_table,
    create_total_energy_comprehensive_table, create_lft_energy_base_statistical_table
)

def example_usage():
    """Example of how to use the new table creation system."""
    
    # Load your data
    df = pd.read_csv("tables/processed_output.csv")
    
    # ========================================================================
    # Method 1: Use the convenience functions (simplest approach)
    # ========================================================================
    print("Method 1: Using convenience functions")
    
    # Create total energy comprehensive table (min/max/mean±std)
    create_total_energy_comprehensive_table(df)
    
    # Create LFT energy delta base statistical table (mean±std only)
    create_lft_energy_base_statistical_table(df)
    
    # ========================================================================
    # Method 2: Use the flexible system for custom configurations
    # ========================================================================
    print("\\nMethod 2: Using flexible system for custom configurations")
    
    # Example 1: Create a min-max table for kinetic energy
    kinetic_energy_var = EnergyVariable(
        column_name="overview_KE_au",
        display_name="Kinetic Energy",
        decimal_places=6,
        use_baseline_normalization=True,
        baseline_type="pdb_min"  # Use PDB minimum as baseline
    )
    create_energy_table(df, kinetic_energy_var, TableType.MIN_MAX)
    
    # Example 2: Create a base statistical table for nuclear-nuclear repulsion
    nuclear_repulsion_var = EnergyVariable(
        column_name="overview_N-N_au",
        display_name="Nuclear-Nuclear Repulsion",
        decimal_places=6,
        use_baseline_normalization=False  # Use raw values, no baseline
    )
    create_energy_table(df, nuclear_repulsion_var, TableType.BASE_STATISTICAL)
    
    # Example 3: Create a comprehensive table for electron-nuclear attraction
    electron_nuclear_var = EnergyVariable(
        column_name="overview_E-N_au",
        display_name="Electron-Nuclear Attraction",
        decimal_places=6,
        use_baseline_normalization=True,
        baseline_type="charge_0_mult_1"  # Use charge=0, multiplicity=1 as reference
    )
    create_energy_table(df, electron_nuclear_var, TableType.COMPREHENSIVE)
    
    # Example 4: Create ground state table for total energy
    total_energy_var = EnergyVariable(
        column_name="total_energies[0]",
        display_name="Total Energy",
        decimal_places=6
    )
    create_ground_state_table(df, total_energy_var)
    
    print("\\nAll tables created successfully!")

def jupyter_notebook_usage():
    """
    This is what you would put in your Jupyter notebook cell to create just 
    the two tables you requested.
    """
    
    # Load your data (assuming 'df' is already loaded in your notebook)
    # df = pd.read_csv("tables/processed_output.csv")  # Uncomment if needed
    
    # Import the functions
    from create_latex_table import (
        create_total_energy_comprehensive_table, 
        create_lft_energy_base_statistical_table
    )
    
    # Create the two requested tables
    print("Creating total energy comprehensive table...")
    create_total_energy_comprehensive_table(df)
    
    print("Creating LFT energy delta base statistical table...")
    create_lft_energy_base_statistical_table(df)
    
    print("Table creation complete!")

if __name__ == "__main__":
    # Run the example
    example_usage()