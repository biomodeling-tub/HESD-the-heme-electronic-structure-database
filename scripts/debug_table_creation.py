#!/usr/bin/env python3
"""
Debug script to investigate why table creation is failing.
"""

import pandas as pd
import numpy as np

# Load the same data as the notebook
print("Loading and processing data...")
df = pd.read_csv("tables/DB.csv")

from preprocessor import DataPreprocessor
preprocessor = DataPreprocessor(df=df)
df, subsets_dict, replacements = preprocessor.process()

print(f"Processed dataframe shape: {df.shape}")
print(f"Columns: {list(df.columns)}")

# Check key columns for table creation
print("\n=== DATA STRUCTURE ANALYSIS ===")

# Check if required columns exist
required_cols = ['total_energies[0]', 'lft_energy_delta', 'file_name', 'charge', 'multiplicity']
for col in required_cols:
    if col in df.columns:
        print(f"✓ {col}: {df[col].dtype}, non-null: {df[col].count()}/{len(df)}")
        if col in ['total_energies[0]', 'lft_energy_delta']:
            print(f"  Range: {df[col].min():.6f} to {df[col].max():.6f}")
    else:
        print(f"✗ {col}: MISSING")

# Check axial ligand columns
print(f"\nAxial ligand columns:")
axial_cols = [col for col in df.columns if 'axial' in col.lower()]
for col in axial_cols:
    print(f"  {col}: {df[col].dtype}")
    if df[col].dtype == 'object':
        print(f"    Unique values: {df[col].unique()}")
    else:
        print(f"    Unique values: {sorted(df[col].unique())}")

# Check charge and multiplicity
print(f"\nCharge values: {sorted(df['charge'].unique())}")
print(f"Multiplicity values: {sorted(df['multiplicity'].unique())}")
print(f"Charge-multiplicity combinations:")
combo_counts = df.groupby(['charge', 'multiplicity']).size()
for (charge, mult), count in combo_counts.items():
    print(f"  q={charge}, m={mult}: {count} samples")

# Test the table creation process step by step
print(f"\n=== TESTING TABLE CREATION PROCESS ===")

# Simulate the prepare_energy_data function for total_energies[0]
print(f"\nTesting total_energies[0] processing:")
energy_col = 'total_energies[0]'

# Filter out missing values
df_energy = df[(df[energy_col] != -1) & df[energy_col].notna()].copy()
print(f"After filtering: {len(df_energy)} samples")

# Check PDB-ID extraction
if 'PDB_ID' not in df_energy.columns:
    df_energy['PDB_ID'] = df_energy['file_name'].str[:4]
print(f"Unique PDB-IDs: {df_energy['PDB_ID'].nunique()}")

# Check axial combinations
if 'axial_combo' not in df_energy.columns:
    df_energy['axial_combo'] = df_energy['axial1'].astype(str) + '-' + df_energy['axial2'].astype(str)
if 'charge_mult' not in df_energy.columns:
    df_energy['charge_mult'] = df_energy['charge'].astype(str) + '-' + df_energy['multiplicity'].astype(str)

print(f"Unique axial combinations: {df_energy['axial_combo'].unique()}")
print(f"Unique charge-mult combinations: {df_energy['charge_mult'].unique()}")

# Check charge=0, multiplicity=1 baseline
charge_0_mult_1_mask = (df_energy['charge'] == 0) & (df_energy['multiplicity'] == 1)
charge_0_mult_1_energies = df_energy.loc[charge_0_mult_1_mask, energy_col]
print(f"\nCharge=0, mult=1 samples: {len(charge_0_mult_1_energies)}")
if len(charge_0_mult_1_energies) > 0:
    baseline_energy = charge_0_mult_1_energies.iloc[0]
    print(f"Baseline energy: {baseline_energy:.6f}")
    
    # Calculate differences
    processed_values = df_energy[energy_col] - baseline_energy
    HARTREE_TO_EV = 27.211386245988
    processed_values = processed_values * HARTREE_TO_EV
    df_energy['processed'] = processed_values
    
    print(f"Processed values range: {processed_values.min():.3f} to {processed_values.max():.3f} eV")
    
    # Group by axial combination and charge-multiplicity
    grouped = df_energy.groupby(['axial_combo', 'charge_mult'])['processed'].agg(['min', 'mean', 'std', 'max', 'count']).reset_index()
    print(f"\nGrouped data shape: {grouped.shape}")
    print("First few groups:")
    print(grouped.head(10))
    
    # Check if pivot works
    try:
        pivot_mean = grouped.pivot(index='axial_combo', columns='charge_mult', values='mean')
        print(f"\nPivot table shape: {pivot_mean.shape}")
        print("Pivot table:")
        print(pivot_mean)
    except Exception as e:
        print(f"Pivot failed: {e}")

else:
    print("ERROR: No charge=0, multiplicity=1 samples found!")

print(f"\n=== TESTING LFT ENERGY DELTA ===")
lft_col = 'lft_energy_delta'
if lft_col in df.columns:
    df_lft = df[(df[lft_col] != -1) & df[lft_col].notna()].copy()
    df_lft = df_lft[df_lft[lft_col] > 0].copy()  # Only positive values
    print(f"LFT samples after filtering: {len(df_lft)}")
    
    if len(df_lft) > 0:
        print(f"LFT range: {df_lft[lft_col].min():.6f} to {df_lft[lft_col].max():.6f}")
else:
    print("LFT energy delta column not found!")