#!/usr/bin/env python3
"""
Debug script to identify why Fe coordination analysis isn't getting color coding.
"""

import pandas as pd
import os
from pathlib import Path

def debug_merge_issue():
    """Debug the merge issue for Fe coordination analysis."""
    
    # Load the Fe coordination data
    fe_coord_path = "/home/pbuser/Desktop/PhD_WORK/heme/prior_analysis/qm_analysis_tables/fe_coordination_analysis.csv"
    fe_coord_df = pd.read_csv(fe_coord_path)
    
    print("Fe Coordination Analysis DataFrame:")
    print(f"Shape: {fe_coord_df.shape}")
    print(f"Columns: {list(fe_coord_df.columns)}")
    print(f"First few filenames: {fe_coord_df['filename'].head().tolist()}")
    print()
    
    # Load the main DB.csv auxiliary data
    db_path = "/home/pbuser/Desktop/PhD_WORK/heme/tables/DB.csv"
    if os.path.exists(db_path):
        with open(db_path, 'r') as f:
            header = f.readline().strip().split(',')
        
        # Find charge/multiplicity/axial columns
        charge_cols = [col for col in header if 'charge' in col.lower()]
        mult_cols = [col for col in header if 'multiplicity' in col.lower()]
        axial_cols = [col for col in header if 'axial' in col.lower()]
        
        print("DB.csv analysis:")
        print(f"Charge columns: {charge_cols}")
        print(f"Multiplicity columns: {mult_cols}")
        print(f"Axial columns: {axial_cols}")
        print()
        
        # Check if the DB.csv has file_name column and some sample entries
        try:
            db_df = pd.read_csv(db_path, usecols=['file_name'] + charge_cols + mult_cols + axial_cols)
            print(f"DB.csv loaded successfully with {len(db_df)} rows")
            print(f"Sample file_name entries: {db_df['file_name'].head().tolist()}")
            print()
        except Exception as e:
            print(f"Error loading DB.csv: {e}")
    
    # Check the single_hemes.csv file which might have better auxiliary data
    single_hemes_path = "/home/pbuser/Desktop/PhD_WORK/heme/tables/single_hemes.csv"
    if os.path.exists(single_hemes_path):
        single_hemes_df = pd.read_csv(single_hemes_path)
        print("Single Hemes DataFrame:")
        print(f"Shape: {single_hemes_df.shape}")
        print(f"Columns: {list(single_hemes_df.columns)}")
        if '# PDB' in single_hemes_df.columns:
            print(f"Sample PDB entries: {single_hemes_df['# PDB'].head().tolist()}")
        print()
    
    # Check the filtered_axial_combinations.csv file
    axial_path = "/home/pbuser/Desktop/PhD_WORK/heme/tables/filtered_axial_combinations.csv"
    if os.path.exists(axial_path):
        axial_df = pd.read_csv(axial_path)
        print("Filtered Axial Combinations DataFrame:")
        print(f"Shape: {axial_df.shape}")
        print(f"Columns: {list(axial_df.columns)}")
        print(f"Sample PDB_ID entries: {axial_df['PDB_ID'].head().tolist()}")
        print()
    
    # Try to match filename patterns
    print("Filename matching analysis:")
    sample_fe_filename = fe_coord_df['filename'].iloc[0]
    sample_fe_base = sample_fe_filename.replace('.json', '.log')
    sample_fe_pdb = fe_coord_df['pdb_id'].iloc[0]
    
    print(f"Sample Fe coordination filename: {sample_fe_filename}")
    print(f"Sample Fe coordination base: {sample_fe_base}")
    print(f"Sample Fe coordination PDB ID: {sample_fe_pdb}")
    
    # See if we can find matches in the auxiliary data
    if os.path.exists(db_path):
        try:
            db_df = pd.read_csv(db_path, usecols=['file_name'])
            matching_files = db_df[db_df['file_name'] == sample_fe_base]
            print(f"Matches in DB.csv for {sample_fe_base}: {len(matching_files)}")
        except Exception as e:
            print(f"Error checking DB.csv matches: {e}")
    
    if os.path.exists(single_hemes_path):
        matching_pdb = single_hemes_df[single_hemes_df['# PDB'] == sample_fe_pdb]
        print(f"Matches in single_hemes.csv for {sample_fe_pdb}: {len(matching_pdb)}")
        if len(matching_pdb) > 0:
            print(f"Matching axial1: {matching_pdb['axial1'].iloc[0]}, axial2: {matching_pdb['axial2'].iloc[0]}")
    
    if os.path.exists(axial_path):
        matching_axial = axial_df[axial_df['PDB_ID'] == sample_fe_pdb]
        print(f"Matches in filtered_axial_combinations.csv for {sample_fe_pdb}: {len(matching_axial)}")
        if len(matching_axial) > 0:
            print(f"Matching axials: {matching_axial['Normalized_Axial1'].iloc[0]}-{matching_axial['Normalized_Axial2'].iloc[0]}")

if __name__ == "__main__":
    debug_merge_issue()