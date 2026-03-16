#!/usr/bin/env python3
"""
Print statistics for HOMO, LUMO, and HOMO-LUMO gap values.

This script:
1. Loads data from processed_output.csv
2. Identifies all HOMO, LUMO, and gap columns
3. Converts values from Hartree to eV
4. Prints mean, min, and max statistics for each column
"""

import pandas as pd
import numpy as np
import re
from pathlib import Path
from repo_paths import resolve_table_input

# Conversion constant (CODATA 2018)
HARTREE_TO_EV = 27.211386245988

def identify_energy_columns(df):
    """
    Identify all HOMO, LUMO, and gap columns in the dataframe.

    Returns:
    --------
    dict : Dictionary with keys 'homo', 'lumo', 'gaps' containing column lists
    """
    columns = {
        'homo': [],
        'lumo': [],
        'gaps': []
    }

    # Patterns for matching columns
    homo_pattern = re.compile(r'^HOMO(-\d+)?$', re.IGNORECASE)
    lumo_pattern = re.compile(r'^LUMO(\+\d+)?$', re.IGNORECASE)
    gap_pattern = re.compile(r'.*gap.*', re.IGNORECASE)

    for col in df.columns:
        if homo_pattern.match(col):
            columns['homo'].append(col)
        elif lumo_pattern.match(col):
            columns['lumo'].append(col)
        elif gap_pattern.match(col):
            # Exclude already-converted eV columns
            if not col.endswith('_eV') and not col.endswith('_ev'):
                columns['gaps'].append(col)

    return columns

def calculate_statistics(df, column):
    """
    Calculate mean, min, max statistics for a column, excluding NaN and -1 values.

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    column : str
        Column name

    Returns:
    --------
    dict : Dictionary with 'mean', 'min', 'max', 'count', 'has_data'
    """
    # Filter out -1 and NaN values
    valid_data = df[column][(df[column] != -1) & df[column].notna()]

    if len(valid_data) == 0:
        return {
            'mean': None,
            'min': None,
            'max': None,
            'count': 0,
            'has_data': False
        }

    return {
        'mean': valid_data.mean(),
        'min': valid_data.min(),
        'max': valid_data.max(),
        'count': len(valid_data),
        'has_data': True
    }

def print_statistics(column_name, stats_ha, stats_ev):
    """
    Print formatted statistics for a column in both Hartree and eV.

    Parameters:
    -----------
    column_name : str
        Name of the column
    stats_ha : dict
        Statistics in Hartree
    stats_ev : dict
        Statistics in eV
    """
    if not stats_ha['has_data']:
        print(f"  {column_name}: No valid data")
        return

    print(f"  {column_name}:")
    print(f"    Hartree: mean = {stats_ha['mean']:>10.6f}, min = {stats_ha['min']:>10.6f}, max = {stats_ha['max']:>10.6f}")
    print(f"    eV:      mean = {stats_ev['mean']:>10.6f}, min = {stats_ev['min']:>10.6f}, max = {stats_ev['max']:>10.6f}")
    print(f"    Valid samples: {stats_ha['count']}")

def main(csv_path=None):
    """
    Main function to load data and print statistics.

    Parameters:
    -----------
    csv_path : str
        Path to the CSV file
    """
    csv_path = str(resolve_table_input('processed_output.csv')) if csv_path is None else csv_path

    # Check if file exists
    if not Path(csv_path).exists():
        print(f"Error: File not found: {csv_path}")
        return

    # Load data
    print(f"Loading data from {csv_path}...")
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows, {len(df.columns)} columns")

    # Identify energy columns
    energy_cols = identify_energy_columns(df)

    print("\n" + "=" * 80)
    print("HOMO, LUMO, AND GAP STATISTICS")
    print("=" * 80)
    print(f"\nConversion factor: 1 Hartree = {HARTREE_TO_EV} eV (CODATA 2018)")

    # Process HOMO columns
    if energy_cols['homo']:
        print("\n" + "-" * 80)
        print("HOMO ORBITALS")
        print("-" * 80)

        # Sort HOMO columns for logical ordering
        homo_sorted = sorted(energy_cols['homo'],
                           key=lambda x: int(re.search(r'-(\d+)', x).group(1)) if '-' in x else 0)

        for col in homo_sorted:
            stats_ha = calculate_statistics(df, col)

            if stats_ha['has_data']:
                # Convert to eV
                valid_data = df[col][(df[col] != -1) & df[col].notna()]
                valid_data_ev = valid_data * HARTREE_TO_EV
                stats_ev = {
                    'mean': valid_data_ev.mean(),
                    'min': valid_data_ev.min(),
                    'max': valid_data_ev.max(),
                    'count': len(valid_data_ev),
                    'has_data': True
                }

                print_statistics(col, stats_ha, stats_ev)
    else:
        print("\nNo HOMO columns found")

    # Process LUMO columns
    if energy_cols['lumo']:
        print("\n" + "-" * 80)
        print("LUMO ORBITALS")
        print("-" * 80)

        # Sort LUMO columns for logical ordering
        lumo_sorted = sorted(energy_cols['lumo'],
                           key=lambda x: int(re.search(r'\+(\d+)', x).group(1)) if '+' in x else 0)

        for col in lumo_sorted:
            stats_ha = calculate_statistics(df, col)

            if stats_ha['has_data']:
                # Convert to eV
                valid_data = df[col][(df[col] != -1) & df[col].notna()]
                valid_data_ev = valid_data * HARTREE_TO_EV
                stats_ev = {
                    'mean': valid_data_ev.mean(),
                    'min': valid_data_ev.min(),
                    'max': valid_data_ev.max(),
                    'count': len(valid_data_ev),
                    'has_data': True
                }

                print_statistics(col, stats_ha, stats_ev)
    else:
        print("\nNo LUMO columns found")

    # Process gap columns
    if energy_cols['gaps']:
        print("\n" + "-" * 80)
        print("HOMO-LUMO GAPS")
        print("-" * 80)

        for col in sorted(energy_cols['gaps']):
            stats_ha = calculate_statistics(df, col)

            if stats_ha['has_data']:
                # Convert to eV
                valid_data = df[col][(df[col] != -1) & df[col].notna()]
                valid_data_ev = valid_data * HARTREE_TO_EV
                stats_ev = {
                    'mean': valid_data_ev.mean(),
                    'min': valid_data_ev.min(),
                    'max': valid_data_ev.max(),
                    'count': len(valid_data_ev),
                    'has_data': True
                }

                print_statistics(col, stats_ha, stats_ev)
    else:
        print("\nNo gap columns found")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total HOMO columns: {len(energy_cols['homo'])}")
    print(f"Total LUMO columns: {len(energy_cols['lumo'])}")
    print(f"Total gap columns: {len(energy_cols['gaps'])}")

    # Collect all columns with data for reference tables
    all_cols = []
    if energy_cols['homo']:
        all_cols.extend([(col, 'HOMO') for col in sorted(energy_cols['homo'],
                        key=lambda x: int(re.search(r'-(\d+)', x).group(1)) if '-' in x else 0)])
    if energy_cols['lumo']:
        all_cols.extend([(col, 'LUMO') for col in sorted(energy_cols['lumo'],
                        key=lambda x: int(re.search(r'\+(\d+)', x).group(1)) if '+' in x else 0)])
    if energy_cols['gaps']:
        all_cols.extend([(col, 'Gap') for col in sorted(energy_cols['gaps'])])

    # Create quick reference table in Hartree
    print("\n" + "-" * 80)
    print("QUICK REFERENCE TABLE (Hartree)")
    print("-" * 80)
    print(f"{'Column':<20} {'Mean (Ha)':>12} {'Min (Ha)':>12} {'Max (Ha)':>12} {'N':>6}")
    print("-" * 80)

    for col, col_type in all_cols:
        stats_ha = calculate_statistics(df, col)
        if stats_ha['has_data']:
            mean_ha = stats_ha['mean']
            min_ha = stats_ha['min']
            max_ha = stats_ha['max']
            count = stats_ha['count']

            print(f"{col:<20} {mean_ha:>12.6f} {min_ha:>12.6f} {max_ha:>12.6f} {count:>6}")

    # Create quick reference table in eV
    print("\n" + "-" * 80)
    print("QUICK REFERENCE TABLE (eV)")
    print("-" * 80)
    print(f"{'Column':<20} {'Mean (eV)':>12} {'Min (eV)':>12} {'Max (eV)':>12} {'N':>6}")
    print("-" * 80)

    for col, col_type in all_cols:
        stats_ha = calculate_statistics(df, col)
        if stats_ha['has_data']:
            valid_data = df[col][(df[col] != -1) & df[col].notna()]
            valid_data_ev = valid_data * HARTREE_TO_EV
            mean_ev = valid_data_ev.mean()
            min_ev = valid_data_ev.min()
            max_ev = valid_data_ev.max()
            count = len(valid_data_ev)

            print(f"{col:<20} {mean_ev:>12.6f} {min_ev:>12.6f} {max_ev:>12.6f} {count:>6}")

    print("=" * 80)

if __name__ == "__main__":
    main()
