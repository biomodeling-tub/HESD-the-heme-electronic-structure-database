#!/usr/bin/env python3
"""
Flexible LaTeX table creation system for energy analysis.
Supports multiple table types and energy variables with configurable baseline normalization.
"""

import pandas as pd
import numpy as np
import os
from enum import Enum
from repo_paths import derived_table_path, resolve_table_input

# Global verbosity setting
verbose = False

class TableType(Enum):
    """Enumeration of available table types."""
    MIN_MAX = "min_max"
    COMPREHENSIVE = "comprehensive"
    BASE_STATISTICAL = "base_statistical"
    BASELINE = "baseline"

class EnergyVariable:
    """Configuration for an energy variable."""
    
    def __init__(self, column_name, display_name, decimal_places=6, 
                 use_baseline_normalization=True, baseline_type="pdb_min"):
        self.column_name = column_name
        self.display_name = display_name
        self.decimal_places = decimal_places
        self.use_baseline_normalization = use_baseline_normalization
        self.baseline_type = baseline_type  # "pdb_min" or "charge_0_mult_1"

def format_scientific_latex(value, decimals=2):
    """Format a number in scientific notation for LaTeX display."""
    if pd.isna(value) or value is None:
        return "--"
    
    # Handle zero separately
    if abs(value) < 1e-15:
        return "0"
    
    # Convert to scientific notation
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10 ** exponent)
    
    # Format based on exponent magnitude
    if -2 <= exponent <= 3:
        # Use regular notation for reasonable range
        return f"{value:.{decimals}f}"
    else:
        # Use scientific notation for very large or very small numbers
        return f"{mantissa:.{decimals}f} \\times 10^{{{exponent}}}"

def prepare_energy_data(df, energy_var, convert_to_ev=True):
    """
    Prepare and normalize energy data for table creation.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing energy data
    energy_var : EnergyVariable
        Energy variable configuration
    convert_to_ev : bool
        If True, convert from atomic units to electron volts
        
    Returns:
    --------
    pd.DataFrame : Processed DataFrame with normalized energy values
    """
    HARTREE_TO_EV = 27.211386245988  # CODATA 2018 value
    
    if verbose:
        print(f"Preparing data for {energy_var.display_name} ({energy_var.column_name})")
    
    # Filter out missing values (-1 and NaN)
    df_energy = df[(df[energy_var.column_name] != -1) & df[energy_var.column_name].notna()].copy()
    
    # For LFT energy delta, only include positive values
    if 'lft_energy_delta' in energy_var.column_name.lower():
        df_energy = df_energy[df_energy[energy_var.column_name] > 0].copy()
        if verbose:
            print(f"LFT Energy Delta: Filtered to include only positive values ({len(df_energy)} remaining)")
    
    if df_energy.empty:
        if verbose:
            print(f"No valid {energy_var.column_name} data found after filtering")
        return None
    
    # Extract PDB-ID from filename (first 4 characters)
    if 'PDB_ID' not in df_energy.columns:
        df_energy['PDB_ID'] = df_energy['file_name'].str[:4]
    
    # Create combinations if not already present
    if 'axial_combo' not in df_energy.columns:
        df_energy['axial_combo'] = df_energy['axial1'].astype(str) + '-' + df_energy['axial2'].astype(str)
    if 'charge_mult' not in df_energy.columns:
        df_energy['charge_mult'] = df_energy['charge'].astype(str) + '-' + df_energy['multiplicity'].astype(str)
    
    # Handle axial combinations - use existing 'axials' column if available
    if verbose:
        print(f"DEBUG: Available columns: {list(df_energy.columns)}")
        print(f"DEBUG: 'axials' in columns: {'axials' in df_energy.columns}")
    
    if 'axials' in df_energy.columns:
        # Use existing pre-computed axials column
        if verbose:
            print(f"DEBUG: axials sample before assignment: {df_energy['axials'].head()}")
        df_energy['axial_combo_decoded'] = df_energy['axials']
        if verbose:
            print(f"Using existing 'axials' column with {df_energy['axial_combo_decoded'].nunique()} unique combinations")
            print(f"DEBUG: axial_combo_decoded sample after assignment: {df_energy['axial_combo_decoded'].head()}")
    elif df_energy['axial1'].dtype in ['int8', 'int16', 'int32', 'int64', 'Int64']:
        # Numeric format - needs decoding
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'}
        }
        
        df_energy['axial1_decoded'] = df_energy['axial1'].map(REPLACEMENTS['axial1'])
        df_energy['axial2_decoded'] = df_energy['axial2'].map(REPLACEMENTS['axial2'])
        df_energy['axial_combo_decoded'] = df_energy['axial1_decoded'] + '-' + df_energy['axial2_decoded']
        if verbose:
            print(f"Decoded numeric axials to {df_energy['axial_combo_decoded'].nunique()} unique combinations")
    else:
        # Text format - create from individual columns
        df_energy['axial1_decoded'] = df_energy['axial1'].astype(str)
        df_energy['axial2_decoded'] = df_energy['axial2'].astype(str)
        df_energy['axial_combo_decoded'] = df_energy['axial1_decoded'] + '-' + df_energy['axial2_decoded']
        if verbose:
            print(f"Created axial combinations from text columns: {df_energy['axial_combo_decoded'].nunique()} unique combinations")
    
    # Apply baseline normalization
    df_energy_processed = df_energy.copy()
    baseline_col = f'{energy_var.column_name}_processed'
    
    if energy_var.use_baseline_normalization:
        if energy_var.baseline_type == "charge_0_mult_1":
            # Use charge=0, multiplicity=1 state as reference PER PDB-ID
            processed_values = df_energy[energy_var.column_name].copy()
            processed_pdbs = []
            skipped_pdbs = []
            
            for pdb_id in df_energy['PDB_ID'].unique():
                pdb_mask = df_energy['PDB_ID'] == pdb_id
                pdb_data = df_energy.loc[pdb_mask]
                
                # Find charge=0, multiplicity=1 energy within this specific PDB-ID
                pdb_charge_0_mult_1_mask = (pdb_data['charge'] == 0) & (pdb_data['multiplicity'] == 1)
                pdb_charge_0_mult_1_energies = pdb_data.loc[pdb_charge_0_mult_1_mask, energy_var.column_name]
                
                if len(pdb_charge_0_mult_1_energies) > 0:
                    # Use PDB-specific charge=0, multiplicity=1 energy as reference
                    pdb_baseline_energy = pdb_charge_0_mult_1_energies.iloc[0]
                    
                    # Calculate differences from this PDB's 0-1 state reference
                    pdb_energies = pdb_data[energy_var.column_name]
                    pdb_gauged_energies = pdb_energies - pdb_baseline_energy
                    processed_values.loc[pdb_mask] = pdb_gauged_energies
                    processed_pdbs.append(pdb_id)
                    
                    if verbose:
                        print(f"PDB {pdb_id}: Using 0-1 state as reference: {pdb_baseline_energy:.6f} a.u.")
                else:
                    # No charge=0, multiplicity=1 state found for this PDB-ID - skip it
                    skipped_pdbs.append(pdb_id)
                    if verbose:
                        print(f"Warning: PDB {pdb_id} has no charge=0, multiplicity=1 state. Skipping.")
            
            if skipped_pdbs:
                # Remove PDB-IDs that don't have charge=0, multiplicity=1 states
                skip_mask = df_energy['PDB_ID'].isin(skipped_pdbs)
                df_energy_processed = df_energy_processed[~skip_mask].copy()
                processed_values = processed_values[~skip_mask]
                if verbose:
                    print(f"Skipped {len(skipped_pdbs)} PDB-IDs without 0-1 states: {skipped_pdbs}")
                    print(f"Processing {len(processed_pdbs)} PDB-IDs with valid 0-1 references: {processed_pdbs}")
            
            if len(processed_pdbs) == 0:
                if verbose:
                    print("Warning: No PDB-IDs have charge=0, multiplicity=1 states. Using raw values.")
                processed_values = df_energy[energy_var.column_name].copy()
        
        elif energy_var.baseline_type == "pdb_min":
            # Use PDB-ID minimum as baseline
            processed_values = df_energy[energy_var.column_name].copy()
            for pdb_id in df_energy['PDB_ID'].unique():
                pdb_mask = df_energy['PDB_ID'] == pdb_id
                pdb_energies = df_energy.loc[pdb_mask, energy_var.column_name]
                
                if len(pdb_energies) > 0:
                    min_energy = pdb_energies.min()
                    processed_values.loc[pdb_mask] = pdb_energies - min_energy
    else:
        # Use raw values without baseline normalization
        processed_values = df_energy[energy_var.column_name].copy()
    
    # Convert to eV if requested
    if convert_to_ev:
        processed_values = processed_values * HARTREE_TO_EV
    
    df_energy_processed[baseline_col] = processed_values
    
    return df_energy_processed

def create_table_data(df_processed, energy_var, table_type):
    """
    Create aggregated data for table generation.
    
    Parameters:
    -----------
    df_processed : pd.DataFrame
        Processed DataFrame with normalized energy values
    energy_var : EnergyVariable
        Energy variable configuration
    table_type : TableType
        Type of table to create
        
    Returns:
    --------
    dict : Dictionary containing pivot tables and metadata
    """
    baseline_col = f'{energy_var.column_name}_processed'
    
    # Group by axial combination and charge-multiplicity
    if table_type == TableType.MIN_MAX:
        grouped = df_processed.groupby(['axial_combo_decoded', 'charge_mult'])[baseline_col].agg(['min', 'max', 'count']).reset_index()
        pivot_data = {
            'min': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='min'),
            'max': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='max'),
            'count': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='count')
        }
    
    elif table_type == TableType.COMPREHENSIVE:
        grouped = df_processed.groupby(['axial_combo_decoded', 'charge_mult'])[baseline_col].agg(['min', 'mean', 'std', 'max', 'count']).reset_index()
        pivot_data = {
            'min': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='min'),
            'mean': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='mean'),
            'std': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='std'),
            'max': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='max'),
            'count': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='count')
        }
    
    elif table_type == TableType.BASE_STATISTICAL:
        grouped = df_processed.groupby(['axial_combo_decoded', 'charge_mult'])[baseline_col].agg(['mean', 'std', 'count']).reset_index()
        pivot_data = {
            'mean': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='mean'),
            'std': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='std'),
            'count': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='count')
        }
    
    elif table_type == TableType.BASELINE:
        grouped = df_processed.groupby(['axial_combo_decoded', 'charge_mult'])[baseline_col].agg(['mean', 'std', 'count']).reset_index()
        pivot_data = {
            'mean': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='mean'),
            'std': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='std'),
            'count': grouped.pivot(index='axial_combo_decoded', columns='charge_mult', values='count')
        }
    
    return pivot_data

def generate_latex_table(pivot_data, energy_var, table_type, convert_to_ev=True):
    """
    Generate LaTeX table from pivot data.
    
    Parameters:
    -----------
    pivot_data : dict
        Dictionary containing pivot tables
    energy_var : EnergyVariable
        Energy variable configuration
    table_type : TableType
        Type of table being created
    convert_to_ev : bool
        Whether energies were converted to eV
        
    Returns:
    --------
    list : List of LaTeX table lines
    """
    latex_table = []
    latex_table.append("\\begin{maybesafe}")
    latex_table.append("    \\begin{table}[htbp]")
    latex_table.append("    \\centering")
    
    # Create caption based on table type
    units = "eV" if convert_to_ev else "a.u."
    
    # Determine baseline description
    if energy_var.use_baseline_normalization:
        if energy_var.baseline_type == "charge_0_mult_1":
            baseline_text = "differences relative to the 0-1 state reference energy"
        else:
            baseline_text = "differences relative to the minimum energy within each PDB-ID"
    else:
        baseline_text = "raw values"
    
    # Table type specific descriptions
    if table_type == TableType.MIN_MAX:
        stats_text = "minimum and maximum values"
    elif table_type == TableType.COMPREHENSIVE:
        stats_text = "comprehensive statistics (min/mean±std/max)"
    elif table_type == TableType.BASE_STATISTICAL:
        stats_text = "mean ± standard deviation"
    else:  # BASELINE
        stats_text = "baseline statistics (mean ± std)"
    
    caption = f"{energy_var.display_name} {stats_text} of {baseline_text} for combinations of axial ligands and different charge ($q$) and spin ($s=5/2m$) models (in {units})"
    latex_table.append(f"    \\caption{{{caption}}}")
    
    # Create label
    label_name = energy_var.display_name.lower().replace(' ', '_').replace('-', '_')
    table_suffix = table_type.value
    latex_table.append(f"    \\label{{tab:{label_name}_{table_suffix}}}")
    latex_table.append("    \\tiny{")
    
    # Get sorted columns
    first_pivot = list(pivot_data.values())[0]
    sorted_columns = sorted(first_pivot.columns)
    n_cols = len(sorted_columns)
    
    # Create column specification based on table type
    if table_type == TableType.COMPREHENSIVE:
        col_spec = "l|l|" + "c" * n_cols
        latex_table.append(f"        \\begin{{tabular}}{{{col_spec}}}")
        latex_table.append("            \\hline")
        
        # Header row
        header = "            Axial Combination & Statistic"
        for col in sorted_columns:
            header += f" & {col}"
        header += " \\\\\\\\"
        latex_table.append(header)
        latex_table.append("            \\hline")
        
        # Data rows for comprehensive table
        for idx in sorted(first_pivot.index):
            # Min row
            row_line = f"            \\multirow{{3}}{{*}}{{{idx.replace('_', '\\\\_')}}} & Min"
            for col in sorted_columns:
                min_val = pivot_data['min'].loc[idx, col] if col in pivot_data['min'].columns else None
                if pd.notna(min_val) and min_val is not None:
                    min_sci = format_scientific_latex(min_val, decimals=energy_var.decimal_places)
                    row_line += f" & ${min_sci}$"
                else:
                    row_line += " & --"
            row_line += " \\\\\\\\"
            latex_table.append(row_line)
            
            # Mean±std row
            row_line = "            & Mean±Std"
            for col in sorted_columns:
                mean_val = pivot_data['mean'].loc[idx, col] if col in pivot_data['mean'].columns else None
                std_val = pivot_data['std'].loc[idx, col] if col in pivot_data['std'].columns else None
                count_val = pivot_data['count'].loc[idx, col] if col in pivot_data['count'].columns else None
                
                if pd.notna(mean_val) and pd.notna(std_val) and mean_val is not None:
                    mean_sci = format_scientific_latex(mean_val, decimals=energy_var.decimal_places)
                    std_sci = format_scientific_latex(std_val, decimals=energy_var.decimal_places)
                    row_line += f" & ${mean_sci} \\pm {std_sci}$ ({int(count_val) if count_val is not None else 0})"
                else:
                    row_line += " & --"
            row_line += " \\\\\\\\"
            latex_table.append(row_line)
            
            # Max row
            row_line = "            & Max"
            for col in sorted_columns:
                max_val = pivot_data['max'].loc[idx, col] if col in pivot_data['max'].columns else None
                if pd.notna(max_val) and max_val is not None:
                    max_sci = format_scientific_latex(max_val, decimals=energy_var.decimal_places)
                    row_line += f" & ${max_sci}$"
                else:
                    row_line += " & --"
            row_line += " \\\\\\\\"
            latex_table.append(row_line)
            latex_table.append("            \\hline")
    
    else:
        # Standard single-row format for other table types
        col_spec = "l|" + "c" * n_cols
        latex_table.append(f"        \\begin{{tabular}}{{{col_spec}}}")
        latex_table.append("            \\hline")
        
        # Header row
        header = "            Axial Combination/q-m"
        for col in sorted_columns:
            header += f" & {col}"
        header += " \\\\\\\\"
        latex_table.append(header)
        latex_table.append("            \\hline")
        
        # Data rows
        for idx in sorted(first_pivot.index):
            row_line = f"            {idx.replace('_', '\\\\_')}"
            
            for col in sorted_columns:
                if table_type == TableType.MIN_MAX:
                    min_val = pivot_data['min'].loc[idx, col] if col in pivot_data['min'].columns else None
                    max_val = pivot_data['max'].loc[idx, col] if col in pivot_data['max'].columns else None
                    count_val = pivot_data['count'].loc[idx, col] if col in pivot_data['count'].columns else None
                    
                    if pd.notna(min_val) and pd.notna(max_val) and min_val is not None:
                        min_sci = format_scientific_latex(min_val, decimals=energy_var.decimal_places)
                        max_sci = format_scientific_latex(max_val, decimals=energy_var.decimal_places)
                        row_line += f" & ${min_sci}$ - ${max_sci}$ ({int(count_val) if count_val is not None else 0})"
                    else:
                        row_line += " & --"
                
                else:  # BASE_STATISTICAL or BASELINE
                    mean_val = pivot_data['mean'].loc[idx, col] if col in pivot_data['mean'].columns else None
                    std_val = pivot_data['std'].loc[idx, col] if col in pivot_data['std'].columns else None
                    count_val = pivot_data['count'].loc[idx, col] if col in pivot_data['count'].columns else None
                    
                    if pd.notna(mean_val) and pd.notna(std_val) and mean_val is not None:
                        mean_sci = format_scientific_latex(mean_val, decimals=energy_var.decimal_places)
                        std_sci = format_scientific_latex(std_val, decimals=energy_var.decimal_places)
                        row_line += f" & ${mean_sci} \\pm {std_sci}$ ({int(count_val) if count_val is not None else 0})"
                    else:
                        row_line += " & --"
            
            row_line += " \\\\\\\\"
            latex_table.append(row_line)
        
        latex_table.append("            \\hline")
    
    latex_table.append("        \\end{tabular}")
    latex_table.append("    }")
    latex_table.append("    \\end{table}")
    latex_table.append("\\end{maybesafe}")
    
    return latex_table

def create_energy_table(df, energy_var, table_type, convert_to_ev=True, save_to_file=True):
    """
    Create a LaTeX table for a specific energy variable and table type.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing energy data
    energy_var : EnergyVariable
        Energy variable configuration
    table_type : TableType
        Type of table to create
    convert_to_ev : bool
        If True, convert from atomic units to electron volts
    save_to_file : bool
        If True, save the table to a file
        
    Returns:
    --------
    list : List of LaTeX table lines
    """
    if verbose:
        print(f"Creating {table_type.value} table for {energy_var.display_name}")
    
    # Prepare data
    df_processed = prepare_energy_data(df, energy_var, convert_to_ev)
    if df_processed is None:
        return None
    
    # Create table data
    pivot_data = create_table_data(df_processed, energy_var, table_type)
    
    # Generate LaTeX table
    latex_table = generate_latex_table(pivot_data, energy_var, table_type, convert_to_ev)
    
    # Save to file if requested
    if save_to_file:
        filename = str(derived_table_path(
            f'{energy_var.column_name.replace("-", "_").replace("[", "_").replace("]", "_")}_{table_type.value}_latex_table.txt'
        ))
        with open(filename, 'w') as f:
            for line in latex_table:
                f.write(line + '\n')
        if verbose:
            print(f"Table saved to '{filename}'")
    
    return latex_table

def create_ground_state_table(df, energy_var, convert_to_ev=True, save_to_file=True):
    """
    Create a LaTeX table showing ground states for total energy variables.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing energy data
    energy_var : EnergyVariable
        Energy variable configuration (should be total energy)
    convert_to_ev : bool
        If True, convert from atomic units to electron volts
    save_to_file : bool
        If True, save the table to a file
        
    Returns:
    --------
    list : List of LaTeX table lines
    """
    HARTREE_TO_EV = 27.211386245988
    
    if verbose:
        print(f"Creating ground state table for {energy_var.display_name}")
    
    # Filter out missing values
    df_energy = df[(df[energy_var.column_name] != -1) & df[energy_var.column_name].notna()].copy()
    
    if df_energy.empty:
        if verbose:
            print(f"No valid {energy_var.column_name} data found")
        return None
    
    # Extract PDB-ID and create charge-multiplicity combination
    if 'PDB_ID' not in df_energy.columns:
        df_energy['PDB_ID'] = df_energy['file_name'].str[:4]
    if 'charge_mult' not in df_energy.columns:
        df_energy['charge_mult'] = df_energy['charge'].astype(str) + '-' + df_energy['multiplicity'].astype(str)
    
    # Find reference energy from 0-1 state
    charge_0_mult_1_mask = (df_energy['charge'] == 0) & (df_energy['multiplicity'] == 1)
    charge_0_mult_1_energies = df_energy.loc[charge_0_mult_1_mask, energy_var.column_name]
    
    if len(charge_0_mult_1_energies) == 0:
        if verbose:
            print("Warning: No charge=0, multiplicity=1 states found. Cannot create ground state table.")
        return None
    
    reference_energy = charge_0_mult_1_energies.iloc[0]
    if verbose:
        print(f"Using 0-1 state as reference: {reference_energy:.6f} a.u.")
    
    # Calculate gauged energies for each PDB-ID
    pdb_results = []
    
    for pdb_id in sorted(df_energy['PDB_ID'].unique()):
        pdb_mask = df_energy['PDB_ID'] == pdb_id
        pdb_data = df_energy.loc[pdb_mask].copy()
        
        # Calculate gauged energies
        gauged_energies = pdb_data[energy_var.column_name] - reference_energy
        if convert_to_ev:
            gauged_energies = gauged_energies * HARTREE_TO_EV
        
        pdb_result = {'PDB_ID': pdb_id}
        min_energy = gauged_energies.min()
        
        for idx in pdb_data.index:
            charge_mult = pdb_data.loc[idx, 'charge_mult']
            gauged_energy = gauged_energies.loc[idx]
            is_ground_state = abs(gauged_energy - min_energy) < 1e-6
            
            pdb_result[charge_mult] = {
                'energy': gauged_energy,
                'is_ground_state': is_ground_state
            }
        
        pdb_results.append(pdb_result)
    
    # Get all unique charge-multiplicity combinations
    all_charge_mult = set()
    for result in pdb_results:
        for key in result.keys():
            if key != 'PDB_ID':
                all_charge_mult.add(key)
    
    sorted_charge_mult = sorted(all_charge_mult)
    
    # Create LaTeX table
    latex_table = []
    latex_table.append("\\begin{maybesafe}")
    latex_table.append("    \\begin{table}[htbp]")
    latex_table.append("    \\centering")
    
    units = "eV" if convert_to_ev else "a.u."
    caption = f"{energy_var.display_name} gauged against 0-1 state reference for each PDB-ID. Ground states (lowest energy) are shown in bold (in {units})"
    latex_table.append(f"    \\caption{{{caption}}}")
    
    label_name = energy_var.display_name.lower().replace(' ', '_').replace('-', '_')
    latex_table.append(f"    \\label{{tab:{label_name}_ground_states}}")
    latex_table.append("    \\tiny{")
    
    n_cols = len(sorted_charge_mult)
    col_spec = "l|" + "c" * n_cols
    latex_table.append(f"        \\begin{{tabular}}{{{col_spec}}}")
    latex_table.append("            \\hline")
    
    # Header row
    header = "            PDB-ID"
    for charge_mult in sorted_charge_mult:
        header += f" & {charge_mult}"
    header += " \\\\\\\\"
    latex_table.append(header)
    latex_table.append("            \\hline")
    
    # Data rows
    for result in pdb_results:
        pdb_id = result['PDB_ID']
        row_line = f"            {pdb_id}"
        
        for charge_mult in sorted_charge_mult:
            if charge_mult in result:
                energy_info = result[charge_mult]
                energy_val = energy_info['energy']
                is_ground_state = energy_info['is_ground_state']
                
                energy_sci = format_scientific_latex(energy_val, decimals=energy_var.decimal_places)
                
                if is_ground_state:
                    row_line += f" & $\\mathbf{{{energy_sci}}}$"
                else:
                    row_line += f" & ${energy_sci}$"
            else:
                row_line += " & --"
        
        row_line += " \\\\\\\\"
        latex_table.append(row_line)
    
    latex_table.append("            \\hline")
    latex_table.append("        \\end{tabular}")
    latex_table.append("    }")
    latex_table.append("    \\end{table}")
    latex_table.append("\\end{maybesafe}")
    
    # Save to file if requested
    if save_to_file:
        filename = str(derived_table_path(
            f'{energy_var.column_name.replace("-", "_").replace("[", "_").replace("]", "_")}_ground_states_table.txt'
        ))
        with open(filename, 'w') as f:
            for line in latex_table:
                f.write(line + '\n')
        if verbose:
            print(f"Ground state table saved to '{filename}'")
    
    return latex_table

# Convenience functions for easy table creation
def create_total_energy_comprehensive_table(df):
    """Create comprehensive table for total energy with baseline normalization."""
    total_energy_var = EnergyVariable(
        column_name="total_energies[0]",
        display_name="Total Energy",
        decimal_places=6,
        use_baseline_normalization=True,
        baseline_type="charge_0_mult_1"
    )
    return create_energy_table(df, total_energy_var, TableType.COMPREHENSIVE)

def create_lft_energy_base_statistical_table(df):
    """Create base statistical table for LFT energy delta without baseline normalization."""
    lft_energy_var = EnergyVariable(
        column_name="lft_energy_delta",
        display_name="LFT Energy Delta",
        decimal_places=8,
        use_baseline_normalization=False
    )
    return create_energy_table(df, lft_energy_var, TableType.BASE_STATISTICAL)

def create_lft_energy_baseline_normalized_table(df):
    """Create base statistical table for LFT energy delta with baseline normalization (charge_0_mult_1)."""
    lft_energy_var = EnergyVariable(
        column_name="lft_energy_delta",
        display_name="LFT Energy Delta",
        decimal_places=8,
        use_baseline_normalization=True,
        baseline_type="charge_0_mult_1"
    )
    return create_energy_table(df, lft_energy_var, TableType.BASE_STATISTICAL)


def create_min_max_range_latex_table(df, column_name, display_name,
                                     decimals=6, convert_to_ev=True,
                                     show_row_references=False, save_to_file=True):
    energy_var = EnergyVariable(
        column_name=column_name,
        display_name=display_name,
        decimal_places=decimals,
        use_baseline_normalization=True,
        baseline_type="pdb_min",
    )
    return create_energy_table(df, energy_var, TableType.MIN_MAX, convert_to_ev=convert_to_ev, save_to_file=save_to_file)


def create_comprehensive_latex_table(df, column_name, display_name,
                                     decimals=6, convert_to_ev=True,
                                     show_row_references=False, save_to_file=True):
    energy_var = EnergyVariable(
        column_name=column_name,
        display_name=display_name,
        decimal_places=decimals,
        use_baseline_normalization=True,
        baseline_type="pdb_min",
    )
    return create_energy_table(df, energy_var, TableType.COMPREHENSIVE, convert_to_ev=convert_to_ev, save_to_file=save_to_file)


def create_pdb_baseline_latex_table(df, column_name, display_name,
                                    decimals=6, convert_to_ev=True,
                                    show_row_references=False, save_to_file=True):
    energy_var = EnergyVariable(
        column_name=column_name,
        display_name=display_name,
        decimal_places=decimals,
        use_baseline_normalization=True,
        baseline_type="charge_0_mult_1",
    )
    return create_energy_table(df, energy_var, TableType.BASELINE, convert_to_ev=convert_to_ev, save_to_file=save_to_file)


def create_pdb_ground_state_table(df, column_name, display_name,
                                  decimals=6, convert_to_ev=True,
                                  save_to_file=True):
    energy_var = EnergyVariable(
        column_name=column_name,
        display_name=display_name,
        decimal_places=decimals,
        use_baseline_normalization=True,
        baseline_type="charge_0_mult_1",
    )
    return create_ground_state_table(df, energy_var, convert_to_ev=convert_to_ev, save_to_file=save_to_file)

if __name__ == "__main__":
    # Load data
    data_files = [
        str(resolve_table_input('processed_output.csv')),
        str(resolve_table_input('DB.csv')),
        str(resolve_table_input('data.csv')),
    ]
    df_plots = None
    
    for filename in data_files:
        if os.path.exists(filename):
            try:
                if verbose:
                    print(f"Loading data from {filename}...")
                df_plots = pd.read_csv(filename)
                if verbose:
                    print(f"Successfully loaded {len(df_plots)} rows from {filename}")
                break
            except Exception as e:
                if verbose:
                    print(f"Error loading {filename}: {e}")
                continue
    
    if df_plots is None:
        print("Error: No data file found.")
        exit(1)
    
    # Create the requested tables
    if verbose:
        print("Creating requested tables...")

    # 1. Total electronic energy table with min, max, mean and std
    create_total_energy_comprehensive_table(df_plots)

    # 2. LFT energy delta table with only mean and std (no baseline normalization)
    create_lft_energy_base_statistical_table(df_plots)

    # 3. LFT energy delta table with baseline normalization (charge_0_mult_1)
    create_lft_energy_baseline_normalized_table(df_plots)

    if verbose:
        print("Table creation complete!")
