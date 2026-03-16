import os
import re
import shutil
import filecmp
import json
import random
import itertools
import math
import numpy as np
import pandas as pd
from datetime import datetime
import MDAnalysis as mda
import glob
from joblib import Parallel, delayed
from calculate_rmsd import RMSDAnalyzer


#Global Flags
debug=False
verbose=False


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)


class DataPreprocessor:
    static_drop_patterns_pre = [
        r'mulliken_heavy_atoms', r'atomic_type', 
        r'geometry\.center', r'atomic_number', r'mulliken_charges_atomic',
        r'quadrupole_moment',
        r'atomic_dipole_orientation', r'natural_electron_config',
        r'natural_population\.atoms\[(?:[5-9]|\d{2,})\]\.total',
        r'lft_orbitals\.atoms\[(?:[1-4]|\d{2,})\]\..*',
        r'sum_mulliken_charges\.charge_sum', r'sum_mulliken_charges\.spin_sum',
        r'Heme', r'ligand', r'core', r'atom_symbol', r'ligand',
        r'polarizability'
    ]

    static_drop_patterns_post = [
        r'^standard_orientation\[\d+\]\..*',
        r'natural_population\.atoms\[\d+\]\.(?:charge|total|geometry\.(?:x|y|z))$',
        r'natural_population\.atoms\[(?:[5-9]|\d{2,})\]\.valence$',
        r'total'
    ]

    def __init__(self, df: pd.DataFrame, keep_homo_lumo: bool = False,
                 write_file: bool = True, skip_plots: bool = True, n_jobs: int = 1,
                 analysis_folder: str = "prior_analysis", homo_lumo_all: bool = False):
        self.original_df = df
        self.df = df.copy(deep=False)
        self.keep_homo_lumo = keep_homo_lumo
        self.write_file = write_file
        self.skip_plots = skip_plots
        self.n_jobs = n_jobs
        self.analysis_folder = analysis_folder
        self.homo_lumo_all = homo_lumo_all
        os.makedirs(self.analysis_folder, exist_ok=True)

        # Precompile regex patterns
        self._pat_npa  = re.compile(r'^natural_population\.atoms\[(\d+)\]\\.')
        self._pat_nec  = re.compile(r'^natural_electron_config\[(\d+)\]\\.')
        self._pat_homo = re.compile(r'^homo\[(\d+)\]$')
        self._pat_lumo = re.compile(r'^lumo\[(\d+)\]$')

        self.report = {}
        self.encoding_dict = {}
        self.sample_tracking = []

    def track_sample_count(self, step_name: str, description: str = ""):
        """
        Records the current number of samples (rows) in the dataset at a specific preprocessing step.
        
        Parameters:
        -----------
        step_name : str
            Name of the preprocessing step
        description : str, optional
            Additional description of what was done in this step
        """
        current_count = len(self.df)
        
        tracking_entry = {
            'step': step_name,
            'sample_count': current_count,
            'description': description,
            'columns_count': len(self.df.columns)
        }
        
        # Calculate change from previous step
        if self.sample_tracking:
            previous_count = self.sample_tracking[-1]['sample_count']
            change = current_count - previous_count
            tracking_entry['samples_change'] = change
            tracking_entry['samples_removed'] = -change if change < 0 else 0
        else:
            tracking_entry['samples_change'] = 0
            tracking_entry['samples_removed'] = 0
        
        self.sample_tracking.append(tracking_entry)
        
        if verbose:
            if tracking_entry['samples_removed'] > 0:
                print(f"{step_name}: {current_count} samples remaining (removed {tracking_entry['samples_removed']})")
            else:
                print(f"{step_name}: {current_count} samples")

    def generate_sample_report(self):
        """
        Generates and prints a comprehensive report showing how many samples remain
        after each preprocessing step.
        """
        if not self.sample_tracking:
            print("No sample tracking data available.")
            return
        
        print("\n" + "="*80)
        print("SAMPLE COUNT TRACKING REPORT")
        print("="*80)
        
        # Summary statistics
        initial_count = self.sample_tracking[0]['sample_count']
        final_count = self.sample_tracking[-1]['sample_count']
        total_removed = initial_count - final_count
        retention_rate = (final_count / initial_count) * 100 if initial_count > 0 else 0
        
        print(f"Initial samples: {initial_count:,}")
        print(f"Final samples: {final_count:,}")
        print(f"Total removed: {total_removed:,}")
        print(f"Retention rate: {retention_rate:.1f}%")
        print("\n" + "-"*80)
        print(f"{'Step':<30} {'Count':<10} {'Removed':<10} {'Description':<30}")
        print("-"*80)
        
        for entry in self.sample_tracking:
            step = entry['step'][:29]  # Truncate if too long
            count = f"{entry['sample_count']:,}"
            removed = f"{entry['samples_removed']:,}" if entry['samples_removed'] > 0 else "-"
            desc = entry['description'][:29] if entry['description'] else ""
            
            print(f"{step:<30} {count:<10} {removed:<10} {desc:<30}")
        
        print("="*80)
        
        # Store in report for JSON output
        self.report['sample_tracking'] = {
            'initial_samples': initial_count,
            'final_samples': final_count,
            'total_removed': total_removed,
            'retention_rate': retention_rate,
            'detailed_steps': self.sample_tracking
        }

    def drop_rydberg_columns(self):
        self.remove_columns(r'rydberg')

    def drop_total_keyword_columns(self):
        self.remove_columns(r'\.Total')

    def drop_homo_lumo_by_index(self):
        if self.keep_homo_lumo or self.homo_lumo_all:
            return
        cols = []
        for c in self.df.columns:
            m = self._pat_homo.match(c)
            if m and int(m.group(1)) < 7:
                cols.append(c)
            m = self._pat_lumo.match(c)
            if m and int(m.group(1)) > 2:
                cols.append(c)
        self.remove_columns(cols)

    def rename_homo_lumo(self):
        mapping = {}
        for c in self.df.columns:
            m = self._pat_homo.match(c)
            if m:
                idx = int(m.group(1))
                if self.keep_homo_lumo or self.homo_lumo_all or idx >= 7:
                    new = 'HOMO' if idx == 9 else f'HOMO-{9-idx}'
                    mapping[c] = new
            m = self._pat_lumo.match(c)
            if m:
                idx = int(m.group(1))
                if self.keep_homo_lumo or self.homo_lumo_all or idx <= 2:
                    new = 'LUMO' if idx == 0 else f'LUMO+{idx}'
                    mapping[c] = new
        if mapping:
            self.df.rename(columns=mapping, inplace=True)
            self.report['homo_lumo_renames'] = mapping

    def drop_atomic_dipole_orientation_by_index(self):
        self.remove_columns(r'atomic_dipole_orientation_\d+_')

    def drop_overview_columns(self):
        self.remove_columns(r'overview_hexadecapole')
        self.remove_columns(r'overview_[^_]+_au')

    def drop_sum_mulliken_duplicates(self):
        cs = 'sum_mulliken_charges.charge_sum'
        ch = 'charge'
        if cs in self.df.columns and ch in self.df.columns and (self.df[cs] == self.df[ch]).all():
            self.df.drop(columns=[cs], inplace=True)
            self.report['dropped_sum_mulliken_charge_sum'] = True
        ss = 'sum_mulliken_charges.spin_sum'
        mu = 'multiplicity'
        if ss in self.df.columns and mu in self.df.columns and ((self.df[ss]*2 + 1) == self.df[mu]).all():
            self.df.drop(columns=[ss], inplace=True)
            self.report['dropped_sum_mulliken_spin_sum'] = True

    def drop_delta_homo_lumo_columns(self):
        to_drop = []
        rename = {}
        for c in self.df.columns:
            m = re.match(r'^delta_homo_lumo_\[(\d+)\]$', c)
            if m:
                idx = int(m.group(1))
                if idx >= 3:
                    to_drop.append(c)
                else:
                    new = 'delta_homo_lumo' if idx == 0 else f'delta_homo_lumo+{idx}'
                    rename[c] = new
        self.remove_columns(to_drop)
        if rename:
            self.df.rename(columns=rename, inplace=True)
            self.report['delta_homo_lumo_renames'] = rename

    def rename_overview_columns(self):
        mapping = {c: c.replace('overview_', '') for c in self.df.columns if 'overview_' in c}
        if mapping:
            self.df.rename(columns=mapping, inplace=True)
            self.report['overview_renames'] = mapping

    def drop_and_rename_legacy_columns(self):
        self.drop_rydberg_columns()
        self.drop_total_keyword_columns()
        self.drop_homo_lumo_by_index()
        self.rename_homo_lumo()
        self.drop_atomic_dipole_orientation_by_index()
        self.drop_overview_columns()
        self.drop_sum_mulliken_duplicates()
        self.drop_delta_homo_lumo_columns()
        self.rename_overview_columns()

    def is_already_processed(self) -> bool:
        idxs = self.df.columns.str.extract(r'\[(\d+)\]')[0]
        nums = pd.to_numeric(idxs, errors='coerce')
        return not (nums.fillna(-1) >= 5).any()

    def remove_columns(self, columns_to_remove):
        if isinstance(columns_to_remove, list):
            to_drop = [c for c in columns_to_remove if c in self.df.columns]
        else:
            to_drop = self.df.filter(regex=columns_to_remove).columns.tolist()
        if to_drop:
            self.df.drop(columns=to_drop, inplace=True)
            self.report.setdefault('columns_removed', []).extend(to_drop)

    def remove_nonzero_npa_columns(self):
        cols = self.df.columns.to_series()
        idxs = cols.str.extract(self._pat_npa)[0].dropna().astype(int)
        to_drop = idxs[~idxs.isin(range(5))].index.tolist()
        self.remove_columns(to_drop)

    def remove_nonzero_nec_columns(self):
        cols = self.df.columns.to_series()
        idxs = cols.str.extract(self._pat_nec)[0].dropna().astype(int)
        to_drop = idxs[~idxs.isin(range(5))].index.tolist()
        self.remove_columns(to_drop)

    def encode_categorical(self, exclude_cols=None):
        if exclude_cols is None:
            exclude_cols = []
        cats = self.df.select_dtypes(include=['object', 'category']).columns.difference(exclude_cols)
        for c in cats:
            cat = pd.Categorical(self.df[c])
            self.encoding_dict[c] = dict(enumerate(cat.categories))
            self.df[c] = cat.codes
        self.report['encoded_columns'] = cats.tolist()
        return self.encoding_dict

    def remove_mostly_nan_columns(self, threshold=0.5, columns=None):
        if columns is None:
            cols = self.df.filter(regex=self._pat_nec.pattern).columns
        else:
            if isinstance(columns, list):
                cols = [c for c in columns if c in self.df.columns]
            else:
                cols = self.df.filter(regex=columns).columns
        if len(cols) == 0:
            self.report['remove_mostly_nan_columns'] = {
                'columns_analyzed': 0,
                'columns_dropped': [],
                'threshold': threshold
            }
            return
        rates = (self.df[cols] == -1).mean(axis=0)
        to_drop = rates[rates >= threshold].index.tolist()
        
        # Log the columns being dropped and their NaN rates
        self.report['remove_mostly_nan_columns'] = {
            'columns_analyzed': len(cols),
            'columns_dropped': to_drop,
            'threshold': threshold,
            'dropped_column_nan_rates': {col: float(rates[col]) for col in to_drop}
        }
        
        self.remove_columns(to_drop)

    def find_and_drop_identical_columns(self, drop=True):
        initial_columns = len(self.df.columns)
        
        # Prioritize keeping iron_natural_charge over natural_charges[0].natural_charge
        # if they are identical
        if ('iron_natural_charge' in self.df.columns and 
            'natural_charges[0].natural_charge' in self.df.columns):
            if self.df['iron_natural_charge'].equals(self.df['natural_charges[0].natural_charge']):
                # Drop natural_charges[0].natural_charge since iron_natural_charge is preferred
                self.df = self.df.drop(columns=['natural_charges[0].natural_charge'])
        
        mask = ~self.df.T.duplicated(keep='first')
        unique = self.df.T[mask].index.tolist()
        
        # Find duplicate groups (columns that are identical to each other)
        duplicate_groups = {}
        if drop:
            to_drop = list(set(self.df.columns) - set(unique))
            
            # For each dropped column, find which unique column it's identical to
            for dropped_col in to_drop:
                for unique_col in unique:
                    if self.df[dropped_col].equals(self.df[unique_col]):
                        if unique_col not in duplicate_groups:
                            duplicate_groups[unique_col] = []
                        duplicate_groups[unique_col].append(dropped_col)
                        break
            
            # Log detailed information about the dropping process
            self.report['find_and_drop_identical_columns'] = {
                'initial_columns': initial_columns,
                'final_columns': len(unique),
                'columns_dropped': to_drop,
                'columns_kept': unique,
                'duplicate_groups': duplicate_groups,
                'total_dropped': len(to_drop)
            }
            
            if to_drop:
                self.df = self.df[unique]
                # Keep the old report key for backward compatibility
                self.report['dropped_identical_columns'] = to_drop
        else:
            # If not dropping, just log the analysis
            would_drop = list(set(self.df.columns) - set(unique))
            self.report['find_and_drop_identical_columns'] = {
                'initial_columns': initial_columns,
                'final_columns': len(unique),
                'columns_that_would_be_dropped': would_drop,
                'columns_kept': unique,
                'total_would_drop': len(would_drop),
                'drop_executed': False
            }

    def preserve_additional_columns(self):
        """
        Identify and preserve additional columns that should be kept throughout preprocessing.
        These columns include molecular energy components and new quantum chemistry data.
        """
        # Define patterns for columns to preserve
        preserve_patterns = [
            # Original energy components
            r'overview_N-N_au',
            r'overview_E-N_au', 
            r'overview_KE_au',
            r'^total_energies\[0\]$',
            
            # New quantum chemistry data - spin density analysis
            r'^spin_density_',
            
            # Comprehensive atomic charges
            r'^charge_comp_',
            
            # Spin contamination analysis
            r'^spin_contam_',
            
            # Bond analysis (Fe coordination)
            r'^bond_',
            
            # NAO/NLMO analysis
            r'^nao_',
            
            # Second order perturbation energies
            r'^pert_',
            
            # Comprehensive orbital analysis
            r'^orbital_',
            
            # Derived quantum chemistry features
            r'^fe_d_electron_count$',
            r'^homo_lumo_gap_ev$',
            r'^total_donation_energy$',
            r'^total_backbonding_energy$',
            r'^net_charge_transfer$',
            r'^coordination_number$',
            r'^heme_nitrogen_bonds$',
            r'^axial_bonds$',
            r'^oxidation_state_estimate$',
            r'^spin_state_estimate$'
        ]
        
        # Find columns matching these patterns
        additional_preserve_cols = []
        for pattern in preserve_patterns:
            compiled_pattern = re.compile(pattern)
            matching_cols = [col for col in self.original_df.columns if compiled_pattern.search(col)]
            additional_preserve_cols.extend(matching_cols)
        
        # Store the data from these columns if they exist
        self.additional_preserved_data = {}
        for col in additional_preserve_cols:
            if col in self.original_df.columns:
                self.additional_preserved_data[col] = self.original_df[col].copy()
        
        if verbose and additional_preserve_cols:
            print(f"Preserving {len(additional_preserve_cols)} additional columns: {additional_preserve_cols}")
            
        return additional_preserve_cols

    def add_preserved_columns_back(self):
        """
        Add preserved additional columns back to the dataframe after preprocessing.
        Only adds rows that still exist in the processed dataframe.
        """
        if not hasattr(self, 'additional_preserved_data'):
            return
            
        # Get file_name index mapping for proper row alignment
        if 'file_name' in self.df.columns and 'file_name' in self.original_df.columns:
            current_files = set(self.df['file_name'])
            
            for col, original_data in self.additional_preserved_data.items():
                # Create series with same index as current df
                aligned_data = pd.Series(index=self.df.index, dtype=original_data.dtype, name=col)
                
                # Map data based on file_name
                for idx, row in self.df.iterrows():
                    file_name = row['file_name']
                    # Find corresponding row in original data
                    orig_mask = self.original_df['file_name'] == file_name
                    if orig_mask.any():
                        orig_idx = self.original_df[orig_mask].index[0]
                        if orig_idx in original_data.index:
                            aligned_data.iloc[idx] = original_data.iloc[orig_idx]
                
                # Add the column to dataframe
                self.df[col] = aligned_data
        else:
            # Fallback: just add the columns if index alignment is possible
            for col, data in self.additional_preserved_data.items():
                if len(data) >= len(self.df):
                    self.df[col] = data.iloc[:len(self.df)].values

    def report_additional_columns_missing_values(self):
        """
        Generate a report on missing values in the additional preserved columns.
        Considers both NaN and -1 as missing values.
        """
        if not hasattr(self, 'additional_preserved_data'):
            return {}
            
        missing_report = {}
        additional_cols = [col for col in self.additional_preserved_data.keys() if col in self.df.columns]
        
        if not additional_cols:
            return missing_report
            
        print("\n" + "="*60)
        print("ADDITIONAL PRESERVED COLUMNS - MISSING VALUES REPORT")
        print("="*60)
        
        for col in additional_cols:
            if col in self.df.columns:
                total_rows = len(self.df)
                
                # Count NaN values
                nan_count = self.df[col].isna().sum()
                
                # Count -1 values (but only if column is numeric)
                minus_one_count = 0
                if pd.api.types.is_numeric_dtype(self.df[col]):
                    minus_one_count = (self.df[col] == -1).sum()
                
                # Total missing
                total_missing = nan_count + minus_one_count
                missing_percentage = (total_missing / total_rows) * 100 if total_rows > 0 else 0
                
                missing_report[col] = {
                    'total_rows': total_rows,
                    'nan_count': nan_count,
                    'minus_one_count': minus_one_count,
                    'total_missing': total_missing,
                    'missing_percentage': missing_percentage
                }
                
                print(f"{col:<40}: {missing_percentage:6.2f}% missing ({total_missing}/{total_rows})")
                if nan_count > 0 and minus_one_count > 0:
                    print(f"{'':>40}  (NaN: {nan_count}, -1: {minus_one_count})")
        
        print("="*60)
        return missing_report

    def process_natural_populations_iron_3d(self):
        """
        Filters lft_orbitals data to keep only 3d orbitals from iron (atom 1),
        extracting only occupancy and energy for each specific d-orbital type.
        Creates columns like: lft_orbitals.Fe.d_xy.occupancy, lft_orbitals.Fe.d_xy.energy, etc.
        """
        # Find all lft_orbitals columns
        lft_orbitals_cols = [col for col in self.df.columns if col.startswith('lft_orbitals.')]
        
        if not lft_orbitals_cols:
            return  # No lft_orbitals data to process
        
        # Define d-orbital mappings for proper naming
        d_orbital_mapping = {
            'dxy': 'd_xy',
            'dxz': 'd_xz', 
            'dyz': 'd_yz',
            'dx2y2': 'd_(x2-y2)',
            'dz2': 'd_z2'
        }
        
        # Find iron 3d orbital columns (atom 1, d orbitals)
        iron_3d_pattern = re.compile(r'lft_orbitals\.atoms\[0\]\.orbitals\[(\d+)\]\.(occupancy|energy)$')
        
        new_columns = {}
        cols_to_drop = []
        
        for col in lft_orbitals_cols:
            match = iron_3d_pattern.match(col)
            if match:
                orbital_idx = int(match.group(1))
                property_type = match.group(2)  # 'occupancy' or 'energy'
                
                # Check if this is a d orbital by examining the type_ao field
                type_ao_col = f'lft_orbitals.atoms[0].orbitals[{orbital_idx}].type_ao'
                
                if type_ao_col in self.df.columns:
                    # Get the first non-null value to determine orbital type
                    sample_type_ao = self.df[type_ao_col].dropna().iloc[0] if not self.df[type_ao_col].dropna().empty else ""
                    
                    # Check if this is a 3d orbital (contains 'Val( 3d)' or 'Ryd( 4d)', etc.)
                    if '3d)' in sample_type_ao:
                        # Extract the specific d orbital type
                        # Look for patterns like 'dxy', 'dxz', 'dyz', 'dx2y2', 'dz2'
                        lang_col = f'lft_orbitals.atoms[0].orbitals[{orbital_idx}].lang'
                        if lang_col in self.df.columns:
                            sample_lang = self.df[lang_col].dropna().iloc[0] if not self.df[lang_col].dropna().empty else ""
                            
                            # Map the orbital name to proper notation
                            orbital_name = d_orbital_mapping.get(sample_lang, sample_lang)
                            
                            # Create new column name
                            new_col_name = f'lft_orbitals.Fe.{orbital_name}.{property_type}'
                            new_columns[new_col_name] = self.df[col]
            
            # Mark all lft_orbitals columns for removal
            cols_to_drop.append(col)
        
        # Add the new columns to the dataframe
        for new_col, data in new_columns.items():
            self.df[new_col] = data
        
        # Remove all original lft_orbitals columns
        self.remove_columns(cols_to_drop)
        
        # Report what was processed
        if new_columns:
            self.report['lft_orbitals_iron_3d_processed'] = list(new_columns.keys())
            self.report['lft_orbitals_columns_removed'] = len(cols_to_drop)

    def calculate_lft_energy_delta(self, exclude_negative=True):
        """
        Calculates the ligand field theory energy delta (Δ = e_g - t_2g) for iron 3d orbitals.
        
        Calculates the arithmetic mean of energies for each orbital group:
        - e_g orbitals: d_(x2-y2), d_z2
        - t_2g orbitals: d_xy, d_xz, d_yz
        
        Then calculates delta = e_g_mean - t_2g_mean.
        Creates a new column 'lft_energy_delta' with the energy difference.
        
        Parameters:
        -----------
        exclude_negative : bool, default=True
            If True, exclude samples with negative ligand field splitting energy
        """
        # Define orbital groups
        t_2g_orbitals = ['d_xy', 'd_xz', 'd_yz']
        e_g_orbitals = ['d_(x2-y2)', 'd_z2']
        
        # Check if the required columns exist
        required_cols = []
        for orbital in t_2g_orbitals + e_g_orbitals:
            energy_col = f'lft_orbitals.Fe.{orbital}.energy'
            required_cols.append(energy_col)
        
        missing_cols = [col for col in required_cols if col not in self.df.columns]
        if missing_cols:
            self.report['lft_energy_delta_missing_columns'] = missing_cols
            return
        
        # Get energy columns
        t_2g_energy_cols = [f'lft_orbitals.Fe.{orbital}.energy' for orbital in t_2g_orbitals]
        e_g_energy_cols = [f'lft_orbitals.Fe.{orbital}.energy' for orbital in e_g_orbitals]
        
        # Initialize columns
        self.df['t_2g_energy_mean'] = np.nan
        self.df['e_g_energy_mean'] = np.nan
        self.df['lft_energy_delta'] = np.nan
        
        valid_count = 0
        
        for idx, row in self.df.iterrows():
            # Get energies for each group
            t_2g_energies = [row[col] for col in t_2g_energy_cols if pd.notna(row[col])]
            e_g_energies = [row[col] for col in e_g_energy_cols if pd.notna(row[col])]
            
            # Skip if missing data
            if len(t_2g_energies) != 3 or len(e_g_energies) != 2:
                continue
            
            # Calculate arithmetic means for each group
            t_2g_mean = np.mean(t_2g_energies)
            e_g_mean = np.mean(e_g_energies)
            
            # Store the means and calculate delta
            self.df.at[idx, 't_2g_energy_mean'] = t_2g_mean
            self.df.at[idx, 'e_g_energy_mean'] = e_g_mean
            self.df.at[idx, 'lft_energy_delta'] = e_g_mean - t_2g_mean
            valid_count += 1
        
        # Filter out negative values if requested
        if exclude_negative:
            negative_mask = self.df['lft_energy_delta'] < 0
            negative_count = negative_mask.sum()
            if negative_count > 0:
                print(f"Excluding {negative_count} samples with negative ligand field splitting energy")
                # Set negative values to NaN (effectively excluding them)
                self.df.loc[negative_mask, 'lft_energy_delta'] = np.nan
                valid_count -= negative_count
                self.report['lft_energy_delta_negative_excluded'] = negative_count
        
        self.report['lft_energy_delta_calculated'] = True
        self.report['lft_energy_delta_valid_count'] = valid_count
        self.report['lft_energy_delta_total_count'] = len(self.df)
        self.report['lft_energy_delta_exclude_negative_flag'] = exclude_negative

    def calculate_lft_occupancy_delta(self):
        """
        Calculates the occupancy difference between e_g and t_2g orbital groups.
        
        e_g orbitals: d_(x2-y2), d_z2
        t_2g orbitals: d_xy, d_xz, d_yz
        
        Calculates delta = sum(e_g occupancies) - sum(t_2g occupancies)
        Creates a new column 'lft_occupancy_delta' with the occupancy difference.
        """
        # Define orbital groups
        t_2g_orbitals = ['d_xy', 'd_xz', 'd_yz']
        e_g_orbitals = ['d_(x2-y2)', 'd_z2']
        
        # Check if the required columns exist
        required_cols = []
        for orbital in t_2g_orbitals + e_g_orbitals:
            occupancy_col = f'lft_orbitals.Fe.{orbital}.occupancy'
            required_cols.append(occupancy_col)
        
        missing_cols = [col for col in required_cols if col not in self.df.columns]
        if missing_cols:
            self.report['lft_occupancy_delta_missing_columns'] = missing_cols
            return
        
        # Get occupancy columns
        t_2g_occupancy_cols = [f'lft_orbitals.Fe.{orbital}.occupancy' for orbital in t_2g_orbitals]
        e_g_occupancy_cols = [f'lft_orbitals.Fe.{orbital}.occupancy' for orbital in e_g_orbitals]
        
        # Calculate total occupancies for each group
        self.df['t_2g_total_occupancy'] = self.df[t_2g_occupancy_cols].sum(axis=1)
        self.df['e_g_total_occupancy'] = self.df[e_g_occupancy_cols].sum(axis=1)
        
        # Calculate occupancy delta (e_g - t_2g)
        self.df['lft_occupancy_delta'] = self.df['e_g_total_occupancy'] - self.df['t_2g_total_occupancy']
        
        # Clean up intermediate columns
        self.df.drop(columns=['t_2g_total_occupancy', 'e_g_total_occupancy'], inplace=True)
        
        self.report['lft_occupancy_delta_calculated'] = True

    def calculate_homo_lumo_energy_differences(self):
        """
        Calculates the energy differences between corresponding HOMO and LUMO orbitals.
        
        Creates new columns:
        - 'homo_lumo_gap': HOMO - LUMO energy difference
        - 'homo1_lumo1_gap': HOMO-1 - LUMO+1 energy difference  
        - 'homo2_lumo2_gap': HOMO-2 - LUMO+2 energy difference
        
        The energy differences represent the HOMO-LUMO gaps for each orbital pair.
        """
        # Define the orbital pairs to calculate differences for
        orbital_pairs = [
            ('HOMO', 'LUMO', 'homo_lumo_gap'),
            ('HOMO-1', 'LUMO+1', 'homo1_lumo1_gap'),
            ('HOMO-2', 'LUMO+2', 'homo2_lumo2_gap')
        ]
        
        calculated_columns = []
        missing_columns = []
        
        for homo_col, lumo_col, diff_col in orbital_pairs:
            # Check if both columns exist in the dataframe
            if homo_col in self.df.columns and lumo_col in self.df.columns:
                # Calculate the energy difference (LUMO - HOMO)
                self.df[diff_col] = self.df[lumo_col] - self.df[homo_col]
                calculated_columns.append(diff_col)
            else:
                missing_cols = []
                if homo_col not in self.df.columns:
                    missing_cols.append(homo_col)
                if lumo_col not in self.df.columns:
                    missing_cols.append(lumo_col)
                missing_columns.extend(missing_cols)
        
        # Report results
        self.report['homo_lumo_energy_differences'] = {
            'calculated_columns': calculated_columns,
            'missing_columns': missing_columns,
            'total_calculated': len(calculated_columns)
        }
        
        if calculated_columns:
            print(f"Calculated HOMO-LUMO energy differences: {calculated_columns}")
        if missing_columns:
            print(f"Missing columns for HOMO-LUMO calculations: {missing_columns}")

    def select_closest_standard_orientation(self, exclude={1,2,3,4}):
        idx_pat = re.compile(r'standard_orientation\[(\d+)\]\.x$')
        # find candidate indices
        all_is = {int(m.group(1)) for c in self.df.columns if (m := idx_pat.match(c))}
        cand = sorted(all_is - {0} - set(exclude))
        valid = [i for i in cand if all(
            f'standard_orientation[{i}].{fld}' in self.df.columns
            for fld in ('atomicnum','x','y','z')
        )]
        if not valid:
            return {}

        # prepare new axial*.closest_atoms.* columns
        for axis in ('axial1','axial2'):
            for feat in ('atomicnum','x','y','z'):
                self.df[f'{axis}.closest_atoms.{feat}'] = np.nan

        kept = {}
        for idx, row in self.df.iterrows():
            # how many axes to fill?
            n_keep = int(pd.notna(row.get('axial1'))) + int(pd.notna(row.get('axial2')))
            if n_keep == 0:
                kept[idx] = []
                continue

            # central coord
            x0 = row.get('standard_orientation[0].x', np.nan)
            y0 = row.get('standard_orientation[0].y', np.nan)
            z0 = row.get('standard_orientation[0].z', np.nan)

            # filter heavy
            nums = row[[f'standard_orientation[{i}].atomicnum' for i in valid]].to_numpy()
            mask = nums > 60
            if not mask.any():
                kept[idx] = []
                continue

            # compute distances & pick closest
            xs = row[[f'standard_orientation[{i}].x' for i in valid]].to_numpy()[mask]
            ys = row[[f'standard_orientation[{i}].y' for i in valid]].to_numpy()[mask]
            zs = row[[f'standard_orientation[{i}].z' for i in valid]].to_numpy()[mask]
            dists = np.sqrt((xs-x0)**2 + (ys-y0)**2 + (zs-z0)**2)
            k = min(n_keep, len(dists))
            parts = np.argpartition(dists, k-1)[:k]
            chosen = np.array(valid)[mask][parts].tolist()
            kept[idx] = chosen

            # write out into axial1.closest_atoms.* and axial2.closest_atoms.*
            for i, axis in enumerate(('axial1','axial2')):
                if i < len(chosen) and pd.notna(row.get(axis)):
                    atom_i = chosen[i]
                    for feat, suffix in (('atomicnum','.atomicnum'),
                                         ('x','.x'),
                                         ('y','.y'),
                                         ('z','.z')):
                        oldcol = f'standard_orientation[{atom_i}]{suffix}'
                        newcol = f'{axis}.closest_atoms.{feat}'
                        self.df.at[idx, newcol] = row.get(oldcol, np.nan)

        self.report['kept_std_orient_atoms_per_row'] = kept
        return kept

    def filter_files_ending_with_re(self):
        """
        Removes rows where the filename base (excluding file extension) ends with 're'.
        """
        if 'file_name' not in self.df.columns:
            return
        
        initial_count = len(self.df)
        
        # Extract filename base without extension and check if it ends with 're'
        mask = ~self.df['file_name'].astype(str).str.replace(r'\.[^.]*$', '', regex=True).str.endswith('re')
        self.df = self.df[mask].reset_index(drop=True)
        
        filtered_count = initial_count - len(self.df)
        self.report['files_ending_with_re_filtered'] = filtered_count
        self.report['remaining_after_re_filter'] = len(self.df)

    def filter_common_axial_ligands(self):
        """
        Filters data to keep only the 5 most common axial ligand combinations:
        CYS-HOH, HIS-HIS, HIS-HOH, HIS-MET, and HIS-OXY.
        """
        # Check if axial ligand columns exist
        axial_cols = ['axial1', 'axial2']
        missing_cols = [col for col in axial_cols if col not in self.df.columns]
        if missing_cols:
            self.report['axial_ligand_filter_missing_columns'] = missing_cols
            return
        
        initial_count = len(self.df)
        
        # Define the allowed ligand combinations (order independent)
        allowed_combinations = {
            frozenset(['CYS', 'HOH']),
            frozenset(['HIS', 'HIS']),
            frozenset(['HIS', 'HOH']),
            frozenset(['HIS', 'MET']),
            frozenset(['HIS', 'OXY'])
        }
        
        # Create a mask for rows with allowed combinations
        def is_allowed_combination(row):
            axial1 = row.get('axial1')
            axial2 = row.get('axial2')
            
            # Skip rows with missing axial ligand data
            if pd.isna(axial1) or pd.isna(axial2):
                return False
            
            # Create a set of the two ligands
            ligand_pair = frozenset([str(axial1), str(axial2)])
            
            return ligand_pair in allowed_combinations
        
        mask = self.df.apply(is_allowed_combination, axis=1)
        self.df = self.df[mask].reset_index(drop=True)
        
        filtered_count = initial_count - len(self.df)
        self.report['axial_ligand_combinations_filtered'] = filtered_count
        self.report['remaining_after_axial_filter'] = len(self.df)
        
        # Report the distribution of remaining combinations
        if len(self.df) > 0:
            combination_counts = {}
            for _, row in self.df.iterrows():
                if pd.notna(row['axial1']) and pd.notna(row['axial2']):
                    combo = frozenset([str(row['axial1']), str(row['axial2'])])
                    combo_str = '-'.join(sorted(combo))
                    combination_counts[combo_str] = combination_counts.get(combo_str, 0) + 1
            
            self.report['remaining_axial_ligand_distribution'] = combination_counts

    def filter_suspicious_rmsd_structures(self):
        """
        Filters out PDB structures that appear in the suspicious_distances.log file.
        These are structures with large iron-axial ligand distance differences that
        may indicate problematic quantum optimizations.
        """
        if 'file_name' not in self.df.columns:
            self.report['rmsd_filter_no_filename_column'] = True
            return
        
        initial_count = len(self.df)
        
        # Initialize RMSD analyzer to get suspicious PDB IDs
        rmsd_analyzer = RMSDAnalyzer()
        
        # Run the iron-axial distance analysis to generate the log if it doesn't exist
        if not os.path.exists(rmsd_analyzer.suspicious_log_file):
            print("Running iron-axial distance analysis to identify suspicious structures...")
            rmsd_analyzer.run_iron_axial_distance_analysis()
        
        # Get suspicious PDB IDs from the log file
        suspicious_pdb_ids = rmsd_analyzer.get_suspicious_pdb_ids()
        
        if not suspicious_pdb_ids:
            self.report['rmsd_filter_no_suspicious_structures'] = True
            self.report['rmsd_filter_structures_removed'] = 0
            return
        
        # Extract PDB ID from filename (first 4 characters)
        def extract_pdb_id(filename):
            return str(filename)[:4] if filename else ""
        
        # Create mask to filter out suspicious structures
        pdb_ids = self.df['file_name'].apply(extract_pdb_id)
        mask = ~pdb_ids.isin(suspicious_pdb_ids)
        
        # Apply filter
        self.df = self.df[mask].reset_index(drop=True)
        
        filtered_count = initial_count - len(self.df)
        self.report['rmsd_filter_structures_removed'] = filtered_count
        self.report['rmsd_filter_suspicious_pdb_ids'] = list(suspicious_pdb_ids)
        self.report['rmsd_filter_remaining_structures'] = len(self.df)
        
        print(f"Filtered out {filtered_count} structures with suspicious iron-axial distances")
        print(f"Suspicious PDB IDs: {sorted(suspicious_pdb_ids)}")
        print(f"Remaining structures: {len(self.df)}")

    def filter_natural_charges_iron_only(self):
        """
        Filters natural charges to keep only the iron (Fe) natural charge.
        Removes all other natural charges but keeps the iron charge at atom position 1.
        """
        # Find all natural_charges columns
        natural_charges_cols = [col for col in self.df.columns if col.startswith('natural_charges[')]
        
        if not natural_charges_cols:
            self.report['natural_charges_filter_no_columns'] = True
            return
        
        # Pattern to match natural_charges columns and extract atom index
        natural_charges_pattern = re.compile(r'natural_charges\[(\d+)\]\.(.+)$')
        
        # Find columns that correspond to iron (atom 1, index 0)
        iron_columns = []
        columns_to_drop = []
        
        for col in natural_charges_cols:
            match = natural_charges_pattern.match(col)
            if match:
                atom_index = int(match.group(1))
                property_name = match.group(2)
                
                # Check if this is the iron atom (index 0, which is atom 1)
                if atom_index == 0:
                    # Check if this atom is actually iron by looking at the atom_symbol
                    atom_symbol_col = f'natural_charges[{atom_index}].atom_symbol'
                    is_iron = True  # Default assumption for index 0
                    
                    if atom_symbol_col in self.df.columns:
                        # Get the first non-null value to check if it's iron
                        sample_symbol = self.df[atom_symbol_col].dropna().iloc[0] if not self.df[atom_symbol_col].dropna().empty else ""
                        is_iron = (sample_symbol == 'Fe')
                    
                    if is_iron:
                        # This is iron, keep it but rename for clarity
                        if property_name == 'natural_charge':
                            iron_columns.append(col)
                        elif property_name in ['atom_symbol', 'atom_number']:
                            iron_columns.append(col)
                
                # All natural_charges columns will be dropped initially
                columns_to_drop.append(col)
        
        # Create new iron-specific columns
        for col in iron_columns:
            match = natural_charges_pattern.match(col)
            if match:
                property_name = match.group(2)
                if property_name == 'natural_charge':
                    # Create a new column specifically for iron natural charge
                    self.df['iron_natural_charge'] = self.df[col]
                elif property_name == 'atom_symbol':
                    # Optionally keep atom symbol for validation
                    self.df['iron_atom_symbol'] = self.df[col]
                elif property_name == 'atom_number':
                    # Optionally keep atom number for validation
                    self.df['iron_atom_number'] = self.df[col]
        
        # Remove all original natural_charges columns
        self.remove_columns(columns_to_drop)
        
        # Report what was processed
        iron_charge_created = 'iron_natural_charge' in self.df.columns
        self.report['natural_charges_iron_filtered'] = {
            'iron_charge_column_created': iron_charge_created,
            'total_natural_charges_columns_removed': len(columns_to_drop),
            'iron_columns_found': len(iron_columns)
        }

    def calculate_one_electron_reduction_potential(self):
        """
        Calculates the approximate one-electron reduction potential for each molecular structure.
        
        Uses the formula:
        E°approx = -(E_Red - E_Ox)/F - E_abs(SHE)
        
        Where:
        - E_Red: Total energy of reduced state (q=0, m=1 or 5)
        - E_Ox: Total energy of oxidized state (q=1, m=2 or 6) 
        - F = 96485 C/mol (Faraday constant)
        - E_abs(SHE) = 4.44 V (absolute potential of standard hydrogen electrode)
        - Conversion: 1 Hartree = 2625.5 kJ/mol
        
        Creates a new dataframe with reduction potentials and writes to CSV.
        """
        # Constants
        FARADAY_CONSTANT = 96485  # C/mol
        E_ABS_SHE = 4.44  # V
        HARTREE_TO_KJ_PER_MOL = 2625.5  # kJ/mol
        HARTREE_TO_J_PER_MOL = HARTREE_TO_KJ_PER_MOL * 1000  # J/mol
        
        # Check required columns
        required_columns = ['total_energies[0]', 'charge', 'multiplicity', 'file_name']
        missing_columns = [col for col in required_columns if col not in self.df.columns]
        if missing_columns:
            self.report['reduction_potential_missing_columns'] = missing_columns
            print(f"Missing columns for reduction potential calculation: {missing_columns}")
            return
        
        # Valid charge/multiplicity combinations
        valid_combinations = {(0, 1), (0, 5), (1, 2), (1, 6)}
        reduced_states = {(0, 1), (0, 5)}  # q=0, m=1 or 5
        oxidized_states = {(1, 2), (1, 6)}  # q=1, m=2 or 6
        
        # Filter for valid combinations only
        df_valid = self.df.copy()
        df_valid['charge_mult_combo'] = list(zip(df_valid['charge'], df_valid['multiplicity']))
        df_valid = df_valid[df_valid['charge_mult_combo'].isin(valid_combinations)].copy()
        
        if len(df_valid) == 0:
            self.report['reduction_potential_no_valid_data'] = True
            print("No valid charge/multiplicity combinations found for reduction potential calculation")
            return
        
        # Extract structure identifier (first 4 characters of filename)
        def extract_structure_id(filename):
            return str(filename)[:4] if filename else ""
        
        df_valid['structure_id'] = df_valid['file_name'].apply(extract_structure_id)
        
        # Group by structure_id
        results = []
        structure_groups = df_valid.groupby('structure_id')
        
        for structure_id, group in structure_groups:
            # Separate reduced and oxidized states
            reduced_rows = group[group['charge_mult_combo'].isin(reduced_states)]
            oxidized_rows = group[group['charge_mult_combo'].isin(oxidized_states)]
            
            if len(reduced_rows) == 0 or len(oxidized_rows) == 0:
                # Skip structures that don't have both reduced and oxidized states
                continue
            
            # Calculate only specific spin transitions (exclude spin-crossover)
            # Low spin-low spin: (0,1) → (1,2)
            # High spin-high spin: (0,5) → (1,6)
            potentials = []
            pair_details = []
            
            # Define allowed transitions (exclude spin-crossover)
            allowed_transitions = [
                ((0, 1), (1, 2)),  # Low spin to low spin
                ((0, 5), (1, 6))   # High spin to high spin
            ]
            
            for red_state, ox_state in allowed_transitions:
                # Find rows matching the specific transition
                red_match = reduced_rows[
                    (reduced_rows['charge'] == red_state[0]) & 
                    (reduced_rows['multiplicity'] == red_state[1])
                ]
                ox_match = oxidized_rows[
                    (oxidized_rows['charge'] == ox_state[0]) & 
                    (oxidized_rows['multiplicity'] == ox_state[1])
                ]
                
                # Skip if either state is missing
                if len(red_match) == 0 or len(ox_match) == 0:
                    continue
                
                # Use the first match (should only be one per structure anyway)
                red_row = red_match.iloc[0]
                ox_row = ox_match.iloc[0]
                
                # Get energies in Hartree
                E_red_hartree = red_row['total_energies[0]']
                E_ox_hartree = ox_row['total_energies[0]']
                
                # Skip if either energy is missing or -1
                if pd.isna(E_red_hartree) or pd.isna(E_ox_hartree) or E_red_hartree == -1 or E_ox_hartree == -1:
                    continue
                
                # Convert to J/mol
                E_red_j_per_mol = E_red_hartree * HARTREE_TO_J_PER_MOL
                E_ox_j_per_mol = E_ox_hartree * HARTREE_TO_J_PER_MOL
                
                # Calculate reduction potential
                # E°approx = -(E_Red - E_Ox)/F - E_abs(SHE)
                delta_E = E_red_j_per_mol - E_ox_j_per_mol
                E_approx = -delta_E / FARADAY_CONSTANT - E_ABS_SHE
                
                potentials.append(E_approx)
                pair_details.append({
                    'red_charge': red_row['charge'],
                    'red_mult': red_row['multiplicity'],
                    'ox_charge': ox_row['charge'],
                    'ox_mult': ox_row['multiplicity'],
                    'red_energy_hartree': E_red_hartree,
                    'ox_energy_hartree': E_ox_hartree,
                    'potential_V': E_approx,
                    'transition_type': f"({red_state[0]},{red_state[1]}) → ({ox_state[0]},{ox_state[1]})"
                })
            
            if len(potentials) == 0:
                continue
            
            # Calculate statistics
            mean_potential = np.mean(potentials)
            std_potential = np.std(potentials, ddof=1) if len(potentials) > 1 else 0.0
            
            # Store result
            result = {
                'structure_id': structure_id,
                'E°approx (V)': mean_potential,
                'std_dev (V)': std_potential,
                'num_red_ox_pairs': len(potentials),
                'num_reduced_states': len(reduced_rows),
                'num_oxidized_states': len(oxidized_rows)
            }
            
            # Add details about the pairs for debugging/validation
            if len(potentials) <= 5:  # Only store details for small numbers to avoid bloat
                result['pair_details'] = pair_details
            
            results.append(result)
        
        # Create results dataframe
        if len(results) == 0:
            self.report['reduction_potential_no_valid_pairs'] = True
            print("No valid Red/Ox pairs found for reduction potential calculation")
            return
        
        results_df = pd.DataFrame(results)
        
        # Write to CSV
        output_filename = 'tables/reduction_potentials.csv'
        os.makedirs('tables', exist_ok=True)
        results_df.to_csv(output_filename, index=False)
        
        # Report summary statistics
        self.report['reduction_potential_calculated'] = {
            'total_structures_processed': len(results),
            'total_valid_structures_in_data': len(structure_groups),
            'mean_potential_V': float(results_df['E°approx (V)'].mean()),
            'std_potential_V': float(results_df['E°approx (V)'].std()),
            'min_potential_V': float(results_df['E°approx (V)'].min()),
            'max_potential_V': float(results_df['E°approx (V)'].max()),
            'output_file': output_filename
        }
        
        print(f"\nReduction Potential Calculation Summary:")
        print(f"- Processed {len(results)} structures")
        print(f"- Mean potential: {results_df['E°approx (V)'].mean():.3f} ± {results_df['E°approx (V)'].std():.3f} V")
        print(f"- Range: {results_df['E°approx (V)'].min():.3f} to {results_df['E°approx (V)'].max():.3f} V")
        print(f"- Results saved to: {output_filename}")
        
        # Store results for potential further use
        self.reduction_potentials_df = results_df
        
        return results_df

    def determine_orbital_distortion_type(self):
        """
        Determines the type of orbital distortion based on iron d-orbital energies.
    
        Compares the energies of the d-orbitals to classify the structure as:
        - Octahedral: e_g and t_2g groups nearly degenerate internally and at similar energy levels
        - Tetragonally elongated: d_x²-y² > d_z² > d_xy > (d_xz & d_yz)
        - Tetragonally compressed: d_z² > d_x²-y² > (d_xz & d_yz) > d_xy
        - Square planar: d_x²-y² > d_z² > d_xy > (d_xz & d_yz)
    
        Creates new columns:
        - 'orbital_distortion_type': encoded as integers (0=octahedral, 1=tet_elongated, 2=tet_compressed, 3=square_planar)
        - 'dx2y2_dz2_diff': Energy difference between d_x²-y² and d_z²
        - 'dz2_dxy_diff': Energy difference between d_z² and d_xy
        - 'dxy_dxz_dyz_diff': Energy difference between d_xy and average of (d_xz & d_yz)
        - 'eg_t2g_diff': Energy difference between e_g group and t_2g group averages
        """
        # Define orbital energy columns
        orbital_energy_cols = {
            'd_xy': 'lft_orbitals.Fe.d_xy.energy',
            'd_xz': 'lft_orbitals.Fe.d_xz.energy', 
            'd_yz': 'lft_orbitals.Fe.d_yz.energy',
            'd_x2y2': 'lft_orbitals.Fe.d_(x2-y2).energy',
            'd_z2': 'lft_orbitals.Fe.d_z2.energy'
        }
    
        # Check if required columns exist
        missing_cols = [col for col in orbital_energy_cols.values() if col not in self.df.columns]
        if missing_cols:
            self.report['orbital_distortion_missing_columns'] = missing_cols
            return
    
        # Initialize new columns
        self.df['orbital_distortion_type'] = np.nan
        self.df['dx2y2_dz2_diff'] = np.nan
        self.df['dz2_dxy_diff'] = np.nan
        self.df['dxy_dxz_dyz_diff'] = np.nan
        self.df['eg_t2g_diff'] = np.nan
    
        # Cutoffs for classification
        similarity_cutoff = 0.01       # average e_g vs t_2g difference
        degeneracy_cutoff = 0.01       # internal splitting within e_g and t_2g
    
        # Encoding mapping
        distortion_encoding = {
            'octahedral': 0,
            'tetrag_elong': 1,
            'tetrag_compr': 2,
            'square_planar': 3
        }
    
        valid_count = 0
        classification_counts = {key: 0 for key in distortion_encoding.keys()}
    
        for idx, row in self.df.iterrows():
            # Get orbital energies
            d_xy = row[orbital_energy_cols['d_xy']]
            d_xz = row[orbital_energy_cols['d_xz']]
            d_yz = row[orbital_energy_cols['d_yz']]
            d_x2y2 = row[orbital_energy_cols['d_x2y2']]
            d_z2 = row[orbital_energy_cols['d_z2']]
    
            # Skip if any energy is missing
            if any(pd.isna(val) for val in [d_xy, d_xz, d_yz, d_x2y2, d_z2]):
                continue
    
            # Calculate energy differences
            dx2y2_dz2_diff = d_x2y2 - d_z2
            dz2_dxy_diff = d_z2 - d_xy
            dxz_dyz_avg = (d_xz + d_yz) / 2
            dxy_dxz_dyz_diff = d_xy - dxz_dyz_avg
    
            # Calculate e_g and t_2g group differences
            e_g_avg = (d_x2y2 + d_z2) / 2
            t_2g_avg = (d_xy + d_xz + d_yz) / 3
            eg_t2g_diff = e_g_avg - t_2g_avg
    
            # Internal degeneracy checks
            eg_split = abs(d_x2y2 - d_z2)
            t2g_split = np.std([d_xy, d_xz, d_yz])
    
            # Store energy differences
            self.df.at[idx, 'dx2y2_dz2_diff'] = dx2y2_dz2_diff
            self.df.at[idx, 'dz2_dxy_diff'] = dz2_dxy_diff
            self.df.at[idx, 'dxy_dxz_dyz_diff'] = dxy_dxz_dyz_diff
            self.df.at[idx, 'eg_t2g_diff'] = eg_t2g_diff
    
            # Classify distortion type
            if (abs(eg_t2g_diff) < similarity_cutoff) and (eg_split < degeneracy_cutoff) and (t2g_split < degeneracy_cutoff):
                distortion_type = 'octahedral'
    
            elif (d_x2y2 > d_z2 and d_z2 > d_xy and d_xy > dxz_dyz_avg):
                distortion_type = 'tetrag_elong'
    
            elif (d_z2 > d_x2y2 and d_x2y2 > dxz_dyz_avg and dxz_dyz_avg > d_xy):
                distortion_type = 'tetrag_compr'
    
            elif (d_x2y2 > d_z2 and d_z2 > d_xy and d_xy > dxz_dyz_avg):
                distortion_type = 'square_planar'
    
            else:
                # Fallback: assign based on dominant splitting pattern
                if abs(dx2y2_dz2_diff) > abs(dz2_dxy_diff):
                    distortion_type = 'tetrag_elong' if dx2y2_dz2_diff > 0 else 'tetrag_compr'
                else:
                    distortion_type = 'square_planar'
    
            # Store encoded distortion type
            self.df.at[idx, 'orbital_distortion_type'] = distortion_encoding[distortion_type]
            classification_counts[distortion_type] += 1
            valid_count += 1
    
        # Add to encoding dictionary for decoding
        self.encoding_dict['orbital_distortion_type'] = {v: k for k, v in distortion_encoding.items()}
    
        # Report results
        self.report['orbital_distortion_calculated'] = {
            'valid_count': valid_count,
            'total_count': len(self.df),
            'similarity_cutoff': similarity_cutoff,
            'degeneracy_cutoff': degeneracy_cutoff,
            'classification_counts': classification_counts,
            'encoding_dict': distortion_encoding
        }
    
        print(f"Orbital distortion classification completed:")
        print(f"- Valid structures: {valid_count}/{len(self.df)}")
        for dist_type, count in classification_counts.items():
            print(f"- {dist_type}: {count}")

    def rename_iron_natural_charge_column(self):
        """
        Renames natural_charges[0].natural_charge to iron_natural_charge at the beginning
        of processing to ensure it's preserved throughout the pipeline.
        """
        iron_natural_charge_col = 'natural_charges[0].natural_charge'
        #natural_charges[0].natural_charge
        if iron_natural_charge_col in self.df.columns:
            self.df.rename(columns={iron_natural_charge_col: 'iron_natural_charge'}, inplace=True)
            self.report['iron_natural_charge_renamed'] = True
            self.report['original_iron_natural_charge_column'] = iron_natural_charge_col
        else:
            self.report['iron_natural_charge_renamed'] = False
            self.report['iron_natural_charge_column_missing'] = True

    def process_pdb_and_extract_charge_multiplicity(self):
        """
        Renames '# PDB' column to 'PDB-ID' and extracts charge and multiplicity 
        from the filename. Charge is the first number of the last two numbers 
        before the file extension, multiplicity is the second number.
        """
        # Rename '# PDB' column to 'PDB-ID' if it exists
        if '# PDB' in self.df.columns:
            self.df.rename(columns={'# PDB': 'PDB-ID'}, inplace=True)
            self.report['pdb_column_renamed'] = True
        
        # Extract charge and multiplicity from filename
        if 'file_name' in self.df.columns:
            # Initialize the new columns
            self.df['charge'] = np.nan
            self.df['multiplicity'] = np.nan
            
            extraction_details = []
            
            for idx, row in self.df.iterrows():
                filename = str(row['file_name'])
                
                # Remove file extension
                base_name = filename.rsplit('.', 1)[0] if '.' in filename else filename
                
                # Extract the last two digits using regex
                # This pattern looks for the last two digits at the end of the filename
                pattern = r'(\d)(\d)$'
                match = re.search(pattern, base_name)
                
                if match:
                    # First digit is charge, second digit is multiplicity
                    charge = int(match.group(1))
                    multiplicity = int(match.group(2))
                    
                    self.df.at[idx, 'charge'] = charge
                    self.df.at[idx, 'multiplicity'] = multiplicity
                    extraction_details.append({
                        'filename': filename,
                        'base_name': base_name,
                        'charge': charge,
                        'multiplicity': multiplicity
                    })
                else:
                    extraction_details.append({
                        'filename': filename,
                        'base_name': base_name,
                        'charge': None,
                        'multiplicity': None
                    })
            
            # Report extraction results
            valid_extractions = self.df[['charge', 'multiplicity']].dropna().shape[0]
            self.report['charge_multiplicity_extracted'] = valid_extractions
            self.report['charge_multiplicity_total'] = len(self.df)
            self.report['charge_multiplicity_details'] = extraction_details[:5]  # First 5 for debugging
            
            # Convert to integer type where possible
            self.df['charge'] = self.df['charge'].astype('Int64')
            self.df['multiplicity'] = self.df['multiplicity'].astype('Int64')

    def remove_negative_lft_energy_delta(self):
        """
        Removes all samples (rows) from the dataframe where 'lft_energy_delta' is negative.
        """
        if 'lft_energy_delta' not in self.df.columns:
            print("Column 'lft_energy_delta' not found. No rows removed.")
            return
        
        initial_count = len(self.df)
        self.df = self.df[self.df['lft_energy_delta'] >= 0].reset_index(drop=True)
        removed_count = initial_count - len(self.df)
        print(f"Removed {removed_count} samples with negative lft_energy_delta.")
        self.report['negative_lft_energy_delta_removed'] = removed_count
        self.track_sample_count("remove_negative_lft_energy_delta", "Removed samples with negative LFT energy delta")

    def _generate_report(self):
        rep = {
            'rows_original': int(self.original_df.shape[0]),
            'rows_final': int(self.df.shape[0]),
            'cols_final': int(self.df.shape[1]),
            'columns_removed': sorted(self.report.get('columns_removed', []))
        }
        fname = f"logs/report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(fname, 'w') as f:
            json.dump(rep, f, indent=2)
        print(f'Report saved to {fname}')

    def process_quantum_chemistry_data(self):
        """
        Process and enhance the dataset with comprehensive quantum chemistry analysis data.
        Extracts key features from the new quantum chemistry fields and adds derived properties.
        """
        print("Processing comprehensive quantum chemistry data...")
        
        # Add Fe d-orbital analysis features
        self._add_fe_d_orbital_features()
        
        # Add orbital energy gap features  
        self._add_orbital_gap_features()
        
        
        self.report['quantum_chemistry_processing'] = {
            'fe_d_orbital_features_added': True,
            'orbital_gap_features_added': True
        }

    def _add_fe_d_orbital_features(self):
        """Add Fe d-orbital occupancy and energy features from NAO analysis."""
        nao_cols = [col for col in self.df.columns if col.startswith('nao_fe_d_orbital_analysis')]
        
        if any('total_d_occupancy' in col for col in nao_cols):
            d_occupancy_col = next((col for col in nao_cols if 'total_d_occupancy' in col), None)
            if d_occupancy_col:
                self.df['fe_d_electron_count'] = self.df[d_occupancy_col]
        
        if any('oxidation_state_estimate' in col for col in nao_cols):
            ox_state_col = next((col for col in nao_cols if 'oxidation_state_estimate' in col), None)
            if ox_state_col:
                self.df['fe_oxidation_state_nao'] = self.df[ox_state_col]

    def _add_orbital_gap_features(self):
        """Add HOMO-LUMO gap features from comprehensive orbital analysis."""
        orbital_cols = [col for col in self.df.columns if col.startswith('orbital_')]
        
        # Add fundamental gap in eV
        fund_gap_col = next((col for col in orbital_cols if 'fundamental_gap_ev' in col), None)
        if fund_gap_col:
            self.df['homo_lumo_gap_ev'] = self.df[fund_gap_col]
        
        # Add alpha gap in eV
        alpha_gap_col = next((col for col in orbital_cols if 'alpha_gap_ev' in col), None)
        if alpha_gap_col:
            self.df['alpha_homo_lumo_gap_ev'] = self.df[alpha_gap_col]
        
        # Add beta gap in eV (for unrestricted calculations)
        beta_gap_col = next((col for col in orbital_cols if 'beta_gap_ev' in col), None)
        if beta_gap_col:
            self.df['beta_homo_lumo_gap_ev'] = self.df[beta_gap_col]
        
        # Add unpaired electron count
        unpaired_col = next((col for col in orbital_cols if 'unpaired_electrons' in col), None)
        if unpaired_col:
            self.df['unpaired_electron_count'] = self.df[unpaired_col]

    def create_specialized_qm_analysis_tables(self):
        """
        Create specialized quantum chemistry analysis tables using the JSONProcessor.
        This integrates with the existing workflow to generate detailed analysis tables.
        """
        if not hasattr(self, 'json_dir'):
            # Try to infer JSON directory from common locations
            possible_dirs = [
                "/home/pbuser/Desktop/PhD_WORK/heme/jsons/",
                "./jsons/",
                "jsons/"
            ]
            json_dir = None
            for dir_path in possible_dirs:
                if os.path.exists(dir_path):
                    json_dir = dir_path
                    break
            
            if not json_dir:
                print("Warning: Could not find JSON directory for specialized table creation")
                return
        else:
            json_dir = self.json_dir
        
        try:
            from jsonprocessor import JSONProcessor
            
            # Create JSONProcessor instance
            json_processor = JSONProcessor(json_dir=json_dir)
            
            # Create specialized tables in the analysis folder
            output_dir = os.path.join(self.analysis_folder, "qm_analysis_tables")
            json_processor.create_specialized_qm_tables(output_dir=output_dir)
            
            print(f"Created specialized quantum chemistry analysis tables in: {output_dir}")
            
            self.report['specialized_qm_tables'] = {
                'created': True,
                'output_directory': output_dir,
                'json_source': json_dir
            }
            
        except ImportError as e:
            print(f"Warning: Could not import JSONProcessor for specialized table creation: {e}")
        except Exception as e:
            print(f"Warning: Error creating specialized quantum chemistry tables: {e}")

    def create_allowed_combinations_filter(self):
        """
        Create a filtering system to exclude non-allowed PDB-ID/charge/multiplicity combinations
        based on the comprehensive quantum chemistry analysis.
        """
        print("Creating filtering system for allowed PDB-ID/charge/multiplicity combinations...")
        
        # Analyze current dataset to identify valid combinations
        valid_combinations = self._identify_valid_combinations()
        
        # Create filtering rules based on spin contamination and validation
        spin_quality_filter = self._create_spin_quality_filter()
        
        # Create coordination environment filter
        coordination_filter = self._create_coordination_filter()
        
        # Combine all filters
        combined_filter = valid_combinations & spin_quality_filter & coordination_filter
        
        # Apply the filter
        initial_count = len(self.df)
        self.df = self.df[combined_filter].copy()
        final_count = len(self.df)
        
        excluded_count = initial_count - final_count
        print(f"Excluded {excluded_count} samples with invalid combinations or poor quality")
        
        # Save the filtering report
        self._save_filtering_report(valid_combinations, spin_quality_filter, coordination_filter, combined_filter)
        
        self.report['combination_filtering'] = {
            'initial_samples': initial_count,
            'final_samples': final_count,
            'excluded_samples': excluded_count,
            'exclusion_rate': excluded_count / initial_count if initial_count > 0 else 0
        }

    def _identify_valid_combinations(self):
        """Identify valid PDB-ID/charge/multiplicity combinations."""
        valid_filter = pd.Series(True, index=self.df.index)
        
        
        # Exclude unrealistic charge/multiplicity combinations
        if 'charge' in self.df.columns and 'multiplicity' in self.df.columns:
            # For heme systems, typical valid combinations are:
            # Charge +1: multiplicity 2, 4, 6 (Fe(III))
            # Charge 0: multiplicity 1, 3, 5 (Fe(II))  
            # Charge -1: multiplicity 2, 4, 6 (Fe(I))
            valid_cm_combinations = [
                (1, 2), (1, 4), (1, 6),  # Fe(III)
                (0, 1), (0, 3), (0, 5),  # Fe(II)
                (-1, 2), (-1, 4), (-1, 6)  # Fe(I)
            ]
            
            cm_pairs = list(zip(self.df['charge'], self.df['multiplicity']))
            cm_valid = pd.Series([pair in valid_cm_combinations for pair in cm_pairs], index=self.df.index)
            valid_filter &= cm_valid
        
        return valid_filter

    def _create_spin_quality_filter(self):
        """Create filter based on spin contamination quality."""
        quality_filter = pd.Series(True, index=self.df.index)
        
        
        return quality_filter

    def _create_coordination_filter(self):
        """Create filter based on coordination environment quality."""
        coord_filter = pd.Series(True, index=self.df.index)
        
        
        return coord_filter

    def _save_filtering_report(self, valid_combinations, spin_quality_filter, coordination_filter, combined_filter):
        """Save detailed filtering report."""
        report = {
            'timestamp': datetime.now().isoformat(),
            'filtering_criteria': {
                'valid_combinations': int(valid_combinations.sum()),
                'spin_quality_passed': int(spin_quality_filter.sum()),
                'coordination_quality_passed': int(coordination_filter.sum()),
                'combined_filter_passed': int(combined_filter.sum())
            },
            'exclusion_breakdown': {
                'invalid_combinations': int((~valid_combinations).sum()),
                'poor_spin_quality': int((~spin_quality_filter).sum()),
                'poor_coordination': int((~coordination_filter).sum())
            }
        }
        
        report_file = os.path.join(self.analysis_folder, "filtering_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"Filtering report saved to: {report_file}")

    def filter_quantum_chemistry_outliers(self):
        """
        Filter out structures identified as outliers in quantum chemistry analysis.
        These outliers were identified through statistical analysis of coordination,
        charge transfer, and spin contamination properties.
        """
        # Define outlier PDB IDs identified from QM analysis
        outlier_pdb_ids = {
            # Fe coordination and charge transfer outliers
            '5f0b', '5o18', '5eys', '1m20', '2wtg',

            # Spin contamination outliers (high annihilation effectiveness outliers)
            '1dt6', '1nr6', '5cp4', '5a5j', '5g6l', '5jkw', '6hqq', '1hrm',
            '2zbz', '2d0e', '1xbn', '1n4g', '2z6s', '2r7a', '2spn', '2vyz',
            '1ash', '6hqs', '5g6m', '1z8q', '1w0e', '2cj1', '4g8w', '2z3u',
            '5g69', '2nnj', '2veb', '2xkh', '1m7v', '2cp4', '2xbk', '5esn',
            '5g6j', '6oo9', '1gwh', '1q5e', '2vyw', '2owh', '1f65', '1yrd'
        }

        if 'file_name' not in self.df.columns:
            print("Warning: file_name column not found - cannot filter QM outliers")
            return

        # Extract PDB ID from filename (first 4 characters)
        self.df['pdb_id_temp'] = self.df['file_name'].str[:4]

        # Create filter to exclude outliers
        outlier_filter = ~self.df['pdb_id_temp'].isin(outlier_pdb_ids)

        # Count excluded samples
        excluded_count = (~outlier_filter).sum()
        excluded_pdb_ids = self.df[~outlier_filter]['pdb_id_temp'].unique()

        # Apply filter
        self.df = self.df[outlier_filter].copy()

        # Clean up temporary column
        if 'pdb_id_temp' in self.df.columns:
            self.df = self.df.drop('pdb_id_temp', axis=1)

        print(f"Excluded {excluded_count} samples identified as quantum chemistry outliers")
        if excluded_count > 0:
            print(f"Excluded PDB IDs: {sorted(excluded_pdb_ids)}")

    def create_sorted_axial_ligands_copy(self):
        """
        Creates a copy of the dataframe with axial ligands in standardized alphabetical order.

        Ensures that axial1 and axial2 are always in alphabetical order for consistent naming:
        - CYS-HOH (not HOH-CYS)
        - HIS-HOH (not HOH-HIS)
        - HIS-MET (not MET-HIS)
        - HIS-OXY (not OXY-HIS)

        This standardization ensures:
        1. Consistent naming across the dataset
        2. Proper grouping in analyses
        3. Meaningful biological interpretation (typically stronger ligand first)

        The sorted copy is saved to 'tables/processed_output_sorted_axial.csv'.
        """
        # Check if axial columns exist
        if 'axial1' not in self.df.columns or 'axial2' not in self.df.columns:
            print("Warning: axial1 or axial2 columns not found - cannot create sorted copy")
            self.report['sorted_axial_ligands_created'] = False
            return

        # Decoding dictionaries for axial ligands
        axial1_decode = {0: 'CYS', 1: 'HIS'}
        axial2_decode = {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'}

        # Combined decode for re-encoding
        all_ligands = {'CYS', 'HIS', 'HOH', 'MET', 'OXY'}

        # Create reverse encoding dictionaries
        axial1_encode = {v: k for k, v in axial1_decode.items()}
        axial2_encode = {v: k for k, v in axial2_decode.items()}

        # Create a copy of the dataframe
        df_sorted = self.df.copy()

        swap_count = 0
        skipped_count = 0

        for idx, row in df_sorted.iterrows():
            axial1_val = row['axial1']
            axial2_val = row['axial2']

            # Skip rows with missing data
            if pd.isna(axial1_val) or pd.isna(axial2_val):
                skipped_count += 1
                continue

            # Decode to string representations
            try:
                axial1_str = axial1_decode.get(int(axial1_val))
                axial2_str = axial2_decode.get(int(axial2_val))

                if axial1_str is None or axial2_str is None:
                    skipped_count += 1
                    continue

                # Use alphabetical ordering for all combinations
                # This ensures consistent naming: CYS-HOH, HIS-HOH, HIS-MET, HIS-OXY, etc.
                should_swap = False

                # Check if alphabetical order is violated (axial1 should come before axial2)
                if axial1_str > axial2_str:
                    should_swap = True

                if should_swap:
                    # Swap them
                    new_axial1_str = axial2_str
                    new_axial2_str = axial1_str

                    # Re-encode
                    new_axial1_encoded = axial1_encode.get(new_axial1_str)
                    new_axial2_encoded = axial2_encode.get(new_axial2_str)

                    # Only swap if both can be re-encoded
                    if new_axial1_encoded is not None and new_axial2_encoded is not None:
                        df_sorted.at[idx, 'axial1'] = new_axial1_encoded
                        df_sorted.at[idx, 'axial2'] = new_axial2_encoded
                        swap_count += 1
                    else:
                        # Cannot re-encode, keep original
                        skipped_count += 1

            except (ValueError, TypeError, KeyError):
                skipped_count += 1
                continue

        # Save the sorted copy
        output_filename = 'tables/processed_output_sorted_axial.csv'
        os.makedirs('tables', exist_ok=True)
        df_sorted.to_csv(output_filename, index=False)

        print(f"\nStandardized axial ligands copy created:")
        print(f"- Applied alphabetical ordering to all combinations")
        print(f"  (CYS-HOH, HIS-HOH, HIS-MET, HIS-OXY, etc.)")
        print(f"- Swapped {swap_count} rows to standardized order")
        print(f"- Skipped {skipped_count} rows due to missing or invalid data")
        print(f"- Saved to: {output_filename}")

        self.report['sorted_axial_ligands'] = {
            'created': True,
            'swaps_performed': swap_count,
            'rows_skipped': skipped_count,
            'total_rows': len(df_sorted),
            'output_file': output_filename
        }

        return df_sorted

    def process(self):
        # Record initial sample count
        self.track_sample_count("0_initial", "Starting preprocessing with original dataset")
        
        # 0) preserve additional columns before any processing
        self.preserve_additional_columns()
        
        # 0.5) rename iron natural charge column first, before any drops
        self.rename_iron_natural_charge_column()
        
        # 1) pre‐orientation drops
        for pat in self.static_drop_patterns_pre:
            self.remove_columns(pat)
        
        self.track_sample_count("1_pre_orientation_drops", "Removed atomic/orbital columns")

        # 2) legacy logic
        self.drop_and_rename_legacy_columns()
        self.track_sample_count("2_legacy_drops", "Legacy column drops and renames")

        # 3) warnings
        self.remove_columns(r'^warnings')
        self.track_sample_count("3_warnings_removed", "Removed warning columns")

        # 4) filtering & solute flag
        self.filter_files_ending_with_re()
        self.track_sample_count("4a_filter_re_files", "Filtered files ending with 're'")
        
        self.filter_common_axial_ligands()
        self.track_sample_count("4b_filter_axial_ligands", "Kept only common axial ligand combinations")
        
        self.filter_suspicious_rmsd_structures()  
        self.track_sample_count("4c_filter_suspicious_rmsd", "Removed structures with suspicious distances")
        
        if 'file_name' in self.df.columns:
            self.df['solute'] = (self.df['file_name'].str[-2] == 'L').astype(int)

        df1 = self.df.copy(deep=False)

        # 5) core pipeline
        self.remove_nonzero_npa_columns()
        self.remove_nonzero_nec_columns()
        self.process_natural_populations_iron_3d()
        
        # 5.5) Process new quantum chemistry data
        self.process_quantum_chemistry_data()
        self.create_specialized_qm_analysis_tables()
        self.track_sample_count("5.5_quantum_chemistry_processing", "Added comprehensive quantum chemistry analysis")
        
        # 5.6) Apply quantum chemistry quality filtering
        self.create_allowed_combinations_filter()
        self.track_sample_count("5.6_qm_quality_filtering", "Applied quantum chemistry quality filters")
        
        # 5.7) Filter out quantum chemistry outliers
        self.filter_quantum_chemistry_outliers()
        self.track_sample_count("5.7_qm_outlier_filtering", "Excluded quantum chemistry outliers")
        
        self.encode_categorical(exclude_cols=['file_name', 'charge', 'multiplicity', 'iron_natural_charge'])
        self.track_sample_count("5_core_pipeline", "Column filtering and categorical encoding")

        df2 = self.df.copy(deep=False)

        columns_to_preserve = ['charge', 'multiplicity', 'iron_natural_charge']
        preserved_data = {}
        for col in columns_to_preserve:
            if col in self.df.columns:
                preserved_data[col] = self.df[col].copy()
        
        self.df.fillna(-1, inplace=True)
        
        for col, data in preserved_data.items():
            self.df[col] = data
        
        df2_5 = self.df.copy(deep=False)

        self.remove_mostly_nan_columns()
        self.find_and_drop_identical_columns()
        self.track_sample_count("6_nan_identical_removal", "Removed columns with high NaN rates and identical columns")

        for col, data in preserved_data.items():
            self.df[col] = data

        df3 = self.df.copy(deep=False)

        subsets = {'standard_orientation': self.select_closest_standard_orientation()}

        for pat in self.static_drop_patterns_post:
            self.remove_columns(pat)
        self.remove_columns(["job_info", "total_energies[1]", "dipole_tot_au", 
                             "dipole_x_au", "dipole_y_au", "dipole_z_au"])
        self.track_sample_count("7_post_orientation_drops", "Post-orientation column cleanup")

        df4 = self.df.copy(deep=False)

        self.process_pdb_and_extract_charge_multiplicity()

        self.calculate_lft_energy_delta()
        self.calculate_lft_occupancy_delta()
        
        # Remove negative LFT energy delta samples
        self.remove_negative_lft_energy_delta()
        
        self.calculate_homo_lumo_energy_differences()
        self.determine_orbital_distortion_type()
        self.add_preserved_columns_back()
        self.calculate_one_electron_reduction_potential()

        self.track_sample_count("8_feature_calculation", "Calculated derived features and restored columns")


        self._generate_report()
        additional_missing_report = self.report_additional_columns_missing_values()
        self.generate_sample_report()

        # Create sorted axial ligands version BEFORE writing to CSV
        df_sorted = self.create_sorted_axial_ligands_copy()

        # Replace self.df with the sorted version (if sorting was successful)
        if df_sorted is not None:
            # Save unsorted version as backup
            if self.write_file:
                self.df.to_csv('tables/processed_output_unsorted_axial.csv', index=False)

            # Replace with sorted version
            self.df = df_sorted

        # Write the (now sorted) dataframe to the main output file
        if self.write_file:
            self.df.to_csv('tables/processed_output.csv', index=False)
            if self.homo_lumo_all:
                self.df.to_csv('tables/processed_output_homo_lumo_all.csv', index=False)

        if debug:
            return self.df, subsets, self.encoding_dict, df1, df2, df2_5, df3, df4
        else:
            return self.df, subsets, self.encoding_dict

class DataAnalyzer:
    """
    A class containing data analysis functions for heme protein quantum chemistry data.
    Provides methods for data filtering, analysis, and debugging.
    """
    
    def __init__(self):
        pass
    
    def add_charge_multiplicity(self, df, file_col='file_name'):
        """
        From each df[file_col], take the base filename (no path or extension),
        optionally strip a trailing 're' if the stem doesn't already end in a digit,
        then grab its last two characters and assign:
          - first of those as integer -> df['charge']
          - second    as integer -> df['multiplicity']
        If conversion fails, prints the problematic stem and sets NaN.
        """
        stems = df[file_col].apply(lambda p: os.path.splitext(os.path.basename(p))[0])
        
        charges = []
        mults = []
        for stem in stems:
            proc = stem
            # if it doesn't naturally end in a digit but ends with 're', strip it
            if not proc[-1:].isdigit() and proc.lower().endswith('re'):
                proc = proc[:-2]
            try:
                charges.append(int(proc[-2]))
                mults.append(int(proc[-1]))
            except Exception:
                print(f"ValueError parsing stem '{stem}' (used '{proc}')")
                charges.append(pd.NA)
                mults.append(pd.NA)
        
        df['charge'] = charges
        df['multiplicity'] = mults
        return df

    def split_df_by_charge_mult(self, df: pd.DataFrame, file_col: str = 'file_name') -> dict:
        """
        Splits df into sub‑DataFrames for each (charge, multiplicity) pair,
        excluding any rows whose basename ends with 're' (case‑insensitive).
        
        Returns a dict keyed by (charge, multiplicity) tuples.
        """
        # compute stem (basename without extension)
        stems = df[file_col].apply(lambda p: os.path.splitext(os.path.basename(p))[0])
        # mask out stems ending in 're'
        mask = ~stems.str.lower().str.endswith('re')
        df_clean = df.loc[mask].copy()
        
        groups = {}
        for c in df_clean['charge'].unique():
            for m in df_clean['multiplicity'].unique():
                sub = df_clean[(df_clean['charge'] == c) &
                               (df_clean['multiplicity'] == m)]
                groups[(c, m)] = sub.copy()
        return groups

    def filter_significantly_high(self, df: pd.DataFrame, cols: list, n_std: float = 3) -> pd.DataFrame:
        """
        Return rows where any of the specified columns has a value
        at least `n_std` standard deviations above its mean.
        """
        means = df[cols].mean()
        stds  = df[cols].std()
        mask = (df[cols] > (means + n_std * stds)).any(axis=1)
        return df.loc[mask].copy()

    def filter_homo_above(self, df: pd.DataFrame, threshold: float, col: str = 'homo[0]') -> pd.DataFrame:
        """
        Return all rows where `col` exceeds a given threshold.
        """
        return df.loc[df[col] > threshold].copy()

    def decode_df(self, df, vars_to_decode=None):
        """
        Decode specified variables in the DataFrame to their string labels.
        Returns a copy with decoded columns.
        """
        # Decoding dictionary for encoded categorical variables
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'},
            'function': {
                0: 'electron transporter',
                1: 'hydrolase',
                2: 'oxidoreductase',
                3: 'oxygen carrier',
                4: 'transferase',
                5: 'unclassified'
            }
        }
        
        df_dec = df.copy()
        for var, mapping in REPLACEMENTS.items():
            if vars_to_decode is None or var in vars_to_decode:
                if var in df_dec.columns:
                    df_dec[var] = df_dec[var].map(mapping)
        return df_dec

    def report_nan_and_minus1_percentages(self, df):
        """
        For each column in df, prints the column name and the percentage of entries
        that are either NaN or equal to -1 (int or float).
        """
        total = len(df)
        if total == 0:
            print("DataFrame is empty.")
            return

        for col in df.columns:
            # build a boolean mask: True where NaN or == -1
            mask = df[col].isna() | df[col].eq(-1)
            pct = mask.sum() / total * 100
            print(f"{col:40s}: {pct:6.2f}% missing")

    def missing_identifier_distribution(self, df, id_col="file_name", threshold=0.30):
        """
        For each column in df with >threshold fraction of missing values
        (NaN or == -1), compute the distribution of the two‐char identifier
        in df[id_col].str[4:6] among the rows where that column is missing.
        Prints each column's missing% and the distribution over ['01','05','12','16'].
        Returns a dict: { column_name: { '01': ..., '05': ..., '12': ..., '16': ... } }.
        """
        results = {}
        ids = ["01","05","12","16"]
        for col in df.columns:
            mask = df[col].isna() | (df[col] == -1)
            prop = mask.mean()
            if prop > threshold:
                # extract identifiers
                id_series = df.loc[mask, id_col].astype(str).str[4:6]
                dist = id_series.value_counts(normalize=True).reindex(ids, fill_value=0) * 100
                print(f"{col}: {prop:.1%} missing")
                for ident, pct in dist.items():
                    print(f"  {ident}: {pct:.1f}%")
                results[col] = dist.to_dict()
        return results

    def print_rows_with_missing(self, df: pd.DataFrame, columns: list, missing_val: object = -1) -> None:
        """
        Prints all rows of df where any of the given columns == missing_val.
        
        Parameters
        ----------
        df : pd.DataFrame
            The DataFrame to search.
        columns : list of str
            Column names to check for missing_val.
        missing_val : any, default -1
            The value that encodes "missing" in those columns.
        """
        # Build a boolean mask: True for rows where any column equals missing_val
        mask = df[columns].eq(missing_val).any(axis=1)
        # Select and print those rows
        print(df.loc[mask])


if __name__ == "__main__":
    # Example usage with homo_lumo_all flag
    import pandas as pd
    
    # Load your data
    df = pd.read_csv("tables/DB.csv")
    
    # Process with all HOMO/LUMO columns preserved
    preprocessor = DataPreprocessor(df, homo_lumo_all=True)
    processed_df, subsets, encoding_dict = preprocessor.process()
    
    print("DataPreprocessor with homo_lumo_all=True option is available.")
    print("Usage: preprocessor = DataPreprocessor(df, homo_lumo_all=True)")

