#!/usr/bin/env python3
"""
Plot energy differences between charge-multiplicity states vs HOMO-LUMO gaps.

This script:
1. Calculates energy differences relative to the 0-1 state (charge=0, mult=1)
2. Uses the 'homo_lumo_gap' column from processed_output.csv (in Hartree)
3. Converts HOMO-LUMO gap from Hartree to eV (× 27.211386245988)
4. Creates scatter plots showing the relationship between energy differences and gaps
5. Supports grouping by axial ligand combinations and PDB-IDs

Note: Uses ONLY homo_lumo_gap column (not homo1_lumo1_gap or homo2_lumo2_gap)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Import centralized color mappings from plots.py
try:
    from plots import get_axial_ligand_colors, get_charge_multiplicity_colors
except ImportError:
    # Fallback: define simple color getters
    def get_axial_ligand_colors(combinations, extended=False):
        colors = ['#9467bd', '#0072B2', '#009E73', '#D55E00', '#F0E442']
        return {combo: colors[i % len(colors)] for i, combo in enumerate(combinations)}

    def get_charge_multiplicity_colors(combinations, extended=False):
        color_map = {
            '0,1': '#0072B2',
            '0_1': '#0072B2',
            '0-1': '#0072B2',
            'q=0,m=1': '#0072B2',
            '0,5': '#009E73',
            '0_5': '#009E73',
            '0-5': '#009E73',
            'q=0,m=5': '#009E73',
            '1,2': '#D55E00',
            '1_2': '#D55E00',
            '1-2': '#D55E00',
            'q=1,m=2': '#D55E00',
            '1,5': '#56B4E9',
            '1_5': '#56B4E9',
            '1-5': '#56B4E9',
            'q=1,m=5': '#56B4E9',
            '1,6': '#595959',
            '1_6': '#595959',
            '1-6': '#595959',
            'q=1,m=6': '#595959',
        }
        fallback_colors = ['#0072B2', '#009E73', '#D55E00', '#595959', '#56B4E9']
        return {
            combo: color_map.get(combo, fallback_colors[i % len(fallback_colors)])
            for i, combo in enumerate(combinations)
        }

# Global settings
HARTREE_TO_EV = 27.211386245988  # CODATA 2018 value
VERBOSE = False

# Plot styling
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

class EnergyGapPlotter:
    """Plot energy differences vs HOMO-LUMO gaps for heme protein structures."""

    def __init__(self, csv_path='tables/processed_output.csv', output_dir='plots'):
        """
        Initialize the plotter.

        Parameters:
        -----------
        csv_path : str
            Path to the processed output CSV file
        output_dir : str
            Directory to save plots
        """
        self.csv_path = csv_path
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Load data
        self.df = pd.read_csv(csv_path)
        if VERBOSE:
            print(f"Loaded {len(self.df)} rows from {csv_path}")

        # Prepare data
        self._prepare_data()

    def _prepare_data(self):
        """Prepare and normalize energy data."""
        # Filter valid total energy values
        self.df_valid = self.df[
            (self.df['total_energies[0]'] != -1) &
            self.df['total_energies[0]'].notna()
        ].copy()

        if VERBOSE:
            print(f"Filtered to {len(self.df_valid)} valid energy rows")

        # Extract PDB-ID from filename
        if 'PDB_ID' not in self.df_valid.columns:
            self.df_valid['PDB_ID'] = self.df_valid['file_name'].str[:4]

        # Create charge-multiplicity label
        if 'charge_mult' not in self.df_valid.columns:
            self.df_valid['charge_mult'] = (
                self.df_valid['charge'].astype(str) + '-' +
                self.df_valid['multiplicity'].astype(str)
            )

        # Handle axial ligand combinations
        self._prepare_axial_ligands()

        # Calculate normalized energies
        self._calculate_normalized_energies()

        # Convert HOMO-LUMO gaps to eV
        self._convert_gaps_to_ev()

    def _prepare_axial_ligands(self):
        """Decode axial ligand combinations."""
        if 'axials' in self.df_valid.columns:
            self.df_valid['axial_combo'] = self.df_valid['axials']
        elif 'axial1' in self.df_valid.columns and 'axial2' in self.df_valid.columns:
            # Check if numeric encoding
            if self.df_valid['axial1'].dtype in ['int64', 'int32', 'int8', 'Int64']:
                AXIAL1_MAP = {0: 'CYS', 1: 'HIS'}
                AXIAL2_MAP = {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'}
                self.df_valid['axial1_decoded'] = self.df_valid['axial1'].map(AXIAL1_MAP)
                self.df_valid['axial2_decoded'] = self.df_valid['axial2'].map(AXIAL2_MAP)
                self.df_valid['axial_combo'] = (
                    self.df_valid['axial1_decoded'] + '-' +
                    self.df_valid['axial2_decoded']
                )
            else:
                self.df_valid['axial_combo'] = (
                    self.df_valid['axial1'].astype(str) + '-' +
                    self.df_valid['axial2'].astype(str)
                )
        else:
            print("Warning: No axial ligand columns found")
            self.df_valid['axial_combo'] = 'Unknown'

    def _calculate_normalized_energies(self):
        """Calculate energy differences relative to 0-1 state per PDB-ID."""
        self.df_valid['energy_diff_eV'] = np.nan

        processed_pdbs = []
        skipped_pdbs = []

        for pdb_id in self.df_valid['PDB_ID'].unique():
            pdb_mask = self.df_valid['PDB_ID'] == pdb_id
            pdb_data = self.df_valid.loc[pdb_mask]

            # Find charge=0, multiplicity=1 reference energy
            ref_mask = (pdb_data['charge'] == 0) & (pdb_data['multiplicity'] == 1)
            ref_energies = pdb_data.loc[ref_mask, 'total_energies[0]']

            if len(ref_energies) > 0:
                ref_energy = ref_energies.iloc[0]

                # Calculate differences in Hartree, then convert to eV
                energy_diffs = (pdb_data['total_energies[0]'] - ref_energy) * HARTREE_TO_EV
                self.df_valid.loc[pdb_mask, 'energy_diff_eV'] = energy_diffs
                processed_pdbs.append(pdb_id)

                if VERBOSE:
                    print(f"PDB {pdb_id}: reference = {ref_energy:.6f} Ha")
            else:
                skipped_pdbs.append(pdb_id)
                if VERBOSE:
                    print(f"Warning: PDB {pdb_id} has no 0-1 state, skipping")

        # Remove PDB-IDs without reference states
        if skipped_pdbs:
            self.df_valid = self.df_valid[
                ~self.df_valid['PDB_ID'].isin(skipped_pdbs)
            ].copy()
            print(f"Excluded {len(skipped_pdbs)} PDB-IDs without 0-1 reference states")

        if VERBOSE:
            print(f"Calculated normalized energies for {len(processed_pdbs)} PDB-IDs")

    def _convert_gaps_to_ev(self):
        """Convert HOMO-LUMO gap from Hartree to eV."""
        # Use only the main homo_lumo_gap column (not homo1_lumo1_gap or homo2_lumo2_gap)
        gap_column = 'homo_lumo_gap'

        if gap_column not in self.df_valid.columns:
            raise ValueError(f"Column '{gap_column}' not found in dataframe")

        # Convert to eV
        self.df_valid['homo_lumo_gap_eV'] = self.df_valid[gap_column] * HARTREE_TO_EV

        if VERBOSE:
            print(f"Converted {gap_column} to homo_lumo_gap_eV")
            print(f"  Original range: {self.df_valid[gap_column].min():.6f} - {self.df_valid[gap_column].max():.6f} Ha")
            print(f"  Converted range: {self.df_valid['homo_lumo_gap_eV'].min():.3f} - {self.df_valid['homo_lumo_gap_eV'].max():.3f} eV")

    def plot_by_charge_multiplicity(self, gap_column='homo_lumo_gap_eV',
                                    figsize=(14, 10), save=True):
        """
        Create scatter plots for each charge-multiplicity state.

        Parameters:
        -----------
        gap_column : str
            Column name for HOMO-LUMO gap (in eV)
        figsize : tuple
            Figure size (width, height)
        save : bool
            Whether to save the figure
        """
        # Get all unique charge-multiplicity states (excluding 0-1 reference)
        states = sorted([s for s in self.df_valid['charge_mult'].unique() if s != '0-1'])

        # Create subplots
        n_states = len(states)
        n_cols = 2
        n_rows = (n_states + 1) // 2

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        axes = axes.flatten() if n_states > 1 else [axes]

        # Color palette for axial ligands
        axial_combos = sorted(self.df_valid['axial_combo'].unique())
        color_map = get_axial_ligand_colors(axial_combos)

        for idx, state in enumerate(states):
            ax = axes[idx]
            state_data = self.df_valid[self.df_valid['charge_mult'] == state].copy()

            if len(state_data) == 0:
                ax.text(0.5, 0.5, f'No data for {state}',
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_xlabel('HOMO-LUMO Gap (eV)')
                ax.set_ylabel('Energy Difference (eV)')
                ax.set_title(f'State {state}')
                continue

            # Plot by axial ligand combination
            for axial in axial_combos:
                axial_data = state_data[state_data['axial_combo'] == axial]
                if len(axial_data) > 0:
                    ax.scatter(axial_data[gap_column],
                             axial_data['energy_diff_eV'],
                             c=[color_map[axial]],
                             label=axial,
                             alpha=0.6,
                             s=50,
                             edgecolors='black',
                             linewidth=0.5)

            # Add reference line at y=0
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1,
                      label='0-1 reference')

            ax.set_xlabel('HOMO-LUMO Gap (eV)')
            ax.set_ylabel('Energy Difference (eV)')
            ax.set_title(f'State q={state.split("-")[0]}, m={state.split("-")[1]} (n={len(state_data)})')
            ax.grid(True, alpha=0.3)

            # Add legend only to first subplot
            if idx == 0:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                         fontsize=8, frameon=True)

        # Remove empty subplots
        for idx in range(n_states, len(axes)):
            fig.delaxes(axes[idx])

        plt.tight_layout()

        if save:
            filename = os.path.join(self.output_dir,
                                   'energy_diff_vs_homo_lumo_gap_by_state.png')
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved plot to {filename}")

        return fig, axes

    def plot_combined_all_states(self, gap_column='homo_lumo_gap_eV',
                                figsize=(12, 8), save=True):
        """
        Create a single plot with all charge-multiplicity states.

        Parameters:
        -----------
        gap_column : str
            Column name for HOMO-LUMO gap (in eV)
        figsize : tuple
            Figure size (width, height)
        save : bool
            Whether to save the figure
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Get all unique states (excluding 0-1)
        states = sorted([s for s in self.df_valid['charge_mult'].unique() if s != '0-1'])

        # Color palette for states
        color_map = get_charge_multiplicity_colors(states)
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

        for idx, state in enumerate(states):
            state_data = self.df_valid[self.df_valid['charge_mult'] == state]
            if len(state_data) > 0:
                ax.scatter(state_data[gap_column],
                         state_data['energy_diff_eV'],
                         c=[color_map[state]],
                         label=f'q={state.split("-")[0]}, m={state.split("-")[1]} (n={len(state_data)})',
                         alpha=0.6,
                         s=60,
                         marker=markers[idx % len(markers)],
                         edgecolors='black',
                         linewidth=0.5)

        # Add reference line
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=2,
                  label='0-1 reference')

        ax.set_xlabel('HOMO-LUMO Gap (eV)', fontsize=12)
        ax.set_ylabel('Energy Difference from 0-1 State (eV)', fontsize=12)
        ax.set_title('Energy Difference vs HOMO-LUMO Gap: All States', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, frameon=True)

        plt.tight_layout()

        if save:
            filename = os.path.join(self.output_dir,
                                   'energy_diff_vs_homo_lumo_gap_all_states.png')
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved plot to {filename}")

        return fig, ax

    def plot_by_axial_ligands(self, gap_column='homo_lumo_gap_eV',
                             figsize=(16, 10), save=True):
        """
        Create separate plots for each axial ligand combination.

        Parameters:
        -----------
        gap_column : str
            Column name for HOMO-LUMO gap (in eV)
        figsize : tuple
            Figure size (width, height)
        save : bool
            Whether to save the figure
        """
        axial_combos = sorted(self.df_valid['axial_combo'].unique())

        n_axials = len(axial_combos)
        n_cols = 3
        n_rows = (n_axials + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        axes = axes.flatten() if n_axials > 1 else [axes]

        # Get all states (excluding 0-1)
        states = sorted([s for s in self.df_valid['charge_mult'].unique() if s != '0-1'])
        color_map = get_charge_multiplicity_colors(states)
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

        for idx, axial in enumerate(axial_combos):
            ax = axes[idx]
            axial_data = self.df_valid[self.df_valid['axial_combo'] == axial]

            if len(axial_data) == 0:
                ax.text(0.5, 0.5, f'No data for {axial}',
                       ha='center', va='center', transform=ax.transAxes)
                continue

            # Plot each state
            for state_idx, state in enumerate(states):
                state_data = axial_data[axial_data['charge_mult'] == state]
                if len(state_data) > 0:
                    ax.scatter(state_data[gap_column],
                             state_data['energy_diff_eV'],
                             c=[color_map[state]],
                             label=f'{state} (n={len(state_data)})',
                             alpha=0.6,
                             s=50,
                             marker=markers[state_idx % len(markers)],
                             edgecolors='black',
                             linewidth=0.5)

            # Add reference line
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1)

            ax.set_xlabel('HOMO-LUMO Gap (eV)', fontsize=10)
            ax.set_ylabel('Energy Difference (eV)', fontsize=10)
            ax.set_title(f'{axial} (n={len(axial_data)})', fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)

            if idx == 0:
                ax.legend(fontsize=7, frameon=True, loc='best')

        # Remove empty subplots
        for idx in range(n_axials, len(axes)):
            fig.delaxes(axes[idx])

        plt.tight_layout()

        if save:
            filename = os.path.join(self.output_dir,
                                   'energy_diff_vs_homo_lumo_gap_by_axial.png')
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved plot to {filename}")

        return fig, axes

    def plot_correlation_heatmap(self, save=True, figsize=(8, 6)):
        """
        Create a correlation plot showing relationship between energy differences and HOMO-LUMO gap.

        Parameters:
        -----------
        save : bool
            Whether to save the figure
        figsize : tuple
            Figure size
        """
        # Use only homo_lumo_gap_eV
        gap_col = 'homo_lumo_gap_eV'
        energy_col = 'energy_diff_eV'

        # Calculate correlations for each state
        states = sorted([s for s in self.df_valid['charge_mult'].unique() if s != '0-1'])

        correlation_data = []

        for state in states:
            state_data = self.df_valid[self.df_valid['charge_mult'] == state]
            if len(state_data) > 10:  # Need enough data points
                corr = state_data[[energy_col, gap_col]].corr().iloc[0, 1]
                n_samples = len(state_data)
                correlation_data.append({
                    'State': state,
                    'Correlation': corr,
                    'N_samples': n_samples
                })

        if len(correlation_data) == 0:
            print("Not enough data for correlation plot")
            return None, None

        corr_df = pd.DataFrame(correlation_data)

        # Create bar plot
        fig, ax = plt.subplots(figsize=figsize)

        colors = ['green' if c > 0 else 'red' for c in corr_df['Correlation']]
        bars = ax.bar(corr_df['State'], corr_df['Correlation'], color=colors, alpha=0.7, edgecolor='black')

        # Add value labels on bars
        for i, (bar, row) in enumerate(zip(bars, corr_df.itertuples())):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{row.Correlation:.3f}\n(n={row.N_samples})',
                   ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=9)

        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax.set_xlabel('Charge-Multiplicity State', fontsize=12)
        ax.set_ylabel('Pearson Correlation', fontsize=12)
        ax.set_title('Correlation: Energy Difference vs HOMO-LUMO Gap\n(homo_lumo_gap)', fontsize=13)
        ax.set_ylim([-1, 1])
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()

        if save:
            filename = os.path.join(self.output_dir,
                                   'energy_gap_correlation_heatmap.png')
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved correlation plot to {filename}")

        return fig, ax

    def generate_all_plots(self):
        """Generate all available plots."""
        print("="*80)
        print("GENERATING ENERGY DIFFERENCE VS HOMO-LUMO GAP PLOTS")
        print("="*80)
        print(f"\nUsing HOMO-LUMO gap column: 'homo_lumo_gap' (converted to eV)")
        print(f"Gap range: {self.df_valid['homo_lumo_gap'].min():.6f} - {self.df_valid['homo_lumo_gap'].max():.6f} Ha")
        print(f"Gap range: {self.df_valid['homo_lumo_gap_eV'].min():.3f} - {self.df_valid['homo_lumo_gap_eV'].max():.3f} eV")

        print("\n1. Plotting by charge-multiplicity states...")
        self.plot_by_charge_multiplicity()

        print("\n2. Plotting all states combined...")
        self.plot_combined_all_states()

        print("\n3. Plotting by axial ligand combinations...")
        self.plot_by_axial_ligands()

        print("\n4. Creating correlation heatmap...")
        self.plot_correlation_heatmap()

        print("\n" + "="*80)
        print("ALL PLOTS COMPLETED")
        print("="*80)

        # Print summary statistics
        self._print_summary_statistics()

    def _print_summary_statistics(self):
        """Print summary statistics about the data."""
        print("\nSUMMARY STATISTICS:")
        print("-" * 80)

        states = sorted([s for s in self.df_valid['charge_mult'].unique() if s != '0-1'])

        print(f"Total structures analyzed: {len(self.df_valid)}")
        print(f"Total PDB-IDs: {self.df_valid['PDB_ID'].nunique()}")
        print(f"Axial combinations: {', '.join(sorted(self.df_valid['axial_combo'].unique()))}")

        print(f"\nCharge-multiplicity states:")
        for state in states:
            count = len(self.df_valid[self.df_valid['charge_mult'] == state])
            print(f"  {state}: {count} structures")

        print(f"\nEnergy difference ranges (eV):")
        for state in states:
            state_data = self.df_valid[self.df_valid['charge_mult'] == state]
            if len(state_data) > 0:
                min_e = state_data['energy_diff_eV'].min()
                max_e = state_data['energy_diff_eV'].max()
                mean_e = state_data['energy_diff_eV'].mean()
                print(f"  {state}: {min_e:.3f} to {max_e:.3f} eV (mean: {mean_e:.3f} eV)")


if __name__ == "__main__":
    # Create plotter instance
    plotter = EnergyGapPlotter(csv_path='tables/processed_output.csv',
                              output_dir='plots')

    # Generate all plots
    plotter.generate_all_plots()

    print("\nPlots saved to 'plots/' directory")
