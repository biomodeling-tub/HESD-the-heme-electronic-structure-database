import os
import re
import shutil
import filecmp
import json
import random
import itertools
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
from datetime import datetime
import MDAnalysis as mda
from scipy.stats import gaussian_kde
import glob
from joblib import Parallel, delayed


#Global Flags
debug=False
verbose=False


class JSONProcessor:
    def __init__(self,
                 json_dir="/home/pbuser/Desktop/PhD_WORK/heme/jsons/",
                 output_csv="/home/pbuser/Desktop/PhD_WORK/heme/DB.csv"):
        self.json_dir = json_dir
        self.output_csv = output_csv

    def process_files(self):
        """
        Process all JSON files in the directory, extract data, and save to a CSV file.
        """
        data_list = []
        for filename in os.listdir(self.json_dir):
            if filename.lower().endswith('.json'):
                json_path = os.path.join(self.json_dir, filename)
                try:
                    with open(json_path, 'r') as f:
                        json_data = json.load(f)
                    flattened_data = self.flatten_json(json_data)
                    data_list.append(flattened_data)
                except (json.JSONDecodeError, OSError) as e:
                    print(f"Error reading file {filename}: {e}")

        df = pd.DataFrame(data_list)

        # If you have a "# PDB" column and want to sort by it plus file_name:
        if "# PDB" in df.columns and "file_name" in df.columns:
            df.sort_values(by=["# PDB", "file_name"], inplace=True)

        # Lift & clean up all the new columns (but keep standard_orientation[...] as-is)
        df = JSONProcessor.modify_dataframe(df)

        df.to_csv(self.output_csv, index=False)
        print(f"Consolidated data saved to {self.output_csv}")

    def flatten_json(self, json_obj, parent_key='', sep='.'):
        """
        Recursively flattens a nested JSON object into a single-level dict.
        Lists become keys like "listName[0].field", "listName[1].field", etc.
        """
        items = []
        for k, v in json_obj.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            if isinstance(v, dict):
                items.extend(self.flatten_json(v, new_key, sep=sep).items())
            elif isinstance(v, list):
                for i, item in enumerate(v):
                    # e.g. key → "standard_orientation[0].center", etc.
                    items.extend(
                        self.flatten_json({f"{new_key}[{i}]": item}, '', sep=sep).items()
                    )
            else:
                items.append((new_key, v))
        return dict(items)

    @staticmethod
    def modify_dataframe(df):
        """
        1) Drop unused columns
        2) Rename homo/lumo columns
        3) Rename dipole, polarizability, quadrupole, overview, atomic dips, Mulliken
        4) Leave standard_orientation[...] columns unchanged
        """
        # 1) Drop columns we never want in the final CSV:
        columns_to_drop = [
            'final_structure', 'termination_status', 'Unnamed: 0',
            'multiheme', 'file_path'
        ]
        
        # Add new quantum chemistry data fields that are too detailed for the main CSV
        # These will be processed separately or kept in specialized analysis tables
        detailed_qm_fields = [
            # Keep validation fields but drop detailed analysis dictionaries
            'nao_nlmo_analysis.natural_atomic_orbitals.all_orbitals',
            'second_order_perturbation_analysis.fe_interactions.fe_as_donor',
            'second_order_perturbation_analysis.fe_interactions.fe_as_acceptor',
            'bond_analysis.bond_orbitals',
            'comprehensive_orbital_analysis.orbital_energies'
        ]
        df = df.drop(columns=columns_to_drop, errors='ignore')

        rename_map = {}

        # 2) HOMO/LUMO renaming:
        for col in df.columns:
            if 'alpha_homo_lumo.homo' in col:
                rename_map[col] = col.replace('alpha_homo_lumo.homo', 'homo')
            if 'alpha_homo_lumo.lumo' in col:
                rename_map[col] = col.replace('alpha_homo_lumo.lumo', 'lumo')

        # 3) dipole, polarizability, quadrupole, overview, atomic dips, Mulliken
        for col in df.columns:
            # molecular dipole moment
            if col.startswith('dipole_moment.'):
                rename_map[col] = col.replace('dipole_moment.', 'dipole_')
            # static polarizability
            if col.startswith('dipole_polarizability.'):
                rename_map[col] = col.replace('dipole_polarizability.', 'polarizability_')
            # first quadrupole tensor
            if col.startswith('quadrupole_moment.'):
                rename_map[col] = col.replace('quadrupole_moment.', 'quadrupole_')
            # overview.* → overview_
            if col.startswith('overview.'):
                rename_map[col] = col.replace('overview.', 'overview_')
            # atomic dipole orientations → atomic_dipole_orientation_[index]_[x/y/z]
            if col.startswith('atomic_dipole_orientation'):
                # e.g. "atomic_dipole_orientation[3].x" → "atomic_dipole_orientation_3_x"
                new = col.replace('[', '_').replace(']', '').replace('.', '_')
                rename_map[col] = new
            # per-atom Mulliken blocks → similar flattening
            if col.startswith('mulliken_charges_atomic'):
                new = col.replace('[', '_').replace(']', '').replace('.', '_')
                rename_map[col] = new

        # 4) New quantum chemistry data fields renaming:
        for col in df.columns:
            # Spin density analysis
            if col.startswith('mulliken_spin_densities.'):
                rename_map[col] = col.replace('mulliken_spin_densities.', 'spin_density_')
            if col.startswith('mulliken_spin_densities_heavy_atoms'):
                new = col.replace('[', '_').replace(']', '').replace('.', '_')
                rename_map[col] = new.replace('mulliken_spin_densities_heavy_atoms', 'spin_density_heavy')
            
            # Comprehensive atomic charges
            if col.startswith('atomic_charges_comprehensive.'):
                # Keep summary fields, rename nested ones
                if '.summary.' in col or '.validation.' in col:
                    rename_map[col] = col.replace('atomic_charges_comprehensive.', 'charge_comp_')
            
            # Spin contamination analysis
            if col.startswith('spin_contamination_analysis.'):
                if any(field in col for field in ['multiplicity', 'theoretical_s_squared', 'contamination_metrics', 'validation']):
                    rename_map[col] = col.replace('spin_contamination_analysis.', 'spin_contam_')
            
            # Bond analysis (coordination focused)
            if col.startswith('bond_analysis.'):
                if any(field in col for field in ['validation', 'fe_coordination.coordination_analysis']):
                    rename_map[col] = col.replace('bond_analysis.', 'bond_')
            
            # NAO/NLMO analysis
            if col.startswith('nao_nlmo_analysis.'):
                if any(field in col for field in ['fe_d_orbital_analysis', 'validation', 'coordination_hybridization.hybridization_summary']):
                    rename_map[col] = col.replace('nao_nlmo_analysis.', 'nao_')
            
            # Second order perturbation energies
            if col.startswith('second_order_perturbation_analysis.'):
                if any(field in col for field in ['interaction_summary', 'charge_transfer_analysis', 'validation']):
                    rename_map[col] = col.replace('second_order_perturbation_analysis.', 'pert_')
            
            # Comprehensive orbital analysis
            if col.startswith('comprehensive_orbital_analysis.'):
                if any(field in col for field in ['homo_lumo_analysis', 'orbital_statistics', 'fe_orbital_contributions', 'validation']):
                    rename_map[col] = col.replace('comprehensive_orbital_analysis.', 'orbital_')

        # 5) Do NOT rename any standard_orientation[...] columns; keep them as-is.

        df = df.rename(columns=rename_map)
        return df

    def create_specialized_qm_tables(self, output_dir="tables"):
        """
        Create specialized tables for detailed quantum chemistry analysis data.
        These tables focus on specific aspects like Fe coordination, orbital analysis, etc.
        """
        import os
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Process all JSON files to extract specialized data
        fe_coordination_data = []
        orbital_analysis_data = []
        charge_transfer_data = []
        spin_contamination_data = []
        
        for filename in os.listdir(self.json_dir):
            if filename.lower().endswith('.json'):
                json_path = os.path.join(self.json_dir, filename)
                try:
                    with open(json_path, 'r') as f:
                        json_data = json.load(f)
                    
                    # Extract Fe coordination summary
                    fe_coord = self._extract_fe_coordination_summary(json_data, filename)
                    if fe_coord:
                        fe_coordination_data.append(fe_coord)
                    
                    # Extract orbital analysis summary
                    orbital = self._extract_orbital_summary(json_data, filename)
                    if orbital:
                        orbital_analysis_data.append(orbital)
                    
                    # Extract charge transfer summary
                    charge_transfer = self._extract_charge_transfer_summary(json_data, filename)
                    if charge_transfer:
                        charge_transfer_data.append(charge_transfer)
                    
                    # Extract spin contamination summary
                    spin_contam = self._extract_spin_contamination_summary(json_data, filename)
                    if spin_contam:
                        spin_contamination_data.append(spin_contam)
                
                except (json.JSONDecodeError, OSError) as e:
                    print(f"Error reading file {filename}: {e}")
        
        # Create DataFrames and save to CSV
        if fe_coordination_data:
            fe_df = pd.DataFrame(fe_coordination_data)
            fe_df.to_csv(os.path.join(output_dir, "fe_coordination_analysis.csv"), index=False)
            print(f"Created Fe coordination analysis table: {len(fe_coordination_data)} entries")
        
        if orbital_analysis_data:
            orbital_df = pd.DataFrame(orbital_analysis_data)
            orbital_df.to_csv(os.path.join(output_dir, "orbital_analysis.csv"), index=False)
            print(f"Created orbital analysis table: {len(orbital_analysis_data)} entries")
        
        if charge_transfer_data:
            ct_df = pd.DataFrame(charge_transfer_data)
            ct_df.to_csv(os.path.join(output_dir, "charge_transfer_analysis.csv"), index=False)
            print(f"Created charge transfer analysis table: {len(charge_transfer_data)} entries")
        
        if spin_contamination_data:
            sc_df = pd.DataFrame(spin_contamination_data)
            sc_df.to_csv(os.path.join(output_dir, "spin_contamination_analysis.csv"), index=False)
            print(f"Created spin contamination analysis table: {len(spin_contamination_data)} entries")

    def _extract_fe_coordination_summary(self, json_data, filename):
        """Extract Fe coordination analysis summary for specialized table."""
        summary = {"filename": filename}
        
        # Basic file info
        summary["pdb_id"] = filename[:4]
        
        # Bond analysis data
        bond_analysis = json_data.get("bond_analysis", {})
        if bond_analysis.get("validation", {}).get("fe_found"):
            coord_analysis = bond_analysis.get("fe_coordination", {}).get("coordination_analysis", {})
            summary.update({
                "fe_atom_index": bond_analysis.get("fe_coordination", {}).get("fe_atom_index"),
                "coordination_number": coord_analysis.get("coordination_number", 0),
                "heme_nitrogen_bonds": coord_analysis.get("heme_nitrogen_bonds", 0),
                "axial_bonds": coord_analysis.get("axial_bonds", 0),
                "oxidation_state_estimate": coord_analysis.get("oxidation_state_estimate"),
                "spin_state_estimate": coord_analysis.get("spin_state_estimate")
            })
        
        # NAO d-orbital analysis
        nao_analysis = json_data.get("nao_nlmo_analysis", {})
        fe_d_analysis = nao_analysis.get("fe_d_orbital_analysis", {})
        if fe_d_analysis:
            summary.update({
                "d_electron_count": fe_d_analysis.get("total_d_occupancy", 0),
                "nao_oxidation_state": fe_d_analysis.get("oxidation_state_estimate"),
                "nao_spin_state": fe_d_analysis.get("spin_state_analysis", {}).get("spin_state")
            })
        
        return summary if len(summary) > 2 else None  # Only return if we have actual data

    def _extract_orbital_summary(self, json_data, filename):
        """Extract orbital analysis summary for specialized table."""
        summary = {"filename": filename, "pdb_id": filename[:4]}
        
        orbital_analysis = json_data.get("comprehensive_orbital_analysis", {})
        if orbital_analysis.get("validation", {}).get("has_orbital_data"):
            homo_lumo = orbital_analysis.get("homo_lumo_analysis", {})
            stats = orbital_analysis.get("orbital_statistics", {})
            
            summary.update({
                "is_unrestricted": orbital_analysis.get("validation", {}).get("is_unrestricted", False),
                "alpha_homo": homo_lumo.get("alpha_homo"),
                "alpha_lumo": homo_lumo.get("alpha_lumo"),
                "alpha_gap_au": homo_lumo.get("alpha_gap"),
                "alpha_gap_ev": homo_lumo.get("alpha_gap") * 27.2114 if homo_lumo.get("alpha_gap") else None,
                "beta_homo": homo_lumo.get("beta_homo"),
                "beta_lumo": homo_lumo.get("beta_lumo"),
                "beta_gap_au": homo_lumo.get("beta_gap"),
                "beta_gap_ev": homo_lumo.get("beta_gap") * 27.2114 if homo_lumo.get("beta_gap") else None,
                "fundamental_gap_au": homo_lumo.get("fundamental_gap"),
                "fundamental_gap_ev": homo_lumo.get("fundamental_gap") * 27.2114 if homo_lumo.get("fundamental_gap") else None,
                "total_electrons": stats.get("total_electrons", 0),
                "unpaired_electrons": stats.get("unpaired_electrons", 0)
            })
        
        return summary if len(summary) > 2 else None

    def _extract_charge_transfer_summary(self, json_data, filename):
        """Extract charge transfer analysis summary for specialized table."""
        summary = {"filename": filename, "pdb_id": filename[:4]}
        
        pert_analysis = json_data.get("second_order_perturbation_analysis", {})
        if pert_analysis.get("validation", {}).get("has_perturbation_data"):
            interaction_summary = pert_analysis.get("interaction_summary", {})
            charge_transfer = pert_analysis.get("charge_transfer_analysis", {})
            
            summary.update({
                "total_interactions": interaction_summary.get("total_interactions", 0),
                "fe_donor_count": interaction_summary.get("fe_donor_count", 0),
                "fe_acceptor_count": interaction_summary.get("fe_acceptor_count", 0),
                "average_interaction_energy": interaction_summary.get("average_energy", 0),
                "ligand_to_fe_donation_energy": charge_transfer.get("total_donation_energy", 0),
                "fe_to_ligand_backbonding_energy": charge_transfer.get("total_backbonding_energy", 0),
                "net_charge_transfer": charge_transfer.get("total_donation_energy", 0) - charge_transfer.get("total_backbonding_energy", 0)
            })
        
        return summary if len(summary) > 2 else None

    def _extract_spin_contamination_summary(self, json_data, filename):
        """Extract spin contamination analysis summary for specialized table."""
        summary = {"filename": filename, "pdb_id": filename[:4]}
        
        spin_analysis = json_data.get("spin_contamination_analysis", {})
        if spin_analysis.get("multiplicity"):
            contamination_metrics = spin_analysis.get("contamination_metrics", {})
            validation = spin_analysis.get("validation", {})
            
            summary.update({
                "multiplicity": spin_analysis.get("multiplicity"),
                "theoretical_s_squared": spin_analysis.get("theoretical_s_squared"),
                "final_s_squared": spin_analysis.get("final_spin", {}).get("S_squared") if spin_analysis.get("final_spin") else None,
                "absolute_contamination": contamination_metrics.get("absolute_contamination"),
                "percent_contamination": contamination_metrics.get("percent_contamination"),
                "contamination_severity": contamination_metrics.get("severity"),
                "annihilation_before": spin_analysis.get("annihilation", {}).get("before") if spin_analysis.get("annihilation") else None,
                "annihilation_after": spin_analysis.get("annihilation", {}).get("after") if spin_analysis.get("annihilation") else None,
                "annihilation_effectiveness": contamination_metrics.get("annihilation_effectiveness"),
                "is_valid_spin_state": validation.get("is_valid_spin_state", True),
                "warnings_count": len(validation.get("warnings", []))
            })
        
        return summary if len(summary) > 2 else None

