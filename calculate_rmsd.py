#!/usr/bin/env python3
"""
Calculate RMSD between PDB and XYZ heme structures.

This script processes all subdirectories in the PDB/ folder and calculates
the RMSD between {subfolder_name}_system_protonated.pdb and {subfolder_name}_g16.xyz
files, focusing only on heavy atoms of the heme porphyrin ring.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import align
import warnings
warnings.filterwarnings('ignore')
from repo_paths import PDB_DIR, derived_table_path, resolve_table_input

# Import centralized color mappings from plots.py
try:
    from plots import (AXIAL_LIGAND_COLOR_MAP, CHARGE_MULTIPLICITY_COLOR_MAP, 
                      COLORBLIND_PALETTE, EXTENDED_COLORBLIND_PALETTE,
                      get_axial_ligand_colors, get_charge_multiplicity_colors)
except ImportError:
    # Fallback if plots.py is not available
    AXIAL_LIGAND_COLOR_MAP = {
        'CYS-HOH': '#9467bd',  # Purple
        'HIS-HIS': '#0072B2',  # Blue
        'HIS-HOH': '#009E73',  # Green
        'HIS-MET': '#D55E00',  # Red
        'HIS-OXY': '#F0E442'   # Yellow
    }
    CHARGE_MULTIPLICITY_COLOR_MAP = {
        '0,1': '#0072B2', '0_1': '#0072B2', '0-1': '#0072B2',
        '0,2': '#009E73', '0_2': '#009E73', '0-2': '#009E73',
        '0,3': '#D55E00', '0_3': '#D55E00', '0-3': '#D55E00',
        '0,5': '#009E73', '0_5': '#009E73', '0-5': '#009E73',
        '1,2': '#D55E00', '1_2': '#D55E00', '1-2': '#D55E00',
        '1,4': '#ff7f0e', '1_4': '#ff7f0e', '1-4': '#ff7f0e',
        '1,6': '#F0E442', '1_6': '#F0E442', '1-6': '#F0E442',
        '-1,2': '#2ca02c', '-1_2': '#2ca02c', '-1-2': '#2ca02c',
        '-1,4': '#d62728', '-1_4': '#d62728', '-1-4': '#d62728',
        '-1,6': '#9467bd', '-1_6': '#9467bd', '-1-6': '#9467bd'
    }
    COLORBLIND_PALETTE = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00']
    
    def get_axial_ligand_colors(combinations, extended=False):
        color_map = {}
        for i, combo in enumerate(combinations):
            if combo in AXIAL_LIGAND_COLOR_MAP:
                color_map[combo] = AXIAL_LIGAND_COLOR_MAP[combo]
            else:
                color_map[combo] = COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)]
        return color_map

verbose=False

def read_xyz_file(xyz_file):
    """
    Read XYZ file and extract atom information.
    
    Args:
        xyz_file (str): Path to XYZ file
        
    Returns:
        tuple: (elements, coordinates) where elements is list of element symbols
               and coordinates is numpy array of coordinates
    """
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    # Parse XYZ format
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    # Extract atom types and coordinates
    elements = []
    coordinates = []
    
    for i in range(2, 2 + n_atoms):
        parts = lines[i].strip().split()
        element = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        
        elements.append(element)
        coordinates.append([x, y, z])
    
    coordinates = np.array(coordinates)
    return elements, coordinates


def get_porphyrin_core_atoms(universe):
    """
    Get selection of porphyrin core atoms (Fe, N, and core carbons).
    
    Args:
        universe (MDAnalysis.Universe): Universe object
        
    Returns:
        MDAnalysis.AtomGroup: Selection of porphyrin core atoms
    """
    # Define porphyrin core atom names based on standard heme nomenclature
    # Fe: Fe center
    # NA, NB, NC, ND: four nitrogen atoms
    # C1A, C1B, C1C, C1D: alpha carbons connecting pyrrole rings  
    # C2A, C2B, C2C, C2D, C3A, C3B, C3C, C3D: beta carbons in pyrrole rings
    # C4A, C4B, C4C, C4D: gamma carbons in pyrrole rings
    # CHA, CHB, CHC, CHD: meso carbons (bridges between pyrrole rings)
    
    core_atoms = [
        "FE", "NA", "NB", "NC", "ND",
        "C1A", "C2A", "C3A", "C4A",
        "C1B", "C2B", "C3B", "C4B", 
        "C1C", "C2C", "C3C", "C4C",
        "C1D", "C2D", "C3D", "C4D",
        "CHA", "CHB", "CHC", "CHD"
    ]
    
    # Create selection string
    selection_str = "resname HEM and (" + " or ".join([f"name {atom}" for atom in core_atoms]) + ")"
    selection = universe.select_atoms(selection_str)
    return selection


def match_atoms_by_order(pdb_universe, xyz_elements, xyz_coords, verbose=False):
    """
    Match atoms between PDB and XYZ by filtering heavy atoms only and matching by order.
    
    Args:
        pdb_universe: MDAnalysis Universe for PDB
        xyz_elements: List of element symbols from XYZ
        xyz_coords: Numpy array of XYZ coordinates
        verbose: Whether to print detailed information
        
    Returns:
        tuple: (pdb_coords, xyz_coords_matched) or (None, None) if failed
    """
    # Get porphyrin core atoms from PDB
    pdb_core_atoms = get_porphyrin_core_atoms(pdb_universe)
    
    if len(pdb_core_atoms) == 0:
        if verbose:
            print("No porphyrin core atoms found in PDB")
        return None, None
    
    # Get expected element sequence from PDB
    pdb_elements = [atom.name for atom in pdb_core_atoms]
    pdb_coords = pdb_core_atoms.positions
    
    # Filter XYZ to get only heavy atoms (non-hydrogen)
    xyz_heavy_indices = []
    xyz_heavy_elements = []
    xyz_heavy_coords = []
    
    for i, element in enumerate(xyz_elements):
        if element.upper() != 'H':  # Exclude hydrogen atoms
            xyz_heavy_indices.append(i)
            xyz_heavy_elements.append(element)
            xyz_heavy_coords.append(xyz_coords[i])
    
    xyz_heavy_coords = np.array(xyz_heavy_coords)
    
    # Reduced verbosity - only print for debugging if needed
    # print(f"PDB core atoms ({len(pdb_core_atoms)}): {pdb_elements}")
    # print(f"XYZ heavy atoms ({len(xyz_heavy_elements)}): {xyz_heavy_elements[:len(pdb_elements)]}")
    
    # Check if we have enough heavy atoms in XYZ
    if len(xyz_heavy_elements) < len(pdb_core_atoms):
        if verbose:
            print(f"XYZ has fewer heavy atoms ({len(xyz_heavy_elements)}) than PDB core atoms ({len(pdb_core_atoms)})")
        return None, None
    
    # Take first N heavy atoms from XYZ that match the number of core atoms
    xyz_matched_coords = xyz_heavy_coords[:len(pdb_core_atoms)]
    xyz_matched_elements = xyz_heavy_elements[:len(pdb_core_atoms)]
    
    # Verify element correspondence (allowing for case differences)
    element_map = {'FE': 'Fe', 'N': 'N', 'C': 'C', 'O': 'O', 'S': 'S'}
    
    mismatches = 0
    for i, (pdb_atom, xyz_element) in enumerate(zip(pdb_elements, xyz_matched_elements)):
        # Extract element from PDB atom name
        if pdb_atom.startswith('FE'):
            expected_element = 'Fe'
        elif pdb_atom.startswith('N'):
            expected_element = 'N'
        elif pdb_atom.startswith('C'):
            expected_element = 'C'
        elif pdb_atom.startswith('O'):
            expected_element = 'O'
        elif pdb_atom.startswith('S'):
            expected_element = 'S'
        else:
            expected_element = pdb_atom[0]
            
        if expected_element.upper() != xyz_element.upper():
            if verbose:
                print(f"Element mismatch at position {i}: PDB {pdb_atom} (expected {expected_element}) vs XYZ {xyz_element}")
            mismatches += 1
    
    if mismatches > len(pdb_core_atoms) * 0.15:  # Allow up to 15% mismatches
        if verbose:
            print(f"Too many element mismatches ({mismatches}/{len(pdb_core_atoms)})")
        return None, None
    
    # Only print mismatches if significant
    if mismatches > 2 and verbose:
        print(f"Warning: {mismatches} element mismatches found, but proceeding with calculation")
    
    return pdb_coords, xyz_matched_coords


def calculate_rmsd_between_structures(pdb_file, xyz_file, verbose=False):
    """
    Calculate RMSD between PDB and XYZ structures using porphyrin core atoms.
    
    Args:
        pdb_file (str): Path to PDB file
        xyz_file (str): Path to XYZ file
        verbose: Whether to print detailed information
        
    Returns:
        float: RMSD value in Angstroms, or None if calculation failed
    """
    try:
        # Load PDB structure
        pdb_universe = mda.Universe(pdb_file)
        
        # Load XYZ structure  
        xyz_elements, xyz_coords = read_xyz_file(xyz_file)
        
        # Match atoms between structures
        pdb_coords, xyz_coords_matched = match_atoms_by_order(pdb_universe, xyz_elements, xyz_coords, verbose=verbose)
        
        if pdb_coords is None or xyz_coords_matched is None:
            return None
        
        # Calculate RMSD using the Kabsch algorithm (optimal superposition)
        # Center both coordinate sets
        pdb_centered = pdb_coords - np.mean(pdb_coords, axis=0)
        xyz_centered = xyz_coords_matched - np.mean(xyz_coords_matched, axis=0)
        
        # Calculate cross-covariance matrix
        H = np.dot(xyz_centered.T, pdb_centered)
        
        # SVD decomposition
        U, S, Vt = np.linalg.svd(H)
        
        # Calculate rotation matrix
        R = np.dot(Vt.T, U.T)
        
        # Ensure proper rotation (not reflection)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = np.dot(Vt.T, U.T)
        
        # Apply rotation to XYZ coordinates
        xyz_rotated = np.dot(xyz_centered, R)
        
        # Calculate RMSD
        diff = pdb_centered - xyz_rotated
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        return rmsd
        
    except Exception as e:
        if verbose:
            print(f"Error calculating RMSD for {pdb_file} and {xyz_file}: {str(e)}")
        return None


def extract_element_from_atom_name(atom_name):
    """
    Extract element symbol from atom name using common PDB naming conventions.
    
    Args:
        atom_name (str): Atom name from PDB file
        
    Returns:
        str: Element symbol
    """
    # Remove whitespace and convert to uppercase
    name = atom_name.strip().upper()
    
    # Special cases for common atom names
    if name.startswith('FE'):
        return 'Fe'
    elif name.startswith('XH') or name == 'HOH':  # Water oxygen
        return 'O'
    elif name.startswith('H') and len(name) > 1:  # Hydrogen atoms
        return 'H'
    elif name.startswith('C'):
        return 'C'
    elif name.startswith('N'):
        return 'N'
    elif name.startswith('O'):
        return 'O'
    elif name.startswith('S'):
        return 'S'
    elif name.startswith('P'):
        return 'P'
    else:
        # Try to extract first one or two characters as element
        if len(name) >= 2 and name[:2] in ['BR', 'CL', 'CA', 'MG', 'ZN', 'FE']:
            return name[:2].capitalize()
        elif len(name) >= 1:
            return name[0]
        else:
            return 'Unknown'


def get_iron_axial_distances(pdb_file, xyz_file, verbose=False):
    """
    Calculate distances between heme Fe and closest atoms from each axial ligand separately.
    Uses proper identification of axial ligands by residue number in PDB and corresponding
    atom positions in XYZ files.
    
    Args:
        pdb_file (str): Path to PDB file
        xyz_file (str): Path to XYZ file
        verbose: Whether to print detailed information
        
    Returns:
        dict: Dictionary with distance information for each axial ligand, or None if failed
    """
    try:
        # Load PDB structure
        pdb_universe = mda.Universe(pdb_file)
        
        # Find iron atom in PDB
        iron_pdb = pdb_universe.select_atoms("name FE")
        if len(iron_pdb) == 0:
            if verbose:
                print(f"No Fe atom found in PDB: {pdb_file}")
            return None
        
        iron_pos_pdb = iron_pdb.positions[0]
        iron_pdb_idx = iron_pdb[0].index  # Get Fe's index in PDB
        
        # Load XYZ structure first to understand total atom count
        _, xyz_coords = read_xyz_file(xyz_file)
        
        # Get all atoms in order from PDB
        all_pdb_atoms = pdb_universe.atoms
        
        # Check atom count consistency
        if len(all_pdb_atoms) != len(xyz_coords):
            if verbose:
                print(f"WARNING: Atom count mismatch - PDB has {len(all_pdb_atoms)} atoms, XYZ has {len(xyz_coords)} atoms")
            with open(self.suspicious_log_file, "a") as f:
                pdb_id = Path(pdb_file).parent.name
                f.write(f"\nWARNING: {pdb_id} - Atom count mismatch (PDB: {len(all_pdb_atoms)}, XYZ: {len(xyz_coords)})\n")
        
        # Find Fe position in XYZ (should be at same index as PDB)
        iron_xyz_idx = iron_pdb_idx
        if iron_xyz_idx >= len(xyz_coords):
            if verbose:
                print(f"Fe index {iron_xyz_idx} exceeds XYZ atom count {len(xyz_coords)}")
            return None
        
        iron_pos_xyz = xyz_coords[iron_xyz_idx]
        
        # Find axial ligand atoms separately for each residue in PDB
        axial_res2 = pdb_universe.select_atoms("resid 2")
        axial_res3 = pdb_universe.select_atoms("resid 3")
        
        if len(axial_res2) == 0 and len(axial_res3) == 0:
            if verbose:
                print(f"No axial ligand atoms found in residues 2 or 3 in PDB: {pdb_file}")
            return None
        
        result = {}
        coord_info = []
        
        # Process residue 2 (axial ligand 1)
        if len(axial_res2) > 0:
            resname_2 = axial_res2[0].resname
            
            # Find closest atom in residue 2
            closest_atom_pdb = None
            closest_distance_pdb = float('inf')
            
            for atom in axial_res2:
                dist_pdb = np.linalg.norm(iron_pos_pdb - atom.position)
                
                if dist_pdb < closest_distance_pdb:
                    closest_distance_pdb = dist_pdb
                    closest_atom_pdb = atom
            
            result['axial1_resname'] = resname_2
            result['axial1_pdb_distance'] = closest_distance_pdb
            
            # Find corresponding atom in XYZ file
            if closest_atom_pdb is not None:
                closest_xyz_idx = closest_atom_pdb.index
                if closest_xyz_idx < len(xyz_coords):
                    closest_xyz_pos = xyz_coords[closest_xyz_idx]
                    dist_xyz = np.linalg.norm(iron_pos_xyz - closest_xyz_pos)
                    result['axial1_xyz_distance'] = dist_xyz
                    
                    # Extract element from atom name
                    element = extract_element_from_atom_name(closest_atom_pdb.name)
                    
                    coord_info.append({
                        'axial': 1,
                        'pdb_atom_idx': closest_atom_pdb.index,
                        'xyz_atom_idx': closest_xyz_idx,
                        'iron_pdb': iron_pos_pdb.copy(),
                        'iron_xyz': iron_pos_xyz.copy(),
                        'axial_pdb': closest_atom_pdb.position.copy(),
                        'axial_xyz': closest_xyz_pos.copy(),
                        'distance_pdb': closest_distance_pdb,
                        'distance_xyz': dist_xyz,
                        'resname': resname_2,
                        'atom_name': closest_atom_pdb.name,
                        'element': element
                    })
        
        # Process residue 3 (axial ligand 2)
        if len(axial_res3) > 0:
            resname_3 = axial_res3[0].resname
            
            # Find closest atom in residue 3
            closest_atom_pdb = None
            closest_distance_pdb = float('inf')
            
            for atom in axial_res3:
                dist_pdb = np.linalg.norm(iron_pos_pdb - atom.position)
                
                if dist_pdb < closest_distance_pdb:
                    closest_distance_pdb = dist_pdb
                    closest_atom_pdb = atom
            
            result['axial2_resname'] = resname_3
            result['axial2_pdb_distance'] = closest_distance_pdb
            
            # Find corresponding atom in XYZ file
            if closest_atom_pdb is not None:
                closest_xyz_idx = closest_atom_pdb.index
                if closest_xyz_idx < len(xyz_coords):
                    closest_xyz_pos = xyz_coords[closest_xyz_idx]
                    dist_xyz = np.linalg.norm(iron_pos_xyz - closest_xyz_pos)
                    result['axial2_xyz_distance'] = dist_xyz
                    
                    # Extract element from atom name
                    element = extract_element_from_atom_name(closest_atom_pdb.name)
                    
                    coord_info.append({
                        'axial': 2,
                        'pdb_atom_idx': closest_atom_pdb.index,
                        'xyz_atom_idx': closest_xyz_idx,
                        'iron_pdb': iron_pos_pdb.copy(),
                        'iron_xyz': iron_pos_xyz.copy(),
                        'axial_pdb': closest_atom_pdb.position.copy(),
                        'axial_xyz': closest_xyz_pos.copy(),
                        'distance_pdb': closest_distance_pdb,
                        'distance_xyz': dist_xyz,
                        'resname': resname_3,
                        'atom_name': closest_atom_pdb.name,
                        'element': element
                    })
        
        # Check for large distance differences and log them
        suspicious_distances = []
        
        for coord in coord_info:
            distance_diff = abs(coord['distance_xyz'] - coord['distance_pdb'])  # Make absolute (positive)
            if distance_diff > 1.5:
                coord['distance_difference'] = distance_diff
                suspicious_distances.append(coord)
        
        # Write suspicious distances to file
        if suspicious_distances:
            pdb_id = Path(pdb_file).parent.name
            with open(self.suspicious_log_file, "a") as f:
                f.write(f"\n{'='*80}\n")
                f.write(f"LARGE DISTANCE DIFFERENCES FOUND FOR {pdb_id}\n")
                f.write(f"PDB File: {pdb_file}\n")
                f.write(f"XYZ File: {xyz_file}\n")
                f.write(f"Total atoms - PDB: {len(all_pdb_atoms)}, XYZ: {len(xyz_coords)}\n")
                f.write(f"{'='*80}\n")
                
                for coord in suspicious_distances:
                    f.write(f"\nAxial Ligand {coord['axial']} ({coord['resname']}):\n")
                    f.write(f"Closest atom: {coord['atom_name']} ({coord['element']}) at PDB index {coord['pdb_atom_idx']}\n")
                    f.write(f"PDB distance: {coord['distance_pdb']:.3f} Å\n")
                    f.write(f"XYZ distance: {coord['distance_xyz']:.3f} Å\n")
                    f.write(f"Distance difference: {coord['distance_difference']:.3f} Å (> 1.5 Å threshold)\n")
                    f.write(f"Iron PDB coords: [{coord['iron_pdb'][0]:.3f}, {coord['iron_pdb'][1]:.3f}, {coord['iron_pdb'][2]:.3f}]\n")
                    f.write(f"Iron XYZ coords: [{coord['iron_xyz'][0]:.3f}, {coord['iron_xyz'][1]:.3f}, {coord['iron_xyz'][2]:.3f}]\n")
                    f.write(f"Axial PDB coords: [{coord['axial_pdb'][0]:.3f}, {coord['axial_pdb'][1]:.3f}, {coord['axial_pdb'][2]:.3f}]\n")
                    f.write(f"Axial XYZ coords: [{coord['axial_xyz'][0]:.3f}, {coord['axial_xyz'][1]:.3f}, {coord['axial_xyz'][2]:.3f}]\n")
                
                f.write(f"\n")
        
        return result
        
    except Exception as e:
        if verbose:
            print(f"Error calculating Fe-axial distances for {pdb_file} and {xyz_file}: {str(e)}")
            import traceback
            traceback.print_exc()
        return None




class RMSDAnalyzer:
    """
    A class for comprehensive RMSD and structural analysis of heme proteins.
    Provides methods for RMSD calculation, iron-axial distance analysis, and structure validation.

    Note: When reading processed_output.csv from the preprocessor, axial ligands are
    pre-sorted alphabetically (axial1 < axial2). This ensures consistent ordering across
    all analyses and visualizations. The sorting maintains combinations like:
    - CYS-HOH (not HOH-CYS)
    - HIS-OXY (not OXY-HIS)
    - HIS-MET (not MET-HIS)
    """
    
    def __init__(self, pdb_base_dir=None, verbose=False, filter_axial_ligands=True):
        self.pdb_base_dir = Path(pdb_base_dir) if pdb_base_dir is not None else PDB_DIR
        self.suspicious_log_file = str(derived_table_path("suspicious_distances.log"))
        self.verbose = verbose
        self.filter_axial_ligands = filter_axial_ligands
        
        # Define the allowed axial ligand combinations (same as in preprocessor.py)
        self.allowed_axial_combinations = {
            frozenset(['CYS', 'HOH']),
            frozenset(['HIS', 'HIS']),
            frozenset(['HIS', 'HOH']),
            frozenset(['HIS', 'MET']),
            frozenset(['HIS', 'OXY'])
        }
        
    def normalize_axial_ligand_name(self, ligand_name):
        """
        Normalize axial ligand names by treating HSD as HIS.
        
        Args:
            ligand_name (str): Original ligand name
            
        Returns:
            str: Normalized ligand name
        """
        if ligand_name is None:
            return None
        
        normalized = str(ligand_name).upper()
        # Treat HSD (protonated histidine) as HIS
        if normalized == 'HSD':
            return 'HIS'
        return normalized
    
    def is_allowed_axial_combination(self, axial1, axial2):
        """
        Check if the given axial ligand combination is in the allowed list.
        Treats HSD as HIS for compatibility.
        
        Args:
            axial1 (str): First axial ligand residue name
            axial2 (str): Second axial ligand residue name
            
        Returns:
            bool: True if combination is allowed, False otherwise
        """
        if not self.filter_axial_ligands:
            return True
            
        if axial1 is None or axial2 is None:
            return False
            
        # Normalize ligand names (HSD -> HIS)
        norm_axial1 = self.normalize_axial_ligand_name(axial1)
        norm_axial2 = self.normalize_axial_ligand_name(axial2)
        
        ligand_pair = frozenset([norm_axial1, norm_axial2])
        return ligand_pair in self.allowed_axial_combinations
    
    def collect_filtered_axial_combinations(self):
        """
        Collect all PDB IDs and their axial ligand combinations that are not part of 
        the allowed combinations and would be filtered out from the iron-plane analysis.
        
        Returns:
            pandas.DataFrame: DataFrame with PDB_ID, Axial1, Axial2, and Normalized_Combination columns
        """
        filtered_results = []
        
        for subdir in self.pdb_base_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name
                
                pdb_file = subdir / f"{pdb_id}_system_protonated.pdb"
                xyz_file = subdir / f"{pdb_id}_g16.xyz"
                
                if pdb_file.exists() and xyz_file.exists():
                    # Get axial ligand information
                    axial_distances = get_iron_axial_distances(str(pdb_file), str(xyz_file), verbose=False)
                    if axial_distances is not None:
                        axial1 = axial_distances.get('axial1_resname', 'UNK')
                        axial2 = axial_distances.get('axial2_resname', 'UNK')
                        
                        # Check if this combination would be filtered out
                        if not self.is_allowed_axial_combination(axial1, axial2):
                            # Normalize the ligand names for reporting
                            norm_axial1 = self.normalize_axial_ligand_name(axial1)
                            norm_axial2 = self.normalize_axial_ligand_name(axial2)
                            norm_combination = '-'.join(sorted([norm_axial1, norm_axial2]))
                            
                            filtered_results.append({
                                'PDB_ID': pdb_id,
                                'Original_Axial1': axial1,
                                'Original_Axial2': axial2,
                                'Normalized_Axial1': norm_axial1,
                                'Normalized_Axial2': norm_axial2,
                                'Normalized_Combination': norm_combination,
                                'Original_Combination': '-'.join(sorted([str(axial1), str(axial2)]))
                            })
        
        return pd.DataFrame(filtered_results)
    
    def save_filtered_axial_combinations_report(self, output_file=None):
        """
        Generate and save a report of PDB IDs and their axial ligand combinations 
        that are filtered out from the iron-plane distance analysis.
        
        Args:
            output_file (str): Output file path for the filtered combinations report
            
        Returns:
            pandas.DataFrame: DataFrame with the filtered combinations
        """
        if output_file is None:
            output_file = str(derived_table_path("filtered_axial_combinations_report.csv"))

        if self.verbose:
            print("Collecting PDB IDs and axial combinations that are filtered out...")
        
        filtered_df = self.collect_filtered_axial_combinations()
        
        if not filtered_df.empty:
            # Save to CSV
            filtered_df.to_csv(output_file, index=False)
            
            if self.verbose:
                print(f"\nFiltered Axial Combinations Report:")
                print(f"Total filtered structures: {len(filtered_df)}")
                print(f"Unique original combinations: {filtered_df['Original_Combination'].nunique()}")
                print(f"Unique normalized combinations: {filtered_df['Normalized_Combination'].nunique()}")
                
                print(f"\nTop 10 filtered combinations:")
                combination_counts = filtered_df['Original_Combination'].value_counts().head(10)
                for combo, count in combination_counts.items():
                    print(f"  {combo}: {count} structures")
                
                print(f"\nReport saved to: {output_file}")
        else:
            if self.verbose:
                print("No structures were filtered out - all have allowed axial ligand combinations")
                
        return filtered_df
        
    def initialize_log(self):
        """Initialize the suspicious distances log file."""
        try:
            with open(self.suspicious_log_file, "w") as f:
                f.write("IRON-AXIAL LIGAND DISTANCE DIFFERENCES LOG\n")
                f.write("=" * 50 + "\n")
                f.write("This file logs cases where iron-axial ligand distance differences > 1.5 Å\n")
                f.write("Large differences may indicate significant structural changes during optimization.\n")
                f.write("=" * 50 + "\n")
        except Exception:
            pass
    
    def get_suspicious_pdb_ids(self):
        """
        Extract PDB IDs from the suspicious_distances.log file.
        
        Returns:
            set: Set of PDB IDs that have suspicious distance differences
        """
        suspicious_ids = set()
        
        if not os.path.exists(self.suspicious_log_file):
            return suspicious_ids
            
        try:
            with open(self.suspicious_log_file, 'r') as f:
                content = f.read()
                
            # Look for lines that contain "LARGE DISTANCE DIFFERENCES FOUND FOR"
            import re
            pattern = r"LARGE DISTANCE DIFFERENCES FOUND FOR (\w+)"
            matches = re.findall(pattern, content)
            suspicious_ids.update(matches)
            
        except Exception as e:
            if self.verbose:
                print(f"Error reading suspicious distances log: {e}")
            
        return suspicious_ids
    
    def process_all_structures(self, verbose=None):
        """
        Process all structures in PDB subdirectories and calculate RMSDs.
        
        Args:
            verbose: Override the instance verbose setting
        
        Returns:
            pandas.DataFrame: DataFrame with PDB IDs and RMSD values
        """
        if verbose is None:
            verbose = self.verbose
        results = []
        
        # Find all subdirectories that contain both required files
        for subdir in self.pdb_base_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name
                
                # Construct expected file paths
                pdb_file = subdir / f"{pdb_id}_system_protonated.pdb"
                xyz_file = subdir / f"{pdb_id}_g16.xyz"
                
                # Check if both files exist
                if pdb_file.exists() and xyz_file.exists():
                    # Check axial ligand filtering if enabled
                    if self.filter_axial_ligands:
                        axial_distances = get_iron_axial_distances(str(pdb_file), str(xyz_file), verbose=False)
                        if axial_distances is not None:
                            axial1 = axial_distances.get('axial1_resname')
                            axial2 = axial_distances.get('axial2_resname')
                            if not self.is_allowed_axial_combination(axial1, axial2):
                                if verbose:
                                    print(f"  Skipping {pdb_id} - axial combination {axial1}-{axial2} not allowed")
                                continue
                    
                    if verbose:
                        print(f"Processing {pdb_id}...")
                    
                    # Calculate RMSD
                    rmsd = calculate_rmsd_between_structures(str(pdb_file), str(xyz_file), verbose=verbose)
                    
                    if rmsd is not None:
                        results.append({
                            'PDB_ID': pdb_id,
                            'RMSD': rmsd
                        })
                        if verbose:
                            print(f"  RMSD: {rmsd:.3f} Å")
                    else:
                        if verbose:
                            print(f"  Failed to calculate RMSD")
                else:
                    if verbose:
                        if not pdb_file.exists():
                            print(f"Missing PDB file: {pdb_file}")
                        if not xyz_file.exists():
                            print(f"Missing XYZ file: {xyz_file}")
        
        # Create DataFrame
        df = pd.DataFrame(results)
        return df

    def standardize_axial_positions(self, distances, verbose=False):
        """
        Standardize axial ligand positions to ensure consistent alphabetical ordering.

        Applies alphabetical sorting to all axial ligand combinations:
        - CYS-HOH (not HOH-CYS)
        - HIS-HIS
        - HIS-HOH (not HOH-HIS)
        - HIS-MET (not MET-HIS)
        - HIS-OXY (not OXY-HIS)

        HSD/HSE/HSP are normalized to HIS for comparison.

        Args:
            distances: Dictionary with axial ligand distance information
            verbose: Whether to print swap information

        Returns:
            Dictionary with potentially swapped axial positions
        """
        if distances is None:
            return None

        axial1_resname = distances.get('axial1_resname')
        axial2_resname = distances.get('axial2_resname')

        # Normalize histidine variants to HIS for comparison
        norm_axial1 = 'HIS' if axial1_resname in ['HSD', 'HSE', 'HSP'] else axial1_resname
        norm_axial2 = 'HIS' if axial2_resname in ['HSD', 'HSE', 'HSP'] else axial2_resname

        # Check if we need to swap (alphabetical order: axial1 should come before axial2)
        if norm_axial1 > norm_axial2:
            if verbose:
                print(f"  Swapping positions: {axial1_resname} (pos1) ↔ {axial2_resname} (pos2) → {axial2_resname} (pos1), {axial1_resname} (pos2)")

            # Create swapped distances dictionary
            swapped = {
                'axial1_resname': distances['axial2_resname'],
                'axial2_resname': distances['axial1_resname'],
            }

            # Swap distance values if they exist
            if 'axial1_pdb_distance' in distances:
                swapped['axial2_pdb_distance'] = distances['axial1_pdb_distance']
            if 'axial1_xyz_distance' in distances:
                swapped['axial2_xyz_distance'] = distances['axial1_xyz_distance']
            if 'axial2_pdb_distance' in distances:
                swapped['axial1_pdb_distance'] = distances['axial2_pdb_distance']
            if 'axial2_xyz_distance' in distances:
                swapped['axial1_xyz_distance'] = distances['axial2_xyz_distance']

            return swapped
        else:
            # No swap needed - already in alphabetical order
            return distances

    def process_iron_axial_distances(self, verbose=None, pdb_ids_filter=None):
        """
        Process all structures and calculate Fe-axial ligand distances for each axial ligand separately.
        
        Args:
            verbose: Override the instance verbose setting
            pdb_ids_filter: List of PDB IDs to process. If None, processes all available.
        
        Returns:
            pandas.DataFrame: DataFrame with PDB IDs, distances, and differences for each axial ligand
        """
        if verbose is None:
            verbose = self.verbose
        results = []
        
        for subdir in self.pdb_base_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name
                
                # Skip if filtering is enabled and this PDB ID is not in the filter list
                if pdb_ids_filter is not None and pdb_id not in pdb_ids_filter:
                    continue
                
                pdb_file = subdir / f"{pdb_id}_system_protonated.pdb"
                xyz_file = subdir / f"{pdb_id}_g16.xyz"
                
                if pdb_file.exists() and xyz_file.exists():
                    if verbose:
                        print(f"Processing Fe-axial distances for {pdb_id}...")
                    
                    distances = get_iron_axial_distances(str(pdb_file), str(xyz_file), verbose=verbose)

                    if distances is not None:
                        # Standardize axial positions (HOH should be in position 2 for HIS-HOH)
                        distances = self.standardize_axial_positions(distances, verbose=verbose)

                        # Check if axial ligand combination is allowed
                        axial1_resname = distances.get('axial1_resname')
                        axial2_resname = distances.get('axial2_resname')
                        if not self.is_allowed_axial_combination(axial1_resname, axial2_resname):
                            if verbose:
                                print(f"  Skipping {pdb_id} - axial combination {axial1_resname}-{axial2_resname} not allowed")
                            continue
                        # Process axial ligand 1 (residue 2)
                        if 'axial1_pdb_distance' in distances and 'axial1_xyz_distance' in distances:
                            difference1 = abs(distances['axial1_xyz_distance'] - distances['axial1_pdb_distance'])  # Make absolute
                            results.append({
                                'PDB_ID': pdb_id,
                                'Axial_Ligand': 1,
                                'Axial_Resname': distances['axial1_resname'],
                                'PDB_Iron_Axial_Distance': distances['axial1_pdb_distance'],
                                'XYZ_Iron_Axial_Distance': distances['axial1_xyz_distance'],
                                'Distance_Difference': difference1
                            })
                            if verbose:
                                print(f"  Axial 1 ({distances['axial1_resname']}) - PDB: {distances['axial1_pdb_distance']:.3f} Å, XYZ: {distances['axial1_xyz_distance']:.3f} Å, Diff: {difference1:.3f} Å")
                        
                        # Process axial ligand 2 (residue 3)
                        if 'axial2_pdb_distance' in distances and 'axial2_xyz_distance' in distances:
                            difference2 = abs(distances['axial2_xyz_distance'] - distances['axial2_pdb_distance'])  # Make absolute
                            results.append({
                                'PDB_ID': pdb_id,
                                'Axial_Ligand': 2,
                                'Axial_Resname': distances['axial2_resname'],
                                'PDB_Iron_Axial_Distance': distances['axial2_pdb_distance'],
                                'XYZ_Iron_Axial_Distance': distances['axial2_xyz_distance'],
                                'Distance_Difference': difference2
                            })
                            if verbose:
                                print(f"  Axial 2 ({distances['axial2_resname']}) - PDB: {distances['axial2_pdb_distance']:.3f} Å, XYZ: {distances['axial2_xyz_distance']:.3f} Å, Diff: {difference2:.3f} Å")
                        
                        if not ('axial1_pdb_distance' in distances or 'axial2_pdb_distance' in distances):
                            if verbose:
                                print(f"  No valid axial distances found")
                    else:
                        if verbose:
                            print(f"  Failed to calculate Fe-axial distances")
        
        return pd.DataFrame(results)
    
    def create_rmsd_histogram(self, df, output_file="plots/rmsd_histogram.png", exclude_suspicious=True, verbose=None):
        """
        Create histogram of RMSD values with statistics, optionally excluding suspicious PDB IDs.
        Bars are colored by axial ligand combinations for the most common 5 combinations.
        
        Args:
            df (pandas.DataFrame): DataFrame with RMSD values
            output_file (str): Output filename for the histogram
            exclude_suspicious (bool): Whether to exclude suspicious PDB IDs from plotting
            verbose: Override the instance verbose setting
        """
        if verbose is None:
            verbose = self.verbose
        if df.empty:
            if verbose:
                print("No data to plot!")
            return
        
        # Filter out suspicious PDB IDs if requested
        if exclude_suspicious:
            suspicious_ids = self.get_suspicious_pdb_ids()
            if suspicious_ids:
                initial_count = len(df)
                df = df[~df['PDB_ID'].isin(suspicious_ids)].copy()
                filtered_count = initial_count - len(df)
                if verbose:
                    print(f"Excluded {filtered_count} suspicious structures from RMSD histogram")
                if filtered_count > 0 and verbose:
                    print(f"Excluded PDB IDs: {sorted(suspicious_ids)}")
        
        if df.empty:
            print("No data remaining after filtering suspicious structures!")
            return
        
        # Get axial ligand data for coloring
        plot_df = df.copy()
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        if 'PDB_ID' in plot_df.columns:
            # Get axial ligand information
            valid_pdb_ids = plot_df['PDB_ID'].unique()
            axial_df = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
            
            if not axial_df.empty:
                # Create axial combinations by pivoting data
                axial_pivot = axial_df.pivot_table(
                    index='PDB_ID', 
                    columns='Axial_Ligand', 
                    values='Axial_Resname', 
                    aggfunc='first'
                ).reset_index()
                
                # Create proper axial combinations (axial1-axial2) with HSD normalization
                if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                    # Normalize ligand names (HSD -> HIS) and create sorted combination
                    # Note: Axial ligands from processed_output.csv are already sorted,
                    # but we sort here for consistency and to handle data from other sources
                    def create_normalized_combo(row):
                        if pd.notna(row[1]) and pd.notna(row[2]):
                            norm_axial1 = self.normalize_axial_ligand_name(row[1])
                            norm_axial2 = self.normalize_axial_ligand_name(row[2])
                            return '-'.join(sorted([norm_axial1, norm_axial2]))
                        else:
                            return 'Unknown-Unknown'
                    axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                elif 1 in axial_pivot.columns:
                    norm_axial1 = axial_pivot[1].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = norm_axial1.astype(str) + '-Unknown'
                elif 2 in axial_pivot.columns:
                    norm_axial2 = axial_pivot[2].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = 'Unknown-' + norm_axial2.astype(str)
                else:
                    axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
                
                axial_ligand_info = axial_pivot[['PDB_ID', 'Axial_Combo']]
                
                # Merge with plot data
                plot_df = plot_df.merge(axial_ligand_info, on='PDB_ID', how='left')
                
                # Check if we have valid axial combination data
                if 'Axial_Combo' in plot_df.columns and plot_df['Axial_Combo'].notna().sum() > 0:
                    use_axial_colors = True
                    # Only use the 5 most common combinations (as defined in allowed_axial_combinations)
                    allowed_combos = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
                    axial_combinations = [combo for combo in allowed_combos 
                                        if combo in plot_df['Axial_Combo'].values]
                    
                    # Use centralized color mapping for axial combinations
                    axial_color_map = get_axial_ligand_colors(axial_combinations)
                    
                    if verbose:
                        print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
                        for combo in axial_combinations:
                            count = (plot_df['Axial_Combo'] == combo).sum()
                            print(f"  {combo}: {count} structures")
        
        # Calculate statistics
        rmsd_values = plot_df['RMSD']
        min_rmsd = rmsd_values.min()
        max_rmsd = rmsd_values.max()
        mean_rmsd = rmsd_values.mean()
        std_rmsd = rmsd_values.std()
        
        if verbose:
            print(f"\nRMSD Statistics:")
            print(f"Minimum: {min_rmsd:.3f} Å")
            print(f"Maximum: {max_rmsd:.3f} Å")
            print(f"Mean: {mean_rmsd:.3f} Å")
            print(f"Standard deviation: {std_rmsd:.3f} Å")
            print(f"Number of structures: {len(rmsd_values)}")
        
        # Create histogram
        plt.figure(figsize=(10, 6))
        
        if use_axial_colors and axial_color_map is not None and axial_combinations is not None:
            # Create stacked histogram by axial combination
            combo_data = {}
            for combo in axial_combinations:
                combo_mask = plot_df['Axial_Combo'] == combo
                combo_values = plot_df.loc[combo_mask, 'RMSD'].dropna()
                combo_data[combo] = combo_values.values
            
            # Create stacked histogram - filter out empty combinations
            valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
            data_arrays = [combo_data[combo] for combo in valid_combos]
            colors = [axial_color_map[combo] for combo in valid_combos]
            labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
            
            # Create stacked histogram
            plt.hist(data_arrays, bins=np.linspace(0, 0.4, 51), color=colors, alpha=0.8, 
                    edgecolor='black', label=labels, stacked=True)
            
            # Add legend
            plt.legend(title="Axial Ligand Combinations", fontsize=12, title_fontsize=14)
            
        else:
            # Fall back to simple histogram
            plt.hist(rmsd_values, bins=np.linspace(0, 0.4, 51), alpha=0.7, color='skyblue', edgecolor='black')
        
        # Formatting
        plt.xlabel('RMSD (Å)', fontsize=16)
        plt.ylabel('Frequency', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.xlim(0, 0.4)
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        if verbose:
            print(f"\nHistogram saved as: {output_file}")
    
    def create_distance_difference_histogram(self, df, output_file="plots/iron-ax-dist_diff.png", exclude_suspicious=True):
        """
        Create separate histograms of Fe-axial ligand distance differences for each axial ligand, 
        where each bar is colored by axial ligand combinations for the most common 5 combinations.
        Optionally excludes suspicious PDB IDs from plotting.
        """
        if df.empty or 'Distance_Difference' not in df.columns:
            if self.verbose:
                print("No distance difference data to plot!")
            return

        # Filter out suspicious structures
        if exclude_suspicious:
            suspicious_ids = self.get_suspicious_pdb_ids()
            if suspicious_ids:
                df = df[~df['PDB_ID'].isin(suspicious_ids)]
                if self.verbose:
                    print(f"Excluded {len(suspicious_ids)} suspicious structures")

        if df.empty:
            print("No data remaining after filtering suspicious structures!")
            return

        # Get axial ligand data for coloring
        plot_df = df.copy()
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        if 'PDB_ID' in plot_df.columns:
            # Get axial ligand information
            valid_pdb_ids = plot_df['PDB_ID'].unique()
            axial_df = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
            
            if not axial_df.empty:
                # Create axial combinations by pivoting data
                axial_pivot = axial_df.pivot_table(
                    index='PDB_ID', 
                    columns='Axial_Ligand', 
                    values='Axial_Resname', 
                    aggfunc='first'
                ).reset_index()
                
                # Create proper axial combinations (axial1-axial2) with HSD normalization
                if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                    # Normalize ligand names (HSD -> HIS) and create sorted combination
                    # Note: Axial ligands from processed_output.csv are already sorted,
                    # but we sort here for consistency and to handle data from other sources
                    def create_normalized_combo(row):
                        if pd.notna(row[1]) and pd.notna(row[2]):
                            norm_axial1 = self.normalize_axial_ligand_name(row[1])
                            norm_axial2 = self.normalize_axial_ligand_name(row[2])
                            return '-'.join(sorted([norm_axial1, norm_axial2]))
                        else:
                            return 'Unknown-Unknown'
                    axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                elif 1 in axial_pivot.columns:
                    norm_axial1 = axial_pivot[1].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = norm_axial1.astype(str) + '-Unknown'
                elif 2 in axial_pivot.columns:
                    norm_axial2 = axial_pivot[2].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = 'Unknown-' + norm_axial2.astype(str)
                else:
                    axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
                
                axial_ligand_info = axial_pivot[['PDB_ID', 'Axial_Combo']]
                
                # Merge with plot data
                plot_df = plot_df.merge(axial_ligand_info, on='PDB_ID', how='left')
                
                # Check if we have valid axial combination data
                if 'Axial_Combo' in plot_df.columns and plot_df['Axial_Combo'].notna().sum() > 0:
                    use_axial_colors = True
                    # Only use the 5 most common combinations (as defined in allowed_axial_combinations)
                    allowed_combos = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
                    axial_combinations = [combo for combo in allowed_combos 
                                        if combo in plot_df['Axial_Combo'].values]
                    
                    # Use centralized color mapping for axial combinations
                    axial_color_map = get_axial_ligand_colors(axial_combinations)
                    
                    if self.verbose:
                        print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
                        for combo in axial_combinations:
                            count = (plot_df['Axial_Combo'] == combo).sum()
                            print(f"  {combo}: {count} structures")

        # Create figure with subplots for each axial ligand
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Function to create histogram for given axial ligand
        def create_single_histogram(ax, axial_ligand_num, title):
            axial_data = plot_df[plot_df['Axial_Ligand'] == axial_ligand_num]
            
            if axial_data.empty:
                #ax.set_title(f"{title} (No Data)", fontsize=16)
                ax.set_xlabel('Fe-Axial Distance Difference (Å)', fontsize=14)
                ax.set_ylabel('Frequency', fontsize=14)
                return
            
            distances = axial_data['Distance_Difference'].dropna()
            
            if use_axial_colors and axial_color_map is not None and axial_combinations is not None:
                # Create stacked histogram by axial combination
                combo_data = {}
                for combo in axial_combinations:
                    combo_mask = axial_data['Axial_Combo'] == combo
                    combo_values = axial_data.loc[combo_mask, 'Distance_Difference'].dropna()
                    combo_data[combo] = combo_values.values
                
                # Create stacked histogram - filter out empty combinations
                valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
                data_arrays = [combo_data[combo] for combo in valid_combos]
                colors = [axial_color_map[combo] for combo in valid_combos]
                labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
                
                # Create stacked histogram
                ax.hist(data_arrays, bins=np.linspace(0, 0.7, 51), color=colors, alpha=0.8, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add dual legends with highlighted ligands
                if ax == ax1:  # Left subplot - highlight first ligand
                    # Create labels with first ligand in bold
                    bold_labels = []
                    for combo in valid_combos:
                        first_ligand, second_ligand = combo.split('-')
                        count = len(combo_data[combo])
                        bold_labels.append(f"$\\mathbf{{{first_ligand}}}$-{second_ligand} (n={count})")
                    
                    ax.legend(bold_labels, title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                             frameon=False, fancybox=False, shadow=False)
                    
                elif ax == ax2:  # Right subplot - highlight second ligand
                    # Create labels with second ligand in bold
                    bold_labels = []
                    for combo in valid_combos:
                        first_ligand, second_ligand = combo.split('-')
                        count = len(combo_data[combo])
                        bold_labels.append(f"{first_ligand}-$\\mathbf{{{second_ligand}}}$ (n={count})")
                    
                    ax.legend(bold_labels, title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                             frameon=False, fancybox=False, shadow=False)
                
            else:
                # Fall back to simple histogram
                ax.hist(distances, bins=np.linspace(0, 0.7, 51), alpha=0.7, color='lightcoral', edgecolor='black')
            
            # Formatting
            ax.set_xlabel('Fe-Axial Distance Difference (Å)', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlim(0, 0.7)
            #ax.set_title(title, fontsize=16)
            ax.grid(True, alpha=0.3)
            
            # Print statistics
            if self.verbose:
                print(f"\n{title} Statistics:")
                print(f"n={len(distances)}, mean={np.mean(distances):.3f}, std={np.std(distances):.3f}")
        
        # Create histograms for each axial ligand
        create_single_histogram(ax1, 1, 'Axial Ligand 1 Distance Differences')
        create_single_histogram(ax2, 2, 'Axial Ligand 2 Distance Differences')

        # Add a) and b) annotations to subplots
        ax1.text(0.05, 0.95, 'a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')
        ax2.text(0.05, 0.95, 'b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')

        # Calculate overall statistics
        all_values = plot_df['Distance_Difference']
        mean_val = all_values.mean()
        std_val = all_values.std()
        min_val = all_values.min()
        max_val = all_values.max()

        if self.verbose:
            print(f"\nOverall Distance Difference Statistics:")
            print(f"Min: {min_val:.3f} Å, Max: {max_val:.3f} Å")
            print(f"Mean: {mean_val:.3f} ± {std_val:.3f} Å")
            print(f"Total measurements: {len(all_values)}")

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()

        if self.verbose:
            print(f"Stacked distance difference histogram saved as: {output_file}")

    def create_actual_distances_histogram(self, df, output_file="plots/iron-ax-dist.png", exclude_suspicious=True):
        """
        Create separate histograms of Fe-axial ligand distances for each axial ligand, 
        where each bar is colored by axial ligand combinations for the most common 5 combinations.
        Optionally excludes suspicious PDB IDs from plotting.
        """
        if df.empty or 'XYZ_Iron_Axial_Distance' not in df.columns:
            if self.verbose:
                print("No distance data to plot!")
            return

        # Filter out suspicious structures
        if exclude_suspicious:
            suspicious_ids = self.get_suspicious_pdb_ids()
            if suspicious_ids:
                df = df[~df['PDB_ID'].isin(suspicious_ids)]
                if self.verbose:
                    print(f"Excluded {len(suspicious_ids)} suspicious structures")

        if df.empty:
            print("No data remaining after filtering suspicious structures!")
            return

        # Get axial ligand data for coloring
        plot_df = df.copy()
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        if 'PDB_ID' in plot_df.columns:
            # Get axial ligand information
            valid_pdb_ids = plot_df['PDB_ID'].unique()
            axial_df = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
            
            if not axial_df.empty:
                # Create axial combinations by pivoting data
                axial_pivot = axial_df.pivot_table(
                    index='PDB_ID', 
                    columns='Axial_Ligand', 
                    values='Axial_Resname', 
                    aggfunc='first'
                ).reset_index()
                
                # Create proper axial combinations (axial1-axial2) with HSD normalization
                if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                    # Normalize ligand names (HSD -> HIS) and create sorted combination
                    def create_normalized_combo(row):
                        if pd.notna(row[1]) and pd.notna(row[2]):
                            norm_axial1 = self.normalize_axial_ligand_name(row[1])
                            norm_axial2 = self.normalize_axial_ligand_name(row[2])
                            if norm_axial1 and norm_axial2:
                                return '-'.join(sorted([norm_axial1, norm_axial2]))
                        return 'Unknown-Unknown'
                    axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                elif 1 in axial_pivot.columns:
                    norm_axial1 = axial_pivot[1].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = norm_axial1.astype(str) + '-Unknown'
                elif 2 in axial_pivot.columns:
                    norm_axial2 = axial_pivot[2].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = 'Unknown-' + norm_axial2.astype(str)
                else:
                    axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
                
                axial_ligand_info = axial_pivot[['PDB_ID', 'Axial_Combo']]
                
                # Merge with plot data
                plot_df = plot_df.merge(axial_ligand_info, on='PDB_ID', how='left')
                
                # Check if we have valid axial combination data
                if 'Axial_Combo' in plot_df.columns and plot_df['Axial_Combo'].notna().sum() > 0:
                    use_axial_colors = True
                    # Only use the 5 most common combinations (as defined in allowed_axial_combinations)
                    allowed_combos = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
                    axial_combinations = [combo for combo in allowed_combos 
                                        if combo in plot_df['Axial_Combo'].values]
                    
                    # Use centralized color mapping for axial combinations
                    axial_color_map = get_axial_ligand_colors(axial_combinations)
                    
                    if self.verbose:
                        print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
                        for combo in axial_combinations:
                            count = (plot_df['Axial_Combo'] == combo).sum()
                            print(f"  {combo}: {count} structures")

        # Create figure with subplots for each axial ligand
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Function to create histogram for given axial ligand
        def create_single_histogram(ax, axial_ligand_num, title):
            axial_data = plot_df[plot_df['Axial_Ligand'] == axial_ligand_num]
            
            if axial_data.empty:
                ax.set_xlabel('Fe-Axial Distance (Å)', fontsize=14)
                ax.set_ylabel('Frequency', fontsize=14)
                return
            
            distances = axial_data['XYZ_Iron_Axial_Distance'].dropna()
            
            if use_axial_colors and axial_color_map is not None and axial_combinations is not None:
                # Create stacked histogram by axial combination
                combo_data = {}
                for combo in axial_combinations:
                    combo_mask = axial_data['Axial_Combo'] == combo
                    combo_values = axial_data.loc[combo_mask, 'XYZ_Iron_Axial_Distance'].dropna()
                    combo_data[combo] = combo_values.values
                
                # Create stacked histogram - filter out empty combinations
                valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
                data_arrays = [combo_data[combo] for combo in valid_combos]
                colors = [axial_color_map[combo] for combo in valid_combos]
                labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
                
                # Create stacked histogram
                ax.hist(data_arrays, bins=np.linspace(1.5, 3.5, 51), color=colors, alpha=0.8, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add dual legends with highlighted ligands
                if ax == ax1:  # Left subplot - highlight first ligand
                    # Create labels with first ligand in bold
                    bold_labels = []
                    for combo in valid_combos:
                        first_ligand, second_ligand = combo.split('-')
                        count = len(combo_data[combo])
                        bold_labels.append(f"$\\mathbf{{{first_ligand}}}$-{second_ligand} (n={count})")
                    
                    ax.legend(bold_labels, title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                             frameon=False, fancybox=False, shadow=False)
                    
                elif ax == ax2:  # Right subplot - highlight second ligand
                    # Create labels with second ligand in bold
                    bold_labels = []
                    for combo in valid_combos:
                        first_ligand, second_ligand = combo.split('-')
                        count = len(combo_data[combo])
                        bold_labels.append(f"{first_ligand}-$\\mathbf{{{second_ligand}}}$ (n={count})")
                    
                    ax.legend(bold_labels, title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                             frameon=False, fancybox=False, shadow=False)
                
            else:
                # Fall back to simple histogram
                ax.hist(distances, bins=np.linspace(1.5, 3.5, 51), alpha=0.7, color='lightcoral', edgecolor='black')
            
            # Formatting
            ax.set_xlabel('Fe-Axial Distance (Å)', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlim(1.5, 3.5)
            ax.grid(True, alpha=0.3)
            
            # Print statistics
            if self.verbose:
                print(f"\n{title} Statistics:")
                print(f"n={len(distances)}, mean={np.mean(distances):.3f}, std={np.std(distances):.3f}")
        
        # Create histograms for each axial ligand
        create_single_histogram(ax1, 1, 'Axial Ligand 1 Distances')
        create_single_histogram(ax2, 2, 'Axial Ligand 2 Distances')

        # Add a) and b) annotations to subplots
        ax1.text(0.05, 0.95, 'a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')
        ax2.text(0.05, 0.95, 'b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')

        # Calculate overall statistics
        all_values = plot_df['XYZ_Iron_Axial_Distance']
        mean_val = all_values.mean()
        std_val = all_values.std()
        min_val = all_values.min()
        max_val = all_values.max()

        if self.verbose:
            print(f"\nOverall Distance Statistics:")
            print(f"Min: {min_val:.3f} Å, Max: {max_val:.3f} Å")
            print(f"Mean: {mean_val:.3f} ± {std_val:.3f} Å")
            print(f"Total measurements: {len(all_values)}")

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()

        if self.verbose:
            print(f"Actual distances histogram saved as: {output_file}")


    
    def run_rmsd_analysis(self):
        """Run RMSD analysis on all structures."""
        df_rmsd = self.process_all_structures()
        if not df_rmsd.empty:
            # Save RMSD results to CSV
            csv_file = str(derived_table_path("rmsd_results.csv"))
            df_rmsd.to_csv(csv_file, index=False)
            if self.verbose:
                print(f"\nRMSD results saved to: {csv_file}")
            
            # Create RMSD histogram (excludes suspicious structures by default)
            self.create_rmsd_histogram(df_rmsd)
            if self.verbose:
                print(f"RMSD analysis: Processed {len(df_rmsd)} structure pairs successfully.")
            return df_rmsd
        else:
            print("No valid structure pairs found for RMSD analysis!")
            return pd.DataFrame()
    
    def run_iron_axial_distance_analysis(self):
        """Run iron-axial ligand distance analysis."""        
        df_distances = self.process_iron_axial_distances()
        if not df_distances.empty:
            # Save distance results to CSV
            distance_csv = str(derived_table_path("iron_axial_distances.csv"))
            df_distances.to_csv(distance_csv, index=False)
            if self.verbose:
                print(f"\nIron-axial distance results saved to: {distance_csv}")
            
            # Create distance difference histogram (excludes suspicious structures by default)
            print("create_distance_difference_histogram called")
            self.create_distance_difference_histogram(df_distances)
            
            # Create actual distances histogram 
            self.create_actual_distances_histogram(df_distances)
            
            if self.verbose:
                print(f"Distance analysis: Processed {len(df_distances)} structure pairs successfully.")
            return df_distances
        else:
            print("No valid structure pairs found for distance analysis! 2")
            return pd.DataFrame()
    
    def calculate_heme_plane_distance(self, pdb_file, xyz_file):
        """
        Calculate the distance of heme Fe from the heme plane for both PDB and XYZ structures.
        The heme plane is defined by the four porphyrin nitrogen atoms.
        
        Args:
            pdb_file (str): Path to PDB file
            xyz_file (str): Path to XYZ file
            
        Returns:
            dict: Dictionary with Fe-plane distances for both structures, or None if failed
        """
        try:
            # Load PDB structure
            pdb_universe = mda.Universe(pdb_file)
            
            # Find Fe and nitrogen atoms in PDB
            iron_pdb = pdb_universe.select_atoms("name FE")
            nitrogens_pdb = pdb_universe.select_atoms("resname HEM and (name NA or name NB or name NC or name ND)")
            
            if len(iron_pdb) == 0 and self.verbose:
                print(f"No Fe atom found in PDB: {pdb_file}")
                return None
            
            if len(nitrogens_pdb) != 4 and self.verbose:
                print(f"Expected 4 nitrogen atoms, found {len(nitrogens_pdb)} in PDB: {pdb_file}")
                return None
            
            iron_pos_pdb = iron_pdb.positions[0]
            nitrogen_positions_pdb = nitrogens_pdb.positions
            
            # Calculate heme plane from nitrogen atoms (PDB)
            plane_normal_pdb, plane_distance_pdb = self._calculate_plane_from_points(nitrogen_positions_pdb)
            iron_plane_distance_pdb = self._point_to_plane_distance(iron_pos_pdb, plane_normal_pdb, plane_distance_pdb)
            
            # Load XYZ structure and match atoms
            xyz_elements, xyz_coords = read_xyz_file(xyz_file)
            
            # Get Fe index from PDB (should be same in XYZ)
            iron_pdb_idx = iron_pdb[0].index
            if iron_pdb_idx >= len(xyz_coords) and self.verbose:
                print(f"Fe index {iron_pdb_idx} exceeds XYZ atom count {len(xyz_coords)}")
                return None
            
            iron_pos_xyz = xyz_coords[iron_pdb_idx]
            
            # Get nitrogen positions from XYZ using PDB indices
            nitrogen_indices = [atom.index for atom in nitrogens_pdb]
            nitrogen_positions_xyz = []
            
            for idx in nitrogen_indices:
                if idx >= len(xyz_coords) and self.verbose:
                    print(f"Nitrogen index {idx} exceeds XYZ atom count {len(xyz_coords)}")
                    return None
                nitrogen_positions_xyz.append(xyz_coords[idx])
            
            nitrogen_positions_xyz = np.array(nitrogen_positions_xyz)
            
            # Calculate heme plane from nitrogen atoms (XYZ)
            plane_normal_xyz, plane_distance_xyz = self._calculate_plane_from_points(nitrogen_positions_xyz)
            iron_plane_distance_xyz = self._point_to_plane_distance(iron_pos_xyz, plane_normal_xyz, plane_distance_xyz)
            
            return {
                'pdb_iron_plane_distance': abs(iron_plane_distance_pdb),
                'xyz_iron_plane_distance': abs(iron_plane_distance_xyz),
                'distance_difference': abs(abs(iron_plane_distance_xyz) - abs(iron_plane_distance_pdb)),
                'pdb_iron_position': iron_pos_pdb.copy(),
                'xyz_iron_position': iron_pos_xyz.copy(),
                'pdb_nitrogen_positions': nitrogen_positions_pdb.copy(),
                'xyz_nitrogen_positions': nitrogen_positions_xyz.copy()
            }
            
        except Exception as e:
            print(f"Error calculating Fe-plane distance for {pdb_file} and {xyz_file}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def _calculate_plane_from_points(self, points):
        """
        Calculate plane normal and distance from a set of points using least squares fitting.
        
        Args:
            points (np.ndarray): Array of 3D points (N x 3)
            
        Returns:
            tuple: (normal_vector, plane_distance) where normal is unit vector and distance is from origin
        """
        centroid = np.mean(points, axis=0)
        centered_points = points - centroid
        _, _, V = np.linalg.svd(centered_points)
        normal = V[-1]  
        normal = normal / np.linalg.norm(normal)
        plane_distance = np.dot(normal, centroid)
        return normal, plane_distance
    
    def _point_to_plane_distance(self, point, plane_normal, plane_distance):
        """
        Calculate the signed distance from a point to a plane.
        
        Args:
            point (np.ndarray): 3D point coordinates
            plane_normal (np.ndarray): Unit normal vector of the plane
            plane_distance (float): Distance from origin to plane along normal
            
        Returns:
            float: Signed distance from point to plane
        """
        return np.dot(plane_normal, point) - plane_distance
    
    def process_iron_plane_distances(self):
        """
        Process all structures and calculate Fe-heme plane distances.
        Also extracts axial ligand combination for each structure.
        
        Returns:
            pandas.DataFrame: DataFrame with PDB IDs, Fe-plane distances, and axial ligand combinations
        """
        results = []
        
        for subdir in self.pdb_base_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name
                
                pdb_file = subdir / f"{pdb_id}_system_protonated.pdb"
                xyz_file = subdir / f"{pdb_id}_g16.xyz"
                
                if pdb_file.exists() and xyz_file.exists():
                    if self.verbose:
                        print(f"Processing Fe-plane distance for {pdb_id}...")
                    
                    plane_data = self.calculate_heme_plane_distance(str(pdb_file), str(xyz_file))
                    
                    # Also get axial ligand information and filter early
                    axial_combination = "Unknown"
                    axial_distances = get_iron_axial_distances(str(pdb_file), str(xyz_file), verbose=False)
                    if axial_distances is not None:
                        # Extract axial ligand names
                        axial1 = axial_distances.get('axial1_resname', 'UNK')
                        axial2 = axial_distances.get('axial2_resname', 'UNK')
                        
                        # Check if axial ligand combination is allowed (skip if not)
                        if not self.is_allowed_axial_combination(axial1, axial2):
                            if self.verbose:
                                print(f"  Skipping {pdb_id} - axial combination {axial1}-{axial2} not allowed")
                            continue
                            
                        # Create sorted combination string using normalized names for consistency
                        # Note: Axial ligands from processed_output.csv are already sorted alphabetically
                        norm_axial1 = self.normalize_axial_ligand_name(axial1)
                        norm_axial2 = self.normalize_axial_ligand_name(axial2)
                        axial_combination = '-'.join(sorted([norm_axial1, norm_axial2]))
                    else:
                        # Skip if we can't determine axial ligand information
                        if self.verbose:
                            print(f"  Skipping {pdb_id} - could not determine axial ligand combination")
                        continue
                    
                    if plane_data is not None:
                        results.append({
                            'PDB_ID': pdb_id,
                            'PDB_Iron_Plane_Distance': plane_data['pdb_iron_plane_distance'],
                            'XYZ_Iron_Plane_Distance': plane_data['xyz_iron_plane_distance'],
                            'Plane_Distance_Difference': plane_data['distance_difference'],
                            'Axial_Combination': axial_combination
                        })
                        if self.verbose:
                            print(f"  PDB: {plane_data['pdb_iron_plane_distance']:.3f} Å, XYZ: {plane_data['xyz_iron_plane_distance']:.3f} Å, Diff: {plane_data['distance_difference']:.3f} Å")
                            print(f"  Axial combination: {axial_combination}")
                    else:
                        if self.verbose:
                            print(f"  Failed to calculate Fe-plane distance")
        
        return pd.DataFrame(results)
    
    def create_iron_plane_histogram(self, df, output_file="plots/iron_plane_distances.png", 
                                   show_charge_mult_breakdown=True, charge_mult_data=None):
        """
        Create histogram of Fe-heme plane distance differences with optional charge-multiplicity breakdown.
        
        Args:
            df (pandas.DataFrame): DataFrame with Fe-plane distance data
            output_file (str): Output filename for the histogram
            show_charge_mult_breakdown (bool): Whether to show breakdown by charge-multiplicity combinations
            charge_mult_data (pandas.DataFrame): DataFrame with PDB_ID, charge, and multiplicity columns
        """
        if df.empty or 'Plane_Distance_Difference' not in df.columns and self.verbose:
            print("No Fe-plane distance data to plot!")
            return
        
        # Prepare working dataframe
        plot_df = df.copy()
        
        # Try to merge with charge-multiplicity data if provided
        if show_charge_mult_breakdown and charge_mult_data is not None:
            if 'PDB_ID' in plot_df.columns and all(col in charge_mult_data.columns for col in ['PDB_ID', 'charge', 'multiplicity']):
                # Merge the datasets
                plot_df = plot_df.merge(charge_mult_data[['PDB_ID', 'charge', 'multiplicity']], 
                                       on='PDB_ID', how='left')
                
                # Create charge-multiplicity combination labels
                plot_df['charge_mult_combo'] = plot_df.apply(
                    lambda row: f"q={int(row['charge'])},m={int(row['multiplicity'])}" 
                    if pd.notna(row['charge']) and pd.notna(row['multiplicity']) else "Unknown", axis=1
                )
                
                # Check if we have valid charge-multiplicity data
                valid_combo_mask = plot_df['charge_mult_combo'] != "Unknown"
                if valid_combo_mask.sum() > 0:
                    show_breakdown = True
                    if self.verbose:
                        print(f"Successfully merged charge-multiplicity data for {valid_combo_mask.sum()}/{len(plot_df)} structures")
                else:
                    show_breakdown = False
                    if self.verbose:
                        print("Warning: No valid charge-multiplicity combinations found after merge")
            else:
                show_breakdown = False
                if self.verbose:
                    print("Warning: Cannot merge charge-multiplicity data - missing required columns")
        else:
            show_breakdown = False
            if not show_charge_mult_breakdown:
                if self.verbose:
                    print("Charge-multiplicity breakdown disabled")
            else:
                print("No charge-multiplicity data provided")
        
        # Calculate overall statistics
        distance_diffs = plot_df['Plane_Distance_Difference']
        min_diff = distance_diffs.min()
        max_diff = distance_diffs.max()
        mean_diff = distance_diffs.mean()
        std_diff = distance_diffs.std()
        
        if self.verbose:
            print(f"\nFe Porphyrin Plane Dist. Difference Statistics:")
            print(f"Minimum: {min_diff:.3f} Å")
            print(f"Maximum: {max_diff:.3f} Å")
            print(f"Mean: {mean_diff:.3f} Å")
            print(f"Standard deviation: {std_diff:.3f} Å")
            print(f"Number of structures: {len(distance_diffs)}")
        
        # Define color scheme for charge-multiplicity combinations
        color_map = {
            'q=0,m=1': '#0072B2',  # Blue
            'q=0,m=5': '#009E73',  # Green  
            'q=1,m=2': '#D55E00',  # Orange
            'q=1,m=6': '#CC79A7'   # Purple
        }
        
        # Create figure
        plt.figure(figsize=(12, 7))
        
        if show_breakdown and 'charge_mult_combo' in plot_df.columns:
            # Get unique combinations and their colors
            unique_combos = sorted([combo for combo in plot_df['charge_mult_combo'].unique() if combo != "Unknown"])
            
            if len(unique_combos) > 0:
                # Prepare data arrays for stacked histogram
                combo_data = {}
                combo_counts = {}
                combo_stats = {}
                
                for combo in unique_combos:
                    combo_mask = plot_df['charge_mult_combo'] == combo
                    combo_data[combo] = plot_df.loc[combo_mask, 'Plane_Distance_Difference'].values
                    combo_counts[combo] = len(combo_data[combo])
                    
                    if len(combo_data[combo]) > 0:
                        combo_stats[combo] = {
                            'mean': np.mean(combo_data[combo]),
                            'std': np.std(combo_data[combo]),
                            'min': np.min(combo_data[combo]),
                            'max': np.max(combo_data[combo])
                        }
                
                # Create stacked histogram
                data_arrays = [combo_data[combo] for combo in unique_combos]
                colors = [color_map.get(combo, '#808080') for combo in unique_combos]
                labels = [f"{combo} (n={combo_counts[combo]})" for combo in unique_combos]
                
                # Create stacked histogram using matplotlib's built-in stacking
                plt.hist(data_arrays, bins=np.linspace(0, 0.12, 51), alpha=0.7, color=colors, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add legend
                plt.legend(loc='upper right', fontsize=14)
                if self.verbose:
                    print(f"\nBreakdown by charge-multiplicity combination:")
                for combo in unique_combos:
                    if combo in combo_stats:
                        stats = combo_stats[combo]
                        if self.verbose:
                            print(f"{combo}: {combo_counts[combo]} structures, Mean: {stats['mean']:.3f}±{stats['std']:.3f} Å")
            else:
                if self.verbose:
                    print("No valid charge-multiplicity combinations found, falling back to simple histogram")
                show_breakdown = False
        
        if not show_breakdown:
            # Fall back to simple histogram
            plt.hist(distance_diffs, bins=np.linspace(0, 0.12, 51), alpha=0.7, color='lightgreen', edgecolor='black')
        
        # Formatting for publication quality
        plt.xlabel('Fe Porphyrin Plane Distance Difference (Å)', fontsize=16)
        plt.ylabel('Frequency', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.xlim(0, 0.12)
        plt.grid(True, alpha=0.3)
        
        # Save plot in publication quality
        plt.tight_layout()
        plt.savefig(output_file, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.show()
        if self.verbose:
            print(f"\nFe-plane distance histogram saved as: {output_file}")

    def create_iron_displacement_barplot_by_charge_mult(self, df, output_file="plots/iron_displacement_charge_mult_bars.png",
                                                        charge_mult_data=None, bin_width=0.01, max_distance=0.12):
        """
        Create bar plot of Fe-heme plane distance differences with bars colored by charge-multiplicity states.
        Each bar represents a distance bin, and is colored proportionally based on the charge-multiplicity
        distribution within that bin.

        Args:
            df (pandas.DataFrame): DataFrame with Fe-plane distance data (must have 'Plane_Distance_Difference' column)
            output_file (str): Output filename for the plot
            charge_mult_data (pandas.DataFrame): DataFrame with PDB_ID, charge, and multiplicity columns
            bin_width (float): Width of each distance bin in Angstroms (default: 0.01)
            max_distance (float): Maximum distance to plot in Angstroms (default: 0.12)
        """
        if df.empty or 'Plane_Distance_Difference' not in df.columns:
            if self.verbose:
                print("No Fe-plane distance data to plot!")
            return

        # Prepare working dataframe
        plot_df = df.copy()

        # Merge with charge-multiplicity data if provided
        if charge_mult_data is not None:
            if 'PDB_ID' in plot_df.columns and all(col in charge_mult_data.columns for col in ['PDB_ID', 'charge', 'multiplicity']):
                # Merge the datasets
                plot_df = plot_df.merge(charge_mult_data[['PDB_ID', 'charge', 'multiplicity']],
                                       on='PDB_ID', how='left')

                # Create charge-multiplicity combination labels
                plot_df['charge_mult_combo'] = plot_df.apply(
                    lambda row: f"q={int(row['charge'])},m={int(row['multiplicity'])}"
                    if pd.notna(row['charge']) and pd.notna(row['multiplicity']) else "Unknown", axis=1
                )

                # Check if we have valid charge-multiplicity data
                valid_combo_mask = plot_df['charge_mult_combo'] != "Unknown"
                if valid_combo_mask.sum() > 0:
                    has_charge_mult = True
                    if self.verbose:
                        print(f"Successfully merged charge-multiplicity data for {valid_combo_mask.sum()}/{len(plot_df)} structures")
                else:
                    has_charge_mult = False
                    if self.verbose:
                        print("Warning: No valid charge-multiplicity combinations found after merge")
            else:
                has_charge_mult = False
                if self.verbose:
                    print("Warning: Cannot merge charge-multiplicity data - missing required columns")
        else:
            has_charge_mult = False
            if self.verbose:
                print("No charge-multiplicity data provided")

        if not has_charge_mult:
            if self.verbose:
                print("Cannot create charge-multiplicity barplot without valid charge-multiplicity data")
            return

        # Create distance bins
        bins = np.arange(0, max_distance + bin_width, bin_width)
        bin_centers = bins[:-1] + bin_width / 2
        bin_labels = [f"{b:.3f}-{b+bin_width:.3f}" for b in bins[:-1]]

        # Assign each structure to a bin
        plot_df['distance_bin'] = pd.cut(plot_df['Plane_Distance_Difference'],
                                         bins=bins, labels=bin_labels, include_lowest=True)

        # Remove structures that fall outside the bin range
        plot_df = plot_df[plot_df['distance_bin'].notna()]

        if len(plot_df) == 0:
            if self.verbose:
                print("No structures fall within the specified distance range")
            return

        # Get unique charge-multiplicity combinations
        unique_combos = sorted([combo for combo in plot_df['charge_mult_combo'].unique() if combo != "Unknown"])

        if len(unique_combos) == 0:
            if self.verbose:
                print("No valid charge-multiplicity combinations found")
            return

        # Define color scheme for charge-multiplicity combinations
        color_map = {
            'q=0,m=1': '#0072B2',  # Blue
            'q=0,m=5': '#009E73',  # Green
            'q=1,m=2': '#D55E00',  # Orange
            'q=1,m=6': '#CC79A7'   # Purple
        }

        # Add fallback colors for any unexpected combinations
        import matplotlib.cm as cm
        extra_colors = cm.tab10(np.linspace(0, 1, 10))
        for i, combo in enumerate(unique_combos):
            if combo not in color_map:
                color_map[combo] = extra_colors[i % 10]

        # Create data structure for stacked bars
        stacked_data = {}
        bin_totals = {}

        for bin_label in bin_labels:
            stacked_data[bin_label] = {}
            subset = plot_df[plot_df['distance_bin'] == bin_label]
            bin_totals[bin_label] = len(subset)

            for combo in unique_combos:
                combo_count = len(subset[subset['charge_mult_combo'] == combo])
                stacked_data[bin_label][combo] = combo_count

        # Filter out empty bins
        non_empty_bins = [bin_label for bin_label in bin_labels if bin_totals.get(bin_label, 0) > 0]

        if len(non_empty_bins) == 0:
            if self.verbose:
                print("No non-empty bins found")
            return

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 7))

        # Create numeric positions for bars
        x_pos = np.arange(len(non_empty_bins))

        # Create stacked bars
        bottom = np.zeros(len(non_empty_bins))
        for combo in unique_combos:
            values = [stacked_data[bin_label][combo] for bin_label in non_empty_bins]

            if sum(values) > 0:  # Only plot if there's data
                bars = ax.bar(x_pos, values, bottom=bottom, width=0.7,
                             color=color_map[combo], label=combo, alpha=0.8, edgecolor='black', linewidth=0.5)
                bottom += np.array(values)

        # Set x-axis labels
        ax.set_xticks(x_pos)
        ax.set_xticklabels(non_empty_bins, rotation=45, ha='right', fontsize=10)

        # Add total count labels on top of bars
        for i, bin_label in enumerate(non_empty_bins):
            total = bin_totals[bin_label]
            ax.annotate(f'{int(total)}',
                       xy=(x_pos[i], total),
                       xytext=(0, 3), textcoords="offset points",
                       ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Formatting
        ax.set_xlabel('Fe Porphyrin Plane Distance Difference (Å)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Count', fontsize=14, fontweight='bold')
        ax.set_title('Fe Displacement from Porphyrin Plane by Charge-Multiplicity State',
                    fontsize=16, fontweight='bold', pad=20)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(loc='upper right', fontsize=12, title='Charge-Multiplicity', title_fontsize=13, framealpha=0.9)
        ax.grid(True, alpha=0.3, axis='y')

        # Set ylim to leave space at the top
        max_val = max(bin_totals.values()) if bin_totals else 1
        ax.set_ylim(0, max_val * 1.2)

        # Print statistics if verbose
        if self.verbose:
            print(f"\nFe Displacement Barplot Statistics:")
            print(f"Total structures: {len(plot_df)}")
            print(f"Number of bins with data: {len(non_empty_bins)}")
            print(f"Bin width: {bin_width} Å")
            print(f"\nBreakdown by charge-multiplicity combination:")
            for combo in unique_combos:
                combo_count = len(plot_df[plot_df['charge_mult_combo'] == combo])
                print(f"  {combo}: {combo_count} structures ({combo_count/len(plot_df)*100:.1f}%)")

        # Save plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.show()

        if self.verbose:
            print(f"\nFe displacement barplot saved as: {output_file}")

    def create_iron_displacement_normalized_by_group(self, df, output_file="plots/iron_displacement_normalized_by_group.png",
                                                      charge_mult_data=None, bin_width=0.01, max_distance=0.12):
        """
        Create normalized bar plot showing what fraction of each charge-multiplicity group
        falls into each displacement bin. This shows the distribution pattern for each group.

        Unlike the stacked barplot which shows absolute counts, this shows relative percentages
        within each charge-multiplicity state, helping identify if certain states have different
        displacement patterns.

        Args:
            df (pandas.DataFrame): DataFrame with Fe-plane distance data (must have 'Plane_Distance_Difference' column)
            output_file (str): Output filename for the plot
            charge_mult_data (pandas.DataFrame): DataFrame with PDB_ID, charge, and multiplicity columns
            bin_width (float): Width of each distance bin in Angstroms (default: 0.01)
            max_distance (float): Maximum distance to plot in Angstroms (default: 0.12)
        """
        if df.empty or 'Plane_Distance_Difference' not in df.columns:
            if self.verbose:
                print("No Fe-plane distance data to plot!")
            return

        # Prepare working dataframe
        plot_df = df.copy()

        # Merge with charge-multiplicity data if provided
        if charge_mult_data is not None:
            if 'PDB_ID' in plot_df.columns and all(col in charge_mult_data.columns for col in ['PDB_ID', 'charge', 'multiplicity']):
                # Merge the datasets
                plot_df = plot_df.merge(charge_mult_data[['PDB_ID', 'charge', 'multiplicity']],
                                       on='PDB_ID', how='left')

                # Create charge-multiplicity combination labels
                plot_df['charge_mult_combo'] = plot_df.apply(
                    lambda row: f"q={int(row['charge'])},m={int(row['multiplicity'])}"
                    if pd.notna(row['charge']) and pd.notna(row['multiplicity']) else "Unknown", axis=1
                )

                # Check if we have valid charge-multiplicity data
                valid_combo_mask = plot_df['charge_mult_combo'] != "Unknown"
                if valid_combo_mask.sum() > 0:
                    has_charge_mult = True
                    if self.verbose:
                        print(f"Successfully merged charge-multiplicity data for {valid_combo_mask.sum()}/{len(plot_df)} structures")
                else:
                    has_charge_mult = False
                    if self.verbose:
                        print("Warning: No valid charge-multiplicity combinations found after merge")
            else:
                has_charge_mult = False
                if self.verbose:
                    print("Warning: Cannot merge charge-multiplicity data - missing required columns")
        else:
            has_charge_mult = False
            if self.verbose:
                print("No charge-multiplicity data provided")

        if not has_charge_mult:
            if self.verbose:
                print("Cannot create normalized plot without valid charge-multiplicity data")
            return

        # Create distance bins
        bins = np.arange(0, max_distance + bin_width, bin_width)
        bin_labels = [f"{b:.3f}-{b+bin_width:.3f}" for b in bins[:-1]]

        # Assign each structure to a bin
        plot_df['distance_bin'] = pd.cut(plot_df['Plane_Distance_Difference'],
                                         bins=bins, labels=bin_labels, include_lowest=True)

        # Remove structures that fall outside the bin range
        plot_df = plot_df[plot_df['distance_bin'].notna()]

        if len(plot_df) == 0:
            if self.verbose:
                print("No structures fall within the specified distance range")
            return

        # Get unique charge-multiplicity combinations
        unique_combos = sorted([combo for combo in plot_df['charge_mult_combo'].unique() if combo != "Unknown"])

        if len(unique_combos) == 0:
            if self.verbose:
                print("No valid charge-multiplicity combinations found")
            return

        # Define color scheme for charge-multiplicity combinations (same as before)
        color_map = {
            'q=0,m=1': '#0072B2',  # Blue
            'q=0,m=5': '#009E73',  # Green
            'q=1,m=2': '#D55E00',  # Orange
            'q=1,m=6': '#CC79A7'   # Purple
        }

        # Add fallback colors for any unexpected combinations
        import matplotlib.cm as cm
        extra_colors = cm.tab10(np.linspace(0, 1, 10))
        for i, combo in enumerate(unique_combos):
            if combo not in color_map:
                color_map[combo] = extra_colors[i % 10]

        # Calculate normalized distributions for each charge-multiplicity group
        group_distributions = {}
        group_totals = {}

        for combo in unique_combos:
            combo_data = plot_df[plot_df['charge_mult_combo'] == combo]
            group_totals[combo] = len(combo_data)
            group_distributions[combo] = {}

            for bin_label in bin_labels:
                count = len(combo_data[combo_data['distance_bin'] == bin_label])
                # Normalize by total count in this group to get percentage
                percentage = (count / group_totals[combo] * 100) if group_totals[combo] > 0 else 0
                group_distributions[combo][bin_label] = percentage

        # Filter out empty bins (bins where no group has any data)
        non_empty_bins = [bin_label for bin_label in bin_labels
                         if any(group_distributions[combo][bin_label] > 0 for combo in unique_combos)]

        if len(non_empty_bins) == 0:
            if self.verbose:
                print("No non-empty bins found")
            return

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 7))

        # Create grouped bars (side-by-side for each charge-multiplicity state)
        x_pos = np.arange(len(non_empty_bins))
        bar_width = 0.8 / len(unique_combos)

        for i, combo in enumerate(unique_combos):
            values = [group_distributions[combo][bin_label] for bin_label in non_empty_bins]
            offset = (i - len(unique_combos)/2 + 0.5) * bar_width

            bars = ax.bar(x_pos + offset, values, bar_width,
                         color=color_map[combo], label=f"{combo} (n={group_totals[combo]})",
                         alpha=0.8, edgecolor='black', linewidth=0.5)

            # Add percentage labels on bars (only if > 5% to avoid clutter)
            for j, (bar, val) in enumerate(zip(bars, values)):
                if val > 5:  # Only label if > 5%
                    height = bar.get_height()
                    ax.annotate(f'{val:.1f}%',
                               xy=(bar.get_x() + bar.get_width()/2, height),
                               xytext=(0, 2), textcoords="offset points",
                               ha='center', va='bottom', fontsize=8, rotation=0)

        # Set x-axis labels
        ax.set_xticks(x_pos)
        ax.set_xticklabels(non_empty_bins, rotation=45, ha='right', fontsize=10)

        # Formatting
        ax.set_xlabel('Fe Porphyrin Plane Distance Difference (Å)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Percentage of Group (%)', fontsize=14, fontweight='bold')
        ax.set_title('Normalized Distribution: Fe Displacement by Charge-Multiplicity State',
                    fontsize=16, fontweight='bold', pad=20)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(loc='upper right', fontsize=11, title='Charge-Multiplicity (total count)',
                 title_fontsize=11, framealpha=0.9)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_ylim(0, 100)  # Percentage scale

        # Print statistics if verbose
        if self.verbose:
            print(f"\nNormalized Distribution Statistics:")
            print(f"Total structures: {len(plot_df)}")
            print(f"Number of bins with data: {len(non_empty_bins)}")
            print(f"Bin width: {bin_width} Å")
            print(f"\nDistribution by charge-multiplicity combination:")
            for combo in unique_combos:
                print(f"\n{combo} ({group_totals[combo]} structures):")
                # Show top 3 bins for this group
                sorted_bins = sorted([(bin_label, group_distributions[combo][bin_label])
                                     for bin_label in non_empty_bins],
                                    key=lambda x: x[1], reverse=True)
                for bin_label, pct in sorted_bins[:3]:
                    if pct > 0:
                        print(f"  {bin_label} Å: {pct:.1f}%")

        # Save plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.show()

        if self.verbose:
            print(f"\nNormalized displacement plot saved as: {output_file}")

    def run_iron_plane_analysis(self, charge_mult_data=None):
        """Run Fe-heme plane distance analysis."""
        df_plane = self.process_iron_plane_distances()
        if not df_plane.empty:
            # Save all Fe-plane distance results to CSV
            plane_csv = str(derived_table_path("iron_plane_distances.csv"))
            df_plane.to_csv(plane_csv, index=False)
            if self.verbose:
                print(f"\nAll Fe-plane distance results saved to: {plane_csv}")
            
            # Create filtered version for plots based on processed_output.csv PDB IDs
            # Note: processed_output.csv contains axial ligands sorted alphabetically (axial1 < axial2)
            try:
                processed_df = pd.read_csv(resolve_table_input("processed_output.csv"))
                if 'file_name' in processed_df.columns:
                    # Extract PDB codes from file names (first 4 characters of log files)
                    processed_pdb_ids = set()
                    for file_name in processed_df['file_name'].unique():
                        if isinstance(file_name, str) and file_name.endswith('.log'):
                            pdb_code = file_name[:4]
                            processed_pdb_ids.add(pdb_code)
                    
                    # Filter Fe-plane data to only include PDB IDs from processed data
                    df_plane_plots = df_plane[df_plane['PDB_ID'].isin(processed_pdb_ids)].copy()
                    
                    # Save filtered data for plots
                    plots_csv = str(derived_table_path("iron_plane_distances_plots.csv"))
                    df_plane_plots.to_csv(plots_csv, index=False)
                    
                    if self.verbose:
                        print(f"Filtered Fe-plane data for plots saved to: {plots_csv}")
                        print(f"Filtered data contains {len(df_plane_plots)} structures (from {len(df_plane)} total)")
                        print(f"PDB IDs in processed data: {len(processed_pdb_ids)}")
                        print(f"PDB IDs with Fe-plane data: {len(df_plane['PDB_ID'].unique())}")
                        print(f"Common PDB IDs: {len(df_plane_plots['PDB_ID'].unique())}")
                    
                    # Use filtered data for plotting analysis
                    plot_data = df_plane_plots
                else:
                    if self.verbose:
                        print("Warning: No 'file_name' column found in processed_output.csv, using all data")
                    plot_data = df_plane
            except FileNotFoundError:
                if self.verbose:
                    print("Warning: processed_output.csv not found, using all Fe-plane data for plots")
                plot_data = df_plane
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Error reading processed_output.csv ({e}), using all data")
                plot_data = df_plane
            
            # Create histogram with charge-multiplicity breakdown if data is provided
            #self.create_iron_plane_histogram(plot_data, charge_mult_data=charge_mult_data)
            
            # Create comprehensive Fe-plane analysis plot
            # Try to get axial ligand data for color coding - only process PDB IDs from filtered data
            valid_pdb_ids = plot_data['PDB_ID'].unique()
            axial_df_filtered = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
            
            if not axial_df_filtered.empty:
                if self.verbose:
                    print(f"Axial data for filtered PDB IDs: {len(axial_df_filtered)} measurements")
                unique_combos = len(axial_df_filtered.groupby('PDB_ID')['Axial_Resname'].apply(lambda x: '-'.join(sorted(x.astype(str)))).unique())
                if self.verbose:
                    print(f"Unique axial combinations: {unique_combos}")
                
                #self.create_comprehensive_iron_plane_analysis(plot_data, 
                #                                            color_by_axial_ligands=True,
                #                                            axial_data=axial_df_filtered)
            else:
                if self.verbose:
                    print("No axial data found for filtered PDB IDs")
                #self.create_comprehensive_iron_plane_analysis(plot_data)
            if self.verbose:
                print(f"Fe-plane analysis: Processed {len(df_plane)} structure pairs successfully.")
            return df_plane
        else:
            if self.verbose:
                print("No valid structure pairs found for Fe-plane analysis!")
            return pd.DataFrame()
    
    def create_absolute_iron_plane_distances_histogram(self, df, output_file="plots/fe-por_dist.png"):
        """
        Create histograms of absolute Fe-heme plane distances,
        where each bar is colored by axial ligand combinations for the most common 5 combinations.
        Creates separate plots for PDB and XYZ distances.
        """
        if df.empty or 'PDB_Iron_Plane_Distance' not in df.columns or 'XYZ_Iron_Plane_Distance' not in df.columns:
            print("No Fe-plane distance data to plot!")
            return

        # Get axial ligand data for coloring
        plot_df = df.copy()
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        if 'PDB_ID' in plot_df.columns:
            # Get axial ligand information
            valid_pdb_ids = plot_df['PDB_ID'].unique()
            axial_df = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
            
            if not axial_df.empty:
                # Create axial combinations by pivoting data
                axial_pivot = axial_df.pivot_table(
                    index='PDB_ID', 
                    columns='Axial_Ligand', 
                    values='Axial_Resname', 
                    aggfunc='first'
                ).reset_index()
                
                # Create proper axial combinations (axial1-axial2) with HSD normalization
                if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                    # Normalize ligand names (HSD -> HIS) and create sorted combination
                    # Note: Axial ligands from processed_output.csv are already sorted,
                    # but we sort here for consistency and to handle data from other sources
                    def create_normalized_combo(row):
                        if pd.notna(row[1]) and pd.notna(row[2]):
                            norm_axial1 = self.normalize_axial_ligand_name(row[1])
                            norm_axial2 = self.normalize_axial_ligand_name(row[2])
                            return '-'.join(sorted([norm_axial1, norm_axial2]))
                        else:
                            return 'Unknown-Unknown'
                    axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                elif 1 in axial_pivot.columns:
                    norm_axial1 = axial_pivot[1].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = norm_axial1.astype(str) + '-Unknown'
                elif 2 in axial_pivot.columns:
                    norm_axial2 = axial_pivot[2].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = 'Unknown-' + norm_axial2.astype(str)
                else:
                    axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
                
                axial_ligand_info = axial_pivot[['PDB_ID', 'Axial_Combo']]
                
                # Merge with plot data
                plot_df = plot_df.merge(axial_ligand_info, on='PDB_ID', how='left')
                
                # Check if we have valid axial combination data
                if 'Axial_Combo' in plot_df.columns and plot_df['Axial_Combo'].notna().sum() > 0:
                    use_axial_colors = True
                    # Only use the 5 most common combinations (as defined in allowed_axial_combinations)
                    allowed_combos = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
                    axial_combinations = [combo for combo in allowed_combos 
                                        if combo in plot_df['Axial_Combo'].values]
                    
                    # Use centralized color mapping for axial combinations
                    axial_color_map = get_axial_ligand_colors(axial_combinations)
                    
                    if self.verbose:
                        print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
                        for combo in axial_combinations:
                            count = (plot_df['Axial_Combo'] == combo).sum()
                            print(f"  {combo}: {count} structures")

        # Create figure with subplots for PDB, XYZ, and comprehensive analysis
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6))
        
        # Function to create histogram for given data type
        def create_single_histogram(ax, distance_col, title, xlim_max=0.35):
            distances = plot_df[distance_col].dropna()
            
            if use_axial_colors and axial_color_map is not None and axial_combinations is not None:
                # Create stacked histogram by axial combination
                combo_data = {}
                for combo in axial_combinations:
                    combo_mask = plot_df['Axial_Combo'] == combo
                    combo_values = plot_df.loc[combo_mask, distance_col].dropna()
                    combo_data[combo] = combo_values.values
                
                # Create stacked histogram - filter out empty combinations
                valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
                data_arrays = [combo_data[combo] for combo in valid_combos]
                colors = [axial_color_map[combo] for combo in valid_combos]
                labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
                
                # Create stacked histogram
                ax.hist(data_arrays, bins=np.linspace(0, xlim_max, 51), color=colors, alpha=0.8, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add single legend only to the first subplot (leftmost)
                if ax == ax1:
                    # Create legend without box/frame for all combinations
                    from matplotlib.patches import Patch
                    legend_elements = [Patch(facecolor=axial_color_map[combo], label=combo) 
                                     for combo in valid_combos]
                    ax.legend(handles=legend_elements, title="Axial Ligand Combinations", 
                             fontsize=10, title_fontsize=12, loc='upper right',
                             frameon=False, fancybox=False, shadow=False,
                             bbox_to_anchor=(1.0, 1.0))
                
            else:
                # Fall back to simple histogram
                ax.hist(distances, bins=np.linspace(0, xlim_max, 51), alpha=0.7, color='lightblue', edgecolor='black')
            
            # Formatting
            ax.set_xlabel('Fe-Porphyrin Plane Distance (Å)', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlim(0, xlim_max)
            #ax.set_title(title, fontsize=16)
            ax.grid(True, alpha=0.3)
            
            # Print statistics
            if self.verbose:
                print(f"\n{title} Statistics:")
                print(f"n={len(distances)}, mean={np.mean(distances):.3f}, std={np.std(distances):.3f}")
        
        # Create histograms
        create_single_histogram(ax1, 'PDB_Iron_Plane_Distance', 'PDB (Crystal) Fe-Plane Distances')
        create_single_histogram(ax2, 'XYZ_Iron_Plane_Distance', 'XYZ (Optimized) Fe-Plane Distances')
        
        # Create third subplot with comprehensive analysis (distance differences)
        if 'Plane_Distance_Difference' in plot_df.columns:
            create_single_histogram(ax3, 'Plane_Distance_Difference', 'Fe-Plane Distance Differences', xlim_max=0.12)
        else:
            ax3.text(0.5, 0.5, 'No distance difference\ndata available', 
                    ha='center', va='center', transform=ax3.transAxes, fontsize=14)
            ax3.set_xlabel('Fe-Porphyrin Plane Distance Difference (Å)', fontsize=14)
            ax3.set_ylabel('Frequency', fontsize=14)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        if self.verbose:
            print(f"\nThree-panel Fe-porphyrin distance histogram saved as: {output_file}")


    def create_comprehensive_iron_plane_analysis(
        self, df, 
        output_file="plots/comprehensive_iron_plane_analysis.png",
        color_by_axial_ligands=False,
        axial_data=None,
        charge_mult_data=None,
        **kwargs
    ):
        """
        Create two-subplot analysis:
        1. Fe-plane distance differences colored by axial ligand combinations
        2. Fe-plane distance differences with charge-multiplicity breakdown
        
        Args:
            df (pandas.DataFrame): DataFrame with Fe-plane distance data
            output_file (str): Output filename for the plot
            color_by_axial_ligands (bool): Whether to color bars by axial ligand combinations
            axial_data (pandas.DataFrame): DataFrame with PDB_ID and axial ligand information
            charge_mult_data (pandas.DataFrame): DataFrame with charge-multiplicity information
        """
        if df.empty or 'Plane_Distance_Difference' not in df.columns:
            if self.verbose:
                print("No Fe-plane distance difference data to plot!")
            return
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Extract distance difference data
        plot_df = df.copy()
        differences = plot_df['Plane_Distance_Difference'].dropna()
        
        if differences.empty:
            if self.verbose:
                print("No valid distance difference data to plot!")
            return
        
        # SUBPLOT 1: Axial ligand colored histogram
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        if color_by_axial_ligands and axial_data is not None:
            # Try to merge with axial ligand data
            if 'PDB_ID' in plot_df.columns:
                # Get axial ligand combinations from axial_data
                # Create proper combinations by grouping axial ligands 1 and 2 separately
                axial_pivot = axial_data.pivot_table(
                    index='PDB_ID', 
                    columns='Axial_Ligand', 
                    values='Axial_Resname', 
                    aggfunc='first'
                ).reset_index()
                
                # Create proper axial combinations (axial1-axial2) with HSD normalization
                if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                    # Normalize ligand names (HSD -> HIS) and create sorted combination
                    # Note: Axial ligands from processed_output.csv are already sorted,
                    # but we sort here for consistency and to handle data from other sources
                    def create_normalized_combo(row):
                        if pd.notna(row[1]) and pd.notna(row[2]):
                            norm_axial1 = self.normalize_axial_ligand_name(row[1])
                            norm_axial2 = self.normalize_axial_ligand_name(row[2])
                            return '-'.join(sorted([norm_axial1, norm_axial2]))
                        else:
                            return 'Unknown-Unknown'
                    axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                elif 1 in axial_pivot.columns:
                    norm_axial1 = axial_pivot[1].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = norm_axial1.astype(str) + '-Unknown'
                elif 2 in axial_pivot.columns:
                    norm_axial2 = axial_pivot[2].apply(self.normalize_axial_ligand_name)
                    axial_pivot['Axial_Combo'] = 'Unknown-' + norm_axial2.astype(str)
                else:
                    axial_pivot['Axial_Combo'] = 'Unknown-Unknown'
                
                axial_ligand_info = axial_pivot[['PDB_ID', 'Axial_Combo']]
                
                # Merge with plot data
                plot_df = plot_df.merge(axial_ligand_info, on='PDB_ID', how='left')
                
                # Check if we have valid axial combination data
                if 'Axial_Combo' in plot_df.columns and plot_df['Axial_Combo'].notna().sum() > 0:
                    use_axial_colors = True
                    axial_combinations = sorted(plot_df['Axial_Combo'].dropna().unique())
                    
                    # Define color map for different axial combinations
                    colors_palette = ['#9467bd', '#0072B2', '#009E73', '#D55E00', '#F0E442',
                                    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
                    axial_color_map = {combo: colors_palette[i % len(colors_palette)] 
                                     for i, combo in enumerate(axial_combinations)}
                    
                    if self.verbose:
                        print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
                        for combo in axial_combinations:
                            count = (plot_df['Axial_Combo'] == combo).sum()
                            print(f"  {combo}: {count} structures")
                else:
                    if self.verbose:
                        print("Warning: Could not merge axial ligand data - using default colors")
            else:
                if self.verbose:
                    print("Warning: No PDB_ID column found - cannot use axial ligand coloring")
        
        # Function to create individual subplot
        def create_subplot(ax, plot_data, title, color_scheme='axial'):
            if color_scheme == 'axial' and use_axial_colors and axial_color_map is not None and axial_combinations is not None:
                # Create stacked histogram by axial combination
                combo_data = {}
                for combo in axial_combinations:
                    combo_mask = plot_data['Axial_Combo'] == combo
                    combo_values = plot_data.loc[combo_mask, 'Plane_Distance_Difference'].dropna()
                    combo_data[combo] = combo_values.values
                
                # Create stacked histogram - filter out empty combinations
                valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
                data_arrays = [combo_data[combo] for combo in valid_combos]
                colors = [axial_color_map[combo] for combo in valid_combos]
                labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
                
                # Define bins based on all data
                bins = np.linspace(0, 0.12, 51)
                
                # Create stacked histogram
                ax.hist(data_arrays, bins=bins, color=colors, alpha=0.8, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add legend without frame
                ax.legend(title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                         frameon=False, fancybox=False, shadow=False)
                
            else:
                # Fall back to simple histogram
                ax.hist(differences, bins=np.linspace(0, 0.12, 51), alpha=0.7, color='lightgreen', edgecolor='black')
            
            # Formatting
            ax.set_xlabel('Fe Porphyrin Plane Distance Difference (Å)', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlim(0, 0.12)
            ax.grid(True, alpha=0.3)
        
        # Create first subplot with axial ligand coloring
        create_subplot(ax1, plot_df, 'Axial Ligand Breakdown', 'axial')
        
        # Add a) annotation to first subplot
        ax1.text(0.05, 0.95, 'a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')
        
        # SUBPLOT 2: Create second subplot with same axial ligand coloring as the first
        create_subplot(ax2, plot_df, 'Axial Ligand Breakdown', 'axial')
        
        # Add b) annotation to second subplot
        ax2.text(0.05, 0.95, 'b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', horizontalalignment='left')
        
        # Calculate and display statistics
        mean_diff = differences.mean()
        std_diff = differences.std()
        
        if self.verbose:
            print(f"\nFe Porphyrin Plane Distance Difference Statistics:")
            print(f"Mean: {mean_diff:.3f} ± {std_diff:.3f} Å")
            print(f"Range: {differences.min():.3f} - {differences.max():.3f} Å")
            print(f"Number of structures: {len(differences)}")
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.show()
        
        # Print summary
        if self.verbose:
            print(f"\nComprehensive Fe-plane analysis saved as: {output_file}")
            if use_axial_colors and axial_combinations is not None:
                print(f"Colored by {len(axial_combinations)} axial ligand combinations")
            print(f"Analyzed {len(differences)} structures")

    
    def run_absolute_iron_plane_analysis(self):
        """Run absolute Fe-heme plane distance analysis."""        
        df_plane = self.process_iron_plane_distances()
        if not df_plane.empty:
            # Create absolute distance histogram
            self.create_absolute_iron_plane_distances_histogram(df_plane)
            if self.verbose:
                print(f"Absolute Fe-plane analysis: Processed {len(df_plane)} structure pairs successfully.")
            return df_plane
        else:
            if self.verbose:
                print("No valid structure pairs found for absolute Fe-plane analysis!")
            return pd.DataFrame()
    
    def run_all_analyses(self, rmsd=True, axial=True, iron_plane=True, absolute_iron_plane=True, distortion=False, charge_mult_data=None):
        """
        Run all analyses based on specified flags.
        
        Args:
            rmsd (bool): Whether to run RMSD analysis
            axial (bool): Whether to run Fe-axial distance analysis
            iron_plane (bool): Whether to run Fe-heme plane distance analysis
            absolute_iron_plane (bool): Whether to run absolute Fe-heme plane distance analysis
            distortion (bool): Whether to run distortion analysis
            charge_mult_data (pandas.DataFrame): DataFrame with PDB_ID, charge, and multiplicity columns for enhanced plotting
            
        Returns:
            dict: Dictionary containing analysis results
        """
        if self.verbose:
            print("Starting comprehensive analysis of heme structures...")
            print("=" * 60)
        
        # Initialize log file
        self.initialize_log()
        
        results = {}
        
        if rmsd:
            results['rmsd'] = self.run_rmsd_analysis()
        
        if axial:
            results['axial_distances'] = self.run_iron_axial_distance_analysis()
        
        if iron_plane:
            results['iron_plane_distances'] = self.run_iron_plane_analysis(charge_mult_data=charge_mult_data)
        
        if absolute_iron_plane:
            results['absolute_iron_plane_distances'] = self.run_absolute_iron_plane_analysis()
        
        if distortion:
            print("Distortion analysis moved to plots.py - use plots.plots_class().create_distortion_modes_histogram()")
        if self.verbose:
            # Final summary
            print("\n" + "=" * 60)
            print("\nAnalysis Complete!")
            print("Generated files:")
            if rmsd:
                print("- rmsd_results.csv")
                print("- rmsd_histogram.png")
            if axial:
                print("- iron_axial_distances.csv") 
                print("- iron_axial_distance_differences.png")
                print("- iron_axial_distance_differences_combined.png")
            if iron_plane:
                print("- iron_plane_distances.csv")
                print("- iron_plane_distances.png")
            if absolute_iron_plane:
                print("- absolute_iron_plane_distances.png")
            if distortion:
                print("- distortion_comparison_histograms.png")
                print("- distortion_comparison_summary.csv")
            print("- suspicious_distances.log (if any distance differences > 1.5 Å found)")
        return results


def main():
    """Main function to run the RMSD analysis and additional analyses."""
    parser = argparse.ArgumentParser(
        description="Run RMSD, Fe-axial, and Fe-plane analyses on prepared heme structures."
    )
    parser.add_argument(
        "--pdb-base-dir",
        default=str(PDB_DIR),
        help="Base directory containing prepared PDB/<pdb_id>/ structure folders.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    parser.add_argument(
        "--include-all-axial-ligands",
        action="store_true",
        help="Disable filtering to the common axial ligand combinations.",
    )
    parser.add_argument("--skip-rmsd", action="store_true", help="Skip RMSD analysis.")
    parser.add_argument("--skip-axial", action="store_true", help="Skip Fe-axial distance analysis.")
    parser.add_argument("--skip-iron-plane", action="store_true", help="Skip Fe-plane distance analysis.")
    parser.add_argument(
        "--skip-absolute-iron-plane",
        action="store_true",
        help="Skip absolute Fe-plane distance analysis.",
    )
    parser.add_argument(
        "--distortion",
        action="store_true",
        help="Request distortion analysis messaging.",
    )
    args = parser.parse_args()

    analyzer = RMSDAnalyzer(
        pdb_base_dir=args.pdb_base_dir,
        verbose=args.verbose,
        filter_axial_ligands=not args.include_all_axial_ligands,
    )
    analyzer.run_all_analyses(
        rmsd=not args.skip_rmsd,
        axial=not args.skip_axial,
        iron_plane=not args.skip_iron_plane,
        absolute_iron_plane=not args.skip_absolute_iron_plane,
        distortion=args.distortion,
    )


if __name__ == "__main__":
    main()
