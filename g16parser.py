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
import argparse
from joblib import Parallel, delayed
from pathlib import Path
from repo_paths import LOGFILES_DIR, JSONS_DIR, canonical_table_path


#Global Flags
debug=False
verbose=False


class g16parser:
    def __init__(self, log_dir=None, json_dir=None, csv_path=None):
        self.log_dir = str(Path(log_dir) if log_dir is not None else LOGFILES_DIR)
        self.json_dir = str(Path(json_dir) if json_dir is not None else JSONS_DIR)
        self.parsed_data = {}
        if csv_path is not None:
            self.csv_data = pd.read_csv(csv_path)
        else:
            self.csv_data = pd.read_csv(canonical_table_path("pyDISH.csv"))
        os.makedirs(self.json_dir, exist_ok=True)

    def removekey_copy(self, dictionary, key):
        r = dict(dictionary)
        del r[key]
        return r

    def removekey(self, key):
        del self.parsed_data[key]

    def parse(self):
        """
        Iterates over all .log files in self.log_dir and calls parse_gaussian_logfile()
        for each, creating a similarly-named .json file in self.json_dir.
        """
        for filename in os.listdir(self.log_dir):
            if filename.lower().endswith('.log'):
                log_path = os.path.join(self.log_dir, filename)
                json_filename = os.path.splitext(filename)[0] + '.json'
                json_path = os.path.join(self.json_dir, json_filename)
                self.parse_gaussian_logfile(log_path, json_path)
        print("Parsing complete. Looking for empty files.")
        self.delete_json_files_without_homo_lumo()
        print("Empty files deleted.")
        return

    def parse_gaussian_logfile(self, log_filepath, output_json_path):
        """
        Reads the Gaussian logfile, extracts information, and writes it to a JSON file.
        """
        if not os.path.isfile(log_filepath):
            raise FileNotFoundError(f"Could not find file: {log_filepath}")
        with open(log_filepath, 'r') as infile:
            log_contents = infile.read()
        with open(log_filepath, 'r') as infile:
            log_lines = infile.readlines()
        self.parsed_data = {
            "file_name": os.path.basename(log_filepath),
            "file_path": log_filepath,
            "warnings": [],
            "energies": {},
            "final_structure": None,
            "termination_status": None,
        }
        # existing extractors
        self.extract_job_info(log_contents)
        self.extract_energy_info(log_contents)
        orbital_energies = self.extract_orbital_energies(log_contents)
        self.parsed_data["energies"] = orbital_energies
        self.extract_homo_lumo()
        self.extract_standard_orientation(log_lines)
        self.extract_geometry(log_lines)
        self.extract_thermodynamics(log_contents)
        self.extract_spectroscopy(log_contents)
        self.extract_population_analysis(log_contents)
        self.extract_warnings(log_contents)
        self.extract_sum_mulliken_charges(log_contents)
        self.extract_mulliken_charges_heavy_atoms(log_contents)
        # new extractors
        self.extract_mulliken_charges_atomic(log_contents)
        self.extract_mulliken_spin_densities(log_contents)
        self.extract_mulliken_spin_densities_heavy_atoms(log_contents)
        self.extract_spin_contamination_analysis(log_contents)
        self.extract_bond_analysis_from_nbo(log_contents)
        self.extract_nao_nlmo_analysis(log_contents)
        self.extract_second_order_perturbation_analysis(log_contents)
        self.extract_comprehensive_orbital_energies(log_contents)
        self.extract_dipole_moment(log_contents)
        self.extract_dipole_polarizability(log_contents)
        self.extract_quadrupole_moment(log_contents)
        self.extract_atomic_dipole_orientation(log_contents)
        self.extract_overview(log_contents)
        # continue existing flow
        self.extract_natural_population_analysis(log_contents)
        self.add_geometry_to_npa()
        self.extract_natural_populations(log_contents)
        self.extract_natural_electron_config(log_contents)
        self.extract_natural_charges(log_contents)
        # Run comprehensive atomic charges analysis after all individual extractions
        self.extract_multiple_atomic_charges(log_contents)
        if "optimized_geometry" in self.parsed_data:
            self.removekey(key="optimized_geometry")
        self.add_csv_information(os.path.basename(log_filepath))
        with open(output_json_path, 'w') as outfile:
            json.dump(self.parsed_data, outfile, indent=2)
        if 'debug' in globals() and debug:
            print(f"Parsed data from {log_filepath} was saved to {output_json_path}")
        return

    def extract_job_info(self, log_contents):
        """Extracts general job information."""
        job_info_pattern = re.compile(r"\\#.*\\n(.*?)\\n\\n", re.DOTALL)
        job_info = job_info_pattern.search(log_contents)
        self.parsed_data["job_info"] = job_info.group(1).strip() if job_info else None

    def extract_energy_info(self, log_contents):
        """Extracts total energy information."""
        energy_pattern = re.compile(r"SCF Done:  E\(\w+\) =\s+([-\d\.]+)")
        energy_match = energy_pattern.findall(log_contents)
        self.parsed_data["total_energies"] = [float(e) for e in energy_match] if energy_match else None

    def extract_orbital_energies(self, log_contents):
        """
        Extracts occupied and virtual orbital energies.
        """
        pattern = re.compile(
            r'^\s*(Alpha|Beta)\s+(occ|virt)\.\s+eigenvalues\s+--\s+(.*)$',
            re.MULTILINE
        )
        energies = {
            "alpha_occ":  [],
            "alpha_virt": [],
            "beta_occ":   [],
            "beta_virt":  []
        }
        current_key = None
        for match in pattern.finditer(log_contents):
            spin_type = match.group(1).lower()  # "alpha" or "beta"
            occ_virt = match.group(2).lower()   # "occ" or "virt"
            floats_str = match.group(3)
            float_values = [float(x) for x in floats_str.split()]
            new_key = f"{spin_type}_{occ_virt}"
            if new_key != current_key:
                energies[new_key] = []
                current_key = new_key
            energies[new_key].extend(float_values)
        return energies

    def extract_homo_lumo(self):
        """
        Extracts HOMO and LUMO values from orbital energies.
        """
        alpha_occ = self.parsed_data.get("energies", {}).get("alpha_occ", [])
        alpha_virt = self.parsed_data.get("energies", {}).get("alpha_virt", [])
        homo = alpha_occ[-10:]
        lumo = alpha_virt[:10]
        self.parsed_data["alpha_homo_lumo"] = {
            "homo": homo,
            "lumo": lumo
        }
        self.parsed_data.pop("energies", None)
        return

    def extract_geometry(self, log_lines):
        """Extracts the optimized geometry of the molecule, including center number,
        atomic number, atomic type, and the X, Y, Z coordinates."""
        geometry_pattern = re.compile(
            r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$'
        )
        geometry_list = []
        expected_center = 1
        for line in log_lines:
            match = geometry_pattern.match(line)
            if match:
                center_number = int(match.group(1))
                if center_number != expected_center:
                    break
                atomic_number = int(match.group(2))
                atomic_type = int(match.group(3))
                x = float(match.group(4))
                y = float(match.group(5))
                z = float(match.group(6))
                geometry_list.append({
                    "center": center_number,
                    "atomic_number": atomic_number,
                    "atomic_type": atomic_type,
                    "x": x,
                    "y": y,
                    "z": z
                })
                expected_center += 1
        if geometry_list:
            self.parsed_data["optimized_geometry"] = geometry_list

    def extract_standard_orientation(self, log_lines):
        """
        Extracts the final Standard orientation coordinates from the Gaussian logfile lines.
        Stores a list of atoms with their center, atomic number, atomic type, and x,y,z.
        """
        std_orientations = []
        capturing = False
        header_count = 0
        temp_list = []
        for i, line in enumerate(log_lines):
            if 'Standard orientation:' in line:
                capturing = True
                header_count = 0
                temp_list = []
                continue
            if capturing:
                # skip separator and header lines
                if header_count < 2:
                    if '----' in line:
                        header_count += 1
                    continue
                # stop at next separator line
                if '----' in line:
                    capturing = False
                    std_orientations = temp_list
                    break
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        center = int(parts[0])
                        atomic_number = int(parts[1])
                        atomic_type = int(parts[2])
                        x = float(parts[3])
                        y = float(parts[4])
                        z = float(parts[5])
                        temp_list.append({
                            "center": center,
                            "atomic_number": atomic_number,
                            "atomic_type": atomic_type,
                            "x": x,
                            "y": y,
                            "z": z
                        })
                    except ValueError:
                        continue
        if std_orientations:
            self.parsed_data["standard_orientation"] = std_orientations
        else:
            self.parsed_data["standard_orientation"] = None

    def extract_thermodynamics(self, log_contents):
        """Extracts thermodynamic properties like ZPE, enthalpy, and entropy."""
        zpe_pattern = re.compile(r"Zero-point correction=\s+([-\d\.]+)")
        enthalpy_pattern = re.compile(r"Thermal correction to Enthalpy=\s+([-\d\.]+)")
        entropy_pattern = re.compile(r"Thermal correction to Gibbs Free Energy=\s+([-\d\.]+)")
        zpe = zpe_pattern.search(log_contents)
        enthalpy = enthalpy_pattern.search(log_contents)
        entropy = entropy_pattern.search(log_contents)
        self.parsed_data["thermodynamics"] = {
            "zero_point_energy": float(zpe.group(1)) if zpe else None,
            "enthalpy": float(enthalpy.group(1)) if enthalpy else None,
            "entropy": float(entropy.group(1)) if entropy else None
        }

    def extract_spectroscopy(self, log_contents):
        """Extracts spectroscopic data such as IR frequencies and intensities."""
        frequency_pattern = re.compile(r"Frequencies --\s+([-\d\.\s]+)")
        intensity_pattern = re.compile(r"IR Inten\(KM/Mole\) --\s+([-\d\.\s]+)")
        frequencies = frequency_pattern.findall(log_contents)
        intensities = intensity_pattern.findall(log_contents)
        self.parsed_data["spectroscopy"] = {
            "frequencies": [float(f) for sublist in frequencies for f in sublist.split()],
            "intensities": [float(i) for sublist in intensities for i in sublist.split()]
        }

    def extract_population_analysis(self, log_contents):
        """Extracts Mulliken charges or other population analysis results."""
        mulliken_pattern = re.compile(r"Mulliken charges and spin densities:\s*\n(?:\s*\d+\s+\w+\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s*\n)+")
        mulliken_match = mulliken_pattern.search(log_contents)
        if mulliken_match:
            charges = [float(line.split()[1]) for line in mulliken_match.group(1).splitlines()]
            self.parsed_data["mulliken_charges"] = charges

    def extract_warnings(self, log_contents):
        """Extracts any warnings present in the log file."""
        pattern = re.compile(r'^\s*Warning:\s*(.*)$', re.MULTILINE)
        warnings = pattern.findall(log_contents)
        self.parsed_data["warnings"] = warnings if warnings else None

    def extract_sum_mulliken_charges(self, log_contents):
        """
        Extracts the sum of Mulliken charges and spin densities from a line like:
        "Sum of Mulliken charges =   1.00000   5.00000"
        """
        pattern = re.compile(r"Sum of Mulliken charges\s*=\s*([-\d.]+)\s+([-\d.]+)")
        match = pattern.search(log_contents)
        if match:
            charge_sum = float(match.group(1))
            spin_sum = float(match.group(2))
            self.parsed_data["sum_mulliken_charges"] = {"charge_sum": charge_sum, "spin_sum": spin_sum}
        else:
            self.parsed_data["sum_mulliken_charges"] = None

    def extract_mulliken_charges_heavy_atoms(self, log_contents):
        """
        Extracts the Mulliken charges and spin densities with hydrogens summed into heavy atoms.
        Expects a block that starts with:
        "Mulliken charges and spin densities with hydrogens summed into heavy atoms:"
        followed by rows such as:
             1  Fe   0.036481   4.255523
             2  N   -0.112201   0.076067
             ...
        """
        header_pattern = re.compile(
            r"Mulliken charges and spin densities with hydrogens summed into heavy atoms:\s*\n",
            re.IGNORECASE
        )
        header_match = header_pattern.search(log_contents)
        if not header_match:
            self.parsed_data["mulliken_heavy_atoms"] = None
            return
        start_index = header_match.end()
        lines = log_contents[start_index:].splitlines()
        row_pattern = re.compile(r"^\s*(\d+)\s+(\w+)\s+([-\d.]+)\s+([-\d.]+)")
        heavy_atoms = []
        for line in lines:
            m = row_pattern.match(line)
            if m:
                index = int(m.group(1))
                element = m.group(2)
                charge = float(m.group(3))
                spin = float(m.group(4))
                heavy_atoms.append({
                    "index": index,
                    "element": element,
                    "mulliken_charge": charge,
                    "spin_density": spin
                })
            else:
                if heavy_atoms:
                    break
        self.parsed_data["mulliken_heavy_atoms"] = heavy_atoms

    def extract_natural_population_analysis(self, log_contents):
        """
        Scans the logfile text for the "Summary of Natural Population Analysis" table
        and extracts each atom's data. If optimized geometry data (with center number,
        atomic number, atomic type, and X, Y, Z coordinates) is available, it is merged into
        each atom's information. The result is stored in self.parsed_data["natural_population"].
        """
        block_pattern = re.compile(
            r'(Summary of Natural Population Analysis:.*?)(?:={5,}|\Z)',
            re.DOTALL
        )
        block_match = block_pattern.search(log_contents)
        if not block_match:
            self.parsed_data["natural_population"] = None
            return
        npa_block = block_match.group(1)
        line_pattern = re.compile(
            r'^\s*([A-Z][a-z]?)\s+(\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s*$',
            re.MULTILINE
        )
        total_pattern = re.compile(
            r'^\s*\*\s*Total\s*\*\s*'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s+'
            r'([+-]?\d+\.\d+)\s*$',
            re.MULTILINE
        )
        atoms_data = []
        for match in line_pattern.finditer(npa_block):
            atom_symbol = match.group(1)
            atom_no = int(match.group(2))
            charge = float(match.group(3))
            core = float(match.group(4))
            valence = float(match.group(5))
            rydberg = float(match.group(6))
            total = float(match.group(7))
            atoms_data.append({
                "atom_symbol": atom_symbol,
                "atom_number": atom_no,
                "charge": charge,
                "core": core,
                "valence": valence,
                "rydberg": rydberg,
                "total": total
            })
        optimized_geometry = self.parsed_data.get("optimized_geometry")
        if optimized_geometry:
            for atom in atoms_data:
                for geom in optimized_geometry:
                    if geom.get("center") == atom["atom_number"]:
                        atom["geometry"] = geom
                        break
        total_match = total_pattern.search(npa_block)
        total_info = None
        if total_match:
            total_info = {
                "charge_sum":  float(total_match.group(1)),
                "core_sum":    float(total_match.group(2)),
                "valence_sum": float(total_match.group(3)),
                "rydberg_sum": float(total_match.group(4)),
                "total_sum":   float(total_match.group(5))
            }
        self.parsed_data["natural_population"] = {
            "atoms": atoms_data,
            "summary": total_info
        }

    def add_geometry_to_npa(self):
        """
        Adds the optimized geometry data to each atom in the natural population analysis.
        Assumes that the order of atoms in 'optimized_geometry' corresponds to the atom numbers.
        For each atom in the natural population data, a new key "geometry" is added
        containing its [x, y, z] coordinates.
        """
        geometry = self.parsed_data.get("optimized_geometry")
        npa_data = self.parsed_data.get("natural_population")
        if geometry is None or npa_data is None or "atoms" not in npa_data:
            return
        for atom in npa_data["atoms"]:
            atom_number = atom.get("atom_number")
            if atom_number is not None and atom_number <= len(geometry):
                atom["geometry"] = geometry[atom_number - 1]

    def extract_natural_populations(self, log_contents):
        """
        Extracts NATURAL POPULATIONS: Natural atomic orbital occupancies data
        for the first five atoms only. Stores data in a format compatible with
        JSON and table processing.
        
        The data structure follows the pattern:
        natural_populations.atoms[i].orbitals[j] = {
            "nao": int,
            "lang": str,
            "type_ao": str,
            "occupancy": float,
            "energy": float
        }
        """
        block_pattern = re.compile(
            r'NATURAL POPULATIONS:\s*Natural atomic orbital occupancies\s*\n'
            r'\s*\n'
            r'\s*NAO\s+Atom\s+No\s+lang\s+Type\(AO\)\s+Occupancy\s+Energy\s*\n'
            r'\s*-+\s*\n'
            r'(.*?)(?=\n\s*\n|\Z)',
            re.DOTALL
        )
        
        block_match = block_pattern.search(log_contents)
        if not block_match:
            self.parsed_data["lft_orbitals"] = None
            return
        
        populations_block = block_match.group(1)
        
        # Pattern to match each orbital line
        line_pattern = re.compile(
            r'^\s*(\d+)\s+([A-Z][a-z]?)\s+(\d+)\s+(\w+)\s+(.+?)\s+'
            r'([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s*$',
            re.MULTILINE
        )
        
        # Dictionary to store data by atom number
        atoms_data = {}
        
        for match in line_pattern.finditer(populations_block):
            nao = int(match.group(1))
            atom_symbol = match.group(2)
            atom_no = int(match.group(3))
            lang = match.group(4)
            type_ao = match.group(5).strip()
            occupancy = float(match.group(6))
            energy = float(match.group(7))
            
            # Only process first 5 atoms
            if atom_no > 5:
                continue
                
            # Initialize atom entry if not exists
            if atom_no not in atoms_data:
                atoms_data[atom_no] = {
                    "atom_symbol": atom_symbol,
                    "atom_number": atom_no,
                    "orbitals": []
                }
            
            # Add orbital data
            orbital_data = {
                "nao": nao,
                "lang": lang,
                "type_ao": type_ao,
                "occupancy": occupancy,
                "energy": energy
            }
            
            atoms_data[atom_no]["orbitals"].append(orbital_data)
        
        # Convert to list format ordered by atom number
        atoms_list = []
        for atom_no in sorted(atoms_data.keys()):
            atoms_list.append(atoms_data[atom_no])
        
        # Store in parsed_data
        self.parsed_data["lft_orbitals"] = {
            "atoms": atoms_list
        }

    def extract_nbo_bond_orbitals(self, log_contents):
        """
        Parses the 'Bond orbital/ Coefficients/ Hybrids' NBO sections from the Gaussian log.
        """
        block_pattern = re.compile(
            r'\(Occupancy\)\s+Bond orbital/ Coefficients/ Hybrids\s*\n[-\s]+\n'
            r'(.*?)\n\s*\n',
            re.DOTALL
        )
        block_match = block_pattern.search(log_contents)
        if not block_match:
            self.parsed_data["nbo_bond_orbitals"] = None
            return
        bond_block = block_match.group(1)
        orbital_start_pattern = re.compile(
            r'^\s*(\d+)\.\s*\(\s*([\d.]+)\)\s+(\S+)\s+(.*)$',
            re.MULTILINE
        )
        starts = list(orbital_start_pattern.finditer(bond_block))
        nbo_data = []
        for i, match in enumerate(starts):
            orb_index = int(match.group(1))
            occupancy = float(match.group(2))
            orb_type = match.group(3)
            atom_line = match.group(4).strip()
            start_pos = match.end()
            if i < len(starts) - 1:
                end_pos = starts[i+1].start()
            else:
                end_pos = len(bond_block)
            contribution_text = bond_block[start_pos:end_pos].strip()
            contrib_pattern = re.compile(
                r'^\s*\(\s*([\d.]+)%\)\s+([-\d.]+)\*([A-Za-z]+)\s+(\d+)?(.*?)(?=\n\s*\(\s*[\d.]+%|\Z)',
                re.MULTILINE | re.DOTALL
            )
            contributions = []
            for c_match in contrib_pattern.finditer(contribution_text):
                percentage = float(c_match.group(1))
                coefficient = float(c_match.group(2))
                atom_sym = c_match.group(3)
                atom_idx = c_match.group(4)
                remainder = c_match.group(5)
                atom_index = int(atom_idx) if atom_idx else None
                orbitals_info = []
                orb_pattern = re.compile(r'([spdfSPDF])\s*(?:([-\d.]+))?\(\s*([-\d.]+)%\)')
                for o_match in orb_pattern.finditer(remainder):
                    orb_label = o_match.group(1).lower()
                    exponent_str = o_match.group(2)
                    orb_percent_str = o_match.group(3)
                    orb_exponent = float(exponent_str) if exponent_str else None
                    orb_percentage = float(orb_percent_str)
                    orbitals_info.append({
                        "orbital": orb_label,
                        "exponent": orb_exponent,
                        "percentage": orb_percentage
                    })
                float_pattern = re.compile(r'[-+]?\d*\.\d+|[-+]?\d+')
                expansions_list = [float(x) for x in float_pattern.findall(remainder)]
                contributions.append({
                    "percentage": percentage,
                    "coefficient": coefficient,
                    "atom_symbol": atom_sym,
                    "atom_index": atom_index,
                    "orbitals": orbitals_info,
                    "coeff_expansion": expansions_list
                })
            nbo_data.append({
                "index": orb_index,
                "occupancy": occupancy,
                "type": orb_type,
                "atom_line": atom_line,
                "contributions": contributions
            })
        self.parsed_data["nbo_bond_orbitals"] = nbo_data

    def extract_nho_directionality(self, log_contents):
        """
        Extracts the 'NHO Directionality and "Bond Bending"' table from the Gaussian logfile.
        """
        def parse_float_or_none(value):
            return None if value == "--" else float(value)
        pattern = re.compile(
            r'^\s*(\d+)\.\s+(\S+)\s+(.*?)\s+'
            r'([-\d.]+|--)\s+([-\d.]+|--)\s+([-\d.]+|--)\s+([-\d.]+|--)'
            r'\s+([-\d.]+|--)\s+([-\d.]+|--)\s+([-\d.]+|--)\s+([-\d.]+|--)$',
            re.MULTILINE
        )
        matches = pattern.findall(log_contents)
        nho_data = []
        for match in matches:
            orb_index_str, bond_type, bond_label = match[0], match[1], match[2]
            col_values = match[3:]
            float_cols = [parse_float_or_none(val) for val in col_values]
            line_of_centers = {"theta": float_cols[0], "phi": float_cols[1]}
            hybrid1 = {"theta": float_cols[2], "phi": float_cols[3], "dev": float_cols[4]}
            hybrid2 = {"theta": float_cols[5], "phi": float_cols[6], "dev": float_cols[7]}
            nho_data.append({
                "index": int(orb_index_str),
                "bond_type": bond_type,
                "bond_label": bond_label,
                "line_of_centers": line_of_centers,
                "hybrid1": hybrid1,
                "hybrid2": hybrid2
            })
        self.parsed_data["nho_directionality"] = nho_data

    def extract_natural_electron_config(self, log_contents):
        """
        Extracts the "Natural Electron Configuration" block from the log.
        It looks for the header starting with "Atom  No" and ending with
        "Natural Electron Configuration", and then captures all text until
        the "NATURAL BOND ORBITAL ANALYSIS:" header (or end-of-file).
        """
        block_start_pattern = re.compile(
            r'Atom\s+No\s+.*?Natural Electron Configuration\s*\n[-\s]+\n'
            r'(.*?)(?=\n\s*NATURAL BOND ORBITAL ANALYSIS:|\Z)',
            re.DOTALL | re.MULTILINE
        )
        block_match = block_start_pattern.search(log_contents)
        if not block_match:
            self.parsed_data["natural_electron_config"] = None
            return
        block_text = block_match.group(1)
        line_pattern = re.compile(r'^\s*([A-Z][a-z]?)\s+(\d+)\s+(.*?)\s*$', re.MULTILINE)
        orbital_pattern = re.compile(r'(\[core\])|(\d+[spdfSPDF]\(\s*[\d.]+\s*\))')
        data_list = []
        for match in line_pattern.finditer(block_text):
            atom_symbol = match.group(1)
            atom_no = int(match.group(2))
            config_str = match.group(3).strip()
            orbitals_found = orbital_pattern.findall(config_str)
            core_present = False
            orbitals_dict = {}
            for core_group, orb_group in orbitals_found:
                if core_group:
                    core_present = True
                elif orb_group:
                    orb_label, occupancy = re.match(r'^(\d+[spdfSPDF])\(\s*([\d.]+)\s*\)$', orb_group).groups()
                    orbitals_dict[orb_label] = float(occupancy)
            data_list.append({
                "atom_symbol": atom_symbol,
                "atom_number": atom_no,
                "core": core_present,
                "orbitals": orbitals_dict
            })
        self.parsed_data["natural_electron_config"] = data_list

    def extract_natural_charges(self, log_contents):
        """
        Extracts natural charges from the NBO analysis section.
        Looks for a table with the pattern:
                        Natural  -----------------------------------------------
            Atom  No    Charge         Core      Valence    Rydberg      Total
         -----------------------------------------------------------------------
             Fe    1   -1.30258      8.99566     5.27320    0.03371    14.30258
        """
        # Pattern to find the natural charge table
        block_pattern = re.compile(
            r'Natural\s+[-]+\s*\n\s*Atom\s+No\s+Charge\s+Core\s+Valence\s+Rydberg\s+Total\s*\n'
            r'\s*[-]+\s*\n'
            r'(.*?)(?=\n\s*[-]+|\n\s*\n|\Z)',
            re.DOTALL
        )
        
        block_match = block_pattern.search(log_contents)
        if not block_match:
            self.parsed_data["natural_charges"] = None
            return
        
        charges_block = block_match.group(1)
        
        # Pattern to match each atom line
        # Example: Fe    1   -1.30258      8.99566     5.27320    0.03371    14.30258
        line_pattern = re.compile(
            r'^\s*([A-Z][a-z]?)\s+(\d+)\s+'
            r'([-+]?\d+\.\d+)\s+'
            r'([-+]?\d+\.\d+)\s+'
            r'([-+]?\d+\.\d+)\s+'
            r'([-+]?\d+\.\d+)\s+'
            r'([-+]?\d+\.\d+)\s*$',
            re.MULTILINE
        )
        
        atoms_data = []
        for match in line_pattern.finditer(charges_block):
            atom_symbol = match.group(1)
            atom_no = int(match.group(2))
            charge = float(match.group(3))
            core = float(match.group(4))
            valence = float(match.group(5))
            rydberg = float(match.group(6))
            total = float(match.group(7))
            
            atoms_data.append({
                "atom_symbol": atom_symbol,
                "atom_number": atom_no,
                "natural_charge": charge,
                "core": core,
                "valence": valence,
                "rydberg": rydberg,
                "total": total
            })
        
        self.parsed_data["natural_charges"] = atoms_data

    def extract_mulliken_charges_atomic(self, log_contents):
        """
        Extracts blocks of per-atom Mulliken charges. Keeps all occurrences indexed,
        defaulting to the last one.
        """
        pattern = re.compile(
            r"Mulliken charges:\s*\n((?:\s*\d+\s+\w+\s+[-+]?\d*\.\d+\s*\n?)+)",
            re.IGNORECASE
        )
        matches = list(pattern.finditer(log_contents))
        if not matches:
            self.parsed_data["mulliken_charges_atomic"] = None
            return
        charges_dict = {}
        for idx, m in enumerate(matches):
            block = m.group(1).strip().splitlines()
            atom_charges = []
            for line in block:
                parts = line.split()
                if len(parts) >= 3:
                    atom_charges.append({
                        "atom_index": int(parts[0]),
                        "element": parts[1],
                        "mulliken_charge": float(parts[2])
                    })
            charges_dict[idx] = atom_charges
        self.parsed_data["mulliken_charges_atomic"] = charges_dict
        # default: last occurrence
        last_idx = max(charges_dict.keys())
        self.parsed_data["mulliken_charges_atomic_default"] = charges_dict[last_idx]

    def extract_dipole_moment(self, log_contents):
        """
        Extracts the molecular dipole moment in atomic units from the log.
        """
        pattern = re.compile(
            r"Electric dipole moment.*?\(au\).*?\n\s*Tot\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*x\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*y\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*z\s*([-\d\.]+)D([+-]\d+)",
            re.DOTALL
        )
        m = pattern.search(log_contents)
        if not m:
            self.parsed_data["dipole_moment"] = None
            return
        def d_to_f(mag, expo):
            return float(mag) * 10**int(expo)
        tot = d_to_f(m.group(1), m.group(2))
        x = d_to_f(m.group(3), m.group(4))
        y = d_to_f(m.group(5), m.group(6))
        z = d_to_f(m.group(7), m.group(8))
        self.parsed_data["dipole_moment"] = {
            "tot_au": tot,
            "x_au": x,
            "y_au": y,
            "z_au": z
        }

    def extract_dipole_polarizability(self, log_contents):
        """
        Extracts the static dipole polarizability (Alpha) in atomic units.
        If multiple orientations are found, prefers dipole orientation over input orientation.
        """
        pattern = re.compile(
            r"Dipole polarizability.*?Alpha\s*\(([^)]+)\).*?Alpha\(0;0\):\s*\n\s*.*?\(au\).*?\n\s*iso\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*aniso\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*xx\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*yx\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*yy\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*zx\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*zy\s*([-\d\.]+)D([+-]\d+)"
            r".*?\n\s*zz\s*([-\d\.]+)D([+-]\d+)",
            re.DOTALL
        )
        
        matches = list(pattern.finditer(log_contents))
        if not matches:
            self.parsed_data["dipole_polarizability"] = None
            return

        def d_to_f(mag, expo):
            return float(mag) * 10**int(expo)

        # Parse all matches and organize by orientation type
        orientations = {}
        for match in matches:
            orientation = match.group(1).strip()
            raw_values = match.groups()[1:]  # Skip the orientation group
            pairs = list(zip(raw_values[0::2], raw_values[1::2]))
            values = [d_to_f(mag, exp) for mag, exp in pairs]
            keys = ["iso", "aniso", "xx", "yx", "yy", "zx", "zy", "zz"]
            orientations[orientation] = dict(zip(keys, values))

        # Store all orientations for completeness
        self.parsed_data["dipole_polarizability_all_orientations"] = orientations
        
        # Choose preferred orientation: dipole > input > first found
        preferred_orientation = None
        if "dipole orientation" in orientations:
            preferred_orientation = "dipole orientation"
        elif "input orientation" in orientations:
            preferred_orientation = "input orientation"
        else:
            # Take the first one found
            preferred_orientation = list(orientations.keys())[0]
        
        self.parsed_data["dipole_polarizability"] = orientations[preferred_orientation]
        self.parsed_data["dipole_polarizability_orientation"] = preferred_orientation

    def extract_quadrupole_moment(self, log_contents):
        """
        Extracts the first quadrupole moment tensor (Debye-Å) from the log.
        """
        pattern = re.compile(
            r"Quadrupole moment.*?\(field-independent basis.*?\)\s*\n\s*XX=\s*([-\d\.]+)"
            r"\s*YY=\s*([-\d\.]+)\s*ZZ=\s*([-\d\.]+)\s*\n\s*XY=\s*([-\d\.]+)"
            r"\s*XZ=\s*([-\d\.]+)\s*YZ=\s*([-\d\.]+)",
            re.DOTALL
        )
        m = pattern.search(log_contents)
        if not m:
            self.parsed_data["quadrupole_moment"] = None
            return
        vals = [float(v) for v in m.groups()]
        self.parsed_data["quadrupole_moment"] = {
            "XX": vals[0], "YY": vals[1], "ZZ": vals[2],
            "XY": vals[3], "XZ": vals[4], "YZ": vals[5]
        }

    def extract_atomic_dipole_orientation(self, log_contents):
        """
        Extracts per-atom dipole orientation vectors (x,y,z) in degrees from the log.
        """
        pattern = re.compile(
            r"Dipole orientation:\s*\n((?:\s*\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s*\n?)+)",
            re.IGNORECASE
        )
        m = pattern.search(log_contents)
        if not m:
            self.parsed_data["atomic_dipole_orientation"] = None
            return
        block = m.group(1).strip().splitlines()
        orientations = []
        for line in block:
            parts = line.split()
            orientations.append({
                "atom_index": int(parts[0]),
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3])
            })
        self.parsed_data["atomic_dipole_orientation"] = orientations

    def extract_overview(self, log_contents):
        """
        Extracts the final overview: spatial extent, net charge, field‐independent
        dipole, quadrupole, traceless quad, octapole, hexadecapole, N-N, E-N, and KE.
        """
        od = {}
        # electronic spatial extent
        m = re.search(r"Electronic spatial extent.*?<R\*\*2>=\s*([-\d\.]+)", log_contents)
        od["spatial_extent_R2_au"] = float(m.group(1)) if m else None
        # net charge
        m = re.search(r"Charge=\s*([-\d\.]+)", log_contents)
        od["net_charge"] = float(m.group(1)) if m else None
        # field‐independent dipole
        m = re.search(
            r"Dipole moment .*?X=\s*([-\d\.]+)\s*Y=\s*([-\d\.]+)\s*Z=\s*([-\d\.]+)\s*Tot=\s*([-\d\.]+)",
            log_contents, re.DOTALL
        )
        if m:
            od["overview_dipole"] = {
                "X": float(m.group(1)),
                "Y": float(m.group(2)),
                "Z": float(m.group(3)),
                "Tot": float(m.group(4))
            }
        else:
            od["overview_dipole"] = None
        # traceless quadrupole
        m = re.search(
            r"Traceless Quadrupole moment.*?\n\s*XX=\s*([-\d\.]+)\s*YY=\s*([-\d\.]+)\s*ZZ=\s*([-\d\.]+)"
            r"\s*\n\s*XY=\s*([-\d\.]+)\s*XZ=\s*([-\d\.]+)\s*YZ=\s*([-\d\.]+)",
            log_contents, re.DOTALL
        )
        if m:
            vals = [float(v) for v in m.groups()]
            od["traceless_quadrupole"] = {
                "XX": vals[0], "YY": vals[1], "ZZ": vals[2],
                "XY": vals[3], "XZ": vals[4], "YZ": vals[5]
            }
        else:
            od["traceless_quadrupole"] = None
        # octapole
        m = re.search(
            r"Octapole moment.*?\n\s*XXX=\s*([-\d\.]+)\s*YYY=\s*([-\d\.]+)\s*ZZZ=\s*([-\d\.]+)\s*XYY=\s*([-\d\.]+)"
            r".*?\n\s*XXY=\s*([-\d\.]+)\s*XXZ=\s*([-\d\.]+)\s*XZZ=\s*([-\d\.]+)\s*YZZ=\s*([-\d\.]+)"
            r".*?\n\s*YYZ=\s*([-\d\.]+)\s*XYZ=\s*([-\d\.]+)",
            log_contents, re.DOTALL
        )
        if m:
            vals = [float(v) for v in m.groups()]
            od["octapole"] = {
                "XXX": vals[0], "YYY": vals[1], "ZZZ": vals[2],
                "XYY": vals[3], "XXY": vals[4], "XXZ": vals[5],
                "XZZ": vals[6], "YZZ": vals[7], "YYZ": vals[8],
                "XYZ": vals[9]
            }
        else:
            od["octapole"] = None
        # hexadecapole
        m = re.search(
            r"Hexadecapole moment.*?\n\s*XXXX=\s*([-\d\.]+)\s*YYYY=\s*([-\d\.]+)\s*ZZZZ=\s*([-\d\.]+)\s*XXXY=\s*([-\d\.]+)"
            r".*?\n\s*XXXZ=\s*([-\d\.]+)\s*YYYX=\s*([-\d\.]+)\s*YYYZ=\s*([-\d\.]+)\s*ZZZX=\s*([-\d\.]+)"
            r".*?\n\s*ZZZY=\s*([-\d\.]+)\s*XXYY=\s*([-\d\.]+)\s*XXZZ=\s*([-\d\.]+)\s*YYZZ=\s*([-\d\.]+)"
            r".*?\n\s*XXYZ=\s*([-\d\.]+)\s*YYXZ=\s*([-\d\.]+)\s*ZZXY=\s*([-\d\.]+)",
            log_contents, re.DOTALL
        )
        if m:
            vals = [float(v) for v in m.groups()]
            od["hexadecapole"] = {
                "XXXX": vals[0], "YYYY": vals[1], "ZZZZ": vals[2],
                "XXXY": vals[3], "XXXZ": vals[4], "YYYX": vals[5],
                "YYYZ": vals[6], "ZZZX": vals[7], "ZZZY": vals[8],
                "XXYY": vals[9], "XXZZ": vals[10], "YYZZ": vals[11],
                "XXYZ": vals[12], "YYXZ": vals[13], "ZZXY": vals[14]
            }
        else:
            od["hexadecapole"] = None
        # N-N, E-N, KE
        m = re.search(
            r"N-N=\s*([-\d\.D+]+)\s*E-N=\s*([-\d\.D+]+)\s*KE=\s*([-\d\.D+]+)",
            log_contents
        )
        if m:
            od["N-N_au"] = float(m.group(1).replace('D','E'))
            od["E-N_au"] = float(m.group(2).replace('D','E'))
            od["KE_au"]  = float(m.group(3).replace('D','E'))
        else:
            od["N-N_au"] = od["E-N_au"] = od["KE_au"] = None
        self.parsed_data["overview"] = od

    def add_csv_information(self, log_filename):
        """
        Adds information from the CSV file to the parsed data.
        """
        identifier = log_filename[:4]
        matching_row = self.csv_data[self.csv_data['# PDB'].astype(str) == identifier]
        if not matching_row.empty:
            row_dict = matching_row.iloc[0].to_dict()
            self.parsed_data.update({key: self._convert_to_serializable(value) for key, value in row_dict.items()})
        else:
            print(f"No matching entry in CSV for identifier: {identifier}")

    def _convert_to_serializable(self, value):
        """
        Ensures the value is JSON serializable.
        """
        if isinstance(value, (np.int64, np.int32)):
            return int(value)
        elif isinstance(value, (np.float64, np.float32)):
            return float(value)
        elif isinstance(value, np.bool_):
            return bool(value)
        else:
            return value

    def delete_json_files_without_homo_lumo(self):
        """
        Deletes incomplete JSON files without HOMO/LUMO data.
        """
        for filename in os.listdir(self.json_dir):
            if filename.lower().endswith(".json"):
                json_path = os.path.join(self.json_dir, filename)
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                    alpha_occ = data.get("alpha_homo_lumo", {}).get("homo", [])
                    alpha_virt = data.get("alpha_homo_lumo", {}).get("lumo", [])
                    if not alpha_occ or not alpha_virt:
                        os.remove(json_path)
                        print(f"Deleted incomplete JSON (no alpha_occ/alpha_virt): {filename}")
                except (json.JSONDecodeError, OSError) as e:
                    print(f"Error reading JSON file {filename}: {e}")
        return

    def extract_mulliken_spin_densities(self, log_contents):
        """
        Extracts Mulliken charges and spin densities from the gaussian logfile.
        Finds sections like:
        'Mulliken charges and spin densities:'
        """
        pattern = re.compile(
            r"Mulliken charges and spin densities:\s*\n"
            r"((?:\s*\d+\s+\w+\s+[-+]?\d*\.\d+\s+[-+]?\d*\.\d+\s*\n?)+)"
            r".*?Sum of Mulliken charges\s*=\s*([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)",
            re.DOTALL | re.IGNORECASE
        )
        
        match = pattern.search(log_contents)
        if not match:
            self.parsed_data["mulliken_spin_densities"] = None
            return
        
        charges_block = match.group(1).strip()
        charge_sum = float(match.group(2))
        spin_sum = float(match.group(3))
        
        atoms_data = []
        line_pattern = re.compile(r"^\s*(\d+)\s+(\w+)\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)")
        
        for line in charges_block.splitlines():
            line_match = line_pattern.match(line)
            if line_match:
                atom_index = int(line_match.group(1))
                element = line_match.group(2)
                mulliken_charge = float(line_match.group(3))
                spin_density = float(line_match.group(4))
                
                atoms_data.append({
                    "atom_index": atom_index,
                    "element": element,
                    "mulliken_charge": mulliken_charge,
                    "spin_density": spin_density
                })
        
        self.parsed_data["mulliken_spin_densities"] = {
            "atoms": atoms_data,
            "charge_sum": charge_sum,
            "spin_sum": spin_sum
        }

    def extract_mulliken_spin_densities_heavy_atoms(self, log_contents):
        """
        Extracts Mulliken charges and spin densities with hydrogens summed into heavy atoms.
        Finds sections like:
        'Mulliken charges and spin densities with hydrogens summed into heavy atoms:'
        """
        pattern = re.compile(
            r"Mulliken charges and spin densities with hydrogens summed into heavy atoms:\s*\n"
            r"((?:\s*\d+\s+\w+\s+[-+]?\d*\.\d+\s+[-+]?\d*\.\d+\s*\n?)+)",
            re.DOTALL | re.IGNORECASE
        )
        
        match = pattern.search(log_contents)
        if not match:
            self.parsed_data["mulliken_spin_densities_heavy_atoms"] = None
            return
        
        charges_block = match.group(1).strip()
        atoms_data = []
        line_pattern = re.compile(r"^\s*(\d+)\s+(\w+)\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)")
        
        for line in charges_block.splitlines():
            line_match = line_pattern.match(line)
            if line_match:
                atom_index = int(line_match.group(1))
                element = line_match.group(2)
                mulliken_charge = float(line_match.group(3))
                spin_density = float(line_match.group(4))
                
                atoms_data.append({
                    "atom_index": atom_index,
                    "element": element,
                    "mulliken_charge": mulliken_charge,
                    "spin_density": spin_density
                })
        
        self.parsed_data["mulliken_spin_densities_heavy_atoms"] = atoms_data

    def extract_multiple_atomic_charges(self, log_contents):
        """
        Unified function to extract and organize all available atomic charges.
        Consolidates different charge types into a single comprehensive structure.
        """
        atomic_charges = {
            "available_types": [],
            "mulliken": {
                "basic": None,
                "with_spin": None,
                "heavy_atoms_only": None,
                "sum_charges": None,
                "sum_spin": None
            },
            "natural": {
                "charges": None,
                "populations": None
            },
            "summary": {
                "total_atoms": 0,
                "available_atoms": {},
                "charge_comparison": {}
            }
        }
        
        # Check what's already extracted and organize it
        if self.parsed_data.get("mulliken_charges_atomic"):
            atomic_charges["available_types"].append("mulliken_basic")
            atomic_charges["mulliken"]["basic"] = self.parsed_data["mulliken_charges_atomic"]
            if isinstance(atomic_charges["mulliken"]["basic"], dict):
                # Get the last (most recent) set of charges
                last_key = max(atomic_charges["mulliken"]["basic"].keys())
                atomic_charges["mulliken"]["basic"] = atomic_charges["mulliken"]["basic"][last_key]
        elif self.parsed_data.get("mulliken_heavy_atoms"):
            # Use heavy atoms data as basic Mulliken charges if available
            atomic_charges["available_types"].append("mulliken_basic")
            atomic_charges["mulliken"]["basic"] = [
                {
                    "atom_index": atom.get("atom_number", i+1),
                    "element": atom.get("element_symbol", atom.get("atom_symbol", "Unknown")),
                    "mulliken_charge": atom.get("mulliken_charge", 0.0)
                }
                for i, atom in enumerate(self.parsed_data["mulliken_heavy_atoms"])
            ]
        
        if self.parsed_data.get("mulliken_spin_densities"):
            atomic_charges["available_types"].append("mulliken_with_spin")
            atomic_charges["mulliken"]["with_spin"] = self.parsed_data["mulliken_spin_densities"]["atoms"]
            atomic_charges["mulliken"]["sum_charges"] = self.parsed_data["mulliken_spin_densities"]["charge_sum"]
            atomic_charges["mulliken"]["sum_spin"] = self.parsed_data["mulliken_spin_densities"]["spin_sum"]
        
        if self.parsed_data.get("mulliken_spin_densities_heavy_atoms"):
            atomic_charges["available_types"].append("mulliken_heavy_atoms")
            atomic_charges["mulliken"]["heavy_atoms_only"] = self.parsed_data["mulliken_spin_densities_heavy_atoms"]
        
        if self.parsed_data.get("natural_charges"):
            atomic_charges["available_types"].append("natural_charges")
            atomic_charges["natural"]["charges"] = self.parsed_data["natural_charges"]
        
        if self.parsed_data.get("natural_population"):
            atomic_charges["available_types"].append("natural_populations")
            atomic_charges["natural"]["populations"] = self.parsed_data["natural_population"]
        
        # Create summary information
        if atomic_charges["mulliken"]["basic"]:
            atomic_charges["summary"]["total_atoms"] = len(atomic_charges["mulliken"]["basic"])
            
            # Create atom index mapping
            for atom_data in atomic_charges["mulliken"]["basic"]:
                atom_idx = atom_data["atom_index"]
                element = atom_data["element"]
                atomic_charges["summary"]["available_atoms"][atom_idx] = element
        elif self.parsed_data.get("natural_charges"):
            atomic_charges["summary"]["total_atoms"] = len(self.parsed_data["natural_charges"])
            
            # Create atom index mapping from natural charges
            for atom_data in self.parsed_data["natural_charges"]:
                atom_idx = atom_data["atom_number"]
                element = atom_data["atom_symbol"]
                atomic_charges["summary"]["available_atoms"][atom_idx] = element
        
        # Compare charge values if multiple types available
        if len(atomic_charges["available_types"]) > 1:
            atomic_charges["summary"]["charge_comparison"] = self._compare_charge_values(atomic_charges)
        
        # Detect calculation type
        if "mulliken_with_spin" in atomic_charges["available_types"]:
            multiplicity = self._extract_multiplicity(log_contents)
            atomic_charges["summary"]["calculation_type"] = "open_shell" if multiplicity > 1 else "closed_shell"
            atomic_charges["summary"]["multiplicity"] = multiplicity
        else:
            atomic_charges["summary"]["calculation_type"] = "closed_shell"
            atomic_charges["summary"]["multiplicity"] = 1
        
        self.parsed_data["atomic_charges_comprehensive"] = atomic_charges

    def _compare_charge_values(self, atomic_charges):
        """
        Compare charge values between different methods for the same atoms.
        """
        comparison = {}
        
        if (atomic_charges["mulliken"]["basic"] and 
            atomic_charges["natural"]["charges"]):
            
            mulliken_dict = {atom["atom_index"]: atom["mulliken_charge"] 
                           for atom in atomic_charges["mulliken"]["basic"]}
            natural_dict = {atom["atom_number"]: atom["natural_charge"] 
                          for atom in atomic_charges["natural"]["charges"]}
            
            for atom_idx in mulliken_dict:
                if atom_idx in natural_dict:
                    comparison[atom_idx] = {
                        "mulliken": mulliken_dict[atom_idx],
                        "natural": natural_dict[atom_idx],
                        "difference": abs(mulliken_dict[atom_idx] - natural_dict[atom_idx])
                    }
        
        return comparison

    def _extract_multiplicity(self, log_contents):
        """
        Extract multiplicity from the log file to determine calculation type.
        """
        mult_pattern = re.compile(r"Multiplicity\s*=\s*(\d+)", re.IGNORECASE)
        match = mult_pattern.search(log_contents)
        return int(match.group(1)) if match else 1

    def extract_spin_contamination_analysis(self, log_contents):
        """
        Comprehensive extraction and analysis of spin contamination data.
        Extracts S**2 values, annihilation data, and calculates contamination metrics.
        """
        spin_analysis = {
            "initial_spin": None,
            "final_spin": None,
            "annihilation": None,
            "multiplicity": None,
            "theoretical_s_squared": None,
            "contamination_metrics": {},
            "validation": {}
        }
        
        # Extract multiplicity
        multiplicity = self._extract_multiplicity(log_contents)
        spin_analysis["multiplicity"] = multiplicity
        
        # Calculate theoretical S**2 value
        s_value = (multiplicity - 1) / 2.0
        theoretical_s_squared = s_value * (s_value + 1)
        spin_analysis["theoretical_s_squared"] = theoretical_s_squared
        
        # Extract initial guess spin values
        initial_pattern = re.compile(
            r"Initial guess <Sx>=\s*([-\d\.]+)\s*<Sy>=\s*([-\d\.]+)\s*<Sz>=\s*([-\d\.]+)\s*<S\*\*2>=\s*([-\d\.]+)\s*S=\s*([-\d\.]+)",
            re.IGNORECASE
        )
        initial_match = initial_pattern.search(log_contents)
        if initial_match:
            spin_analysis["initial_spin"] = {
                "Sx": float(initial_match.group(1)),
                "Sy": float(initial_match.group(2)),
                "Sz": float(initial_match.group(3)),
                "S_squared": float(initial_match.group(4)),
                "S": float(initial_match.group(5))
            }
        
        # Extract final SCF spin values (multiple possible formats)
        final_patterns = [
            re.compile(r"<Sx>=\s*([-\d\.]+)\s*<Sy>=\s*([-\d\.]+)\s*<Sz>=\s*([-\d\.]+)\s*<S\*\*2>=\s*([-\d\.]+)\s*S=\s*([-\d\.]+)"),
            re.compile(r"S\*\*2=\s*([-\d\.]+)")
        ]
        
        # Try comprehensive pattern first
        for pattern in final_patterns:
            matches = list(pattern.finditer(log_contents))
            if matches:
                final_match = matches[-1]  # Get the last occurrence (final SCF result)
                if len(final_match.groups()) == 5:
                    spin_analysis["final_spin"] = {
                        "Sx": float(final_match.group(1)),
                        "Sy": float(final_match.group(2)),
                        "Sz": float(final_match.group(3)),
                        "S_squared": float(final_match.group(4)),
                        "S": float(final_match.group(5))
                    }
                elif len(final_match.groups()) == 1:
                    # Just S**2 value
                    s_squared = float(final_match.group(1))
                    s_value_calc = (-1 + (1 + 4*s_squared)**0.5) / 2
                    spin_analysis["final_spin"] = {
                        "Sx": None,
                        "Sy": None,
                        "Sz": None,
                        "S_squared": s_squared,
                        "S": s_value_calc
                    }
                break
        
        # Extract annihilation data
        annihilation_pattern = re.compile(
            r"Annihilation of the first spin contaminant:\s*\n\s*S\*\*2 before annihilation\s+([-\d\.]+),\s*after\s+([-\d\.]+)",
            re.IGNORECASE | re.DOTALL
        )
        annihilation_match = annihilation_pattern.search(log_contents)
        if annihilation_match:
            spin_analysis["annihilation"] = {
                "before": float(annihilation_match.group(1)),
                "after": float(annihilation_match.group(2))
            }
        
        # Calculate contamination metrics
        if spin_analysis["final_spin"] and spin_analysis["final_spin"]["S_squared"] is not None:
            final_s_squared = spin_analysis["final_spin"]["S_squared"]
            contamination = abs(final_s_squared - theoretical_s_squared)
            contamination_percent = (contamination / theoretical_s_squared * 100) if theoretical_s_squared > 0 else 0
            
            spin_analysis["contamination_metrics"] = {
                "absolute_contamination": contamination,
                "percent_contamination": contamination_percent,
                "severity": self._classify_contamination_severity(contamination_percent)
            }
            
            # Add annihilation effectiveness if available
            if spin_analysis["annihilation"]:
                annihilation_improvement = abs(spin_analysis["annihilation"]["before"] - theoretical_s_squared) - abs(spin_analysis["annihilation"]["after"] - theoretical_s_squared)
                spin_analysis["contamination_metrics"]["annihilation_improvement"] = annihilation_improvement
                spin_analysis["contamination_metrics"]["annihilation_effectiveness"] = annihilation_improvement / abs(spin_analysis["annihilation"]["before"] - theoretical_s_squared) * 100 if abs(spin_analysis["annihilation"]["before"] - theoretical_s_squared) > 0 else 0
        
        # Validation checks
        spin_analysis["validation"] = self._validate_spin_state(spin_analysis)
        
        self.parsed_data["spin_contamination_analysis"] = spin_analysis

    def _classify_contamination_severity(self, percent_contamination):
        """
        Classify the severity of spin contamination based on percentage.
        """
        if percent_contamination < 1.0:
            return "negligible"
        elif percent_contamination < 5.0:
            return "minor"
        elif percent_contamination < 10.0:
            return "moderate"
        elif percent_contamination < 20.0:
            return "significant"
        else:
            return "severe"

    def _validate_spin_state(self, spin_analysis):
        """
        Validate the spin state and provide recommendations.
        """
        validation = {
            "is_valid_spin_state": True,
            "warnings": [],
            "recommendations": []
        }
        
        multiplicity = spin_analysis.get("multiplicity", 1)
        theoretical_s_squared = spin_analysis.get("theoretical_s_squared", 0)
        contamination_metrics = spin_analysis.get("contamination_metrics", {})
        
        # Check for closed shell (multiplicity = 1)
        if multiplicity == 1:
            if spin_analysis.get("final_spin") and spin_analysis["final_spin"]["S_squared"] is not None:
                if spin_analysis["final_spin"]["S_squared"] > 0.001:
                    validation["is_valid_spin_state"] = False
                    validation["warnings"].append("Significant S**2 value in claimed singlet state")
                    validation["recommendations"].append("Check for broken symmetry or use different initial guess")
        
        # Check contamination severity
        severity = contamination_metrics.get("severity", "negligible")
        if severity in ["moderate", "significant", "severe"]:
            validation["warnings"].append(f"Spin contamination is {severity}")
            if severity == "severe":
                validation["is_valid_spin_state"] = False
                validation["recommendations"].append("Consider different method or basis set")
                validation["recommendations"].append("Check geometry optimization convergence")
        
        # Check annihilation effectiveness
        if contamination_metrics.get("annihilation_effectiveness", 0) < 50:
            validation["warnings"].append("Spin annihilation not very effective")
            validation["recommendations"].append("Consider stability analysis")
        
        # Check for overflow cases (marked as ******* in Gaussian)
        if spin_analysis.get("final_spin") and spin_analysis["final_spin"].get("S") and spin_analysis["final_spin"]["S"] > 10:
            validation["is_valid_spin_state"] = False
            validation["warnings"].append("Extreme spin contamination detected (overflow)")
            validation["recommendations"].append("SCF convergence likely failed - recalculate with different settings")
        
        return validation

    def extract_bond_analysis_from_nbo(self, log_contents):
        """
        Extract bond-related information from NBO analysis for coordination analysis.
        While explicit Wiberg bond indices may not be available, this extracts:
        - Fe coordination information from geometry and NBO data
        - Bond orbital descriptions containing Fe
        - Natural population analysis for coordination assessment
        """
        bond_analysis = {
            "fe_coordination": {
                "fe_atom_index": None,
                "fe_bonds": [],
                "coordination_analysis": {}
            },
            "bond_orbitals": [],
            "wiberg_indices": {
                "available": False,
                "estimated_fe_coordination": {},
                "method": "estimated_from_nbo"
            },
            "validation": {
                "has_nbo_analysis": False,
                "fe_found": False,
                "coordination_number": 0
            }
        }
        
        # Check if NBO analysis is present
        nbo_present = "NBO" in log_contents or "Natural Bond Orbital" in log_contents
        bond_analysis["validation"]["has_nbo_analysis"] = nbo_present
        
        if not nbo_present:
            self.parsed_data["bond_analysis"] = bond_analysis
            return
        
        # Find Fe atom index from geometry or natural population analysis
        fe_atom_index = self._find_fe_atom_index()
        if fe_atom_index:
            bond_analysis["fe_coordination"]["fe_atom_index"] = fe_atom_index
            bond_analysis["validation"]["fe_found"] = True
        
        # Extract bond orbital descriptions containing Fe
        bond_orbitals = self._extract_fe_bond_orbitals(log_contents, fe_atom_index)
        bond_analysis["bond_orbitals"] = bond_orbitals
        
        # Estimate coordination from available data
        coordination_analysis = self._analyze_fe_coordination(log_contents, fe_atom_index, bond_orbitals)
        bond_analysis["fe_coordination"]["coordination_analysis"] = coordination_analysis
        
        # Look for explicit Wiberg bond indices (though unlikely to be found)
        wiberg_data = self._extract_wiberg_indices(log_contents)
        if wiberg_data:
            bond_analysis["wiberg_indices"]["available"] = True
            bond_analysis["wiberg_indices"]["data"] = wiberg_data
            bond_analysis["wiberg_indices"]["method"] = "explicit_from_log"
        else:
            # Estimate bond strengths from NBO data
            estimated_bonds = self._estimate_bond_indices_from_nbo(log_contents, fe_atom_index)
            bond_analysis["wiberg_indices"]["estimated_fe_coordination"] = estimated_bonds
        
        # Calculate coordination metrics
        bond_analysis["validation"]["coordination_number"] = len(bond_analysis["bond_orbitals"])
        
        self.parsed_data["bond_analysis"] = bond_analysis

    def _find_fe_atom_index(self):
        """
        Find the Fe atom index from already extracted data.
        """
        # Check natural charges data
        if self.parsed_data.get("natural_charges"):
            for atom in self.parsed_data["natural_charges"]:
                if atom["atom_symbol"].upper() == "FE":
                    return atom["atom_number"]
        
        # Check mulliken heavy atoms data  
        if self.parsed_data.get("mulliken_heavy_atoms"):
            for i, atom in enumerate(self.parsed_data["mulliken_heavy_atoms"]):
                if atom.get("element_symbol", "").upper() == "FE" or atom.get("atom_symbol", "").upper() == "FE":
                    return atom.get("atom_number", i+1)
        
        # Check standard orientation data
        if self.parsed_data.get("standard_orientation"):
            for atom in self.parsed_data["standard_orientation"]:
                if atom.get("atomic_number") == 26:  # Fe atomic number
                    return atom.get("center")  # Use "center" field for atom index
        
        return None

    def _extract_fe_bond_orbitals(self, log_contents, fe_atom_index):
        """
        Extract Fe coordination information from NBO analysis.
        Since direct Fe-N bonds are rare in NBO analysis, this looks for:
        1. Direct covalent bonds involving Fe (if any)
        2. Second-order perturbation interactions with Fe
        """
        if not fe_atom_index:
            return []
        
        bond_orbitals = []
        
        # Pattern 1: Look for direct Fe bonds (rare in heme systems)
        direct_bond_pattern = re.compile(
            rf"BD\s*\(\s*\d+\)\s*(Fe\s*{fe_atom_index}\s*-\s*(\w+)\s*(\d+)|(\w+)\s*(\d+)\s*-\s*Fe\s*{fe_atom_index})",
            re.IGNORECASE
        )
        
        for match in direct_bond_pattern.finditer(log_contents):
            if match.group(2) and match.group(3):  # Fe-X pattern
                partner_element = match.group(2)
                partner_index = int(match.group(3))
            elif match.group(4) and match.group(5):  # X-Fe pattern
                partner_element = match.group(4)
                partner_index = int(match.group(5))
            else:
                continue
                
            bond_orbital = {
                "partner_element": partner_element,
                "partner_atom_index": partner_index,
                "interaction_type": "direct_bond",
                "bond_type": "covalent"
            }
            bond_orbitals.append(bond_orbital)
        
        # Pattern 2: Look for second-order perturbation interactions involving Fe
        # Format: "number. BD (1) N 2 - C 6 / number. LP (6)Fe 1"
        coordination_pattern = re.compile(
            rf"(\d+)\.\s*BD\s*\(\s*\d+\)\s*(\w+)\s*(\d+)\s*-\s*(\w+)\s*(\d+)\s*/.*?LP\*?\s*\(\s*\d+\)Fe\s*{fe_atom_index}",
            re.IGNORECASE | re.DOTALL
        )
        
        existing_partners = set()
        for match in coordination_pattern.finditer(log_contents):
            # This represents a coordination interaction where a bond orbital
            # (like N-C) interacts with Fe lone pairs
            bond_element1 = match.group(2)
            bond_index1 = int(match.group(3))
            bond_element2 = match.group(4)
            bond_index2 = int(match.group(5))
            
            # Check if either atom in the bond is coordinating to Fe
            # N atoms are most likely to be coordinating
            coord_candidates = []
            if bond_element1.upper() in ['N', 'O', 'S']:  # Typical coordinating atoms
                coord_candidates.append((bond_element1, bond_index1))
            if bond_element2.upper() in ['N', 'O', 'S']:
                coord_candidates.append((bond_element2, bond_index2))
            
            for element, index in coord_candidates:
                partner_key = (element, index)
                if partner_key not in existing_partners:
                    existing_partners.add(partner_key)
                    bond_orbital = {
                        "partner_element": element,
                        "partner_atom_index": index,
                        "interaction_type": "second_order_perturbation",
                        "bond_type": "coordination",
                        "via_bond": f"{bond_element1}{bond_index1}-{bond_element2}{bond_index2}"
                    }
                    bond_orbitals.append(bond_orbital)
        
        # Pattern 3: Look for Fe lone pairs interacting with other atoms
        fe_lp_pattern = re.compile(
            rf"(\d+)\.\s*LP\*?\s*\(\s*\d+\)Fe\s*{fe_atom_index}\s*/.*?BD\*?\s*\(\s*\d+\)\s*(\w+)\s*(\d+)\s*-\s*(\w+)\s*(\d+)",
            re.IGNORECASE | re.DOTALL
        )
        
        for match in fe_lp_pattern.finditer(log_contents):
            # Fe lone pair donating to or accepting from other bonds
            bond_element1 = match.group(2)
            bond_index1 = int(match.group(3))
            bond_element2 = match.group(4)
            bond_index2 = int(match.group(5))
            
            # Add coordinating atoms that aren't already included
            for element, index in [(bond_element1, bond_index1), (bond_element2, bond_index2)]:
                if element.upper() in ['N', 'O', 'S', 'C'] and (element, index) not in existing_partners:
                    existing_partners.add((element, index))
                    bond_orbital = {
                        "partner_element": element,
                        "partner_atom_index": index,
                        "interaction_type": "fe_lone_pair_interaction",
                        "bond_type": "coordination",
                        "via_bond": f"{bond_element1}{bond_index1}-{bond_element2}{bond_index2}"
                    }
                    bond_orbitals.append(bond_orbital)
        
        return bond_orbitals

    def _analyze_fe_coordination(self, log_contents, fe_atom_index, bond_orbitals):
        """
        Analyze Fe coordination environment from available data.
        """
        analysis = {
            "coordination_number": 0,
            "heme_nitrogen_bonds": 0,
            "axial_bonds": 0,
            "bond_partners": {},
            "oxidation_state_estimate": None,
            "spin_state_estimate": None
        }
        
        if not fe_atom_index:
            return analysis
        
        # Count bonds by element type (using bond_orbitals parameter)
        bond_partners = {}
        for bond in bond_orbitals:
            element = bond["partner_element"].upper()
            if element not in bond_partners:
                bond_partners[element] = []
            bond_partners[element].append(bond["partner_atom_index"])
        
        analysis["bond_partners"] = bond_partners
        analysis["coordination_number"] = sum(len(indices) for indices in bond_partners.values())
        
        # Estimate heme vs axial bonds (N atoms likely in heme plane)
        if "N" in bond_partners:
            analysis["heme_nitrogen_bonds"] = len(bond_partners["N"])
        
        # Other elements likely axial
        analysis["axial_bonds"] = analysis["coordination_number"] - analysis["heme_nitrogen_bonds"]
        
        # Estimate oxidation state from d-electron count if available
        if self.parsed_data.get("natural_electron_config"):
            for atom_config in self.parsed_data["natural_electron_config"]:
                if atom_config.get("atom_number") == fe_atom_index:
                    d_occupancy = atom_config.get("d_occupancy", 0)
                    # Fe(II) typically has 6 d-electrons, Fe(III) has 5
                    if d_occupancy > 5.5:
                        analysis["oxidation_state_estimate"] = 2
                        analysis["spin_state_estimate"] = "high_spin" if d_occupancy > 6.5 else "low_spin"
                    elif d_occupancy > 4.5:
                        analysis["oxidation_state_estimate"] = 3
                        analysis["spin_state_estimate"] = "high_spin" if d_occupancy > 5.5 else "low_spin"
                    break
        
        return analysis

    def _extract_wiberg_indices(self, log_contents):
        """
        Look for explicit Wiberg bond index matrices (unlikely to be present).
        """
        wiberg_pattern = re.compile(
            r"Wiberg bond index matrix.*?\n(.*?)(?=\n\s*\n|\n\s*[A-Z]|\Z)",
            re.DOTALL | re.IGNORECASE
        )
        
        match = wiberg_pattern.search(log_contents)
        if match:
            # Parse the matrix data
            matrix_text = match.group(1)
            # Implementation would depend on exact format
            return {"matrix_text": matrix_text}
        
        return None

    def _estimate_bond_indices_from_nbo(self, log_contents, fe_atom_index):
        """
        Estimate bond strengths from NBO orbital occupancies and contributions.
        """
        estimated_bonds = {}
        
        if not fe_atom_index or not self.parsed_data.get("bond_analysis", {}).get("bond_orbitals"):
            return estimated_bonds
        
        for bond in self.parsed_data["bond_analysis"]["bond_orbitals"]:
            if bond.get("occupancy") and bond.get("fe_contribution_percent"):
                # Simple estimation based on occupancy and orbital overlap
                # Real Wiberg indices would require more sophisticated calculation
                partner_key = f"{bond['partner_element']}{bond['partner_atom_index']}"
                
                # Rough estimation: bond strength correlates with occupancy and Fe contribution
                occupancy = bond["occupancy"] / 100.0  # Convert percentage
                fe_contribution = bond["fe_contribution_percent"] / 100.0
                
                # Simplified bond index estimate
                estimated_index = occupancy * fe_contribution * 2.0  # Rough scaling
                estimated_bonds[partner_key] = {
                    "partner_element": bond["partner_element"],
                    "partner_atom_index": bond["partner_atom_index"],
                    "estimated_wiberg_index": estimated_index,
                    "basis": "nbo_occupancy_contribution"
                }
        
        return estimated_bonds

    def extract_nao_nlmo_analysis(self, log_contents):
        """
        Extract Natural Atomic Orbital (NAO) and Natural Localized Molecular Orbital (NLMO) 
        composition data for detailed orbital analysis, particularly focused on Fe coordination.
        """
        nao_analysis = {
            "natural_atomic_orbitals": {
                "fe_orbitals": [],
                "coordinating_atom_orbitals": [],
                "all_orbitals": []
            },
            "fe_d_orbital_analysis": {
                "individual_d_orbitals": {},
                "total_d_occupancy": 0.0,
                "d_orbital_ordering": [],
                "oxidation_state_estimate": None,
                "spin_state_analysis": {}
            },
            "coordination_hybridization": {
                "fe_ligand_bonds": [],
                "hybridization_summary": {}
            },
            "orbital_energies": {
                "fe_valence_orbitals": [],
                "homo_lumo_contributions": {}
            },
            "validation": {
                "has_nao_data": False,
                "fe_found": False,
                "total_orbitals_analyzed": 0
            }
        }
        
        # Check if NAO analysis is present
        nao_present = "NATURAL POPULATIONS" in log_contents and "Natural atomic orbital occupancies" in log_contents
        nao_analysis["validation"]["has_nao_data"] = nao_present
        
        if not nao_present:
            self.parsed_data["nao_nlmo_analysis"] = nao_analysis
            return
        
        # Extract Natural Atomic Orbital data
        nao_orbitals = self._extract_nao_orbital_data(log_contents)
        nao_analysis["natural_atomic_orbitals"]["all_orbitals"] = nao_orbitals
        
        # Find Fe atom and analyze its orbitals
        fe_atom_index = self._find_fe_atom_index()
        if fe_atom_index:
            nao_analysis["validation"]["fe_found"] = True
            
            # Extract Fe-specific orbital data
            fe_orbitals = [orbital for orbital in nao_orbitals if orbital.get("atom_number") == fe_atom_index]
            nao_analysis["natural_atomic_orbitals"]["fe_orbitals"] = fe_orbitals
            
            # Analyze Fe d-orbitals specifically
            d_orbital_analysis = self._analyze_fe_d_orbitals(fe_orbitals)
            nao_analysis["fe_d_orbital_analysis"] = d_orbital_analysis
            
            # Extract coordinating atom orbitals (N, S, O)
            coord_orbitals = [orbital for orbital in nao_orbitals 
                            if orbital.get("element", "").upper() in ["N", "S", "O"]]
            nao_analysis["natural_atomic_orbitals"]["coordinating_atom_orbitals"] = coord_orbitals
        
        # Extract coordination hybridization from NBO data
        hybridization_data = self._extract_coordination_hybridization(log_contents, fe_atom_index)
        nao_analysis["coordination_hybridization"] = hybridization_data
        
        # Extract orbital energy information
        orbital_energies = self._extract_orbital_energy_analysis(fe_orbitals if fe_atom_index else [])
        nao_analysis["orbital_energies"] = orbital_energies
        
        nao_analysis["validation"]["total_orbitals_analyzed"] = len(nao_orbitals)
        
        self.parsed_data["nao_nlmo_analysis"] = nao_analysis

    def _extract_nao_orbital_data(self, log_contents):
        """
        Extract detailed Natural Atomic Orbital data from the NAO section.
        """
        orbitals = []
        
        # Pattern for NAO orbital data
        # NAO  Atom  No  lang   Type(AO)    Occupancy      Energy
        nao_pattern = re.compile(
            r"^\s*(\d+)\s+(\w+)\s+(\d+)\s+(\w+)\s+(\w+\([^)]+\))\s+([\d\.]+)\s+([-\d\.]+)",
            re.MULTILINE
        )
        
        # Find the NAO section
        nao_section_start = log_contents.find("NATURAL POPULATIONS:  Natural atomic orbital occupancies")
        if nao_section_start == -1:
            return orbitals
        
        # Extract from NAO section
        nao_section = log_contents[nao_section_start:nao_section_start + 50000]  # Reasonable section size
        
        for match in nao_pattern.finditer(nao_section):
            orbital_data = {
                "nao_number": int(match.group(1)),
                "element": match.group(2),
                "atom_number": int(match.group(3)),
                "orbital_type": match.group(4),  # S, P, D, etc.
                "ao_description": match.group(5),  # Cor(1S), Val(3d), etc.
                "occupancy": float(match.group(6)),
                "energy": float(match.group(7))
            }
            
            # Parse additional orbital information
            ao_desc = orbital_data["ao_description"]
            if "Cor" in ao_desc:
                orbital_data["orbital_class"] = "core"
            elif "Val" in ao_desc:
                orbital_data["orbital_class"] = "valence"
            elif "Ryd" in ao_desc:
                orbital_data["orbital_class"] = "rydberg"
            else:
                orbital_data["orbital_class"] = "unknown"
            
            # Extract specific d-orbital types for transition metals
            if orbital_data["orbital_type"].lower().startswith("d"):
                orbital_data["d_orbital_type"] = orbital_data["orbital_type"]
            
            orbitals.append(orbital_data)
        
        return orbitals

    def _analyze_fe_d_orbitals(self, fe_orbitals):
        """
        Detailed analysis of Fe d-orbital occupancies and energies.
        """
        d_analysis = {
            "individual_d_orbitals": {},
            "total_d_occupancy": 0.0,
            "d_orbital_ordering": [],
            "oxidation_state_estimate": None,
            "spin_state_analysis": {}
        }
        
        # Extract d-orbital information
        d_orbitals = [orbital for orbital in fe_orbitals 
                     if orbital.get("orbital_class") == "valence" and 
                        orbital.get("orbital_type", "").lower().startswith("d")]
        
        for d_orbital in d_orbitals:
            d_type = d_orbital.get("d_orbital_type", d_orbital.get("orbital_type", ""))
            d_analysis["individual_d_orbitals"][d_type] = {
                "occupancy": d_orbital["occupancy"],
                "energy": d_orbital["energy"],
                "description": d_orbital["ao_description"]
            }
            d_analysis["total_d_occupancy"] += d_orbital["occupancy"]
        
        # Order d-orbitals by energy
        d_orbital_energies = [(d_type, data["energy"]) for d_type, data in d_analysis["individual_d_orbitals"].items()]
        d_orbital_energies.sort(key=lambda x: x[1])  # Sort by energy
        d_analysis["d_orbital_ordering"] = [orbital[0] for orbital in d_orbital_energies]
        
        # Estimate oxidation state from d-electron count
        total_d_electrons = d_analysis["total_d_occupancy"]
        if 5.5 < total_d_electrons <= 6.5:
            d_analysis["oxidation_state_estimate"] = 2  # Fe(II)
            d_analysis["spin_state_analysis"]["expected_d_electrons"] = 6
            d_analysis["spin_state_analysis"]["spin_state"] = "high_spin" if total_d_electrons > 6.2 else "low_spin"
        elif 4.5 < total_d_electrons <= 5.5:
            d_analysis["oxidation_state_estimate"] = 3  # Fe(III)
            d_analysis["spin_state_analysis"]["expected_d_electrons"] = 5
            d_analysis["spin_state_analysis"]["spin_state"] = "high_spin" if total_d_electrons > 5.2 else "low_spin"
        elif 3.5 < total_d_electrons <= 4.5:
            d_analysis["oxidation_state_estimate"] = 4  # Fe(IV)
            d_analysis["spin_state_analysis"]["expected_d_electrons"] = 4
        
        d_analysis["spin_state_analysis"]["d_electron_count"] = total_d_electrons
        
        return d_analysis

    def _extract_coordination_hybridization(self, log_contents, fe_atom_index):
        """
        Extract hybridization data for Fe coordination bonds from NBO analysis.
        """
        hybridization_data = {
            "fe_ligand_bonds": [],
            "hybridization_summary": {}
        }
        
        if not fe_atom_index:
            return hybridization_data
        
        # Pattern for bond orbital hybridization
        # ( 58.68%)   0.7660* N   2 s( 34.56%)p 1.88( 65.09%)d 0.01(  0.34%)
        hybridization_pattern = re.compile(
            rf"\(\s*([\d\.]+)%\)\s*([\d\.]+)\*\s*(\w+)\s*({fe_atom_index}|\d+)\s*s\(\s*([\d\.]+)%\)p\s*([\d\.]+)\(\s*([\d\.]+)%\)(?:d\s*([\d\.]+)\(\s*([\d\.]+)%\))?",
            re.IGNORECASE
        )
        
        # Look for Fe hybridization in bond descriptions
        for match in hybridization_pattern.finditer(log_contents):
            atom_number = int(match.group(4))
            if atom_number == fe_atom_index:
                element = match.group(3)
                hybrid_data = {
                    "element": element,
                    "atom_number": atom_number,
                    "contribution_percent": float(match.group(1)),
                    "coefficient": float(match.group(2)),
                    "s_percent": float(match.group(5)),
                    "p_character": float(match.group(6)),
                    "p_percent": float(match.group(7))
                }
                
                # Add d-character if present
                if match.group(8) and match.group(9):
                    hybrid_data["d_character"] = float(match.group(8))
                    hybrid_data["d_percent"] = float(match.group(9))
                
                hybridization_data["fe_ligand_bonds"].append(hybrid_data)
        
        # Summary statistics
        if hybridization_data["fe_ligand_bonds"]:
            total_bonds = len(hybridization_data["fe_ligand_bonds"])
            avg_s = sum(bond.get("s_percent", 0) for bond in hybridization_data["fe_ligand_bonds"]) / total_bonds
            avg_p = sum(bond.get("p_percent", 0) for bond in hybridization_data["fe_ligand_bonds"]) / total_bonds
            avg_d = sum(bond.get("d_percent", 0) for bond in hybridization_data["fe_ligand_bonds"]) / total_bonds
            
            hybridization_data["hybridization_summary"] = {
                "number_of_bonds": total_bonds,
                "average_s_character": avg_s,
                "average_p_character": avg_p,
                "average_d_character": avg_d,
                "hybridization_type": self._classify_hybridization(avg_s, avg_p, avg_d)
            }
        
        return hybridization_data

    def _extract_orbital_energy_analysis(self, fe_orbitals):
        """
        Extract and analyze orbital energy information for Fe.
        """
        energy_analysis = {
            "fe_valence_orbitals": [],
            "homo_lumo_contributions": {}
        }
        
        # Get valence orbitals and their energies
        valence_orbitals = [orbital for orbital in fe_orbitals 
                          if orbital.get("orbital_class") == "valence"]
        
        # Sort by energy
        valence_orbitals.sort(key=lambda x: x["energy"], reverse=True)  # Highest energy first
        
        energy_analysis["fe_valence_orbitals"] = [
            {
                "orbital_type": orbital["orbital_type"],
                "energy": orbital["energy"],
                "occupancy": orbital["occupancy"],
                "description": orbital["ao_description"]
            }
            for orbital in valence_orbitals
        ]
        
        return energy_analysis

    def _classify_hybridization(self, s_percent, p_percent, d_percent):
        """
        Classify hybridization type based on s, p, d character percentages.
        """
        if s_percent > 40:
            return "sp3-like"
        elif 25 < s_percent <= 40:
            return "sp2-like"
        elif 15 < s_percent <= 25:
            return "sp-like"
        elif d_percent > 20:
            return "d-hybridized"
        else:
            return "p-dominant"

    def extract_second_order_perturbation_analysis(self, log_contents):
        """
        Extract second order perturbation theory analysis from NBO for understanding
        donor-acceptor interactions, particularly Fe coordination and charge transfer.
        """
        perturbation_analysis = {
            "fe_interactions": {
                "fe_as_donor": [],
                "fe_as_acceptor": [],
                "high_energy_interactions": [],
                "coordination_interactions": []
            },
            "interaction_summary": {
                "total_interactions": 0,
                "fe_donor_count": 0,
                "fe_acceptor_count": 0,
                "average_energy": 0.0,
                "energy_threshold": 0.5  # kcal/mol
            },
            "charge_transfer_analysis": {
                "ligand_to_fe_donation": [],
                "fe_to_ligand_backbonding": [],
                "total_donation_energy": 0.0,
                "total_backbonding_energy": 0.0
            },
            "validation": {
                "has_perturbation_data": False,
                "fe_found": False,
                "section_count": 0
            }
        }
        
        # Check if second order perturbation analysis is present
        perturbation_present = "Second Order Perturbation Theory Analysis" in log_contents
        perturbation_analysis["validation"]["has_perturbation_data"] = perturbation_present
        
        if not perturbation_present:
            self.parsed_data["second_order_perturbation_analysis"] = perturbation_analysis
            return
        
        # Find Fe atom index
        fe_atom_index = self._find_fe_atom_index()
        if fe_atom_index:
            perturbation_analysis["validation"]["fe_found"] = True
        
        # Extract perturbation data from all sections
        perturbation_sections = self._find_perturbation_sections(log_contents)
        perturbation_analysis["validation"]["section_count"] = len(perturbation_sections)
        
        all_interactions = []
        for section in perturbation_sections:
            interactions = self._parse_perturbation_section(section, fe_atom_index)
            all_interactions.extend(interactions)
        
        # Classify and analyze interactions
        if fe_atom_index:
            self._classify_fe_interactions(all_interactions, fe_atom_index, perturbation_analysis)
            self._analyze_charge_transfer(all_interactions, fe_atom_index, perturbation_analysis)
        
        # Calculate summary statistics
        perturbation_analysis["interaction_summary"]["total_interactions"] = len(all_interactions)
        if all_interactions:
            energies = [interaction["energy_kcal"] for interaction in all_interactions]
            perturbation_analysis["interaction_summary"]["average_energy"] = sum(energies) / len(energies)
        
        self.parsed_data["second_order_perturbation_analysis"] = perturbation_analysis

    def _find_perturbation_sections(self, log_contents):
        """
        Find all second order perturbation theory analysis sections in the log.
        """
        sections = []
        
        # Pattern to find section starts
        section_pattern = re.compile(
            r"Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis.*?"
            r"(?=Second Order Perturbation Theory Analysis|Natural Bond Orbitals|NATURAL POPULATIONS|\Z)",
            re.DOTALL | re.IGNORECASE
        )
        
        for match in section_pattern.finditer(log_contents):
            sections.append(match.group(0))
        
        return sections

    def _parse_perturbation_section(self, section_text, fe_atom_index):
        """
        Parse individual perturbation interactions from a section.
        """
        interactions = []
        
        # Pattern for perturbation lines
        # Format: "146. LP (   1)Fe   1                /***. BD*(   1) N   2 - C   6            1.15    0.70    0.025"
        perturbation_pattern = re.compile(
            r"^\s*(\d+)\.\s+(\w+\*?)\s*\([^)]*\)\s*(\w+)\s+(\d+)\s+(/[^/]*)\s+([0-9.]+)\s+([0-9.-]+)\s+([0-9.-]+)",
            re.MULTILINE
        )
        
        for match in perturbation_pattern.finditer(section_text):
            try:
                donor_number = int(match.group(1))
                donor_type = match.group(2)
                donor_element = match.group(3)
                donor_atom = int(match.group(4))
                acceptor_description = match.group(5).strip()
                energy_kcal = float(match.group(6))
                energy_diff_au = float(match.group(7))
                fock_element = float(match.group(8))
                
                # Parse acceptor information from description
                acceptor_info = self._parse_acceptor_description(acceptor_description)
                
                interaction = {
                    "donor": {
                        "orbital_number": donor_number,
                        "orbital_type": donor_type,
                        "element": donor_element,
                        "atom_number": donor_atom
                    },
                    "acceptor": acceptor_info,
                    "energy_kcal": energy_kcal,
                    "energy_diff_au": energy_diff_au,
                    "fock_element": fock_element,
                    "involves_fe": (donor_atom == fe_atom_index) or 
                                  (acceptor_info.get("atom_number") == fe_atom_index if acceptor_info else False)
                }
                
                interactions.append(interaction)
                
            except (ValueError, IndexError) as e:
                # Skip malformed lines
                continue
        
        return interactions

    def _parse_acceptor_description(self, acceptor_desc):
        """
        Parse acceptor orbital description from the perturbation line.
        """
        acceptor_info = {
            "orbital_number": None,
            "orbital_type": None,
            "element": None,
            "atom_number": None,
            "bond_description": None
        }
        
        # Pattern for acceptor like "/151. LP*(   6)Fe   1" or "/***. BD*(   1) N   2 - C   6"
        acceptor_pattern = re.compile(r"/(\*{3}|\d+)\.\s+(\w+\*?)\s*\([^)]*\)\s*(\w+)\s*(\d+)(?:\s*-\s*(\w+)\s*(\d+))?")
        
        match = acceptor_pattern.search(acceptor_desc)
        if match:
            acceptor_info["orbital_number"] = match.group(1) if match.group(1) != "***" else None
            acceptor_info["orbital_type"] = match.group(2)
            acceptor_info["element"] = match.group(3)
            acceptor_info["atom_number"] = int(match.group(4))
            
            # For bond orbitals, capture the bond description
            if match.group(5) and match.group(6):
                acceptor_info["bond_description"] = f"{match.group(3)}{match.group(4)}-{match.group(5)}{match.group(6)}"
        
        return acceptor_info

    def _classify_fe_interactions(self, interactions, fe_atom_index, analysis):
        """
        Classify Fe-involving interactions by type and significance.
        """
        fe_as_donor = []
        fe_as_acceptor = []
        high_energy = []
        coordination_interactions = []
        
        for interaction in interactions:
            if not interaction["involves_fe"]:
                continue
            
            # Classify by Fe role
            if interaction["donor"]["atom_number"] == fe_atom_index:
                fe_as_donor.append(interaction)
                analysis["interaction_summary"]["fe_donor_count"] += 1
            
            if interaction["acceptor"] and interaction["acceptor"].get("atom_number") == fe_atom_index:
                fe_as_acceptor.append(interaction)
                analysis["interaction_summary"]["fe_acceptor_count"] += 1
            
            # Identify high energy interactions (>10 kcal/mol)
            if interaction["energy_kcal"] > 10.0:
                high_energy.append(interaction)
            
            # Identify coordination interactions (N, O, S to Fe)
            donor_element = interaction["donor"]["element"].upper()
            acceptor_element = interaction["acceptor"].get("element", "") if interaction["acceptor"] else ""
            acceptor_element = acceptor_element.upper() if acceptor_element else ""
            
            if ((donor_element in ["N", "O", "S"] and interaction["acceptor"] and interaction["acceptor"].get("atom_number") == fe_atom_index) or
                (acceptor_element in ["N", "O", "S"] and interaction["donor"]["atom_number"] == fe_atom_index)):
                coordination_interactions.append(interaction)
        
        analysis["fe_interactions"]["fe_as_donor"] = fe_as_donor
        analysis["fe_interactions"]["fe_as_acceptor"] = fe_as_acceptor
        analysis["fe_interactions"]["high_energy_interactions"] = high_energy
        analysis["fe_interactions"]["coordination_interactions"] = coordination_interactions

    def _analyze_charge_transfer(self, interactions, fe_atom_index, analysis):
        """
        Analyze charge transfer patterns involving Fe.
        """
        ligand_to_fe = []
        fe_to_ligand = []
        
        for interaction in interactions:
            if not interaction["involves_fe"]:
                continue
            
            donor_element = interaction["donor"]["element"].upper()
            acceptor_element = interaction["acceptor"].get("element", "") if interaction["acceptor"] else ""
            acceptor_element = acceptor_element.upper() if acceptor_element else ""
            
            # Ligand to Fe donation (typical coordination)
            if (donor_element in ["N", "O", "S", "C"] and 
                interaction["acceptor"] and interaction["acceptor"].get("atom_number") == fe_atom_index):
                ligand_to_fe.append(interaction)
            
            # Fe to ligand backbonding
            if (interaction["donor"]["atom_number"] == fe_atom_index and
                acceptor_element in ["N", "O", "S", "C"]):
                fe_to_ligand.append(interaction)
        
        # Calculate total energies
        total_donation = sum(interaction["energy_kcal"] for interaction in ligand_to_fe)
        total_backbonding = sum(interaction["energy_kcal"] for interaction in fe_to_ligand)
        
        analysis["charge_transfer_analysis"]["ligand_to_fe_donation"] = ligand_to_fe
        analysis["charge_transfer_analysis"]["fe_to_ligand_backbonding"] = fe_to_ligand
        analysis["charge_transfer_analysis"]["total_donation_energy"] = total_donation
        analysis["charge_transfer_analysis"]["total_backbonding_energy"] = total_backbonding

    def extract_comprehensive_orbital_energies(self, log_contents):
        """
        Enhanced orbital energy extraction with alpha/beta sorting, gap analysis,
        and orbital characterization for detailed electronic structure analysis.
        """
        orbital_analysis = {
            "orbital_energies": {
                "alpha_occupied": [],
                "alpha_virtual": [],
                "beta_occupied": [],
                "beta_virtual": []
            },
            "homo_lumo_analysis": {
                "alpha_homo": None,
                "alpha_lumo": None,
                "beta_homo": None,
                "beta_lumo": None,
                "alpha_gap": None,
                "beta_gap": None,
                "fundamental_gap": None
            },
            "orbital_statistics": {
                "total_electrons": 0,
                "alpha_electrons": 0,
                "beta_electrons": 0,
                "unpaired_electrons": 0,
                "highest_occupied_energy": None,
                "lowest_virtual_energy": None
            },
            "fe_orbital_contributions": {
                "fe_d_in_homo": False,
                "fe_d_in_lumo": False,
                "fe_d_character_analysis": {}
            },
            "validation": {
                "has_orbital_data": False,
                "is_unrestricted": False,
                "orbital_count": 0
            }
        }
        
        # Extract orbital energies using existing function
        orbital_energies = self.extract_orbital_energies(log_contents)
        
        if not orbital_energies or not any(orbital_energies.values()):
            self.parsed_data["comprehensive_orbital_analysis"] = orbital_analysis
            return
        
        orbital_analysis["validation"]["has_orbital_data"] = True
        orbital_analysis["orbital_energies"] = orbital_energies
        
        # Determine if calculation is unrestricted (has beta orbitals)
        is_unrestricted = bool(orbital_energies.get("beta_occ") or orbital_energies.get("beta_virt"))
        orbital_analysis["validation"]["is_unrestricted"] = is_unrestricted
        
        # Analyze HOMO/LUMO
        homo_lumo = self._analyze_homo_lumo_comprehensive(orbital_energies)
        orbital_analysis["homo_lumo_analysis"] = homo_lumo
        
        # Calculate orbital statistics
        orbital_stats = self._calculate_orbital_statistics(orbital_energies)
        orbital_analysis["orbital_statistics"] = orbital_stats
        
        # Count total orbitals
        total_orbitals = sum(len(orbs) for orbs in orbital_energies.values())
        orbital_analysis["validation"]["orbital_count"] = total_orbitals
        
        # Analyze Fe d-orbital contributions (if Fe is present and we have NAO data)
        fe_contributions = self._analyze_fe_orbital_contributions(orbital_energies)
        orbital_analysis["fe_orbital_contributions"] = fe_contributions
        
        self.parsed_data["comprehensive_orbital_analysis"] = orbital_analysis

    def _analyze_homo_lumo_comprehensive(self, orbital_energies):
        """
        Comprehensive HOMO/LUMO analysis for both alpha and beta spins.
        """
        homo_lumo = {
            "alpha_homo": None,
            "alpha_lumo": None,
            "beta_homo": None,
            "beta_lumo": None,
            "alpha_gap": None,
            "beta_gap": None,
            "fundamental_gap": None
        }
        
        # Alpha orbitals
        alpha_occ = orbital_energies.get("alpha_occ", [])
        alpha_virt = orbital_energies.get("alpha_virt", [])
        
        if alpha_occ:
            homo_lumo["alpha_homo"] = alpha_occ[-1]  # Highest occupied
        if alpha_virt:
            homo_lumo["alpha_lumo"] = alpha_virt[0]   # Lowest virtual
        
        if homo_lumo["alpha_homo"] is not None and homo_lumo["alpha_lumo"] is not None:
            homo_lumo["alpha_gap"] = homo_lumo["alpha_lumo"] - homo_lumo["alpha_homo"]
        
        # Beta orbitals (for unrestricted calculations)
        beta_occ = orbital_energies.get("beta_occ", [])
        beta_virt = orbital_energies.get("beta_virt", [])
        
        if beta_occ:
            homo_lumo["beta_homo"] = beta_occ[-1]
        if beta_virt:
            homo_lumo["beta_lumo"] = beta_virt[0]
        
        if homo_lumo["beta_homo"] is not None and homo_lumo["beta_lumo"] is not None:
            homo_lumo["beta_gap"] = homo_lumo["beta_lumo"] - homo_lumo["beta_homo"]
        
        # Fundamental gap (smallest gap)
        gaps = [gap for gap in [homo_lumo["alpha_gap"], homo_lumo["beta_gap"]] if gap is not None]
        if gaps:
            homo_lumo["fundamental_gap"] = min(gaps)
        
        return homo_lumo

    def _calculate_orbital_statistics(self, orbital_energies):
        """
        Calculate comprehensive orbital occupation statistics.
        """
        stats = {
            "total_electrons": 0,
            "alpha_electrons": 0,
            "beta_electrons": 0,
            "unpaired_electrons": 0,
            "highest_occupied_energy": None,
            "lowest_virtual_energy": None
        }
        
        # Count electrons
        alpha_occ = orbital_energies.get("alpha_occ", [])
        beta_occ = orbital_energies.get("beta_occ", [])
        
        stats["alpha_electrons"] = len(alpha_occ)
        stats["beta_electrons"] = len(beta_occ) if beta_occ else len(alpha_occ)  # Restricted case
        stats["total_electrons"] = stats["alpha_electrons"] + stats["beta_electrons"]
        stats["unpaired_electrons"] = abs(stats["alpha_electrons"] - stats["beta_electrons"])
        
        # Find energy extremes
        all_occupied = []
        all_virtual = []
        
        if alpha_occ:
            all_occupied.extend(alpha_occ)
        if beta_occ:
            all_occupied.extend(beta_occ)
        
        alpha_virt = orbital_energies.get("alpha_virt", [])
        beta_virt = orbital_energies.get("beta_virt", [])
        
        if alpha_virt:
            all_virtual.extend(alpha_virt)
        if beta_virt:
            all_virtual.extend(beta_virt)
        
        if all_occupied:
            stats["highest_occupied_energy"] = max(all_occupied)
        if all_virtual:
            stats["lowest_virtual_energy"] = min(all_virtual)
        
        return stats

    def _analyze_fe_orbital_contributions(self, orbital_energies):
        """
        Analyze Fe d-orbital contributions to frontier orbitals.
        """
        fe_contributions = {
            "fe_d_in_homo": False,
            "fe_d_in_lumo": False,
            "fe_d_character_analysis": {}
        }
        
        # This analysis would require orbital composition data which isn't directly
        # available from the orbital energies alone. It would need MO coefficients
        # or population analysis. For now, we provide the structure for future enhancement.
        
        # Check if we have NAO data that might inform about Fe d-orbital energies
        if hasattr(self, 'parsed_data') and self.parsed_data.get("nao_nlmo_analysis"):
            nao_data = self.parsed_data["nao_nlmo_analysis"]
            if nao_data.get("fe_d_orbital_analysis"):
                fe_d_analysis = nao_data["fe_d_orbital_analysis"]
                fe_contributions["fe_d_character_analysis"] = {
                    "d_orbital_count": len(fe_d_analysis.get("individual_d_orbitals", {})),
                    "total_d_occupancy": fe_d_analysis.get("total_d_occupancy", 0),
                    "oxidation_state_estimate": fe_d_analysis.get("oxidation_state_estimate"),
                    "d_orbital_energies": [
                        orbital_data.get("energy") 
                        for orbital_data in fe_d_analysis.get("individual_d_orbitals", {}).values()
                        if orbital_data.get("energy") is not None
                    ]
                }
        
        return fe_contributions


def main():
    parser = argparse.ArgumentParser(
        description="Parse Gaussian log files into JSON records."
    )
    parser.add_argument(
        "--log-dir",
        default=str(LOGFILES_DIR),
        help="Directory containing Gaussian .log files.",
    )
    parser.add_argument(
        "--json-dir",
        default=str(JSONS_DIR),
        help="Directory where parsed JSON files will be written.",
    )
    parser.add_argument(
        "--csv-path",
        default=str(canonical_table_path("pyDISH.csv")),
        help="Metadata CSV used to enrich parsed records.",
    )
    parser.add_argument(
        "--log-file",
        help="Optional single Gaussian log file to parse instead of the entire directory.",
    )
    parser.add_argument(
        "--output-json",
        help="Optional output JSON path when --log-file is used.",
    )
    args = parser.parse_args()

    parser_instance = g16parser(
        log_dir=args.log_dir,
        json_dir=args.json_dir,
        csv_path=args.csv_path,
    )

    if args.log_file:
        output_json = args.output_json
        if output_json is None:
            output_json = str(Path(args.json_dir) / f"{Path(args.log_file).stem}.json")
        parser_instance.parse_gaussian_logfile(args.log_file, output_json)
        return

    parser_instance.parse()


if __name__ == "__main__":
    main()
