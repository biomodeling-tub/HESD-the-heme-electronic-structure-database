import numpy as np
import pandas as pd
import pdb_tools
import MDAnalysis as mda
import os
import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    import Bio
    from Bio import PDB
    from Bio.PDB.PDBParser import PDBParser
warnings.filterwarnings("ignore", category=UserWarning)


class ops():
    def convert_ent_to_pdb(file, output):
        cleanup=True
        parser = PDBParser()
        structure = parser.get_structure('structure', file)
        io = Bio.PDB.PDBIO()
        io.set_structure(structure)
        io.save(output)
        if cleanup:
            os.remove(file)
            if os.path.exists(file):
                raise FileExistsError(f"File {file} exists but should have been deleted. ")
        if not os.path.exists(output):
            raise FileNotFoundError(f"File {output} not found, but should exist. ")
        return output

    def HETATM_to_ATM(file, output):
        verbose=False
        lines = pdb_tools.files.read_file(pdb_file=file)
        lines_ = ["CRYST1\n"]
        if verbose:
            print(lines)
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            pdb_tools.line_operations.exchange_HETATM_ATOM(line_dict=line_dict)
            if line_dict["atom_name"].strip() == "FE":
                line_dict["elem_symb"] = "FE"
                line_dict["segment"] = "    "
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        lines_.append("")
        if verbose:
            print(lines_)
        pdb_tools.files.write_file(file=output, lines=lines_)
        return

    def pdb_to_xyz(file, output):
        with open(file, 'r') as pdb_file, open(output, 'w') as xyz_file:
            lines = pdb_file.readlines()
            num_atoms = 0
            atom_lines = []
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    num_atoms += 1
                    atom_lines.append(line)
            xyz_file.write(str(num_atoms) + '\n')
            xyz_file.write('Converted from PDB to XYZ\n')
            for line in atom_lines:
                tokens = line.split()
                atom_symbol = tokens[2]
                atom_symbol = atom_symbol[0]
                if atom_symbol == "F":
                    atom_symbol = "Fe"
                if atom_symbol == "X":
                    atom_symbol = "O"
                #print(atom_symbol)
                x, y, z = float(tokens[5]), float(tokens[6]), float(tokens[7])
                xyz_file.write(f'{atom_symbol} {x:.6f} {y:.6f} {z:.6f}\n')
        return

    def concatenate_xyz(input1, input2, output):
        # Read the content of both XYZ files
        with open(input1, 'r') as xyz_file1, open(input2, 'r') as xyz_file2:
            lines1 = xyz_file1.readlines()
            lines2 = xyz_file2.readlines()
        # Calculate the total number of atoms
        num_atoms1 = int(lines1[0])
        num_atoms2 = int(lines2[0])
        total_atoms = num_atoms1 + num_atoms2
        # Concatenate the XYZ content and update the total atom count in the first line
        concatenated_lines = [f"{total_atoms}\n"] + lines1[1:] + lines2[2:]
        # Write the concatenated content to the output XYZ file
        with open(output, 'w') as output_file:
            output_file.writelines(concatenated_lines)
        return

    def format_xyz_file(pdb_id):
        file, output = f"PDB/{pdb_id}.xyz", f"PDB/{pdb_id}.xyz"
        with open(file, 'r') as f:
            lines = f.readlines()
        atoms, coordinates = [], []
        for line in lines[2:]:
            parts = line.split()
            atoms.append(parts[0])
            coordinates.append([float(coord) for coord in parts[1:]])
        max_widths = [max(len(str(coord)) for coord in column) for column in zip(*coordinates)]
        with open(output, 'w') as f:
            f.write(f'{len(atoms)}\n')
            f.write(f'Heme structure file of {file} \n')
            for atom, coords in zip(atoms, coordinates):
                coord_str = ' '.join(f'{coord:{width}.6f}' for coord, width in zip(coords, max_widths))
                f.write(f'{atom:<2s} {coord_str}\n')
        return

    def xyz_to_pdb(file, output):
        with open(file, 'r') as xyz_file, open(output, 'w') as pdb_file:
            lines_ = []
            pdb_file.write("HEADER    XYZ to PDB Conversion\n")
            residue_number = 1
            atom_number = 1
            for i, line in enumerate(xyz_file):
                if i < 2:
                    continue  # Skip the first two lines of .xyz-file header
                parts = line.split()
                atom_symbol, x, y, z = parts
                pdb_line = f"ATOM     {atom_number:2}  {atom_symbol:2}  HEM H{residue_number:2}      {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}  1.00  0.00           {atom_symbol:4}\n"
                pdb_file.write(pdb_line)
                atom_number += 1
                if atom_number > 25:
                    break
            pdb_file.write("END\n")
        return

    def is_xyz_empty(pdb_id):
        verbose=False
        try:
            with open(f"PDB/{pdb_id}_system_protonated.xyz", 'r') as pdb:
                not any(line.strip() for line in pdb)
                if verbose:
                    print(f"File 'PDB/{pdb_id}_system_protonated.xyz' is empty.")
                ops.pdb_to_xyz(file=f"PDB/{pdb_id}_system_protonated.pdb", output=f"PDB/{pdb_id}_system_protonated.xyz")
                if verbose:
                    print(f"File 'PDB/{pdb_id}_system_protonated.pdb' converted to File 'PDB/{pdb_id}_system_protonated.xyz'.")
                return
        except FileNotFoundError:
            print(f"Error: PDB file 'PDB/{pdb_id}_system_protonated.xyz' not found.")
            return
        except Exception as e:
            print(f"Error while checking if PDB file is empty: {e}")
            return

    def check_pdb_file(pdb_id):
        verbose = True
        try:
            with open(f"PDB/{pdb_id}_system_protonated.pdb", 'r') as pdb:
                pdb_content = pdb.read()
                if verbose:
                    if "9999.00" in pdb_content:
                        print("Charmm protonation is not successful")
                    else:
                        print("Charmm protonation is successful")
        except FileNotFoundError:
            print(f"Error: PDB file PDB/{pdb_id}_system_protonated.pdb not found.")
        except Exception as e:
            print(f"Error while reading PDB file: {e}")
        return