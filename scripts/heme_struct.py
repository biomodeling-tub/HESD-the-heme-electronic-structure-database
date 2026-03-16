import numpy as np
import pandas as pd
import matplotlib as mpl
import pdb_tools
import subprocess
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.distances import distance_array
import os, sys, shutil
from contextlib import redirect_stdout, redirect_stderr
import orient_heme
from collections import Counter
from file_ops import ops
from MDAnalysis.analysis import align, rms
import warnings
import io
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    import Bio
    from Bio import PDB
    from Bio.PDB.PDBParser import PDBParser
warnings.filterwarnings("ignore", category=UserWarning)

test=True
debug=False

class PDBPrepWarning(Warning):
    """Custom class for specific warnings that may happen during execution of the PDBPrep class"""
    pass

class PDBPrep:
    def read_pyDISH(self):
        data = pd.read_csv('tables/pyDISH.csv')
        data = data.drop('Unnamed: 0', axis=1)
        if test: data = data.iloc[44:]
        print("Data loaded from pyDISH.csv with shape:", data.shape)
        print("First few rows of data:")
        print(data.head())
        return data

    def read_pdb(self, pdb_id):
        filename = PDB.PDBList().retrieve_pdb_file(pdb_id, pdir='PDB', file_format='pdb')
        filename = ops.convert_ent_to_pdb(file=f"PDB/pdb{pdb_id}.ent", output=f"PDB/{pdb_id}.pdb")
        filename = f"PDB/{pdb_id}.pdb"
        universe  = mda.Universe(filename)
        heme = self.isolate_heme(output=f"PDB/{pdb_id}_only_heme.pdb", pdb_id=pdb_id)
        if self.debug_read_pdb:
            print("Heme atoms isolated by self.read_pdb: ")
            print(heme.atoms)
        ops.HETATM_to_ATM(file=f"PDB/{pdb_id}_only_heme.pdb", output=f"PDB/{pdb_id}_heme_atom.pdb")
        self.isolate_sides(pdb_id)
        return universe, heme

    def isolate_heme(self, output, pdb_id):
        if not os.path.exists(f"PDB/{pdb_id}.pdb"):
            raise FileNotFoundError(f"File PDB/{pdb_id}.pdb not found but should exist. ")
        heme_name_list = ["HEM", "HEA", "HEB", "HEC", "HEO"]
        lines, lines_ = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}.pdb"), []
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            if line_dict["resname"].strip() in heme_name_list:
                line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
                lines_.append(line_)
        pdb_tools.files.write_file(file=output, lines=lines_)
        heme = mda.Universe(f"PDB/{pdb_id}_only_heme.pdb")
        return heme

    def isolate_sides(self, pdb_id):
        lines_ = []
        if not os.path.exists(f"PDB/{pdb_id}_heme_atom.pdb"):
            raise FileNotFoundError(f"File 'PDB/{pdb_id}_heme_atom.pdb' not found, but should exist. ")
        lines = pdb_tools.files.read_file(f"PDB/{pdb_id}_heme_atom.pdb")
        sidechain_atomtypes = ["CMA", "CMB", "CMC", "CMD", "CAB", "CBB", "CAC", "CBC", "CAA", "CAD"]
        # Add CAA and CAD atom types from propionates to save the orientation, even if we cut props down to methyl groups.
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            if line_dict["atom_name"].strip() in sidechain_atomtypes:
                line_dict["chainID"] = "H"
                line_dict["resi_no"] = "   1"
                #if line_dict["atom_name"].strip() in ["CAA", "CAD"]: -> remember to change rtf file so that the formerly propionates are recognized as CT3 methyls.
                line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
                lines_.append(line_)
        pdb_tools.files.write_file(file=f"PDB/{pdb_id}_heme_sides.pdb", lines=lines_)
        if not os.path.exists(f"PDB/{pdb_id}_heme_sides.pdb"):
            raise FileNotFoundError(f"File 'PDB/{pdb_id}_heme_sides.pdb' was not found, but should exist. ")
        return

    def get_heme_normal(self, heme):
        v1 = heme.atoms.positions[1] - heme.atoms.positions[0]
        v2 = heme.atoms.positions[2] - heme.atoms.positions[0]
        plane_normal = np.cross(v1, v2)
        plane_normal /= np.linalg.norm(plane_normal)
        return plane_normal

    def center_of_mass_(self, residue):
        com = residue.atoms.center_of_mass()
        return com

    def calculate_angle(self, pos1, pos2, pos3):
        vector1 = pos1 - pos2
        vector3 = pos3 - pos2
        dot_product = np.dot(vector1, vector3)
        magnitude1 = np.linalg.norm(vector1)
        magnitude3 = np.linalg.norm(vector3)
        cos_theta = dot_product / (magnitude1 * magnitude3)
        theta_rad = np.arccos(cos_theta)
        theta_deg = np.degrees(theta_rad)
        return theta_deg

    def find_closest_angle_to_180_deg(self, angles):
        angle_diffs = np.abs(np.array(angles) - 180.0)
        closest_index = np.argmin(angle_diffs)
        closest_angle = angles[closest_index]
        return closest_angle, closest_index

    def find_closest_angles_to_180_deg(self, angles):
        angle_diffs = np.abs(np.array(angles) - 180.0)
        closest_indices = np.argsort(angle_diffs)[:2]
        closest_angles = [angles[i] for i in closest_indices]
        return closest_angles, closest_indices

    def create_dummy_axials(self, heme_iron, plane_normal):
        #do not init atom, just give back position, because the atom does not exist and we can treat them as only positions
        dum_pos_n = heme_iron.positions[0] + plane_normal
        dum_neg_n = heme_iron.positions[0] - plane_normal
        if self.debug_get_axials_geometric:
            print("dummy position +n: ", dum_pos_n)
            print("dummy_position -n: ", dum_neg_n)
        return dum_pos_n, dum_neg_n

    def find_single_axial(self, pdb_id, axials, axial1_pydish, axial2_pydish, heme_iron, dum_pos_n, dum_neg_n):
        print("start single axial")
        if self.debug_get_axials_geometric:
            print(f"structure {pdb_id} has only one axial residue in the pyDISH database. In the table are: {axial1_pydish} and {axial2_pydish}. ")
        coms = []
        for residue in axials.residues:
            com = self.center_of_mass_(residue)
            coms.append(com)
        if self.debug_get_axials_geometric:
            print("axials.residues: ", axials.residues)
            print("centers of mass: ", coms)
            print("dum_pos_n: ", dum_pos_n, "\ndum_pos_n type:", type(dum_pos_n))
            print("dum_neg_n: ", dum_neg_n, "\ndum_neg_n type:", type(dum_neg_n))
            #print("heme_iron.atoms.positions: ", heme_iron.atoms.positions, "\nheme_iron.atoms.positions type: ", type(heme_iron.atoms.positions))
            print("heme_iron.atoms.positions[0]: ", heme_iron.atoms.positions[0], "\nheme_iron.atoms.positions[0] type: ", type(heme_iron.atoms.positions[0]))
        for com in coms:
            dum_pos_angle = self.calculate_angle(pos1=com, pos2=heme_iron.atoms.positions[0], pos3=dum_pos_n)
            dum_neg_angle = self.calculate_angle(pos1=com, pos2=heme_iron.atoms.positions[0], pos3=dum_neg_n)
        angles = [dum_pos_angle, dum_neg_angle]
        closest_angle, closest_index = self.find_closest_angle_to_180_deg(angles)
        axials = axials.residues[closest_index]
        if self.debug_get_axials_geometric:
            print("closest angle to 180 degrees: ", closest_angle)
            print("closest_index is: ", closest_index)
            print(axials)
        return axials

    def add_dummy_to_axials(self, axials, dum_pos_n, dum_neg_n):
        axials_dum, dummies = axials, mda.Universe.empty(n_atoms=2, n_residues=2, n_segments=2, atom_resindex=[0,1], residue_segindex=[0,1], trajectory=True, velocities=False, forces=False)
        if self.debug_get_axials_geometric:
            #print("len(dummies)  =             ", len(dummies))
            print("dummies       =             ", dummies)
            print("dummies.atoms =             ", dummies.atoms)
        dummies.add_TopologyAttr("name", ["DUP", "DUN"])
        if self.debug_get_axials_geometric:
            print("dummies.atoms.names =      ", dummies.atoms.names)
        dummies.add_TopologyAttr("resname", ["DUP", "DUN"])
        if self.debug_get_axials_geometric:
            print("dummies.residues =         ", dummies.residues)
        dummies.atoms.positions = [dum_pos_n, dum_neg_n]
        if self.debug_get_axials_geometric:
            print("dummies.atoms.positions =  ", dummies.atoms.positions)
        if self.debug_get_axials_geometric:
            print("dummies.atoms    = ", dummies.atoms)
            print("axials_dum.atoms = ", axials_dum.atoms)
        axials_dum = mda.Merge(axials_dum.atoms, dummies.atoms)
        resids, residue_map = [], {}
        for atom in axials_dum.atoms:
            residue_name = atom.resname
            if residue_name not in residue_map:
                residue_map[residue_name] = len(residue_map) + 1
            resids.append(len(residue_map))
            #print(residue_map)
            #print(resids)
        axials_dum.add_TopologyAttr("resid", [1,2,3,4])
        #axials_dum.add_TopologyAttr("resid", resids)
        if self.debug_get_axials_geometric:
            print("merged Universe of axial choices and dummies: ", axials_dum)
            print("merged Universe atoms:                        ", axials_dum.atoms)
            print("merged Universe atoms positions:              ", axials_dum.atoms.positions)
            print("merged Universe residues:                     ", axials_dum.residues)
            print("merged Universe segments:                     ", axials_dum.segments)
        return axials_dum

    def distance_warnings(self, shortest_distances):
        d12, d13, d14, d23, d24, d34 = shortest_distances[0], shortest_distances[1], shortest_distances[2], shortest_distances[3], shortest_distances[4], shortest_distances[5]
        if self.debug_get_axials_geometric:
            print(d12[2])
        if d12[2] < 4:
            warnings.warn("d(ax1-ax2) is smaller than 4 Angstrom. Please check structure for correctness. ", PDBPrepWarning)
        if d34[2] < 1:
            warnings.warn("d(dum1-dum2) is smaller than 1Angstrom. Dummy atoms are within the VdW radius of the iron, the axial choice may be wrong. ", PDBPrepWarning)
        if d13[2] and d23[2] <1.5:
            warnings.warn("d(ax1-dum1) and d(ax2-dum1) are smaller than 1.5 Angstrom. This may indicate several serious issues with the choice of axials.", PDBPrepWarning)
        if d14[2] and d24[2] <1.5:
            warnings.warn("d(ax1-dum2) and d(ax2-dum2) are smaller than 1.5 Angstrom. This may indicate several serious issues with the choice of axials.", PDBPrepWarning)
        return

    def calculate_shortest_distances(self, axials, dum_pos_n, dum_neg_n):
        """
        self.calculate shortest distances takes from the self instance of the class the debug option self.debug_get_axials_geometric
        and the mda.Universe type object containing the axials and the dummy atoms: axials_dum.
        distances calculated are (in order): [d(ax1-ax2), d(ax1-dum1), d(ax1-dum2), d(ax2-dum1), d(ax2-dum2), d(dum1-dum2)].
        """
        axials_dum = self.add_dummy_to_axials(axials, dum_pos_n, dum_neg_n)
        if self.debug_get_axials_geometric:
            print("merged Universe of axial choices and dummies: ", axials_dum)
            print("merged Universe atoms:                        ", axials_dum.atoms)
            print("merged Universe atoms positions:              ", axials_dum.atoms.positions)
            print("merged Universe residues:                     ", axials_dum.residues)
            print("merged Universe segments:                     ", axials_dum.segments)
            for i in [1,2,3,4]:
                print(axials_dum.select_atoms(f"resid {i}"))
        shortest_distances, group_selections = [], ["resid 1", "resid 2", "resname DUP", "resname DUN"]
        for i in range(len(group_selections)):
            for j in range(i + 1, len(group_selections)):
                group1 = axials_dum.select_atoms(group_selections[i])
                if self.debug_get_axials_geometric:
                    print("group1:     ", group1)
                group2 = axials_dum.select_atoms(group_selections[j])
                if self.debug_get_axials_geometric:
                    print("group2:     ", group2)
                distances = distance_array(group1.positions, group2.positions)
                if self.debug_get_axials_geometric:
                    print("distances:     ", distances)
                shortest_distance = distances.min()
                if self.debug_get_axials_geometric:
                    print(f"shortest distance between residue {group1.resids[0]} and residue {group2.resids[0]} =    ", shortest_distance)
                shortest_distances.append((group_selections[i], group_selections[j], shortest_distance))
                if self.debug_get_axials_geometric:
                    print("Shortest distances between the groups:     ", shortest_distances)
        self.distance_warnings(shortest_distances)
        return shortest_distances

    def find_double_axial(self, axials, heme_iron, dum_pos_n, dum_neg_n):
        coms, angles = [], []
        for residue in axials.residues:
            com = self.center_of_mass_(residue)
            coms.append(com)
        if self.debug_get_axials_geometric:
            print("axials.residues: ", axials.residues)
            print("centers of mass: ", coms)
        for com in coms:
            dum_pos_angle = self.calculate_angle(pos1=com, pos2=heme_iron.atoms.positions[0], pos3=dum_pos_n)
            dum_neg_angle = self.calculate_angle(pos1=com, pos2=heme_iron.atoms.positions[0], pos3=dum_neg_n)
            angles.append(dum_pos_angle)
            angles.append(dum_neg_angle)
        if self.debug_get_axials_geometric:
            print("angles are: ", angles)
        closest_angles, closest_indices = self.find_closest_angles_to_180_deg(angles)
        # check if the indices are one odd and one even, to know that they are on opposite sides of the heme-plane,
        # the Assertion error will be raised if the following expression is not true.
        assert np.logical_xor(closest_indices % 2 == 1, closest_indices % 2 == 0).any(), "We have found two axials on the same side to be the ones closest to 180 degrees, based on position in the list.\n All even entries should be axial1, all odd entries should be axial2."
        if self.debug_get_axials_geometric:
            print("closest angles to 180 degrees: ", closest_angles, "\ntype of closest_angles is: ", type(closest_angles))
            print("closest_indices are: ", closest_indices, "\ntype of closest_indices is: ", type(closest_indices))
        axials = axials.residues[closest_indices//2]
        distances = self.calculate_shortest_distances(axials, dum_pos_n, dum_neg_n)
        print(distances)
        print("axials.resnames: ", axials.resnames)
        if self.debug_get_axials_geometric:
            print(axials)
        return axials

    def get_axials_geometric(self, pdb_id, universe, heme, axials, axial1_pydish, axial2_pydish):
        heme_iron = universe.select_atoms("type FE and (resname HEM or resname HEC or resname HEA or resname HEB or resname HEO)")
        plane_normal = self.get_heme_normal(heme)
        dum_pos_n, dum_neg_n = self.create_dummy_axials(heme_iron, plane_normal)
        if self.debug_get_axials_geometric:
            print("heme plane normal vector: ", plane_normal)
        #if somehow there is no axial in the pyDISH table
        if isinstance(axial1_pydish, float) and isinstance(axial2_pydish, float):
            raise ValueError(f"Both axials values are NaN which is unlikely, the entry {pdb_id} needs to be manually checked. ")
        print("passed check for no-axials, moving on to check for single axial. ")
        #if there is one axial in the pyDISH table and we want to find out which one it is
        # we filtered first by axial name, yet if this proves to still be non-unique, we can use this
        if isinstance(axial1_pydish, float) or isinstance(axial2_pydish, float):
            if self.debug_get_axials_geometric:
                print("pydish axials: ", (axial1_pydish, axial2_pydish))
                print("found axials (one info per atom due to mda-resnames method going through all atoms): ", (axials[:].resnames))
            if len(axials.residues) == 1 and axials[0].resname in (axial1_pydish, axial2_pydish):
                #only one axial found and it has the right residue name
                if self.debug_get_axials_geometric:
                    print("returning from axial advanced search because axial residue satisfies conditions: \n only one axial and axial resname equal to the pydish axial resname.")
                return
            else:
                #only one axial in the pydish table and more than one axial found with the mda atom-selection
                axials = self.find_single_axial(pdb_id, axials, axial1_pydish, axial2_pydish, heme_iron, dum_pos_n, dum_neg_n)
        #if there are two axials in the pyDISH table and we are not sure which axials are the axials
        else:
            print("passed check for single axial, moving on to check for double axials. ")
            if self.debug_get_axials_geometric:
                print("pydish axials: ",(axial1_pydish, axial2_pydish))
                print("found axials (one info per atom due to mda-resnames method): ",(axials[:].resnames))
                print("found axials with MDA: ", axials.residues)
            if len(axials.residues) == 2 and (axial1_pydish, axial2_pydish) == (axials[0].resname, axials[1].resname):
                # if two axials are found and the axials have the same residue name as the residues in the pyDISH table
                if self.debug_get_axials_geometric:
                    print("returning from axial advanced search because axial residues satisfy conditions: \n only two axials and axial resnames equal to the pydish axial resnames.")
                return
            else:
                # if more than two axials are found for a structure in which we expect two axials
                axials = self.find_double_axial(axials, heme_iron, dum_pos_n, dum_neg_n)
        return axials, plane_normal, dum_pos_n, dum_neg_n, heme_iron 

    def get_axial(self, pdb_id, universe, heme, axial1_pydish, axial2_pydish):
        self.verbose, self.cleanup = False, True
        axials = universe.select_atoms("same resid as around 4 (type FE and (resname HEM or resname HEC or resname HEA or resname HEB or resname HEO))")
        axials = axials.select_atoms("not (resname HEM or resname HEC or resname HEA or resname HEB or resname HEO)")
        axials = axials.select_atoms(f"resname {axial1_pydish} or resname {axial2_pydish}")
        axials = axials.select_atoms("not backbone")
        if self.debug_get_axial:
            print(f"Number of residues in axials: {axials.n_residues}.")
        axials, plane_normal, dum_pos_n, dum_neg_n, heme_iron = self.get_axials_geometric(pdb_id, universe, heme, axials, axial1_pydish, axial2_pydish)
        try:
            axials.atoms.write(f"PDB/{pdb_id}_axials.pdb")
            print(f"wrote file PDB/{pdb_id}_axials.pdb")
        except IndexError:
            print("IndexError when writing axials with MDAnalysis, writing empty file instead. ")
            with mda.Writer(f"PDB/{pdb_id}_axials.pdb", bonds=None, n_atoms=0, universe=axials) as pdb_writer:
                pdb_writer.write(axials)
        return axials

    def __try_clean(self, file):
        try:
            os.remove(file)
        except FileNotFoundError:
            print(f"could not remove file {file}, FileNotFound. ")
            pass
        return

    def move_files(self, pdb_id):
        file_list, source_directory, destination_directory = os.listdir("PDB/"), "PDB/", f"PDB/{pdb_id}"
        if os.path.exists(destination_directory):
            shutil.rmtree(destination_directory)
        os.mkdir(destination_directory)
        for filename in file_list:
            if pdb_id in filename:
                source_path = os.path.join(source_directory, filename)
                destination_path = os.path.join(destination_directory, filename)
                if os.path.abspath(source_path) != os.path.abspath(destination_path):
                    try:
                        shutil.move(source_path, destination_path)
                        print(f"Moved {filename} from {source_path} to {destination_path}")
                    except Exception as e:
                        print(f"Error moving {filename}: {e}")
                else:
                    print(f"Skipping move for {filename}: Source and destination paths are the same.")

    def cleaner(self, pdb_id):
        self.__try_clean(file=f"PDB/{pdb_id}.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_only_heme.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_atom.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_prad.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_axials.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_sorted.xyz")
        self.__try_clean(file=f"PDB/{pdb_id}_heme.xyz")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_fixed.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_from_wts.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_married.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_sides.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_types.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_AXI1.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_AXI2.pdb")
        self.__try_clean(file=f"PDB/{pdb_id}_heme_HEME.pdb")
        #self.__try_clean(file=f"PDB/{pdb_id}_charmm.inp")
        #self.__try_clean(file=f"PDB/{pdb_id}_system.pdb")
        self.__try_clean(file=f"rmsfit_{pdb_id}_g16.xyz")
        self.__try_clean(file=f"rmsfit_{pdb_id}_system_protonated.xyz")
        self.move_files(pdb_id)
        return

    def fix_heme(self, pdb_id):
        lines_ = []
        ops.xyz_to_pdb(file = f"PDB/{pdb_id}_heme_sorted.xyz" ,output = f"PDB/{pdb_id}_heme_from_wts.pdb")
        lines = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_heme_from_wts.pdb")
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            line_dict["resi_no"] = "   1"
            #the following if condition is problematically hacky, and needs to be changed depending on the updates of pdb_tools.py by Nereu! -> open up Issue in pdb_tools.py repo. Also, pdbparser package???
            if line_dict["atom_name"].strip() == "FE":
                line_dict["elem_symb"] = "FE"
                line_dict["segment"] = "   "
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        pdb_tools.files.write_file(lines=lines_, file=f"PDB/{pdb_id}_heme_fixed.pdb")
        return

    def fix_heme_atomtypes(self, pdb_id):
        lines, lines_ = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_heme_fixed.pdb"), []
        atomtype_dict = {
            1  : " FE ",
            2  : " CHD",
            3  : " C4C",
            4  : " C3C",
            5  : " NC ",
            6  : " C2C",
            7  : " C1C",
            8  : " CHC",
            9  : " C4B",
            10 : " C3B",
            11 : " NB ",
            12 : " C2B",
            13 : " C1B",
            14 : " CHB",
            15 : " C4A",
            16 : " C3A",
            17 : " NA ",
            18 : " C2A",
            19 : " C1A",
            20 : " CHA",
            21 : " C4D",
            22 : " C3D",
            23 : " ND ",
            24 : " C2D",
            25 : " C1D"
        }
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            line_dict["atom_name"] = atomtype_dict[int(line_dict["serial_no"].strip())]
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        pdb_tools.files.write_file(file=f"PDB/{pdb_id}_heme_types.pdb", lines=lines_)
        return

    def marry_heme(self, pdb_id):
        pdb_tools.operations.fuse_segments(pdb_files=[f"PDB/{pdb_id}_heme_types.pdb", f"PDB/{pdb_id}_heme_sides.pdb"], pdb_output=f"PDB/{pdb_id}_heme_married.pdb")
        lines = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_heme_married.pdb")
        lines_ = []
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            if line_dict["resname"].strip() in {"HEC", "HEA", "HEO", "HEB"}:
                line_dict["resname"] = "HEM "
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        pdb_tools.files.write_file(lines=lines_, file=f"PDB/{pdb_id}_heme_married.pdb")
        return

    def __find_resnames(self, pdb_id, lines):
        residues, seen = [], set()
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            # Create a unique combination of resname and resi_no as a tuple
            res_tuple = (line_dict["resname"], line_dict["resi_no"])
            if self.debug_finders:
                print("res_tuple =     ", res_tuple)
            residues.append(res_tuple)
        seen_add = seen.add
        if self.debug_finders:
            print("seen_add =      ", seen_add)
        # Filter for unique combinations
        residues = [k for k in residues if not (k in seen or seen_add(k))]
        if self.debug_finders:
            print(f"Unique combinations of resnames and resi_no's in file PDB/{pdb_id}_system.pdb before final corrections:", residues)
            print("len(residues) =      ", len(residues))
        # Here we need again a distinction between the three potential numbers of axials: 0,1,2.
        # This is easier than writing a function that deals with all three cases well without a case-distinction.
        # len(residues must be no(axials)+1 due to heme being counted in the search for all unique residues.)
        axial_list = []
        if len(residues) == 1:
            raise ValueError(f"no axials found in PDB/{pdb_id}_system.pdb.")
        if len(residues) == 2:
            axial1_ = residues[1]
            axial1 = axial1_[0]
            axial_list.append(axial1)
            if self.debug_finders:
                print(f"Found one axial {axial1} of datatype {type(axial1)}.")
        if len(residues) == 3:
            axial1_ = residues[1]
            axial1 = axial1_[0]
            axial2_ = residues[2]
            axial2 = axial2_[0]
            axial_list.append(axial1)
            axial_list.append(axial2)
            if self.debug_finders:
                print(f"Found one axial {axial1} of datatype {type(axial1)} and\none axial {axial2} of datatype {type(axial2)}")
        resi_no = len(residues)
        return axial_list, residues, axial1, axial2, resi_no

    def __find_segments(self, pdb_id, lines):
        segments, seen = [], set()
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            segments.append(line_dict["segment"])
        seen_add = seen.add
        segments = [k for k in segments if not (k in seen or seen_add(k))]
        if self.debug_finders:
            print(f"segments in file PDB/{pdb_id}_system.pdb after final corrections: ")
            print(segments)
        try:
            axial1_segname, axial2_segname = segments[1], segments[2]
        except IndexError:
            axial1_segname, axial2_segname = segments[1], None
        return axial1_segname, axial2_segname, segments

    def resi_no_resname(self, line_dict, resi_nos_unique, resi_no):
        line_dict["chainID"], line_dict["segment"], line_dict["atom"] = "H", "HEME", "ATOM  "
        if line_dict["resname"] in ["HEA", "HEB", "HEC", "HEO"]:
            line_dict["resname"] = "HEM"
        if line_dict["resname"] == "HEM":
            line_dict["resi_no"] = "   1"
        if self.debug_files:
            print(line_dict)
        if line_dict["resi_no"].strip() == resi_nos_unique[1]:
            line_dict["resi_no"] = "   2"
        # add here to check for double axial structures, because only they have a second axial (third residue).
        if resi_no == 3:
            if line_dict["resi_no"].strip() == resi_nos_unique[2]:
                line_dict["resi_no"] = "   3"
        return line_dict

    def HIS_HSD(self, line_dict, axial1, axial2, resi_no):
        if line_dict["resname"] == axial1 == "HIS ":
            line_dict["resname"] = "HSD "
        if resi_no == 3:
            if line_dict["resname"] == axial2 == "HIS ":
                line_dict["resname"] = "HSD "
        return line_dict

    def XH2_OH2(self, line_dict, axial1, axial2, resi_no):
        if (line_dict["atom_name"].strip() == "O") and (line_dict["resname"] == "HOH "):
            line_dict["atom_name"] = " XH2"
        if (line_dict["atom_name"].strip() == "O") and (line_dict["resname"].strip() == axial1 == "O"):
            line_dict["atom_name"] = " XH2"
        if resi_no == 3:
            if (line_dict["atom_name"].strip() == "O") and (line_dict["resname"].strip() == axial2 == "O"):
                line_dict["atom_name"] = " XH2"
        if line_dict["resname"] == axial1 == "O   ":
            line_dict["resname"] = "HOH "
        if resi_no == 3:
            if line_dict["resname"] == axial2 == "O   ":
                line_dict["resname"] = "HOH "
        return line_dict

    def fix_fe_atomtype(self, line_dict):
        if line_dict["atom_name"].strip() == "FE":
                line_dict["elem_symb"] = "FE"
        return line_dict

    def assign_segments(self, pdb_id):
        lines, lines_ = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_system.pdb"), []
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            if line_dict["atom_name"].strip() == "FE":
                line_dict["elem_symb"] = "FE"
            if line_dict["resi_no"].strip() == "2":
                line_dict["segment"] = "AXI1"
            if line_dict["resi_no"].strip() == "3":
                line_dict["segment"] = "AXI2"
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        pdb_tools.files.write_file(file=f"PDB/{pdb_id}_system.pdb", lines=lines_)
        lines = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_system.pdb")
        return lines

    def marry_pdbs(self, pdb_id):
        pdb_tools.operations.fuse_segments(pdb_files=[f"PDB/{pdb_id}_heme_married.pdb", f"PDB/{pdb_id}_axials.pdb"], pdb_output=f"PDB/{pdb_id}_system.pdb")
        lines = pdb_tools.files.read_file(pdb_file=f"PDB/{pdb_id}_system.pdb")
        lines_ = []
        axial_list, residues, axial1, axial2, resi_no = self.__find_resnames(pdb_id, lines)
        #check for resi_nos so that we do not mix up two axials into one when they have the same resname
        resi_nos_unique = []
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            if line_dict["resi_no"].strip() not in resi_nos_unique:
                resi_nos_unique.append(line_dict["resi_no"].strip())
        self.debug_files = False
        if self.debug_files:
            print("resi_nos_unique: ", resi_nos_unique)
        for line in lines:
            line_dict = pdb_tools.line_operations.read_pdb_line(line=line)
            line_dict = self.resi_no_resname(line_dict, resi_nos_unique, resi_no)
            line_dict = self.HIS_HSD(line_dict, axial1, axial2, resi_no)
            line_dict = self.XH2_OH2(line_dict, axial1, axial2, resi_no)
            line_dict = self.fix_fe_atomtype(line_dict)
            line_ = pdb_tools.line_operations.create_line(line_dict=line_dict)
            if self.debug_files:
                print(line_)
            lines_.append(line_)
        pdb_tools.files.write_file(file=f"PDB/{pdb_id}_system.pdb", lines=lines_)
        pdb_tools.operations.renumber(pdb_file=f"PDB/{pdb_id}_system.pdb", pdb_file_output=f"PDB/{pdb_id}_system.pdb")
        lines = self.assign_segments(pdb_id)
        axial1_segname, axial2_segname, segments = self.__find_segments(pdb_id=pdb_id, lines=lines)
        return axial_list, residues, axial1, axial2, resi_no, axial1_segname, axial2_segname, segments

    def split_system(self, pdb_id):
        try:
            pdb_tools.operations.split_segments(pdb_file=f"PDB/{pdb_id}_system.pdb", segnames=["HEME", "AXI1", "AXI2"], pdb_id=pdb_id)
        except IndexError:
            pdb_tools.operations.split_segments(pdb_file=f"PDB/{pdb_id}_system.pdb", segnames=["HEME", "AXI1"], pdb_id=pdb_id)
        return

    def __add_patches_ax1(self, axial1):
        patches = []
        if axial1 == "HIS ":
            patches.append("\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "HOH ":
            patches.append("\nPATCH POSH HEME 1 AXI1 2 SETUP\n")
        if axial1 == "CYS ":
            patches.append("\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "MET ":
            patches.append("\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "OXY ":
            patches.append("\n")
        if axial1 == "TYR ":
            patches.append("\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "GLU ":
            patches.append("\nPATCH GLUP AXI1 2 SETUP\nPATCH HIS2 AXI1 2 SETUP\n PATCH DELH AXI1 2 SETUP\n")
        if axial1 == "PER ":
            patches.append("\n")
        if axial1 == "LYS ":
            patches.append("\nPATCH DELK AXI1 2 SETUP\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "ARG ":
            patches.append("\nPATCH DELR AXI1 2 SETUP\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        if axial1 == "TRP ":
            patches.append("\nPATCH HIS2 AXI1 2 SETUP\nPATCH DELH AXI1 2 SETUP\n")
        patches.append("\nHBUILD SELE SEGID AXI1 .and. HYDROGEN END \nCOOR PRINT SELE SEGID AXI1 .and. .not. INIT END\n")
        if self.debug_add_patches_ax1:
            print(patches)
        return patches

    def __add_patches_ax2(self, axial2):
        patches = []
        if axial2 == "HIS ":
            patches.append("\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "HOH ":
            patches.append("\nPATCH POSH HEME 1 AXI2 3 SETUP\n")
        if axial2 == "CYS ":
            patches.append("\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "MET ":
            patches.append("\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "OXY ":
            patches.append("\n")
        if axial2 == "TYR ":
            patches.append("\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "GLU ":
            patches.append("\nPATCH GLUP AXI2 3 SETUP\nPATCH HIS2 AXI2 3 SETUP\n PATCH DELH AXI2 3 SETUP\n")
        if axial2 == "O   ":
            patches.append("\n")
        if axial2 == "PER ":
            patches.append("\n")
        if axial2 == "LYS ":
            patches.append("\nPATCH DELK AXI2 3 SETUP\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "ARG ":
            patches.append("\nPATCH DELR AXI2 3 SETUP\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        if axial2 == "TRP ":
            patches.append("\nPATCH HIS2 AXI2 3 SETUP\nPATCH DELH AXI2 3 SETUP\n")
        patches.append("\nHBUILD SELE SEGID AXI2 .and. HYDROGEN END \nCOOR PRINT SELE SEGID AXI2 .and. .not. INIT END\n")
        if self.debug_add_patches_ax2:
            print(patches)
        return patches

    def create_charmm_file(self, pdb_id, axial1, axial2, axial1_segname, axial2_segname, resi_no):
        if self.debug_create_charmm_file:
            try:
                print("pdb_id =         ", pdb_id, "\naxial1 =         ", axial1, "\naxial2 =         ", axial2, "\naxial1_segname = ", axial1_segname, "\naxial2_segname = ", axial2_segname)
                print("printed without TypeError, so both axial1 and axial2 are known")
            except TypeError:
                print("TypeError printing axial1 and axial2, continuing with only printing axial1 info.")
                print("pdb_id =         ", pdb_id, "\naxial1 =         ", axial1, "\naxial1_segname = ", axial1_segname)
            except AttributeError:
                print("TypeError printing axial1 and axial2, continuing with only printing axial1 info.")
                print("pdb_id =         ", pdb_id, "\naxial1 =         ", axial1, "\naxial1_segname = ", axial1_segname)
        with open("textfiles/charmm_header.inp") as header_file:
            header = header_file.readlines()
        with open("textfiles/charmm_heme.inp") as heme_file:
            heme = heme_file.readlines()
        with open("textfiles/charmm_non_water_axial1.inp") as non_water_axial_file:
            non_water_axial1 = non_water_axial_file.readlines()
        with open("textfiles/charmm_water_axial1.inp") as water_axial_file:
            water_axial1 = water_axial_file.readlines()
        with open("textfiles/charmm_non_water_axial2.inp") as non_water_axial_file:
            non_water_axial2 = non_water_axial_file.readlines()
        with open("textfiles/charmm_water_axial2.inp") as water_axial_file:
            water_axial2 = water_axial_file.readlines()
        with open("textfiles/charmm_ending.inp") as ending_file:
            ending = ending_file.readlines()
        charmm_file_strings = []
        charmm_file_strings.append(header)
        charmm_file_strings.append(heme)
        # axial1
        if axial1 != "HOH ":
            charmm_file_strings.append(non_water_axial1)
        if axial1 == "HOH ":
            charmm_file_strings.append(water_axial1)
        patches_ax1 = self.__add_patches_ax1(axial1)
        charmm_file_strings.append(patches_ax1)
        # axial2
        if resi_no == 3: # again only add this part if there are two axials (= three overall residues)
            if axial2 == None:
                pass
            if axial2 != "HOH " and axial2 is not None:
                charmm_file_strings.append(non_water_axial2)
                patches_ax2 = self.__add_patches_ax2(axial2)
                charmm_file_strings.append(patches_ax2)
            if axial2 == "HOH ":
                charmm_file_strings.append(water_axial2)
                patches_ax2 = self.__add_patches_ax2(axial2)
                charmm_file_strings.append(patches_ax2)
        charmm_file_strings.append(ending)
        #convert list of lists to list of strings with list comprehension
        if self.debug_create_charmm_file:
            print(charmm_file_strings)
        charmm_file_strings = [element for sublist in charmm_file_strings for element in sublist]
        if self.debug_create_charmm_file:
            print(charmm_file_strings)
        #fill in information because you can not turn into f-strings and fill in the info automatically
        if resi_no == 3:
            charmm_file_strings = [f_string.format(pdb_id=pdb_id, axial1=axial1, axial2=axial2, axial1_segname=axial1_segname, axial2_segname=axial2_segname) for f_string in charmm_file_strings]
            with open(f"PDB/{pdb_id}_charmm.inp", "w") as writefile:
                writefile.writelines(charmm_file_strings)
        elif resi_no == 2:
            charmm_file_strings = [f_string.format(pdb_id=pdb_id, axial1=axial1, axial1_segname=axial1_segname) for f_string in charmm_file_strings]
            with open(f"PDB/{pdb_id}_charmm.inp", "w") as writefile:
                writefile.writelines(charmm_file_strings)
        else:
            raise AttributeError(f"Structure {pdb_id} has no axials at CHARMM stage, something went wrong. ")
        return

    def check_charmm_logfile(self, pdb_id):
        multiline_string = """                            /---------\\
                           /           \\
                          /             \\
                         /               \\
                         !  XXXX   XXXX  !
                         !  XXXX   XXXX  !
                         !  XXX     XXX  !
                         !       X       !
                          --\   XXX   /--
                           ! !  XXX  ! !
                           ! !       ! !
                           ! I I I I I !
                           !  I I I I  !
                            \         /
                             --     --
                               \---/
                        XXX             XXX
                       XXXX             XXXX
                       XXXXX           XXXXX
                          XXX         XXX
                            XXX     XXX
                               XXXXX
                              XXX XXX
                            XXX     XXX
                          XXX         XXX
                       XXXXX           XXXXX
                       XXXX             XXXX
                        XXX             XXX       """
        try:
            with open(f"PDB/{pdb_id}_charmm.out", 'r') as file:
                file_content = file.read()
                if multiline_string in file_content:
                    return print("CHARMM run completed without the skull. ")
        except FileNotFoundError:
            return print(f"CHARMM run completed with the skull. Please recheck the structure, with help from the logfile PDB/{pdb_id}/{pdb_id}.out")

    def run_charmm(self, pdb_id):
        cmd = ['charmm', f"<PDB/{pdb_id}_charmm.inp", f">PDB/{pdb_id}_charmm.out"]
        #cmd = ['/Users/jejo/Desktop/charmm/charmm', f"</Users/jejo/Desktop/Work/heme/PDB/{pdb_id}_charmm.inp", f">/Users/jejo/Desktop/Work/heme/PDB/{pdb_id}_charmm.out"]
        process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #process = subprocess.Popen(' '.join(cmd), shell=True, executable="/usr/local/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(' '.join(cmd))
        stdout, stderr = process.communicate()
        if self.debug_create_charmm_file:
            if process.returncode == 0:
                print("CHARMM run completed successfully")
            else:
                print(f"CHARMM run failed with return code {process.returncode}")
                print("Standard Output:")
                print(stdout.decode())
                print("Standard Error:")
                print(stderr.decode())
        self.check_charmm_logfile(pdb_id)
        return

    def create_xtb_file(self, pdb_id, axial1, axial2, axial1_segname, axial2_segname, resi_no):
        if self.debug_create_xtb_file:
            try:
                print("pdb_id =         " + pdb_id)
                print("axial1 =         " + axial1)
                print("axial1_segname = " + axial1_segname)
                if resi_no == 3:
                    print("axial2 =         " + axial2)
                    print("axial2_segname = " + axial2_segname)
            except TypeError:
                print("pdb_id =         " + pdb_id)
                print("axial1 =         " + axial1)
                print("axial1_segname = " + axial1_segname)
        with open("textfiles/xtb_heme.inp") as heme_file:
            heme = heme_file.readlines()
        with open("textfiles/xtb_water_axial1.inp") as water_axial1_file:
            water_axial1 = water_axial1_file.readlines()
        with open("textfiles/xtb_ending.inp") as ending_file:
            ending = ending_file.readlines()
        xtb_file_strings = []
        xtb_file_strings.append(heme)
        if axial1 == "HOH ":
            xtb_file_strings.append(water_axial1)
        if axial1 != "HOH ":
            pass
        if resi_no == 3:
            if axial2 != "HOH ":
                pass
            if axial2 == "HOH ":
                with open(f"PDB/{pdb_id}_system_protonated.xyz", 'r') as input_file:
                    n = int(input_file.readline().strip())
                n_minus_2 = n - 2 # Calculate n-2, n-1, n
                n_minus_1 = n - 1
                new_text = f"\n$fix\n    atoms: {n_minus_2}, {n_minus_1}, {n} "
                if self.debug_create_xtb_file:
                    print(new_text)
                xtb_file_strings.append(new_text)
        xtb_file_strings.append(ending)
        if self.debug_create_xtb_file:
            print(xtb_file_strings)
        xtb_file_strings = [element for sublist in xtb_file_strings for element in sublist]
        if self.debug_create_xtb_file:
            print(xtb_file_strings)
        if resi_no == 3:
            pdb_id, axial1, axial2, axial1_segname, axial2_segname = pdb_id, axial1, axial2, axial1_segname, axial2_segname
            xtb_file_strings = [f_string.format(pdb_id=pdb_id, axial1=axial1, axial2=axial2, axial1_segname=axial1_segname, axial2_segname=axial2_segname) for f_string in xtb_file_strings]
        elif resi_no == 2:
            pdb_id, axial1, axial1_segname = pdb_id, axial1, axial1_segname
            xtb_file_strings = [f_string.format(pdb_id=pdb_id, axial1=axial1, axial1_segname=axial1_segname) for f_string in xtb_file_strings]
        with open(f"PDB/{pdb_id}_xtb.inp", "w") as writefile:
            writefile.writelines(xtb_file_strings)
        return

    def run_xtb(self, pdb_id):
        os.chdir(f'PDB/')
        cmd = [f'xtb {pdb_id}_system_protonated.xyz --opt --input {pdb_id}_xtb.inp >{pdb_id}_xtb.out']
        print(cmd)
        process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            print("XTB run completed successfully")
        else:
            with open("heme_struct.err", "a") as error_file:
                with redirect_stderr(error_file):
                    print(f"XTB run failed with return code {process.returncode}")
                    print("Standard Output:")
                    print(stdout.decode())
                    print("Standard Error:")
                    print(stderr.decode())
            cmd_crude = [f'xtb {pdb_id}_system_protonated.xyz --opt crude --input {pdb_id}_xtb.inp']
            print(cmd_crude)
            process = subprocess.Popen(' '.join(cmd_crude), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                print("XTB crude run completed successfully")
                cmd_fron_xtbopt = [f'xtb xtbopt.xyz --opt --input {pdb_id}_xtb.inp >{pdb_id}_xtb.out']
                print(cmd_fron_xtbopt)
                process = subprocess.Popen(' '.join(cmd_fron_xtbopt), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                if process.returncode == 0:
                    print("XTB run completed successfully")
                else:
                    with open("heme_struct.err", "a") as error_file:
                        with redirect_stderr(error_file):
                            print(f"XTB run failed with return code {process.returncode}")
                            print("Standard Output:")
                            print(stdout.decode())
                            print("Standard Error:")
                            print(stderr.decode())
                    with open(f"xtb_history.txt", 'a') as xtb_history_file:
                        xtb_history_file.write(f"{pdb_id} normal\n")
            else:
                with open("heme_struct.err", "a") as error_file:
                    with redirect_stderr(error_file):
                        print(f"XTB run failed with return code {process.returncode}")
                        print("Standard Output:")
                        print(stdout.decode())
                        print("Standard Error:")
                        print(stderr.decode())
                with open(f"xtb_history.txt", 'a') as xtb_history_file:
                    xtb_history_file.write(f"{pdb_id} crude\n")
        current_directory = os.getcwd()
        parent_directory = os.path.dirname(current_directory)
        os.chdir(parent_directory)
        return

    def extract_structure(self, pdb_id):
        start_marker = " final structure:"
        end_marker = "Bond Distances (Angstroems)"
        try:
            with open(f'PDB/{pdb_id}_xtb.out', 'r') as infile:
                lines = infile.readlines()
            try:
                start_index = next(i for i, line in enumerate(lines) if start_marker in line)
            except StopIteration:
                raise ValueError(f"Start marker '{start_marker}' not found in the file.")
            try:
                end_index = next(i for i, line in enumerate(lines) if end_marker in line)
            except StopIteration:
                raise ValueError(f"End marker '{end_marker}' not found in the file.")
            extracted_lines = lines[start_index + 2:end_index]
            with open(f'PDB/{pdb_id}_g16.xyz', 'w') as outfile:
                outfile.writelines(extracted_lines)
            if self.verbose:
                print(f"Extraction completed. Results written to 'PDB/{pdb_id}_g16.xyz'")
        except FileNotFoundError:
            print(f"Error: File f'PDB/{pdb_id}_xtb.out' not found.")
        except Exception as e:
            print(f"Error during extraction: {e}")

    def compare_axials(self, axial1, axial2, axial1_pydish, axial2_pydish, resi_no):
        try:
            axial1_stripped = str(axial1).strip()
            axial1_pydish_stripped = str(axial1_pydish).strip()
            if resi_no == 3:
                axial2_stripped = str(axial2).strip()
                axial2_pydish_stripped = str(axial2_pydish).strip()
                if axial1_stripped == axial1_pydish_stripped and axial2_stripped == axial2_pydish_stripped:
                    print(f"pydish axials ({axial1_pydish_stripped}, {axial2_pydish_stripped}) are equal to axial1 and axial2 ({axial1_stripped}, {axial2_stripped}) as stripped strings")
                else:
                    raise ValueError("Axials are different than they should be. ")
        except ValueError:
            try:
                # Swap the order of comparison
                if axial1_stripped == axial2_pydish_stripped and axial2_stripped == axial1_pydish_stripped:
                    print(f"pydish axials ({axial1_pydish_stripped}, {axial2_pydish_stripped}) are equal to axial1 and axial2 ({axial1_stripped}, {axial2_stripped}) as stripped strings")
                else:
                    raise ValueError("Axials are different than they should be. ")
            except ValueError:
                print(f"pydish axials ({axial1_pydish_stripped}, {axial2_pydish_stripped}) are NOT equal to axial1 and axial2 ({axial1_stripped}, {axial2_stripped}) as stripped strings")

    def calculate_rmsd(self, pdb_id):
        pdb_directory = f"PDB/{pdb_id}"
        if not os.path.exists(pdb_directory) or not os.listdir(pdb_directory):
            pdb_results_df = pd.DataFrame({'PDBID': [pdb_id],
                                        'RMSD_all': [0],
                                        'RMSD_heme': [0],
                                        'RMSD_ax': [0]})
        else:
            try:
                existing_rmsd_df = pd.read_csv('rmsd.csv')
            except FileNotFoundError:
                existing_rmsd_df = pd.DataFrame(columns=['PDBID', 'RMSD_all', 'RMSD_heme', 'RMSD_ax'])
            all_results = [] 
            xtb = mda.Universe(f"PDB/{pdb_id}/{pdb_id}_g16.xyz")
            charmm = mda.Universe(f"PDB/{pdb_id}/{pdb_id}_system_protonated.xyz")
            align.AlignTraj(xtb, charmm).run()
            rmsd = rms.RMSD(xtb, charmm, select='all', groupselections=['all'], ref_frame=0).run()
            pdb_id_column = [pdb_id] * len(rmsd.results.rmsd[:, 2])
            rmsd_values = rmsd.results.rmsd[:, 2]
            if len(xtb.select_atoms('all')) != len(charmm.select_atoms('all')):
                    raise ValueError(f"Number of atoms in xtb ({len(xtb.select_atoms('all'))}) "
                                    f"and charmm ({len(charmm.select_atoms('all'))}) selections must be the same.")
            xtb_select = xtb.select_atoms('all')[:64]
            charmm_select = charmm.select_atoms('all')[:64]
            rmsd = rms.RMSD(xtb_select, charmm_select).run()
            rmsd_values_heme = rmsd.rmsd[:, 2]
            xtb_select_ax = xtb.select_atoms('all')[64:]
            charmm_select_ax = charmm.select_atoms('all')[64:]
            rmsd = rms.RMSD(xtb_select_ax, charmm_select_ax).run()
            rmsd_values_ax = rmsd.rmsd[:, 2]
            pdb_results_df = pd.DataFrame({'PDBID': pdb_id_column, 'RMSD_all': rmsd_values, 'RMSD_heme': rmsd_values_heme, 'RMSD_ax': rmsd_values_ax})
            all_results.append(pdb_results_df)
            final_rmsd_df = pd.concat([existing_rmsd_df] + all_results, ignore_index=True)
            final_rmsd_df.to_csv('rmsd.csv', index=False)
        return

    def __init__(self, verbose=False, cleanup=True):
        self.verbose = verbose
        self.cleanup = cleanup
        self.debug_flags = {
            "debug": False, 
            "debug_add_patches_ax1": False,
            "debug_add_patches_ax2": False,
            "debug_create_charmm_file": False,
            "debug__init__": False,
            "debug_get_axial": False,
            "debug_read_pdb": False,
            "debug_create_xtb_file": False,
            "debug_get_axials_geometric": False,
            "debug_files": False,
            "debug_finders": False
        }
        self.debug_add_patches_ax1 = False
        self.debug_add_patches_ax2 = False
        self.debug_create_charmm_file = False
        self.debug__init__ = False
        self.debug_get_axial = False
        self.debug_read_pdb = False
        self.debug_create_xtb_file = False
        self.debug_create_charmm_file = False
        self.debug_get_axials_geometric = False
        self.debug_files = False
        self.debug_finders = False

        self.residue_names = []
        self.data = self.read_pyDISH()

    #@staticmethod
    def process_row(self, row, verbose, cleanup, debug_flags):
        pdb_id = row['# PDB']
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        try:
            with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
                if debug_flags.get("debug"):
                    print(f"Processing {pdb_id}...")
                axial1 = row['axial1']
                axial2 = row['axial2']
                if debug_flags.get("debug"):
                    print(f"axial1: {axial1}, axial2: {axial2}")
                if debug: print(f"Processing {pdb_id}...")
                axial1 = row['axial1']
                axial2 = row['axial2']
                if debug: print(f"axial1: {axial1}, axial2: {axial2}")
                if debug: print(pdb_id)
                axial1_pydish, axial2_pydish = row['axial1'], row['axial2']
                if debug:
                    print("axial1_pydish: ", axial1_pydish, "     Value:", repr(axial1_pydish), "     Data type:", type(axial1_pydish))
                    print("axial2_pydish: ", axial2_pydish, "     Value:", repr(axial2_pydish), "     Data type:", type(axial2_pydish))
                universe, heme = self.read_pdb(pdb_id)
                molecule = orient_heme.porphyr(pdb_file=f"PDB/{pdb_id}_heme_atom.pdb", pdb_id=pdb_id)
                axials = self.get_axial(pdb_id, universe, heme, axial1_pydish, axial2_pydish)
                if debug: print("#### self.get_axial complete")
                self.fix_heme(pdb_id)
                self.fix_heme_atomtypes(pdb_id)
                self.marry_heme(pdb_id)
                if debug: print("#### self.(heme-stuff) complete")
                axial_list, residues, axial1, axial2, resi_no, axial1_segname, axial2_segname, segments = self.marry_pdbs(pdb_id)
                if debug: print("#### self.marry_pdbs complete")
                self.split_system(pdb_id)
                if debug: print("#### self.split_system complete")
                self.create_charmm_file(pdb_id, axial1, axial2, axial1_segname, axial2_segname, resi_no)
                if debug: print("#### self.create_charmm_file complete")
                self.run_charmm(pdb_id)
                if debug: print("#### self.run_charmm complete")
                ops.check_pdb_file(pdb_id=pdb_id)
                ops.is_xyz_empty(pdb_id=pdb_id)
                if debug: print("#### file checks complete")
                self.create_xtb_file(pdb_id, axial1, axial2, axial1_segname, axial2_segname, resi_no)
                if debug: print("#### self.create_xtb_file complete")
                self.run_xtb(pdb_id)
                self.extract_structure(pdb_id)
                self.compare_axials(axial1, axial2, axial1_pydish, axial2_pydish, resi_no)
                if cleanup:
                    self.cleaner(pdb_id)
                self.calculate_rmsd()
            result = "Success"
        except Exception as e:
            with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
                print(f"Something went wrong with {pdb_id}: {e}")
            result = "Failed"
        return pdb_id, stdout_buffer.getvalue(), stderr_buffer.getvalue(), result

def run_parallel_processing(pdbprep_instance, data):
    results = []
    with ProcessPoolExecutor(max_workers=1) as executor:
        futures = {
            executor.submit(
                pdbprep_instance.process_row,
                row,
                pdbprep_instance.verbose,
                pdbprep_instance.cleanup,
                pdbprep_instance.debug_flags
            ): row['# PDB']
            for _, row in data.iterrows()
        }
        for future in as_completed(futures):
            results.append(future.result())
    return results

if __name__ == "__main__":
    # Instantiate PDBPrep; data is now stored as pdbprep.data.
    pdbprep = PDBPrep()
    data = pdbprep.data
    results = run_parallel_processing(pdbprep, data)
    
    # Dump logs after parallel processing
    with open("heme_struct.log", "w") as log_file, open("heme_struct.err", "w") as err_file:
        for pdb_id, stdout_log, stderr_log, status in results:
            log_file.write(f"Log for {pdb_id} (Status: {status}):\n{stdout_log}\n")
            err_file.write(f"Error log for {pdb_id} (Status: {status}):\n{stderr_log}\n")


# Automation: Sauerstoffe in gaussian