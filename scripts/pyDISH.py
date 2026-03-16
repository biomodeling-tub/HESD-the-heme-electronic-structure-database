import numpy as np
import pandas as pd
import matplotlib as mpl
import pdb_tools
import subprocess
import MDAnalysis as mda
import os
import orient_heme
#from pdbparser.pdbparser import pdbparser
import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    import Bio
    from Bio import PDB
    from Bio.PDB.PDBParser import PDBParser
warnings.filterwarnings("ignore", category=UserWarning)
from file_ops import ops

class pydish():
    def add_std_dev_column(self):
        grouped = self.multihemes.groupby("# PDB")
        std_values = grouped[self.value_column].transform('std')
        self.multihemes["std_"+self.value_column] = std_values
        return

    def trim_pyDISH(self):
        self.data[['axial1', 'axial2']] = self.data['ligand'].str.split('-', expand=True)
        self.data['multiheme'] = self.data['# PDB'].duplicated(keep=False).astype(int)
        self.multihemes = self.data[self.data["multiheme"]==1]
        self.single_hemes = self.data[self.data["multiheme"]==0]
        value_columns = ["pocket_volume", "saddling", "ruffling", "doming", "breathing", "rotation", "wavx", "wavy", "pro", "nstr", "mstr", "transx", "transy"]
        for self.value_column in value_columns:
            self.add_std_dev_column()
        columns_pre = self.multihemes.columns
        order_multi = ["# PDB", "function", "Heme", "multiheme", "ligand", "axial1", "axial2", "pocket_volume", "std_pocket_volume", "saddling", "std_saddling", "ruffling", "std_ruffling", "doming", "std_doming", "breathing", "std_breathing", "rotation", "std_rotation", "wavx", "std_wavx", "wavy", "std_wavy", "pro", "std_pro", "nstr", "std_nstr", "mstr", "std_mstr", "transx", "std_transx", "transy", "std_transy"]
        self.multihemes = self.multihemes[order_multi]
        columns_post = self.multihemes.columns
        #print([k for k in columns_pre if k not in columns_post])
        #print(self.multihemes[:5])
        order_single = ["# PDB", "function", "Heme", "multiheme", "ligand", "axial1", "axial2", "pocket_volume", "saddling", "ruffling", "doming", "breathing", "rotation", "wavx", "wavy", "pro", "nstr", "mstr", "transx", "transy"]
        self.single_hemes = self.single_hemes[order_single]
        ax_lig = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "O", "None", "", None, "OXY", "PER"]
        self.filtered_single_hemes = self.single_hemes[(self.single_hemes["axial1"].isin(ax_lig)) & (self.single_hemes["axial2"].isin(ax_lig))]
        self.data.to_csv("pyDISH.csv")
        self.multihemes.to_csv("multihemes.csv")
        self.filtered_single_hemes.to_csv("single_hemes_filtered.csv")
        return