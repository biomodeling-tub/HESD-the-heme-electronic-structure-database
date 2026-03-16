import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import Bio
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
import pdb_tools
import MDAnalysis as mda
import os
import orient_heme
import re

# Define custom color map
colors = ['orange', 'blue', 'red', 'green']
cmap = LinearSegmentedColormap.from_list('CustomColors', colors, N=len(colors))

# Define the directory containing the log files
paths = ["Benchmark/Benchmark_logs/Konrad/", "Benchmark/Benchmark_logs/Functional/B3LYP/", "Benchmark/Benchmark_logs/Functional/BP86/", "Benchmark/Benchmark_logs/Propionate_side_chains/"]
file_names = ["1aw305.log", "1cyo05.log", "1ehb05.log", "1mz405.log", "1ycc05.log", "2pcb05.log"]

def read_homo_lumo(file_names, paths):
    fill_dict_occ, fill_dict_unocc = {}, {}

    for path in paths:
        for file_name in file_names:
            with open(os.path.join(path, file_name), 'r', encoding='utf-8') as file:
                lines = file.readlines()

            debug = False
            # Read HOMOs from the logfile
            occupied = [k for k in lines if k.startswith(" Alpha  occ. eigenvalues -- ")]
            if debug:
                occupied = occupied[:10]
            occupied = [k.replace(" Alpha  occ. eigenvalues -- ", "") for k in occupied]
            occupied = [k.strip() for k in occupied if k.strip()]
            occupied = [float(num) for k in occupied for num in k.split()]
            occupied_short = occupied[-10:]
            name = file_name.replace("01.log", "")
            short_path = path[-7:]
            key = short_path+name
            fill_dict_occ[key] = occupied_short

            # Read LUMOs from the logfile
            unoccupied = [k for k in lines if k.startswith(" Alpha virt. eigenvalues -- ")]
            if debug:
                unoccupied = unoccupied[:10]
            unoccupied = [k.replace(" Alpha virt. eigenvalues -- ", "") for k in unoccupied]
            unoccupied = [k.strip() for k in unoccupied if k.strip()]
            unoccupied = [float(num) for k in unoccupied for num in k.split()]
            unoccupied_short = unoccupied[:10]
            fill_dict_unocc[path+name] = unoccupied_short

        HOMOs, LUMOs = pd.DataFrame(fill_dict_occ), pd.DataFrame(fill_dict_unocc)
    return HOMOs, LUMOs

if __name__=="__main__":
    HOMOs, _ = read_homo_lumo(file_names=file_names, paths=paths)

    # Group columns by the first 4 characters of their names (4 groups)
    groups = HOMOs.columns.str[:4]
    grouped_columns = {group: [] for group in groups.unique()}

    for column in HOMOs.columns:
        group = column[:4]
        grouped_columns[group].append(column)

    plt.figure(figsize=(10, 6))
    ax = plt.gca()  # Get the current Axes object

    for i, (group, columns) in enumerate(grouped_columns.items()):
        group_data = HOMOs[columns]
        color = cmap(i)  # Use the custom color map
        group_data.plot(ax=ax, color=color, label='_nolegend_')
    tick_labels = ['HOMO', 'OMO-1', 'OMO-2', 'OMO-3', 'OMO-4', 'OMO-5', 'OMO-6', 'OMO-7', 'OMO-8', 'OMO-9']
    plt.xticks(np.arange(10), tick_labels)
    plt.xlabel('Orbital')
    plt.ylabel('Energy [au]')
    plt.title('HOMO of 05')
    plt.xticks(np.arange(0, 10, step=1))

    # Remove the legend
    ax.get_legend().remove()

    # Add text to the plot
    plt.text(0.1, 0.9, 'blue = b3lyp, red = bp86, green = pbe1pbe no propchains, orange = pbe1pbe', transform=ax.transAxes, fontsize=12)

    plt.tight_layout()
    plt.show()