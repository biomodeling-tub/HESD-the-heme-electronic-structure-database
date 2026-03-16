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
from natsort import natsorted


# Define the directory containing the log files
paths = ["Benchmark/Benchmark_logs/Functional/B3LYP/", "Benchmark/Benchmark_logs/Functional/BP86/", "Benchmark/Benchmark_logs/Propionate_side_chains/"]
file_names = ["1aw301.log", "1cyo01.log", "1ehb01.log", "1mz401.log", "1ycc01.log", "2pcb01.log"]
#debug = True
# Create an empty dictionary to store file names
#for path in pathes:
#occupied_short, unoccupied_short = {}, {}
def read_homo_lumo(file_names, paths):
    fill_dict_occ, fill_dict_unocc = {}, {}

    for path in paths:
        for file_name in file_names:
            with open(path+file_name, 'r', encoding='utf-8') as file:
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
            #print(f"HOMOs from {file_name}")
            #print(occupied_short)
            #print(f"HOMOs from corresp. dict")
            #print(fill_dict_occ[file_name])


            # Read LUMOs from the logfile
            unoccupied = [k for k in lines if k.startswith(" Alpha virt. eigenvalues -- ")]
            if debug:
                unoccupied = unoccupied[:10]
            unoccupied = [k.replace(" Alpha virt. eigenvalues -- ", "") for k in unoccupied]
            unoccupied = [k.strip() for k in unoccupied if k.strip()]
            unoccupied = [float(num) for k in unoccupied for num in k.split()]
            unoccupied_short = unoccupied[:10]
            fill_dict_unocc[path+name] = unoccupied_short
            #(f"LUMOs from {file_name}")
            #print(unoccupied_short)

            # Write the data to the output file
            #output_file.write(f"File: {file_name}\n")
            #output_file.write("HOMOs:\n")
            #output_file.write(",".join(map(str, occupied_short)) + "\n")
            #output_file.write("LUMOs:\n")
            #output_file.write(",".join(map(str, unoccupied_short)) + "\n")

        HOMOs, LUMOs = pd.DataFrame(fill_dict_occ), pd.DataFrame(fill_dict_unocc)
    return HOMOs, LUMOs

if __name__=="__main__":
        HOMOs, LUMOs = read_homo_lumo(file_names=file_names, paths=paths)


        HOMOs.plot()

        # Customize the plot (optional)
        plt.xlabel('Index')
        plt.ylabel('HOMO Values')
        plt.title('HOMO Values by Path')
        plt.legend()

        # Show the plot
        plt.show()
#plt.savefig('sample_plot.png')

# Print a message indicating that the data has been saved
#print(f"Data saved to {output_file_path}")

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
            short_path = path[-7:]
            key = short_path+name
            unoccupied_short = unoccupied[:10]
            fill_dict_unocc[key] = unoccupied_short

        HOMOs, LUMOs = pd.DataFrame(fill_dict_occ), pd.DataFrame(fill_dict_unocc)
    return HOMOs, LUMOs

if __name__=="__main__":
    _, LUMOs = read_homo_lumo(file_names=file_names, paths=paths)  # Read LUMOs
    HOMOs, _ = read_homo_lumo(file_names=file_names, paths=paths)  # Read LUMOs
    # Group columns by the first 4 characters of their names (3 groups)
    #groups = LUMOs.columns.str[:5]
    #unique_groups = groups.unique()

    #plt.figure(figsize=(10, 6))
    #ax = plt.gca()  # Get the current Axes object

    #for i, group in enumerate(unique_groups):
    #    group_columns = [column for column in LUMOs.columns if column.startswith(group)]
    #    group_data = LUMOs[group_columns]
    #    color = cmap(i)  # Use the custom color map
    #    group_data.plot(ax=ax, color=color, label='_nolegend_')  # Add '_nolegend_' to prevent labels

    #plt.xlabel('Orbital Number')
    #plt.ylabel('LUMO Energy')
    #plt.title('LUMO of 16')
    #plt.xticks(np.arange(0, 10, step=1))

    # Remove the legend
    #ax.get_legend().remove()

    # Add text to the plot
    #plt.text(0.1, 0.9, 'blue = b3lyp, red = bp86, green = no propchains, orange = pbe1pbe', transform=ax.transAxes, fontsize=12)

    #plt.tight_layout()
    #plt.show()


    # Calculate the difference between HOMO and LUMO for each group
    homo_minus_lumo = HOMOs - LUMOs

    # Group columns by the first 4 characters of their names (3 groups)
    groups = homo_minus_lumo.columns.str[:5]
    unique_groups = groups.unique()

    plt.figure(figsize=(10, 6))
    ax = plt.gca()  # Get the current Axes object

    for i, group in enumerate(unique_groups):
        group_columns = [column for column in homo_minus_lumo.columns if column.startswith(group)]
        group_data = homo_minus_lumo[group_columns]
        color = cmap(i)  # Use the custom color map
        group_data.plot(ax=ax, color=color, label='_nolegend_')  # Add '_nolegend_' to prevent labels

    tick_labels = ['HOMO', 'HOMO-1', 'HOMO-2', 'HOMO-3', 'HOMO-4', 'HOMO-5', 'HOMO-6', 'HOMO-7', 'HOMO-8', 'HOMO-9']
    plt.xticks(np.arange(10), tick_labels)
    plt.xlabel('Orbital')
    plt.ylabel('HOMO - LUMO Energy')
    plt.title('HOMO - LUMO of 01')
    plt.xticks(np.arange(0, 10, step=1))

    # Remove the legend
    ax.get_legend().remove()

    plt.tight_layout()
    plt.show()


#read LUMOs from logfile
#unoccupied = [k for k in lines if k.startswith(" Alpha virt. eigenvalues -- ")]
#if debug:
#    unoccupied = unoccupied[:10]
#unoccupied = [k.replace(" Alpha virt. eigenvalues -- ", "") for k in unoccupied]
#unoccupied = [k.strip() for k in unoccupied if k.strip()]
#unoccupied = [float(num) for k in unoccupied for num in k.split()]
#unoccupied = natsorted(unoccupied)
#unoccupied_short = unoccupied[:10]
#print("LUMOs")
#print(unoccupied_short)
