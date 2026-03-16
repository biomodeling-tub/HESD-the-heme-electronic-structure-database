import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, re, subprocess, sys, json, glob
from molmod import *
from molmod.io import FCHKFile
from molmod.io.xyz import XYZReader, XYZFile
import requests
import molmod
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

#print("""
#|__database
#|    |__logfiles
#|    |  |__*
#|    |
#|    |__logfilessplit
#|    |__pdb
#|       |__from_rcsb_old
#|       |__prepared
#|
#|__heme.ipynb
#|__split_traj_com.tcl
#|__plots
#|__tables
#""")


erep_pattern = re.compile("nuclear repulsion energy")
dipole_pattern = "Dipole moment (field-independent basis, Debye)"
edisp_pattern = re.compile("Nuclear repulsion after empirical dispersion term")
homo_pattern = re.compile("Alpha  occ. eigenvalues")
polarizability_ex_pattern = re.compile("Dipole polarizability, Alpha")


def TCL_Skript():
    """ change the directory and starts the tcl-script
        the TCL-script takes from a pdb-file that contains the whole protein only the heme of chain A
    theres no gurantee it works, sometimes it fails. In this case you create a pdb-file of the heme with vmd and put it in 'pre_prepared' """
    os.chdir("database/pdb")
    os.system('vmd -dispdev text -e split_traj_com.tcl ')
    os.system('mv *.pdb from_rcsb_old')
    os.system('cp pre_prepared/*.pdb prepared/')
    os.chdir("../..")


class porphyr:

    #defautldihedlist orders all dihderal angles with N or the methyl between the pyrrols to the same order.
    defaultdihedlist = [[0, 4, 2, 1], [0, 4, 2, 3], [0, 4, 6, 5], [0, 4, 6, 19], [0, 10, 8, 7], [0, 10, 8, 9], [0, 10, 12, 1], [0, 10, 12, 11], [0, 16, 14, 13], [0, 16, 14, 15], [0, 16, 18, 7], [0, 16, 18, 17], [0, 22, 20, 19], [0, 22, 20, 21], [0, 22, 24, 13], [0, 22, 24, 23], [1, 2, 3, 5], [1, 2, 4, 0], [1, 2, 4, 6], [1, 12, 10, 0], [1, 12, 10, 8], [1, 12, 11, 9], [2, 1, 12, 10], [2, 1, 12, 11], [2, 4, 0, 10], [2, 4, 0, 16], [2, 4, 0, 22], [2, 4, 6, 19], [3, 2, 1, 12], [3, 2, 4, 0], [3, 5, 6, 19], [4, 0, 10, 8], [4, 0, 10, 12], [4, 0, 16, 14], [4, 0, 16, 18], [4, 0, 22, 20], [4, 0, 22, 24], [4, 2, 1, 12], [4, 6, 19, 20], [5, 3, 2, 1], [5, 6, 4, 0], [5, 6, 19, 20], [6, 4, 0, 10], [6, 4, 0, 16], [6, 4, 0, 22], [6, 4, 2, 1], [6, 19, 20, 21], [6, 19, 20, 22], [7, 8, 9, 11], [7, 8, 10, 0], [7, 8, 10, 12], [7, 18, 16, 0], [7, 18, 16, 14], [7, 18, 17, 15], [8, 7, 18, 16], [8, 7, 18, 17], [8, 10, 0, 4], [8, 10, 0, 16], [8, 10, 0, 22], [8, 10, 12, 1], [9, 8, 7, 18], [9, 8, 10, 0], [9, 11, 12, 1], [10, 0, 4, 2], [10, 0, 4, 6], [10, 0, 16, 14], [10, 0, 16, 18], [10, 0, 22, 20], [10, 0, 22, 24], [10, 8, 7, 18], [10, 12, 1, 2], [11, 9, 8, 7], [11, 12, 1, 2], [11, 12, 10, 0], [12, 1, 2, 3], [12, 1, 2, 4], [12, 10, 0, 4], [12, 10, 0, 16], [12, 10, 0, 22], [12, 10, 8, 7], [13, 14, 15, 17], [13, 14, 16, 0], [13, 14, 16, 18], [13, 24, 22, 0], [13, 24, 22, 20], [13, 24, 23, 21], [14, 13, 24, 22], [14, 13, 24, 23], [14, 16, 0, 4], [14, 16, 0, 10], [14, 16, 0, 22], [14, 16, 18, 7], [15, 14, 13, 24], [15, 14, 16, 0], [15, 17, 18, 7], [16, 0, 4, 2], [16, 0, 4, 6], [16, 0, 10, 8], [16, 0, 10, 12], [16, 0, 22, 20], [16, 0, 22, 24], [16, 14, 13, 24], [16, 18, 7, 8], [17, 15, 14, 13], [17, 18, 7, 8], [17, 18, 16, 0], [18, 7, 8, 9], [18, 7, 8, 10], [18, 16, 0, 4], [18, 16, 0, 10], [18, 16, 0, 22], [18, 16, 14, 13], [19, 6, 4, 0], [19, 6, 4, 2], [19, 6, 5, 3], [19, 20, 21, 23], [19, 20, 22, 0], [19, 20, 22, 24], [20, 19, 6, 4], [20, 19, 6, 5], [20, 22, 0, 4], [20, 22, 0, 10], [20, 22, 0, 16], [20, 22, 24, 13], [21, 20, 19, 6], [21, 20, 22, 0], [21, 23, 24, 13], [22, 0, 4, 2], [22, 0, 4, 6], [22, 0, 10, 8], [22, 0, 10, 12], [22, 0, 16, 14], [22, 0, 16, 18], [22, 20, 19, 6], [22, 24, 13, 14], [23, 21, 20, 19], [23, 24, 13, 14], [23, 24, 22, 0], [24, 13, 14, 15], [24, 13, 14, 16], [24, 22, 0, 4], [24, 22, 0, 10], [24, 22, 0, 16], [24, 22, 20, 19]]
    calldict = {}

    # Atomnamen für die entsprechenden dihedralen Winkel
    # atom names for the dihedral angle
    saddling_compass = [["C1_SW", "N_S" ,"N_N","C1_NW"],["C1_SO", "N_S","N_N", "C1_NO"],["C1_WS", "NW", "NO" ,"C1_OS"],["C1_WN", "NW", "NO" ,"C1_ON"]]
    ruffling_compass = [["C3_SW", "C1_SW", "C1_WS", "C3_WS"],["C3_SO", "C1_SO", "C1_OS", "C3_OS"],["C3_WN", "C1_WN", "C1_NW", "C3_NW"],["C3_NO", "C1_NO", "C1_ON", "C3_ON"]]

    def namema(self,liste):
        return [self.calldict[i] for i in liste]

    def getorders(self,lis):
        return [self.order.index(i) for i in lis]

    def getorder(self,elem):
        return self.order.index(elem)

    def sortorder(self,lis):
        return sorted(lis, key = self.getorder)

    #____________________________________________________________________________________ possible tools for debugging  ____________________________________________________

    def get_symbol_per_an_list(self,list_of_indices):
        oz = {1:"H",6:"C",7:"N",8:"O",16:"S",26:"FE"}
        return[oz[i] for i in list_of_indices]

    def get_symbol_per_an_single(self,index):
        oz = {1:"H",6:"C",7:"N",8:"O",16:"S",26:"FE"}
        return oz[index]

    def get_compassid(self,listofindizies):
        return [self.compassdict[i] for i in listofindizies]

    def compassorder(self):
        return( [self.get_compassid([self.order[i] for i in dfd ]  ) for dfd in self.defaultdihedlist] )

    def get_compassname(self,listofindizies):
        strlist = [self.compassdict[i] for i in listofindizies]
        name = strlist[0]
        for n  in strlist[1:] :
            name = name +"_"+n
        return name

    def compassordername(self):
        return( [self.get_compassname([self.order[i] for i in dfd ]  ) for dfd in self.defaultdihedlist] )

    def neighbourlist(self):
        for i in self.order:
            print(i," - ",list(self.mol.graph.neighbors[i]))

    def porphyrcompass(self, Fe, NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ,WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS, SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO , OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON  ):

        """ shows the oriented heme in Ascii, printing the names of the atoms into the respective positions"""

        print(f"          {C3_NW}---{C3_NO}				    	")
        print(f"             /        \                         ")
        print(f"        {C1_NW}       {C1_NO}					")
        print(f"       /      \      /     \                    ")
        print(f"    {WC2N}      {N_N}       {NC2O}				")
        print(f"     |            |    		    |   		    ")
        print(f"  __{C1_WN}       |       {C1_ON}__				")
        print(f"{C3_WN}    \      |      /     {C3_ON}			")
        print(f"|      {N_W}----{Fe}----{N_O}      |  			")
        print(f"{C3_WS}__  /      |      \   __{C3_OS}			")
        print(f"    {C1_WS}       |       {C1_OS}-				")
        print(f"     |            |    		    |			    ")
        print(f"    {SC2W}      {N_S}       {OC2S}			    ")
        print(f"       \       /     \       /                  ")
        print(f"        {C1_SW}       {C1_SO}					")
        print(f"          \               /                     ")
        print(f"          {C3_SW}---{C3_SO}						")



    def porphycompassnumber(self):
        """print the heme by names"""
        self.porphycompass(self.FE, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW ,self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON  )

    def porphyrcompassname(self):
        """print the heme by indizes"""
        self.porphycompass(*self.get_compassid([self.FE, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW , self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON]))

    #____________________________________________________________________________________ possible tools for debugging  _____________________________________________________

    def set_ordered_xyz(self):
        lines, line = [], ""
        lines.append(f"25\n")
        lines.append(f'Heme structure file of {self.pdb_file.replace("PDB/", "").replace("_heme_ATOM.pdb","")} \n')
        for i in self.order:
            line += str(self.get_symbol_per_an_single(self.mol.graph.numbers[i]))
            line += str(self.mol.coordinates[i]/1.8897259886)
            line = line.replace(" ", "  ")
            line = line.replace("[", "       ")
            line = line.replace("]", "\n")
            lines.append(line)
            line = ""
        file = open(f"PDB/{self.pdb_id}_heme_sorted.xyz", "w").writelines(lines)
        if self.verbose:
            print(lines)
        return

    def format_xyz_file(self):
        self.input_file, self.output_file = f"PDB/{self.pdb_id}_heme_sorted.xyz", f"PDB/{self.pdb_id}_heme.xyz"
        with open(self.input_file, 'r') as f:
            lines = f.readlines()
        atoms, coordinates = [], []
        for line in lines[2:]:
            parts = line.split()
            atoms.append(parts[0])
            coordinates.append([float(coord) for coord in parts[1:]])
        max_widths = [max(len(str(coord)) for coord in column) for column in zip(*coordinates)]
        with open(self.output_file, 'w') as f:
            f.write(f'{len(atoms)}\n')
            f.write(f'Heme structure file of {self.pdb_file.replace("PDB/", "").replace("_heme_ATOM.pdb","")} \n')
            for atom, coords in zip(atoms, coordinates):
                coord_str = ' '.join(f'{coord:{width}.6f}' for coord, width in zip(coords, max_widths))
                f.write(f'{atom:<2s} {coord_str}\n')
        return

    # Here is the part where the lists of dihedrals are made and get_dihed_per_list get a list of Index list of four indizes and calculate the dihdral angle.


    def get_saddling_compass(self): # just the give the name of the viewed atoms in the saddling-dihedrals
        return ["C1_SW_N_S_N_N_C1_NW" ,"C1_SO_N_S_N_N_C1_NO", "C1_WS_NW_NO_C1_OS", "C1_WN_NW_NO_C1_ON"]

    def get_ruffling_compass(self): # just the give the name of the viewed atoms in the ruffling-dihedrals
        return  ["C3_SW_C1_SW_C1_WS_C3_WS","C3_SO_C1_SO_C1_OS_C3_OS","C3_WN_C1_WN_C1_NW_C3_NW","C3_NO_C1_NO_C1_ON_C3_ON"]

    def get_dihed_per_list(self,givenlist):
        dihedral = []
        for i1, i2, i3, i4 in givenlist:
            dihedral.append([[i1, i2, i3,i4],self.namema([i1, i2, i3,i4]), dihed_angle(self.mol.coordinates[[i1, i2, i3,i4]])[0]/deg , [self.order.index(i) for i in [i1,i2,i3,i4]]   ])
        return dihedral

    # get_dihed, get_dihed_ruffling, get_dihed_saddling are the function, that you use later for getting the dihedrals. if you want to figure out more dihedrals, create a list of the lists with four indizies and give it to get_dihed_per_list

    def get_dihed(self):
        dihedlist = [[self.order[i] for i in [io1, io2, io3, io4]]  for io1, io2, io3, io4 in self.defaultdihedlist]
        return self.get_dihed_per_list(dihedlist)

    def get_dihed_ruffling(self):
        return self.get_dihed_per_list(self.ruffling)

    def get_dihed_saddling(self):
        return self.get_dihed_per_list(self.saddling)



    def set_porphyr(self):
        """ set_porphyr creates a Molecule, of the porpyhrinring and Fe. heme is orientated that { R2: north, acid & methyl : west, acid &R3: south, R1: east}
        set_porphyr erstellt ein Molecule, bestend aus den Porpyhrinring und Fe. Das Häm wird dabei anhand der Seitenketten so orientiert, dass Pyrrol mit { R2:Norden, Prop.Säure & Methyl : Westen, Prop.Säure&R3: Süden, R1: Osten}  """
        for i1 in range(self.mol.size):                 # loop over all atoms
            n = list(self.mol.graph.neighbors[i1])      # n are the neigbours of the current atom
            if self.verbose == True:
                print(f"molecular size is {self.mol.size}")
                print(f"neighborlist is: {n}")
                #print(n)
                print(f"mol.numbers[i1] (=neighbour atoms) currently processed are {self.mol.numbers[i1]}")
            if self.mol.numbers[i1] == 26:              # finding FE based on atom number 26 -> n are the axial Ligands and N
                NDict = {}
                for N in n:
                    C1 = list(self.mol.graph.neighbors[N])  # C1 = alpha-C of Pyrrol (neigbour of N)
                    C1Dict ={}
                    for nn in n:
                        try:
                            C1.remove(nn)               # remove N if a N is neigbour of another N
                        except:
                            2
                    for c1 in C1:
                        C2 = list(self.mol.graph.neighbors[c1])
                        C1Dict.update({c1:C2})
                        for nn in n:
                            try:
                                C2.remove(nn)        # only want to have beta-C and not N
                            except:
                                2
                    NDict.update({N:C1Dict})
                    #NDict = {neibour of Fe:{neibour of the neigbours of FE:{neibour of the neighbours of the neibours of Fe}...} ...}
                    # it's supposed to be = {N:{C1:{C2,C2}   }} but there's still the axial ligands who have to be discovered
                allof = [[[c2 for c2 in c1 ]for c1 in N.values()] for N in NDict.values()]      # allof is a list of lists of lists of the indizies
                c =  []
                for a in allof:
                    for b in a:
                        c = c + b   #                           c is one list of indizies
                self.por_index = []  # self.por_index shall contain only the indizes of atoms, of the porphyrinring
                for k,N in zip (NDict.keys(),NDict.values()):
                    for c1 in N.values():
                        for c2 in c1:
                            cc=0
                            if c.count(c2) == 2:   #  We want only the atoms of the porphyrinring and not of the rest-groups. The atoms of the porphyrinring appear at least to times in NDict because they have at least to neibours, that are in the Dict. This way we can find out if the watched atom is in the porphyrinring or in a rest-group
                                cc+=1
                            if cc > 0:
                                self.por_index.append(k)  #
                                for zz in c1:
                                    self.por_index.append(zz)
                                for zz in list(N.keys()):
                                    self.por_index.append( zz )

        self.graph  = self.mol.graph.get_subgraph(self.por_index)  # self.graph is the Molecule; only the prophyrin and FE
        N = []
        C1 = []
        C2 = []
        C3 = []
        for i in set(self.por_index):
            if self.graph.numbers[i] == 7 :
                N.append(i)
            elif 7 in  self.graph.numbers[ list(self.graph.neighbors[i])] and self.graph.numbers[i]!=26 :
                C1.append(i)
        for i in set(self.por_index):
            if i not in N and i not in C1:
                if all( i in C1 for i in list(self.graph.neighbors[i])):
                    C2.append(i)
                elif self.graph.numbers[i]!=26:
                    C3.append(i)
        self.calldict = {}
        for i in set(self.por_index):
            if i in N:
                self.calldict[i] = "N"
            if i in C1:
                self.calldict[i] = "C1"
            if i in C2:
                self.calldict[i] = "C2"
            if i in C3:
                self.calldict[i] = "C3"
            if self.graph.numbers[i]==26:
                self.calldict[i] = "Fe"
        N1  = N[0]                      # now we have all atoms of the Heme and know if it is N, alpha-C (C1), beta-C (C2) or a methyl-C between the pyrrosl (C3). now we want to find out the orientation of the indices using the rest groups.
        C11a,C11b = 1,2
        C11a,C11b = [i for i in  list(self.mol.graph.neighbors[N1]) if i in C1]     # We start, by numbering the pyrrols from 1 to 4 with the clock, starting by a random N
        C1C2 = [i for i in list(self.mol.graph.neighbors[C11b]) if i in C2][0]
        C4C1 = [i for i in list(self.mol.graph.neighbors[C11a]) if i in C2][0]
        C31b = [i for i in list(self.mol.graph.neighbors[C11b]) if i in C3][0]
        C31a = [i for i in list(self.mol.graph.neighbors[C11a]) if i in C3][0]
        C12a = [i for i in list(self.mol.graph.neighbors[C1C2]) if i in C1 and i!=C11b][0]
        N2 = [i for i in list(self.mol.graph.neighbors[C12a]) if i in N][0]
        C32a = [i for i in list(self.mol.graph.neighbors[C12a]) if i in C3 ][0]
        C32b = [i for i in list(self.mol.graph.neighbors[C32a]) if i in C3][0]
        C12b = [i for i in list(self.mol.graph.neighbors[C32b]) if i in C1][0]
        C12a = [i for i in list(self.mol.graph.neighbors[C32a]) if i in C1][0]
        C2C3 = [i for i in list(self.mol.graph.neighbors[C12b]) if i in C2][0]
        C13a = [i for i in list(self.mol.graph.neighbors[C2C3]) if i in C1 and i != C12b]  [0]
        C33a = [i for i in list(self.mol.graph.neighbors[C13a]) if i in C3][0]
        C33b = [i for i in list(self.mol.graph.neighbors[C33a]) if i in C3][0]
        C13b = [i for i in list(self.mol.graph.neighbors[C33b]) if i in C1][0]
        C13a = [i for i in list(self.mol.graph.neighbors[C33a]) if i in C1][0]
        N3 = [i for i in list(self.mol.graph.neighbors[C13a]) if i in N][0]
        C3C4 = [i for i in list(self.mol.graph.neighbors[C13b]) if i in C2][0]
        C14a = [i for i in list(self.mol.graph.neighbors[C3C4]) if i in C1 and i != C13b][0]
        C34a = [i for i in list(self.mol.graph.neighbors[C14a]) if i in C3 ][0]
        C34b = [i for i in list(self.mol.graph.neighbors[C34a]) if i in C3][0]
        C14b = [i for i in list(self.mol.graph.neighbors[C34b]) if i in C1][0]
        C14a = [i for i in list(self.mol.graph.neighbors[C34a]) if i in C1][0]
        N4 = [i for i in list(self.mol.graph.neighbors[C14a]) if i in N][0]
        Fe = [i for i in self.calldict.keys() if self.calldict[i] == "Fe"][0]
        if not all([     all(  i in C1 for i in [C11a,C11b,C12b,C12a]),          all(  i in C2 for i in [C1C2,C4C1]),          all(  i in C3 for i in [C31a,C31b,C32a,C32b]),          all(  i in N for i in [N1,N2])  ]):
            raise ValueError('Error in assignment of atoms.')
        methyldict={}
        for i in C3:
            c = [c for c  in list(self.mol.graph.neighbors[i]) if c not in C3 and c not in C1][0]
            if all([i != 6 for i in self.graph.numbers[[v for v in list(self.mol.graph.neighbors[c]) if v != i]  ]]):
                methyldict[i] = True
            else:
                methyldict[i] = False
            list(self.graph.neighbors[i])
        pairs = [(C31b, C32a), (C32b, C33a), (C33b, C34a), (C34b, C31a)]
        npairs = [N1, N2, N3, N4]
        self.g1 = [C4C1, C11a, C31a, N1, C31b, C11b]
        self.g2 = [C1C2, C12a, C32a, N2, C32b, C12b]
        self.g3 = [C2C3, C13a, C33a, N3, C33b, C13b]
        self.g4 = [C3C4, C14a, C34a, N4, C34b, C14b]
        for p in pairs:                                             # Here we define if the pyrrol binds an acid or another restgroup, thats how we want to orient the Heme
            if methyldict[p[0]] and methyldict[p[1]]:
                methyl = list(p)
                NS = npairs[pairs.index(p)]
                NN = npairs[pairs.index(p)-2]
            if not methyldict[p[0]] and not methyldict[p[1]]:
                acid = list(p)
                NW = npairs[pairs.index(p)]
                NO  = npairs[pairs.index(p)-2]
        for i,g in enumerate( [self.g1,self.g2,self.g3,self.g4]):
            if any(m in g for m in methyl):
                if any(a in g for a in acid):
                    SG = g
                else:
                    OG = g
            else:
                if any(a in g for a in acid):
                    WG = g
                else:
                    NG = g
        WC2N, C1_NW, C3_NW, N_N, C3_NO, C1_NO = NG
        NC2O, C1_ON, C3_ON, N_O, C3_OS, C1_OS = OG
        OC2S, C1_SO, C3_SO, N_S, C3_SW, C1_SW = SG
        SC2W, C1_WS, C3_WS, N_W, C3_WN, C1_WN = WG
        if not (C1_OS in list(self.mol.graph.neighbors[OC2S]) and C1_SO in list(self.mol.graph.neighbors[OC2S])):    # It is not knwon if the heme is mirrored, if it is we mirror the indizes
            NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW = NG
            WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS = OG
            SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO = SG
            OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON = WG

            WC2N, C1_NW, C3_NW, N_N, C3_NO, C1_NO = NG
            SC2W, C1_WS, C3_WS, N_W, C3_WN, C1_WN = OG
            OC2S, C1_SO, C3_SO, N_S, C3_SW, C1_SW = SG
            NC2O, C1_ON, C3_ON, N_O, C3_OS, C1_OS = WG
        else:
            print("no Change")
        self.Fe, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW , self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON  =Fe, NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ,WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS, SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO , OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON
        self.importantdihed = [Fe, NC2O, OC2S, SC2W, WC2N]
        self.order  = [Fe] + NG + OG + SG + WG
        self.order  = [Fe] + [ NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ] + [WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS] + [SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO ] + [OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON]
        self.saddling = [[C1_SW, N_S ,N_N,C1_NW],[C1_SO, N_S,N_N, C1_NO],[C1_WS, N_W, N_O ,C1_OS],[C1_WN, N_W, N_O ,C1_ON]]
        self.ruffling = [[C3_SW, C1_SW, C1_WS, C3_WS],[C3_SO, C1_SO, C1_OS, C3_OS],[C3_WN, C1_WN, C1_NW, C3_NW],[C3_NO, C1_NO, C1_ON, C3_ON]]
        self.compassdict = { Fe:"FE",   NC2O:"NC2O", C1_NO:"C1_NO",C3_NO:"C3_NO", N_N:"N_N", C3_NW:"C3_NW", C1_NW:"C1_NW",             OC2S:"OC2S", C1_OS:"C1_OS", C3_OS:"C3_OS", N_O:"N_O", C3_ON:"C3_ON", C1_ON:"C1_ON",             SC2W:"SC2W", C1_SW:"C1_SW", C3_SW:"C3_SW", N_S:"N_S", C3_SO:"C3_SO", C1_SO:"C1_SO",            WC2N:"WC2N", C1_WN:"C1_WN", C3_WN:"C3_WN", N_W:"N_W", C3_WS:"C3_WS", C1_WS:"C1_WS" }
        if self.verbose == True:
            molmod.io.pdb.dump_pdb(filename=f"PDB/{self.pdb_id}_dump.pdb", molecule=self.mol)
            print(self.order)
            print("molmod internal coordinates:")
            for i in self.order:
                print(self.mol.coordinates[i])
            print("molmod coordinates turned back into pdb coordinates:")
            for i in self.order:
                print(self.mol.coordinates[i]/1.8897259886)
        self.set_ordered_xyz()
        self.format_xyz_file()
        return

    # INIT

    def __init__(self, pdb_file, pdb_id):
        """porpyr(file.pdb) creates a Molecule.heme where the heme is oriented based on the restgroups"""
        self.pdb_file, self.pdb_id = pdb_file, pdb_id
        self.verbose = False
        if self.verbose:
            print(f"porphyr.__init__() reading molecule from file {self.pdb_file}")
        self.mol = Molecule.from_file(self.pdb_file)
        if self.verbose:
            print(f"coordinates read from {self.pdb_file}")
            molmod.io.pdb.dump_pdb(self.pdb_file, molecule=self.mol)
        if self.verbose:
            print(self.mol.coordinates)
            print("converting molecule to molecule.graph data-structure")
            print(f"molecular system size before graph building: {self.mol.size}")
        self.mol.set_default_graph()
        if self.verbose:
            print(f"molecular system size after graph building: {self.mol.size}")
            print("orienting the molecule based on side-chain positions")
        molecule = self.set_porphyr()
        if self.verbose:
            print(f"molecule oriented and output printed to file {self.output_file}")


class dihedpdb:
    """dihedpdb creates tables of the wanted dihedrals by looping through the pdb """
    df = pd.DataFrame()

    def __init__(self, **kwargs):
        """read_keep = False; if True the old table will be added by the new dihedrals, else a new Table will be created """
        if kwargs.get("read_keep"):
            try:
                self.df = pd.read_csv("tables/Dihedral.csv")                        # all dihedrals containig N or the C between the Pyrrols
                self.df = self.df.rename(columns={"Unnamed: 0": "pdb"})
                self.df = self.df.set_index("pdb")
            except:
                self.df = pd.DataFrame()
                print("except_dihed")

            try:
                self.df_saddling = pd.read_csv("tables/Saddling.csv")
                self.df_saddling = self.df_saddling.rename(columns={"Unnamed: 0": "pdb"}) # all saddling dihedrals
                self.df_saddling = self.df_saddling.set_index("pdb")
            except:
                self.df_saddling = pd.DataFrame()
                print("except_saddling")

            try:
                self.df_ruffling = pd.read_csv("tables/Ruffling.csv")
                self.df_ruffling = self.df_ruffling.rename(columns={"Unnamed: 0": "pdb"}) # all ruffling dihedrals
                self.df_ruffling = self.df_ruffling.set_index("pdb")
            except:
                self.df_ruffling = pd.DataFrame()
                print("except ruffling")

        else:
            self.df = pd.DataFrame()
            self.df_ruffling = pd.DataFrame()
            self.df_saddling = pd.DataFrame()

        knowndihed = list(self.df.index)                        #knwondihed is to make the the skript avoiding to calculate a pdb twice
        knownsaddling = list(self.df_saddling.index)
        knownruffling = list(self.df_ruffling.index)

        for k in glob.glob("database/pdb/prepared/*.pdb"):
            knowndihed = list(self.df.index)                        #knwondihed is to make the the skript avoiding to calculate a pdb twice
            knownsaddling = list(self.df_saddling.index)
            knownruffling = list(self.df_ruffling.index)
            if not (k[k.find("/prepared/")+10:-4].upper() in knowndihed) or not (k[k.find("/prepared/")+10:-4].upper() in knownsaddling) or not (k[k.find("/prepared/")+10:-4].upper() in knownruffling):
                try:
                    a = porphyr(k)
                    if not (k[k.find("/prepared/")+10:-4].upper() in knowndihed):
                        try:
                            dihed = a.get_dihed()
                            self.df = self.df.append(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.compassordername())) #)   )
                        except:
                            print(f"problems with {k} dihed" )

                    if not (k[k.find("/prepared/")+10:-4].upper() in knownsaddling):
                        try:
                            saddling = a.get_dihed_saddling()
                        # print(pd.DataFrame( [[i[2] for i in saddling]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_saddling_compass()))
                            self.df_saddling = self.df_saddling.append(pd.DataFrame( [[i[2] for i in saddling]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_saddling_compass()))
                        except:
                            print(f"problems with {k} saddling" )

                    if not (k[k.find("/prepared/")+10:-4].upper() in knownruffling):
                        try:
                            dihed = a.get_dihed_ruffling()
                            self.df_ruffling = self.df_ruffling.append(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_ruffling_compass()))

                        except:
                            print(f"problems with {k} ruffling" )
                except:
                    print("problems with ",k)

        self.df.index = self.df.index.str.upper()
        self.df.to_csv("tables/Dihedral.csv")
        self.df_saddling.index = self.df_saddling.index.str.upper()
        self.df_saddling.to_csv("tables/Saddling.csv")
        self.df_ruffling.index = self.df_ruffling.index.str.upper()
        self.df_ruffling.to_csv("tables/Ruffling.csv")


class prepare_gaussian_logs:
    """preaper_gaussian_log splits logfiles with two calculations (chlor and nbo), if nbo is the first one If nbo is the only one, nothing happens
    the following conditionas are important

    exept of 01 ox state and multiplicity should be in the end of the name of the log file: 1a6g05.log

    If the jobs of chloro and nbo splitted in different jobs it has to be the following ending of names

    "Link1.log":"_01_nbo.log",
    "Link2.log":"_01_chloro.log",
    "05Link3.log": "_05_nbo.log",
    "05Link4.log":"_05_chloro.log",
    "12Link5.log": "_12_nbo.log",
    "12Link6.log":"_12_chloro.log",
    "16Link7.log": "_16_nbo.log",
    "16Link8.log":"_16_chloro.log"

    all log file should be in a file, with thats called with the pdb for example:
    ├── 1a6g
    │   ├── 1a6g05.log
    │   ├── 1a6g12.log
    │   ├── 1a6g16Link7.log
    │   ├── 1a6g16Link8.log
    │   └── 1a6g.log


    """
    finaldirectory = "database/logfilessplit/"

    foulder_of_foulders = "database/logfiles/"


    def defaultsplit(self,foulder,suffix):
        with open(foulder) as f:
            lines = f.readlines()
        indizesofend = [i for i,line  in enumerate(lines) if "Normal termination of Gaussian 16" in line]
        name=foulder[foulder.rindex('/')+1:][:-4]
        if (name[-2:]) == "05":
            ox = "05"
            pdb = name[:-2]
        elif (name[-2:]) == "12":
            ox = "12"
            pdb = name[:-2]
        elif (name[-2:]) == "16":
            ox = "16"
            pdb = name[:-2]
        elif (name[-2:]) == "01":
            ox = "01"
            pdb = name[:-2]
        else:
            ox = "01"
            pdb = name

        if (len(indizesofend)==2):
            nbo = lines[:indizesofend[0]+1]
            chloro = lines[indizesofend[0]+1:]
        #  with open( self.finaldirectory +  name+"_nbo.log", 'w') as fp:
            with open( self.finaldirectory +  pdb+"_"+ ox +"_nbo" + ".log", 'w') as fp:
                for line in nbo:
                    fp.write(line)
    #            with open( self.finaldirectory +  name+"_chloro.log", 'w') as fp:
            with open( self.finaldirectory +  pdb+"_"+ ox +"_chloro" + ".log", 'w') as fp:
                for line in chloro:
                    fp.write(line)
        else:
            with open( self.finaldirectory +  pdb+"_"+ ox +"_nbo" + ".log", 'w') as fp:
                for line in lines:
                    fp.write(line)
        #print("database/logfiles/"+name)

    def core(self):
        for k in glob.glob(self.foulder_of_foulders+"*"):
            j = k[k.rindex("/")+1:]
            l  = glob.glob(k+"/*")
            Linksuffix = {                                                  # its in portant, the calculations folows this notions.
                "Link1.log":"_01_nbo.log",
                "Link2.log":"_01_chloro.log",
                "05Link3.log": "_05_nbo.log",
                "05Link4.log":"_05_chloro.log",
                "12Link5.log": "_12_nbo.log",
                "12Link6.log":"_12_chloro.log",
                "16Link7.log": "_16_nbo.log",
                "16Link8.log":"_16_chloro.log"
            }
            for L in l:
                L2 = L[L.rindex("/")+1:]
                try:
                    sfx = Linksuffix[ [i for i in list(Linksuffix.keys()) if i in L2][0]]
                    os.system( "cp {} {}/{}".format(L,self.finaldirectory,j+sfx)     )
                except:
                    if L[L.rindex("/")+1:-4] == k[k.rindex("/")+1:]:
                        print(L,j+"_01.log")
                        self.defaultsplit(L,j+"_01.log")
                    else:
                        self.defaultsplit(L,L2)
    def __init__(self):
        return


class helpers():
    def coordtomatrix(coordinates):
        Matrix=[[int(float(i)) if int(float(i))==float(i) else float(i) for i in j.split() ] for j in coordinates]
        df1 =pd.DataFrame(Matrix)
        df = df1.rename({0 : 'CenterNumber', 1: 'AtomicNumber', 2:'AtomicType', 3: 'X' ,4:'Y', 5:'Z'  },axis=1)
        return  df

    def dist(tupel1, tupel2):
        a = [(tupel1[i]-tupel2[i])**2 for i in range(0,3)]
        b=np.sqrt(a[0]+a[1]+a[2])
        return b

    def get_geom(streams): # extracts the geometry from the compressed stream - input orientation!
        geom = []
        for item in streams[-1][16:]:
            if item == "":
                break
            geom.append([item.split(",")[0],float(item.split(",")[-3]),float(item.split(",")[-2]),float(item.split(",")[-1])])
        return(geom)

    def dicttodelete(dictofvalues, number=0):
        arbeitsdict={}
        arbeitsdict.update({
                            'pdb':dictofvalues["pdb"],#[0],
                            'Ox':dictofvalues["Ox"],#[0],
                            'spin':dictofvalues["spin"],#[0],
                            'method':dictofvalues["method"],#[0],
                            "e":dictofvalues["e"][number] ,
                            "edisp":dictofvalues["edisp"],
                            "homo":dictofvalues["homo"][0],
                            "lumo":dictofvalues["homo"][1],
                            "chem_pot":dictofvalues["homo"][2],
                            "diff":dictofvalues["homo"][3],
                            "elekphil":dictofvalues["homo"][4],
                            "dipole" : dictofvalues["dipole"]  ,
                            "qpole1" : dictofvalues["qpole"][0] ,
                            "qpole2" : dictofvalues["qpole"][1],
                            "qpole3" : dictofvalues["qpole"][2],
                            "qpole4" : dictofvalues["qpole"][3],
                            "polar-iso":  dictofvalues["polar-iso"],
                            "polar-aniso":  dictofvalues["polar-aniso"]})
        return arbeitsdict

    def delete():
        def a(pdb):
            dfwork = dfges.loc[[pdb]].copy()
            dfbela = dfwork.copy()
            pdb_I = [i+"_"+c*"I" for i in  dfwork.index for c in range(1,len(dfexp.loc[pdb])+1)]
            for c in range(1,len(dfexp.loc[pdb])):
                dfwork = dfwork.append(dfwork)
            #dfbela
            dfwork["pdb_I"]= pdb_I
            return dfwork


class Physical_quantity:
    """it is necessary, that ne logfile contents only one job and the nation ist pdb_ox_solution.

    Physical_quantity take the wanted measurments from the log-file with get-functions and creates a dictionary
    """

    #_________________________________get and set functions_______________________________________________________________________________


    # some othe get functions return an error message


    def get_planeangle(self):   # not used
        streams = self.outstreams
        atoms=["2","5","8","12","13","19"]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]  #you have to define a list of atoms if you want to use, these numbers are random
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        planeangleout = ""
        error = ""
        if len(atoms)%6 != 0:
            error = str(len(atoms)) + " numbers given. Plane angles require sets of 6 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for plane angle "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+";"+int(len(atoms)/6)*";")

        for i in range(int(len(atoms)/6)):
            a = geom[int(atoms[6*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[6*i+1])-1][:4]
            c = geom[int(atoms[6*i+2])-1][:4]
            d = geom[int(atoms[6*i+3])-1][:4]
            e = geom[int(atoms[6*i+4])-1][:4]
            f = geom[int(atoms[6*i+5])-1][:4]

            ab = np.array([a[1]-b[1],a[2]-b[2],a[3]-b[3]]) # Vectors
            bc = np.array([b[1]-c[1],b[2]-c[2],b[3]-c[3]])
            de = np.array([d[1]-e[1],d[2]-e[2],d[3]-e[3]])
            ef = np.array([e[1]-f[1],e[2]-f[2],e[3]-f[3]])

            n1 = np.cross(ab,bc) # Normal vectors
            n2 = np.cross(de,ef)

            planeangle = round(np.degrees(np.arccos(np.dot(n1,n2) / (np.linalg.norm(n1)*np.linalg.norm(n2)))),3)
            planeangle = min(abs(planeangle),abs(180-planeangle))
            planeangleout += str(a[0]+atoms[4*i+0]+" " + b[0]+atoms[4*i+1]+" " + c[0]+atoms[4*i+2]+" " + d[0]+atoms[4*i+3]) + ";" + str(planeangle) + ";"
        return(planeangleout,error)

    def set_planeangle(self):   # not used
        self.planeangle ,self.error["planeangle"]= self.get_planeangle()

    def get_dihedrals(self):    # not used
        streams = self.outstreams
        atoms = ["2","4","6","8"]
        atoms = [1,2,3,4,5,6,7,8,9,10,11, 12]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        dihedralout = ""
        error = ""
        # check input
        if len(atoms)%4 != 0:
            error = str(len(atoms)) + " numbers given. Dihedrals require sets of 4 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for dihedrals "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+";"+int(len(atoms)/4)*";")
        for i in range(int(len(atoms)/4)):
            a = geom[int(atoms[4*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[4*i+1])-1][:4]
            c = geom[int(atoms[4*i+2])-1][:4]
            d = geom[int(atoms[4*i+3])-1][:4]
            ab = np.array([a[1]-b[1],a[2]-b[2],a[3]-b[3]]) # Vectors
            bc = np.array([b[1]-c[1],b[2]-c[2],b[3]-c[3]])
            cd = np.array([c[1]-d[1],c[2]-d[2],c[3]-d[3]])
            n1 = np.cross(ab,bc) # Normal vectors
            n2 = np.cross(bc,cd) - bc
            dihedral = round(np.degrees(np.arccos(np.dot(n1,n2) / (np.linalg.norm(n1)*np.linalg.norm(n2)))),3)
            dihedralout += str(a[0]+atoms[4*i+0]+" " + b[0]+atoms[4*i+1]+" " + c[0]+atoms[4*i+2]+" " + d[0]+atoms[4*i+3]) + ";" + str(dihedral) + ";"
        return(dihedralout,error)

    def set_dihedrals(self):    # not used
        self.dihedral ,self.error["dihedral"]= self.get_dihedrals()

    def get_angles(self):       # not used
        streams = self.outstreams
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        anglesout = ""
        error = ""
        atoms = [1,2,3,4,5,6,7,8,9,10,11, 12]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]
        if len(atoms)%3 != 0:
            error = str(len(atoms)) + " numbers given. Angles require sets of 3 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for angles "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+int(len(atoms)/3)*";")
        for i in range(int(len(atoms)/3)):
            a = geom[int(atoms[3*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[3*i+1])-1][:4]
            c = geom[int(atoms[3*i+2])-1][:4]
            ba = np.array(a[1:]) - np.array(b[1:])
            bc = np.array(c[1:]) - np.array(b[1:])
            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle = np.arccos(cosine_angle)
            anglesout += str(a[0]+atoms[3*i+0]+" " + b[0]+atoms[3*i+1]+" " + c[0]+atoms[3*i+2]) + ";" +str(round(np.degrees(angle),3)) + ";"
        return(anglesout,error)

    def set_angles(self):       # not used
        self.angles ,self.error["angle"]= self.get_angles()

    def get_outstreams(self): # gets the compressed stream information at the end of a Gaussian job
        streams = []
        starts,ends = [],[]
        error = "failed or incomplete job" # default unless "normal termination" is in file
        loglines=self.lines
        for i in range(len(loglines)):
            if "1\\1\\" in loglines[i]:
                starts.append(i)
            if "@" in loglines[i]:
                ends.append(i)
            if "Normal termination" in loglines[i]:
                error = ""
        if len(starts) != len(ends) or len(starts) == 0: #probably redundant
            error = "failed or incomplete job"
            return(streams,error)
        for i in range(len(starts)):
            tmp = ""
            for j in range(starts[i],ends[i]+1,1):
                tmp = tmp + loglines[j][1:-1]
            streams.append(tmp.split("\\"))
        return(streams,error)

    def set_outstreams(self):
        self.outstreams=self.get_outstreams()[0]

    def get_distances(self):
        df=self.df.copy()
        distlist = [helper.dist((df["X"][i],df["Y"][i],df["Z"][i]),((df["X"][j],df["Y"][j],df["Z"][j]))) for i in range(0,df.index.stop)  for j in range(0,df.index.stop) ]
        distmatrix = np.array(distlist).reshape(self.df.index.stop,self.df.index.stop)
        distdict = pd.DataFrame(distmatrix)
        return distdict

    def set_distances(self):
        self.distances = self.get_distances()

    def get_edisp(self):
        filecont=self.lines
        error = "no dispersion correction found in this file"
        disp = 0
        for i in range(len(filecont)-1,1,-1):
            if edisp_pattern.search(filecont[i]):
                e_rep_disp = float(filecont[i].split()[-2]) # hartree
                disp = 1
            if erep_pattern.search(filecont[i]) and disp == 1:
                e_rep_nodisp = float(filecont[i].split()[-2]) # hartree
                e_disp = (e_rep_disp - e_rep_nodisp) * 627.50947
                return(e_disp,"")
        return(None,error)

    def set_edisp(self):
        self.edisp, self.error["edisp"]=self.get_edisp()

    def get_polarizability(self): #iso aniso
        filecont = self.lines
        for i in range(len(filecont)-1,1,-1):
            if polarizability_ex_pattern.search(filecont[i]):
                alpha_iso = float(filecont[i+4].split()[1].replace("D","E"))
                alpha_aniso = float(filecont[i+4].split()[2].replace("D","E"))
                #return(str(alpha_iso) + ";" + str(alpha_aniso) + ";","")
                return ({"iso":alpha_iso, "aniso":alpha_aniso},"")
        error = "no polarizability information found;"
        alpha_aniso, alpha_iso = None,None
        return({"iso":alpha_iso, "aniso":alpha_aniso}   ,error)

    def set_polarizability(self):
        self.polar, self.error["polar"] = self.get_polarizability()

    def get_quadrupole(self):
        for item in self.outstreams[-1]:
            if "Quadrupole" in item:
                q = item[11:].split(",")
                q = [float(i) for i in q]
                q_comps = np.array(([q[0],q[3],q[4]],[q[3],q[1],q[5]],[q[4],q[5],q[2]]))
                q_diag = np.linalg.eig(q_comps)[0]
                q_ampl = np.linalg.norm(q_diag)
                results = (np.max(q_diag),-(np.max(q_diag)+np.min(q_diag)),np.min(q_diag),q_ampl)
                return(results,"")
        return((None,None,None,None),"no quadrupole information found")

    def set_quadrupole(self):
        self.qpol, self.error["qpole"] = self.get_quadrupole()

    def get_homolumo(self): # homo,lumo energies and derived values of last job in file
        loglines=self.lines
        error = ""
        for line in loglines[::-1]:
            if  homo_pattern.search(line):
                homo = float(str.split(line)[-1])
                lumo = float(str.split(loglines[loglines.index(line)+1])[4])
                mu =  (homo+lumo)/2 # chemical potential / negative of molecular electronegativity
                eta = lumo-homo     # hardness/softness
                omega = str(round(mu**2/(2*eta),5))  # electrophilicity index
                return((float(homo),float(lumo) ,float(round(mu,5)) ,float(round(eta,5)) ,float(omega)),error)
                return()
        error = "no orbital information found;;"
        return ([None]*5,error )

    def set_homolumo(self):
        self.homo, self.error["homo"] = self.get_homolumo()

    def get_dipole(self,no=0):   #no muss nocheinmal überprüft werden
        filecont = self.lines
        if no != -1:
            for i in range(len(filecont)-1,0,-1):
                if dipole_pattern in filecont[i]:
                    dipole = str.split(filecont[i+1])[-1]
                    return(float(dipole),"")
        if no == -1:
            for i in range(len(filecont)-1):
                if dipole_pattern in filecont[i]:
                    no += 1
                    dipole = str.split(filecont[i+1])[-1]
                    if no == 2: return(float(dipole),"")
        error = "no dipole information found;"
        return (None, error    )

    def set_dipole(self):
        self.dipole, self.error["dipol"] = self.get_dipole()


    def get_e_hf(self,no):
        stream = self.outstreams[no]
        e_hf = ""
        for item in stream:
            if "HF=" in item:
                e_hf = float(item[3:] )
        return(e_hf,"")

    def set_e(self):
        l = np.array([self.get_e_hf(i) for i in range(0,len(self.outstreams)) ])
        self.e = [float(i) for i in  list(l[:,0])]
        self.error["e"] = list(l[:,1])


    #_________________________________get and set functions_______________________________________________________________________________















    def get_dictionary(self):    # returns the dict of measurments
        return self.quantitiydict

    def get_error(self):   #  returns the dict of erromessages
        return self.error

    def get_bothdicts(self): # returns both
        return self.quantitiydict , self.error

    def get_name_ox_method(self):
        return self.name.split("_")

    def set_name_ox_method(self):   # identify name, multiplicity and soloution by the name of the file
        self.pdb, oxspin, self.method = self.get_name_ox_method()
        self.spin,self.ox = {"01":[1,0], "05":[5,0] ,"12":[2,1] ,"16":[6,1]}[oxspin]
        self.error["pdb"], self.error["ox"], self.error["method"] = ["","",""]

    def get_DataFrame(self):
        index1 = self.lines.index(' Number     Number       Type             X           Y           Z\n' )+2
        index2 = index1 + int([i for i in self.lines if "NAtoms=" in i][0].split()[1])
        coordinates=self.lines[index1:index2]
        df = helpers.coordtomatrix(coordinates)
        return df

    def set_DataFrame(self):
        self.df = self.get_DataFrame()




    def __init__(self, pwd , name = "1f1f"):
        """ Physical_uantity(pwd, name) """
        filename = name +".log"
        filename = pwd
        self.name = filename[filename.rindex("/")+1:][:-4]
        with open(filename) as f:
            self.lines = f.readlines()
        self.error={}
    #_______________________set variables_________________________________________________________________
        try:
            self.set_DataFrame()
        except:
            self.df=pd.DataFrame({None:[None]})

        self.set_outstreams()
        self.set_e()
        self.set_edisp()
        self.set_dipole()
        self.set_homolumo()
        self.set_quadrupole()
        #self.set_distances()
        self.set_polarizability()
        self.set_name_ox_method()
    #  self.set_planeangle()
    # self.set_dihedrals()
        #self.set_angles()

    #_______________________set variables_________________________________________________________________
    #_______________________Make the dict_________________________________________________________________
        self.quantitiydict = {}
        csvdict={}
        # self.quantitiydict["Coordinates"] = self.df.to_numpy()
    # self.quantitiydict["Distances"] = self.distances.to_numpy()
        self.quantitiydict["pdb"] = self.pdb.upper()
        self.quantitiydict["Ox"] = self.ox
        self.quantitiydict["spin"] = self.spin
        self.quantitiydict["method"] = self.method
        self.quantitiydict["e"] = [np.mean(self.e)]
        self.quantitiydict["edisp"]   = self.edisp
        self.quantitiydict["homo"] = self.homo
        self.quantitiydict["dipole"] = self.dipole
        self.quantitiydict["qpole"] = self.qpol
        self.quantitiydict["polar-iso"] = self.polar["iso"]
        self.quantitiydict["polar-aniso"] = self.polar["aniso"]
        self.quantitiydict["method"] = self.method


class onecsv:
    """lops over all files in logfilessplit, takes the values with Physical_Quanitity and creates a csv-table"""
    dictoferrors = {}

    csv = "tables/calculated.csv"

    def get_csv(self):
        for k in self.dictoflogs.keys():
            dictofvalues=self.dictoflogs[k]
            for e,E in enumerate( dictofvalues["e"]):
                if (len(dictofvalues["e"]) ==1):
                    titel = k
                else:
                    titel = k+"_"+{0:"a",1:"b",2:"c",3:"d",4:"e"}[e]
                number = e
                got_answer = False

                arbeitsdict= helpers.dicttodelete(dictofvalues,number)
                data=list(arbeitsdict.values())
                #print(self.df.to_string())
                #print(f"\n\n{data}")
                self.df[titel]=data
        return

    def make_dict(self):    # takes the measuremnts using Physical_quantity()
        for k in glob.glob('database/logfilessplit/*.log'):
            i=k[k.rindex('/')+1:][:-4]
            if k in self.df.columns:
                print(k ,"schon erhalten")
            else:
                try:
                    self.dictoflogs[i], self.dictoferrors[i] = Physical_quantity(k).get_bothdicts()
                except:
                    print( i, "failed")



    def save_csv(self):
        dft = self.df.T
        header_row = dft.iloc[0]
        df_calc = pd.DataFrame(dft.values[1:], columns=header_row)
        if df_calc.iloc[0,0]==df_calc.columns[0]:
            try:
                df_calc = df_calc.drop([0])
            except:
                2
        df_calc.set_index("pdb")
        df_calc.to_csv("tables/calculated.csv")

    def __init__(self):
        self.dictoflogs = {}



        newindex = ['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","diff","elekphil","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"]
        self.df =pd.DataFrame({ "example":['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","diff","elekphil","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"] }, index=newindex) #creates a new Dataframe
        self.make_dict()
        self.get_csv()
        self.save_csv()
        return


