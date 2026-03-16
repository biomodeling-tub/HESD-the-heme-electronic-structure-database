import sys, os
import subprocess as proc 
import textwrap
import inspect
import shutil
from datetime import datetime




class automation():
    
    def __init__(self):
        """ 
        Class that creates the necessary files for running jobs automatically except for the .xyz coordinate files. 
        """
        spin_states=["01", "05", "12", "16"]
        pdb_list = ['1ac8',
 '1aej',
 '1aek',
 '1aem',
 '1aen',
 '1aeo',
 '1aeq',
 '1aet',
 '1aeu',
 '1aom',
 '1arp',
 '1arx',
 '1ash']
        automation.__create_runtime_folders(pdb_list=pdb_list, spin_states=spin_states)
        for pdb_id in pdb_list:
            if not os.path.exists(f"../PDB/{pdb_id}/{pdb_id}_g16.xyz"):
                print(f"{pdb_id}_g16.xyz does not exist. Skipping {pdb_id}.")
            else:
                automation.__create_gaussian_scripts_all(pdb_id=pdb_id)
        return

    def __create_runtime_folders(pdb_list, spin_states):
        """
        Checks if folder names exist, if not, creates folders.  
        """
        for pdb_id in pdb_list:

            for spin_state in spin_states:
                if not os.path.exists(f"{pdb_id}"):
                    os.mkdir(f"{pdb_id}")
                if not os.path.exists(f"{pdb_id}/{pdb_id}{spin_state}/"):
                    os.mkdir(f"{pdb_id}/{pdb_id}{spin_state}/")
        return

    def __create_run_script_01(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                    
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}01.com

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}01/"):
            os.mkdir(f"{pdb_id}/{pdb_id}01/")
        run=open(f"{pdb_id}/{pdb_id}01/g16_{pdb_id}01.sh", "w").write(run01)
        return
    
    def __create_run_script_Link2(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                   
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}_Link2.com

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}01/"):
            os.mkdir(f"{pdb_id}/{pdb_id}01/")
        run=open(f"{pdb_id}/{pdb_id}01/g16_{pdb_id}_Link2.sh", "w").write(run01)
        return
    
    def __create_run_script_05(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                   
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}05.com

        """

        run05=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}05/"):
            os.mkdir(f"{pdb_id}/{pdb_id}05/")
        run=open(f"{pdb_id}/{pdb_id}05/g16_{pdb_id}05.sh", "w").write(run05)
        return
    
    def __create_run_script_Link4(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                     
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}_Link4.com

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}05/"):
            os.mkdir(f"{pdb_id}/{pdb_id}05/")
        run=open(f"{pdb_id}/{pdb_id}05/g16_{pdb_id}_Link4.sh", "w").write(run01)
        return
    
    def __create_run_script_12(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                    
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}12.com

        """

        run12=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}12/"):
            os.mkdir(f"{pdb_id}/{pdb_id}12/")
        run=open(f"{pdb_id}/{pdb_id}12/g16_{pdb_id}12.sh", "w").write(run12)
        return
    
    def __create_run_script_Link6(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                   
        #SBATCH --mem=16G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}_Link6.com

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}12/"):
            os.mkdir(f"{pdb_id}/{pdb_id}12/")
        run=open(f"{pdb_id}/{pdb_id}12/g16_{pdb_id}_Link6.sh", "w").write(run01)
        return
    
    def __create_run_script_16(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                  
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}16.com

        """

        run16=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}16/"):
            os.mkdir(f"{pdb_id}/{pdb_id}16/")
        run=open(f"{pdb_id}/{pdb_id}16/g16_{pdb_id}16.sh", "w").write(run16)
        return
    
    def __create_run_script_Link8(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100          
        #SBATCH --nodes=1                  
        #SBATCH --mem=32G                     
        #SBATCH --ntasks=72             
        #SBATCH --gres=gpu:4             

        module load cuda/11.8
        module load gaussian/16.C02

        g16 {pdb_id}_Link8.com

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}16/"):
            os.mkdir(f"{pdb_id}/{pdb_id}16/")
        run=open(f"{pdb_id}/{pdb_id}16/g16_{pdb_id}_Link8.sh", "w").write(run01)
        return

    def __create_gaussian_scripts_01(pdb_id):
        link0=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %chk={pdb_id}_II_O.chk
        # polar def2tzvp apfd int=(grid=ultrafine)

        {pdb_id}_II_S 

        0 1

        """

        link1=f"""

        --Link1--
        %mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_O.chk
        %chk={pdb_id}_II_Sn.chk
        # apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read pop=nbo scf=qc

        {pdb_id}_II_S_nbo

        0 1

        """

        link2=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_Sn.chk
        %chk={pdb_id}_II_S.chk
        # apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read scrf=(smd,solvent=Chloroform) 

        {pdb_id}_II_S_solv 

        0 1 

        """
        with open(f"../PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords = '\n'.join(coords)
        com01=inspect.cleandoc(link0)+"\n"+coords+"\n\n"+inspect.cleandoc(link1)+"\n\n\n"
        lin2=inspect.cleandoc(link2)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}01/{pdb_id}01.com", "w").write(com01)
        lin=open(f"{pdb_id}/{pdb_id}01/{pdb_id}_Link2.com", "w").write(lin2)
        return

    def __create_gaussian_scripts_05(pdb_id):
        link3=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_Sre.chk
        %chk={pdb_id}_II_Qn.chk
        # apfd/def2TZVP geom=check scf=maxcycle=999 geom=check pop=nbo

        {pdb_id}_II_Q

        0 5

        """

        link4=f"""%mem=8GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_Qn.chk
        %chk={pdb_id}_II_Q.chk
        # apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

        {pdb_id}_II_Q_solv

        0 5

        """

        com05=inspect.cleandoc(link3)+"\n\n\n"
        lin4=inspect.cleandoc(link4)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}05/{pdb_id}05.com", "w").write(com05)
        lin=open(f"{pdb_id}/{pdb_id}05/{pdb_id}_Link4.com", "w").write(lin4)
        return

    def __create_gaussian_scripts_12(pdb_id):
        link5=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_Sre.chk
        %chk={pdb_id}_III_Dn.chk
        # apfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

        {pdb_id}_III_D

        1 2

        """

        link6=f"""%mem=8GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_Dn.chk
        %chk={pdb_id}_III_D.chk
        # apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

        {pdb_id}_III_D_solv

        1 2

        """

        com12=inspect.cleandoc(link5)+"\n\n\n"
        lin6=inspect.cleandoc(link6)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}12/{pdb_id}12.com", "w").write(com12)
        lin=open(f"{pdb_id}/{pdb_id}12/{pdb_id}_Link6.com", "w").write(lin6)
        return
    
    def __create_gaussian_scripts_16(pdb_id):
        link7=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_Dre.chk
        %chk={pdb_id}_III_Hn.chk
        # apfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

        {pdb_id}_III_H

        1 6

        """

        link8=f"""%mem=8GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_Hn.chk
        %chk={pdb_id}_III_H.chk
        # apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

        {pdb_id}_III_H_solv

        1 6

        """

        com16=inspect.cleandoc(link7)+"\n\n\n"
        lin8=inspect.cleandoc(link8)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}16/{pdb_id}16.com", "w").write(com16)
        lin=open(f"{pdb_id}/{pdb_id}16/{pdb_id}_Link8.com", "w").write(lin8)
        return
    
    def __create_rerun_gaussian_scripts_05(pdb_id):
        link3re = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_II_Q

        0 5

        C N O H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        link3re1 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_II_Q

        0 5

        C N H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """

        link3re2 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_II_Q

        0 5

        C N H S 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """
        link3re3 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_II_Q

        0 5

        C N H S O 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """

        with open(f"../PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords_str = '\n'.join(coords)
        link_to_use = (link3re3 if any("O" in line and "S" in line for line in coords) 
               else link3re2 if any("S" in line and "O" not in line for line in coords) 
               else link3re1 if any("O" in line and "S" not in line for line in coords) 
               else link3re)
        com05re=inspect.cleandoc(link_to_use)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}05/{pdb_id}05re.com", "w").write(com05re)
        return
    
    def __create_rerun_gaussian_scripts_12(pdb_id):
        link5re = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_D

        1 2

        C N O H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        link5re1 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_D

        1 2

        C N H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 
        link5re2 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_D

        1 2

        C N H S 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """

        link5re3 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_S.chk
        %chk={pdb_id}_II_Sre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_D

        1 2

        C N H O S 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """

        with open(f"../PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords_str = '\n'.join(coords)
        link_to_use = (link5re3 if any("O" in line and "S" in line for line in coords) 
               else link5re2 if any("S" in line and "O" not in line for line in coords) 
               else link5re1 if any("O" in line and "S" not in line for line in coords) 
               else link5re)
        com12re=inspect.cleandoc(link_to_use)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}12/{pdb_id}12re.com", "w").write(com12re)
        return

    def __create_rerun_gaussian_scripts_16(pdb_id):
        link7re = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_D.chk
        %chk={pdb_id}_III_Dre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_H

        1 6

        C N O H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        link7re1 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_D.chk
        %chk={pdb_id}_III_Dre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_H

        1 6

        C N H 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        link7re2 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_D.chk
        %chk={pdb_id}_III_Dre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_H

        1 6

        C N H S 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        link7re3 = f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_III_D.chk
        %chk={pdb_id}_III_Dre.chk
        # apfd/gen scf=maxcycle=999 geom=check

        {pdb_id}_III_H

        1 6

        C N H S O 0
        6-31G*
        ****
        Fe 0
        TZVP
        ****


        """ 

        with open(f"../PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords_str = '\n'.join(coords)
        link_to_use = (link7re3 if any("O" in line and "S" in line for line in coords) 
               else link7re2 if any("S" in line and "O" not in line for line in coords) 
               else link7re1 if any("O" in line and "S" not in line for line in coords) 
               else link7re)
        com16re=inspect.cleandoc(link_to_use)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}16/{pdb_id}16re.com", "w").write(com16re)
        return

    def __create_gaussian_scripts_all(pdb_id):
        automation.__create_gaussian_scripts_01(pdb_id=pdb_id)
        automation.__create_gaussian_scripts_05(pdb_id=pdb_id)
        automation.__create_gaussian_scripts_12(pdb_id=pdb_id)
        automation.__create_gaussian_scripts_16(pdb_id=pdb_id)
        automation.__create_rerun_gaussian_scripts_05(pdb_id=pdb_id)
        automation.__create_rerun_gaussian_scripts_12(pdb_id=pdb_id)
        automation.__create_rerun_gaussian_scripts_16(pdb_id=pdb_id)
        automation.__create_run_script_01(pdb_id=pdb_id)
        automation.__create_run_script_05(pdb_id=pdb_id)
        automation.__create_run_script_12(pdb_id=pdb_id)
        automation.__create_run_script_16(pdb_id=pdb_id)
        automation.__create_run_script_Link2(pdb_id=pdb_id)
        automation.__create_run_script_Link4(pdb_id=pdb_id)
        automation.__create_run_script_Link6(pdb_id=pdb_id)
        automation.__create_run_script_Link8(pdb_id=pdb_id)
        return
        
inst = automation()
