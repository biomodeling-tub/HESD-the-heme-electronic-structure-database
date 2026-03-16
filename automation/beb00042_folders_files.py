import sys, os
import subprocess as proc
import textwrap
import inspect
import shutil
from datetime import datetime
from pdb_list import pdb_list


pdb_list = ['102m','1a2f','1a6k','1a6m','1a6n','1aa4','1ac8','1aej','1aek','1aem','1aen','1aeo','1aeq','1aet','1aeu','1aom','1arp',
    '1arx','1ash','1bek','1bem','1beq','1bgp','1bva','1bz6','1c52','1c75','1cc5','1ccb','1cce','1ccj','1cck','1ccl','1ccr',
    '1cgn','1cgo','1ch1','1ch2','1ch5','1ch7','1cie','1cif','1cih','1cik','1ck6','1cmt','1cmu','1co6','1co8','1co9','1cot',
    '1cpf','1cpg','1cpo','1cpq','1cpr','1cpt','1cpw','1crh','1crj','1cry','1csu','1csv','1csw','1csx','1ctj','1ctz','1cxc',
    '1cxy','1cyf','1cyo','1dcc','1dp6','1drm','1dt1','1dt6','1dti','1dve','1dws','1dwt','1dxd','1dy7','1e29','1e93','1ebe',
    '1eca','1ecd','1ehb','1ehe','1ehf','1ehg','1ehk','1es1','1eup','1f63','1f65','1fcb','1fcs','1gdv','1geb','1ged','1gjm',
    '1gwe','1gwf','1gwh','1gwo','1gwt','1gwu','1gza','1gzb','1h55','1h58','1h5a','1h5c','1h5d','1h5f','1h5g','1h5h','1h5i',
    '1h5k','1h5l','1h5m','1h6n','1h7k','1hrc','1hrm','1hsr','1hsy','1hzu','1hzv','1io3','1irv','1iwj','1iwk','1iyn','1j77',
    '1jdr','1jf3','1jf4','1jp8','1jp9','1kb0','1kfr','1krj','1kxm','1kzm','1lco','1ldc','1lgf','1lhs','1lqx','1lr6','1ls9',
    '1lsw','1ltd','1ltw','1lwl','1m20','1m2i','1m2m','1m59','1m7v','1m85','1mbs','1mgn','1mkq','1mlh','1mlk','1mll','1mln',
    '1mlo','1mlr','1mn2','1mnh','1mob','1mod','1mti','1myt','1mz4','1n4g','1n5u','1n6b','1nr6','1ntf','1nz5','1obm','1ofj',
    '1ofk','1phc','1q5d','1q5e','1qo4','1r9o','1rap','1raq','1re9','1rf9','1tqn','1u13','1u75','1u7s','1u9u','1ubb','1ue8',
    '1ulw','1uvy','1vag','1vgi','1vxa','1vxb','1vxd','1w0e','1wla','1x3k','1x46','1xbn','1xch','1xj3','1y5l','1yea','1yiq',
    '1yma','1ymb','1yq3','1yq4','1yrd','1ywd','1yyg','1yzp','1yzr','1z8p','1z8q','1zby','1zbz','1zoy','1zvi','2acp','2al0',
    '2all','2amm','2an2','2anz','2as3','2as6','2b0z','2b4z','2bcn','2bgv','2bh4','2bh5','2bk9','2bli','2c2c','2cah','2ccp',
    '2ce0','2cep','2civ','2ciw','2cix','2cj0','2cj1','2cl4','2cmn','2cp4','2cpo','2cyp','2d09','2d0e','2d0s','2eup','2euq',
    '2eur','2eus','2fr7','2g5g','2gdm','2ggn','2ghe','2hbg','2hi4','2hz1','2hz2','2ia8','2icv','2iiz','2ivf','2j0p','2j18',
    '2j19','2jjn','2jjo','2mbw','2mge','2mgg','2mgl','2mgm','2nni','2nnj','2nrl','2nrm','2nwb','2o86','2oh8','2oha','2ohb',
    '2owh','2pt3','2q8q','2q9f','2qrb','2r5l','2r79','2r7a','2rbt','2rbu','2rbv','2rbw','2rbx','2rbz','2rc1','2spm','2spn',
    '2uuq','2v1h','2v1j','2v1k','2v2e','2v7i','2v7j','2veb','2vka','2vku','2vlx','2vly','2vlz','2vn0','2vyw','2vyz','2w6w',
    '2w6x','2w6y','2wgy','2wm5','2wtg','2x08','2x5l','2xbk','2xc3','2xif','2xj6','2xj8','2xkg','2xkh','2xld','2xlh','2xm0',
    '2xm4','2xvx','2xvz','2y5a','2y6a','2y6b','2yca','2ycc','2ycg','2ykz','2yli','2yqb','2z3u','2z6f','2z6s','2z6t','2zby',
    '2zbz','2zon','2zoo','2zwj','2zwt','2zwu','351c','3a1l','3a2g','3a5a','3abb','3aeb','3awm','3awp','3b47','3buj','3c2c',
    '3ccp','3ccx','3cp5','3cqv','3cu4','3cv8','3dam','3dan','3dbm','3dmi','3dp5','3e2n','3e2o','3e4y','3e5l','3ecl','3eh4',
    '3eh5','3emm','3eqm','3erh','3exb','3fkg','3fm1','3fm4','3fmu','3g5f','3gck','3h58','3hb6','3hen','3i6n','3iqb','3iw0',
    '3k9z','3kw4','3l1m','3l61','3l62','3l63','3m23','3m25','3m26','3m28','3m2a','3m2c','3m2e','3m2h','3m2i','3m3a','3m3b',
    '3m5q','3m8m','3nv6','3o89','3oia','3ovu','3p6n','3p6o','3p6q','3p6t','3p6v','3p6w','3p6x','3pc2','3pc3','3pm0','3py4',
    '3qf1','3qjs','3qjv','3qns','3qzx','3qzz','3rgk','3rwl','3s8f','3s8g','3sfd','3sfe','3tm9','3u3e','3ua1','3uas','3uba',
    '3vno','3vnw','3voo','3vp5','3vrm','3vtj','3vu3','3vxi','3vxj','3vym','3wc8','3wec','3wfw','3wxo','3x32','3x33','3x35',
    '3zbm','3zcg','3zch','3ziy','451c','4a6z','4a71','4a7m','4am4','4b4y','4blk','4bll','4bln','4blx','4bly','4bm0','4bm3',
    '4bm4','4c44','4ccx','4cda','4cdy','4cjo','4cp4','4cpp','4czo','4czp','4czq','4czr','4d3m','4d3n','4d3o','4d3t','4d3u',
    '4d3v','4d7h','4d7i','4d7j','4dc7','4dy9','4e2p','4eic','4eid','4f68','4f6d','4f6i','4faa','4fcs','4fdq','4fef','4fia',
    '4fwx','4fwy','4fwz','4g05','4g46','4g47','4g70','4g71','4g7l','4g7q','4g8u','4g8w','4gp5','4h07','4h0b','4h1n','4h8p',
    '4h8q','4hhs','4hka','4hpc','4ict','4ips','4iq9','4jcg','4jm6','4jm8','4jmb','4jms','4jmt','4jmw','4jmz','4jn0','4jpt',
    '4jpu','4jqk','4kpa','4ksz','4ktj','4ktk','4ktl','4kvk','4l1y','4l49','4l4b','4l54','4lpi','4lxj','4mbn','4msf','4njb',
    '4ns2','4nva','4nvb','4nvc','4nvd','4nve','4nvf','4nvg','4nvi','4nvj','4nvl','4nvm','4nvn','4nvo','4nxc','4o4t','4o4z',
    '4o6t','4o6u','4oek','4of9','4oq7','4oqr','4pnj','4pq6','4pqb','4pqc','4qau','4qjq','4rm4','4tpn','4tpo','4tvf','4tx3',
    '4tyx','4u9e','4uax','4ubs','4ug5','4ug6','4ug7','4ug9','4uga','4ugc','4ugd','4uge','4ugf','4ugg','4ugj','4ugk','4ugl',
    '4ugm','4ugo','4ugp','4ugq','4ugr','4ugs','4ugt','4ugx','4ugy','4w7o','4wch','4wgz','4wpz','4wq0','4xmc','4xmd','4xv6',
    '4xv7','4xv8','4y4s','4yxd','4yzr','4z5q','4zf6','4zfb','4zid','5a1p','5a5i','5a5j','5abn','5abo','5al9','5aog','5b4z',
    '5b51','5b72','5b82','5b85','5ccp','5cje','5cn4','5cn5','5cn6','5cn7','5cn9','5cne','5cnf','5cng','5cp4','5cpp','5d5r',
    '5d6m','5eet','5ejt','5ejx','5eoh','5esn','5eyj','5eys','5f0b','5f2a','5fne','5fvo','5g65','5g66','5g68','5g69','5g6a',
    '5g6b','5g6c','5g6e','5g6f','5g6g','5g6i','5g6j','5g6k','5g6l','5g6m','5g6n','5g6p','5g6q','5gnl','5gxg','5hav','5hlq',
    '5ibi','5ibj','5ik1','5ikd','5jkv','5jkw','5jl6','5jp7','5jqr','5jsl','5kdb','5kke','5l1o','5l1p','5l1r','5l1s','5l1t',
    '5l1u','5l1v','5lft','5lth','5m0n','5m3s','5mfa','5mjc','5mxy','5mxz','5n1t','5ncb','5nvi','5nws','5o17','5o18','5o27',
    '5obo','5ocb','5ocf','5omr','5op9','5oxu','5oy2','5t6q','5u5u','5u5v','5u5w','5u5x','5u5y','5u5z','5u60','5u61','5u6u',
    '5u6w','5uec','5ug2','5utc','5uvb','5vcc','5vcg','5via','5vju','5vnu','5vws','5w0c','5wjk','5wv3','5x23','5x24','5x7e',
    '5xkv','5xkw','5xnt','5xw2','5xxi','5xzi','5yce','5ych','5ylw','5ym3','5yqa','5ysm','5z25','5z7e','5zeo','5zm9','5zze',
    '5zzg','6a3k','6a4y','6a7j','6bdd','6bde','6bmt','6cp4','6cpp','6cuk','6d45','6dyi','6ekw','6ekx','6f0b','6f0c','6ff5',
    '6fyj','6gd6','6gd7','6geq','6gii','6gk6','6gmf','6hqk','6hqn','6hqo','6hqq','6hqs','6hqt','6ied','6j95','6jsa','6jt6',
    '6k9i','6k9j','6krc','6l8j','6lco','6m4q','6m4s','6m7e','6m7l','6m8f','6mi0','6mjm','6mx5','6myo','6mys','6myt','6myu',
    '6nzx','6o0a','6onq','6oo9','6oob','6oow','6pqd','6prr','6prs','6pt7','6qq0','6qq1','6qq2','6r1q','6r3w','6ro8','6rqd',
    '6rqe','6rsj','6t0g','6tet','6tev','6u3k','6u87','6u97','6unn','6upg','6vdq','6vjx','6vxv','6vz6','6wj6','6wzc','6x8x',
    '6xa3','6xaj','6xcx','6xk6','6xk7','6xmc','6y1t','6y2y','6yci','6ycl','6yco','6ycp','6zmq','7c74','7ccf','7ccp','7cen',
    '7cez','7cka','7cl7','7cl8','7cl9','7d52','7d5i','7de5','7di3','7dlh','7dls','7dmr','7kcs','7ks8','7l3y','7mhy','7n14',
    '7nhp','7nhq','7nqn','7oo5','7ow9','7red','7reh','7rl2','7sa3']


class automation():
    def __init__(self):
        """
        Class that creates the necessary files for running jobs automatically except for the .xyz coordinate files.
        """
        spin_states=["01", "05", "12", "16"]
        automation.__create_runtime_folders(pdb_list=pdb_list, spin_states=spin_states)
        for pdb_id in pdb_list:
            if not os.path.exists(f"PDB/{pdb_id}/{pdb_id}_g16.xyz"):
                print(f"{pdb_id}_g16.xyz does not exist. Skipping {pdb_id}.")
            else:
                automation.__create_gaussian_scripts_all(pdb_id=pdb_id)
        return

    def __create_runtime_folders(pdb_list, spin_states):
        """Checks if folder names exist, if not, creates folders."""
        for pdb_id in pdb_list:
            for spin_state in spin_states:
                if not os.path.exists(f"{pdb_id}"):
                    os.mkdir(f"{pdb_id}")
                if not os.path.exists(f"{pdb_id}/{pdb_id}{spin_state}/"):
                    os.mkdir(f"{pdb_id}/{pdb_id}{spin_state}/")
        return

    def error_warning_job(pdb_id):
        sh="""#!/bin/bash
        #SBATCH --time=00:01:00
        #SBATCH --mail-type=begin # send email when job begins
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=standard96:test
        #SBATCH --nodes=1
        #SBATCH --ntasks=1

        sleep(5)

        """
        error_warning_file=inspect.cleandoc(sh)+"\n"
        spin_states=["01", "05", "12", "16"]
        for spin_state in spin_states:
            if not os.path.exists(f"{pdb_id}/{pdb_id}{spin_state}/"):
                os.mkdir(f"{pdb_id}/{pdb_id}{spin_state}/")
            run=open(f"{pdb_id}/{pdb_id}{spin_state}/{pdb_id}{spin_state}error.sh", "w").write(error_warning_file)
        return run

    def __create_run_script_01(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}01.com
        cd ../
        rsync -avz {pdb_id}01/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "01"
        sbatch g16_{pdb_id}L2.sh

        """
        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}01/"):
            os.mkdir(f"{pdb_id}/{pdb_id}01/")
        run=open(f"{pdb_id}/{pdb_id}01/g16_{pdb_id}01.sh", "w").write(run01)
        return run

    def __create_run_script_Link2(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}_L2.com
        cd ../
        rsync -avz {pdb_id}01/ /../../scratch/projects/beb00042/{pdb_id}/
        rsync -avz {pdb_id}01/{pdb_id}_II_S.chk {pdb_id}05/
        rsync -avz {pdb_id}01/{pdb_id}_II_S.chk {pdb_id}12/
        python ../beb00042_agent.py "{pdb_id}" "01" "L2"
        rm -r {pdb_id}01/
        sbatch {pdb_id}05/g16_{pdb_id}05re.sh
        sbatch {pdb_id}12/g16_{pdb_id}12re.sh

        """
        runL2=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}01/"):
            os.mkdir(f"{pdb_id}/{pdb_id}01/")
        run=open(f"{pdb_id}/{pdb_id}01/g16_{pdb_id}L2.sh", "w").write(runL2)
        return run

    def __create_run_script_05re(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=2:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}05re.com
        cd ../
        rsync -avz {pdb_id}05/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "05" "re"
        sbatch {pdb_id}05/g16_{pdb_id}05.sh

        """

        run05re=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}05/"):
            os.mkdir(f"{pdb_id}/{pdb_id}05/")
        run=open(f"{pdb_id}/{pdb_id}05/g16_{pdb_id}05re.sh", "w").write(run05re)
        return run

    def __create_run_script_05(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}05.com
        cd ../
        rsync -avz {pdb_id}05/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "05"
        sbatch {pdb_id}05/g16_{pdb_id}L4.sh

        """

        run05=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}05/"):
            os.mkdir(f"{pdb_id}/{pdb_id}05/")
        run=open(f"{pdb_id}/{pdb_id}05/g16_{pdb_id}05.sh", "w").write(run05)
        return run

    def __create_run_script_Link4(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09


        g16 {pdb_id}05_L4.com
        python ../beb00042_chk_agent.py "{pdb_id}" "05" "L4"
        cd ../
        rsync -avz {pdb_id}05/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "05" "L4"
        rm -r {pdb_id}05/

        """

        runL4=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}05/"):
            os.mkdir(f"{pdb_id}/{pdb_id}05/")
        run=open(f"{pdb_id}/{pdb_id}05/g16_{pdb_id}L4.sh", "w").write(runL4)
        return run

    def __create_run_script_12re(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=2:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}12re.com
        cd ../
        rsync -avz {pdb_id}12/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "12" "re"
        sbatch {pdb_id}12/g16_{pdb_id}12.sh

        """
        run12re=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}12/"):
            os.mkdir(f"{pdb_id}/{pdb_id}12/")
        run=open(f"{pdb_id}/{pdb_id}12/g16_{pdb_id}12re.sh", "w").write(run12re)
        return

    def __create_run_script_12(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}12.com
        cd ../
        rsync -avz {pdb_id}12/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "12"
        sbatch {pdb_id}12/g16_{pdb_id}L6.sh

        """
        run12=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}12/"):
            os.mkdir(f"{pdb_id}/{pdb_id}12/")
        run=open(f"{pdb_id}/{pdb_id}12/g16_{pdb_id}12.sh", "w").write(run12)
        return

    def __create_run_script_Link6(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=16G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}12.com
        python ../beb00042_chk_agent.py "{pdb_id}" "12" "L6"
        cd ../
        rsync -avz {pdb_id}12/ /../../scratch/projects/beb00042/{pdb_id}/
        cp {pdb_id}12/{pdb_id}_III_D.chk {pdb_id}16/
        python ../beb00042_agent.py "{pdb_id}" "12" "L6"
        rm -r {pdb_id}12/
        sbatch {pdb_id}16/g16_{pdb_id}16re.sh

        """

        runL6=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}12/"):
            os.mkdir(f"{pdb_id}/{pdb_id}12/")
        run=open(f"{pdb_id}/{pdb_id}12/g16_{pdb_id}L6.sh", "w").write(runL6)
        return run

    def __create_run_script_16re(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=2:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}16re.com
        cd ../
        rsync -avz {pdb_id}16/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "16" "re"
        sbatch {pdb_id}16/g16_{pdb_id}16.sh

        """

        run16re=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}16/"):
            os.mkdir(f"{pdb_id}/{pdb_id}16/")
        run=open(f"{pdb_id}/{pdb_id}16/g16_{pdb_id}16re.sh", "w").write(run16re)
        return

    def __create_run_script_16(pdb_id):
        g16=f"""#!/bin/bash
        #SBATCH --time=4:00:00
        #SBATCH --mail-type=end # send email when job ends
        #SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
        #SBATCH --partition=gpu-a100
        #SBATCH --nodes=1
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09

        g16 {pdb_id}16.com
        cd ../
        rsync -avz {pdb_id}16/ /../../scratch/projects/beb00042/{pdb_id}/
        python ../beb00042_agent.py "{pdb_id}" "16"
        sbatch {pdb_id}16/g16_{pdb_id}L8.sh

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
        #SBATCH -A beb00042
        #SBATCH --mem=32G
        #SBATCH --ntasks=72
        #SBATCH --gres=gpu:4

        module load cuda/11.8
        module load gaussian/16.C02
        module load anaconda3/2023.09


        g16 {pdb_id}16L8.com
        python beb00042_chk_agent.py "{pdb_id}" "16" "L8"
        cd ../
        cp -r {pdb_id}16/ /../../scratch/projects/beb00042/{pdb_id}/
        cd ../
        python beb00042_agent.py "{pdb_id}" "16" "L8"
        rm -r {pdb_id}/

        """

        run01=inspect.cleandoc(g16)+"\n"
        if not os.path.exists(f"{pdb_id}/{pdb_id}16/"):
            os.mkdir(f"{pdb_id}/{pdb_id}16/")
        run=open(f"{pdb_id}/{pdb_id}16/g16_{pdb_id}L8.sh", "w").write(run01)
        return run

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
        with open(f"PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords = '\n'.join(coords)
        com01=inspect.cleandoc(link0)+"\n"+coords+"\n\n"+inspect.cleandoc(link1)+"\n\n\n"
        link2=inspect.cleandoc(link2)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}01/{pdb_id}01.com", "w").write(com01)
        lin=open(f"{pdb_id}/{pdb_id}01/{pdb_id}L2.com", "w").write(link2)
        return

    def __create_gaussian_scripts_05(pdb_id):
        link3=f"""%mem=16GB
        %CPU=0-71
        %gpucpu=0,1,2,3=0,1,2,3
        %oldchk={pdb_id}_II_Sre.chk
        %chk={pdb_id}_II_Qn.chk
        # apfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

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
        link4=inspect.cleandoc(link4)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}05/{pdb_id}05.com", "w").write(com05)
        lin=open(f"{pdb_id}/{pdb_id}05/{pdb_id}L4.com", "w").write(link4)
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
        link6=inspect.cleandoc(link6)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}12/{pdb_id}12.com", "w").write(com12)
        lin=open(f"{pdb_id}/{pdb_id}12/{pdb_id}L6.com", "w").write(link6)
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
        link8=inspect.cleandoc(link8)+"\n\n\n"
        com=open(f"{pdb_id}/{pdb_id}16/{pdb_id}16.com", "w").write(com16)
        lin=open(f"{pdb_id}/{pdb_id}16/{pdb_id}L8.com", "w").write(link8)
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

        with open(f"PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        #coords_str = '\n'.join(coords)
        link_to_use = (link3re3 if any("O" in line and "S" in line for line in coords)
               else link3re2 if any("S" in line and "O" not in line for line in coords)
               else link3re if any("O" in line and "S" not in line for line in coords)
               else link3re1)
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

        with open(f"PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        #coords_str = '\n'.join(coords)
        link_to_use = (link5re3 if any("O" in line and "S" in line for line in coords)
               else link5re2 if any("S" in line and "O" not in line for line in coords)
               else link5re if any("O" in line and "S" not in line for line in coords)
               else link5re1)
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

        with open(f"PDB/{pdb_id}/{pdb_id}_g16.xyz", "r") as file:
            coords = file.read().splitlines()[2:]
        coords_str = '\n'.join(coords)
        link_to_use = (link7re3 if any("O" in line and "S" in line for line in coords)
               else link7re2 if any("S" in line and "O" not in line for line in coords)
               else link7re if any("O" in line and "S" not in line for line in coords)
               else link7re1)
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
        automation.__create_run_script_Link2(pdb_id=pdb_id)

        automation.__create_run_script_05re(pdb_id=pdb_id)
        automation.__create_run_script_05(pdb_id=pdb_id)
        automation.__create_run_script_Link4(pdb_id=pdb_id)

        automation.__create_run_script_12re(pdb_id=pdb_id)
        automation.__create_run_script_12(pdb_id=pdb_id)
        automation.__create_run_script_Link6(pdb_id=pdb_id)

        automation.__create_run_script_16re(pdb_id=pdb_id)
        automation.__create_run_script_16(pdb_id=pdb_id)
        automation.__create_run_script_Link8(pdb_id=pdb_id)

        automation.error_warning_job(pdb_id=pdb_id)
        return

if __name__=="__main__":
    inst = automation()
