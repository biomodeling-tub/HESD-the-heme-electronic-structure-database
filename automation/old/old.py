def __create_jobscript(pdb_id, spin_state):
    jobtext=f"""
    #!/bin/bash --login
    #$ -cwd
    #$ -pe mp 4
    #$ -N gaussian{pdb_id}{spin_state}
    #$ -o gaussian.out
    #$ -j y
    #$ -l h_rt=14400
    #$ -m ea
    #$ -M elizaveta.zhartovska@campus.tu-berlin.de  

    module load g16

    declare -x GAUSS_SCRDIR="${{SCRATCHDIR}}"

    #echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"
    #echo "Path is: $PATH"
    #echo "LdPath is: $LD_LIBRARY_PATH"
    #echo "turbodir is: $TURBODIR"
    # which g03
    #echo "g16rootis: $g16root"
    #echo "GAUSS_SCRDIR is: $GAUSS_SCRDIR"
    #echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"

    STOREDIR=$1
    COMNAME=$2

    if [ ! -s ${{STOREDIR}}/${{COMNAME}} ]; then
    echo "Error: input file ${{STOREDIR}}/${{COMNAME}} does not exist!"
    exit 1
    fi

    SCRATCHDIR=$(mktemp -d /scratch/zhartovska_gaussian.XXXXXXX)
    chmod go+rwx ${{SCRATCHDIR}}
    cd ${{SCRATCHDIR}}
    cp ${{STOREDIR}}/*.chk ${{SCRATCHDIR}}/
    cp ${{STOREDIR}}/*.com ${{SCRATCHDIR}}/

    g16 ${{STOREDIR}}/${{COMNAME}}
    sleep 5

    mv -f ${{SCRATCHDIR}}/*.chk ${{STOREDIR}}/

    ### cp -a ${{SCRATCHDIR}}  ${{STOREDIR}}/${{COMNAME}}_scratch

    ls -al ${{SCRATCHDIR}}

    cd ${{STOREDIR}}
    #formchk *.chk
    rm -r ${{SCRATCHDIR}}
    """
    jobtext=inspect.cleandoc(jobtext)
    automation.__write_jobscript(pdb_id=pdb_id, spin_state=spin_state, jobtext=jobtext)
    return

def __create_jobscript_re(pdb_id, spin_state):
    jobtext=f"""
    #!/bin/bash --login
    #$ -cwd
    #$ -pe mp 4
    #$ -N gaussian{pdb_id}{spin_state}re
    #$ -o gaussian.out
    #$ -j y
    #$ -l h_rt=14400
    #$ -m ea
    #$ -M elizaveta.zhartovska@campus.tu-berlin.de  

    module load g16

    declare -x GAUSS_SCRDIR="${{SCRATCHDIR}}"

    #echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"
    #echo "Path is: $PATH"
    #echo "LdPath is: $LD_LIBRARY_PATH"
    #echo "turbodir is: $TURBODIR"
    # which g03
    #echo "g16rootis: $g16root"
    #echo "GAUSS_SCRDIR is: $GAUSS_SCRDIR"
    #echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"

    STOREDIR=$1
    COMNAME=$2

    if [ ! -s ${{STOREDIR}}/${{COMNAME}} ]; then
    echo "Error: input file ${{STOREDIR}}/${{COMNAME}} does not exist!"
    exit 1
    fi

    SCRATCHDIR=$(mktemp -d /scratch/zhartovska_gaussian.XXXXXXX)
    chmod go+rwx ${{SCRATCHDIR}}
    cd ${{SCRATCHDIR}}
    cp ${{STOREDIR}}/*.chk ${{SCRATCHDIR}}/
    cp ${{STOREDIR}}/*.com ${{SCRATCHDIR}}/

    g16 ${{STOREDIR}}/${{COMNAME}}
    sleep 5

    mv -f ${{SCRATCHDIR}}/*.chk ${{STOREDIR}}/

    ### cp -a ${{SCRATCHDIR}}  ${{STOREDIR}}/${{COMNAME}}_scratch

    ls -al ${{SCRATCHDIR}}

    cd ${{STOREDIR}}
    #formchk *.chk
    rm -r ${{SCRATCHDIR}}
    """ 
    jobtext=inspect.cleandoc(jobtext)
    automation.__write_jobscript_re(pdb_id=pdb_id, spin_state=spin_state, jobtext=jobtext)
    return


def __write_jobscript(pdb_id, spin_state, jobtext):
    if not os.path.exists(f"{pdb_id}/{pdb_id}{spin_state}/"):
        os.mkdir(f"{pdb_id}/{pdb_id}{spin_state}/")
    bash=open(f"{pdb_id}/{pdb_id}{spin_state}/run_gaussian.sh", "w").write(jobtext)
    return

def __write_jobscript_re(pdb_id, spin_state, jobtext):
    if not os.path.exists(f"{pdb_id}/{pdb_id}{spin_state}re/"):
        os.mkdir(f"{pdb_id}/{pdb_id}{spin_state}re/")
    bash=open(f"{pdb_id}/{pdb_id}{spin_state}re/{pdb_id}{spin_state}re.job", "w").write(jobtext)
    return


def __create_jobscripts_all(pdb_ids, spin_states):
    for pdb_id in pdb_ids:
        for spin_state in spin_states:
            automation.__create_jobscript(pdb_id=pdb_id, spin_state=spin_state)
    for pdb_id in pdb_ids:
        for spin_state in spin_states[1:]:
            automation.__create_jobscript_re(pdb_id=pdb_id, spin_state=spin_state)
    return

