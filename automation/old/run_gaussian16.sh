#!/bin/bash --login
#$ -cwd
#$ -pe mp 4
#$ -N gaussian
#$ -o gaussian.out
#$ -j y
#$ -l h_rt=604800
#$ -m ea
#$ -M elizaveta.zhartovska@campus.tu-berlin.de  

module load g16

declare -x GAUSS_SCRDIR="${SCRATCHDIR}"

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

if [ ! -s ${STOREDIR}/${COMNAME} ]; then
 echo "Error: input file ${STOREDIR}/${COMNAME} does not exist!"
 exit 1
fi

SCRATCHDIR=$(mktemp -d /scratch/zhartovska_gaussian.XXXXXXX)
chmod go+rwx ${SCRATCHDIR}
cd ${SCRATCHDIR}
cp ${STOREDIR}/*.chk ${SCRATCHDIR}/
cp ${STOREDIR}/*.com ${SCRATCHDIR}/

g16 ${STOREDIR}/${COMNAME}
sleep 5

mv -f ${SCRATCHDIR}/*.chk ${STOREDIR}/

### cp -a ${SCRATCHDIR}  ${STOREDIR}/${COMNAME}_scratch

ls -al ${SCRATCHDIR}

cd ${STOREDIR}
#formchk *.chk
rm -r ${SCRATCHDIR}

