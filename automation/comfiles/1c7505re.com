--Link3--
%nprocs=12
%mem=16GB
%oldchk=1c75_II_S.chk
%chk=1c75_II_Q.chk
# PBE1PBE/gen scf=maxcycle=999 geom=check

1c75_II_Q

0 5

C N O H 0
6-31G*
****
Fe 0
TZVP
****




