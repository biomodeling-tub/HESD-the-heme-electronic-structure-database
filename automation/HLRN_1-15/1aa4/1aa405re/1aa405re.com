--Link3--
%nprocs=12
%mem=16GB
%oldchk=1aa4_II_S.chk
%chk=1aa4_II_Sre.chk
# apfd/gen scf=maxcycle=999 geom=check

1aa4_II_Q

0 5

C N O H 0
6-31G*
****
Fe 0
TZVP
****
