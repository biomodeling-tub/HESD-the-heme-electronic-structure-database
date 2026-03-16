--Link5--
%nprocs=12
%mem=16GB
%oldchk=1aed_II_S.chk
%chk=1aed_II_Sre.chk
# apfd/gen scf=maxcycle=999 geom=check

1aed_III_D

1 2

C N O H 0
6-31G*
****
Fe 0
TZVP
****
