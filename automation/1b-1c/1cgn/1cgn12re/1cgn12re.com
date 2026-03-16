%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cgn_II_S.chk
%chk=1cgn_II_Sre.chk
# apfd/gen scf=maxcycle=999 geom=check

1cgn_III_D

1 2

C N H 0
6-31G*
****
Fe 0
TZVP
****
