%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cc5_II_S.chk
%chk=1cc5_II_Sre.chk
# apfd/gen scf=maxcycle=999 geom=check

1cc5_II_Q

0 5

C N H 0
6-31G*
****
Fe 0
TZVP
****
