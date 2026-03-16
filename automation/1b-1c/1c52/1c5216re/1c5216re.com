%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1c52_III_D.chk
%chk=1c52_III_Dre.chk
# apfd/gen scf=maxcycle=999 geom=check

1c52_III_H

1 6

C N H 0
6-31G*
****
Fe 0
TZVP
****
