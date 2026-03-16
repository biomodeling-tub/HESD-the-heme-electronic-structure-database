%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=2gsm_III_D.chk
%chk=2gsm_III_Dre.chk
# apfd/gen scf=maxcycle=999 geom=check

2gsm_III_H

1 6

C N O H 0
6-31G*
****
Fe 0
TZVP
****


