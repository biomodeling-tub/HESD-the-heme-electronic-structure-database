%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=2gsm_II_S.chk
%chk=2gsm_II_Sre.chk
# apfd/gen scf=maxcycle=999 geom=check

2gsm_III_D

1 2

C N O H 0
6-31G*
****
Fe 0
TZVP
****


