%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cot_III_D.chk
%chk=1cot_III_Hn.chk
# apfd/def2TZVP scf=maxcycle=999 geom=check pop=nbo

1cot_III_H

1 6


