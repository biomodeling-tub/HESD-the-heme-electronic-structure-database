%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1bem_III_D.chk
%chk=1bem_III_Hn.chk
# apfd/def2TZVP scf=maxcycle=999 geom=check pop=nbo

1bem_III_H

1 6


