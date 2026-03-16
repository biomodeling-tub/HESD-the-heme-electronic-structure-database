%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=2gsm_III_Dre.chk
%chk=2gsm_III_Hn.chk
# apfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

2gsm_III_H

1 6


