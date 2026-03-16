%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1ch1_II_S.chk
%chk=1ch1_III_Dn.chk
# apfd/def2TZVP scf=maxcycle=999 geom=check pop=nbo

1ch1_III_D

1 2


