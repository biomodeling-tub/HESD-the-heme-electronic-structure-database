%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cpf_II_S.chk
%chk=1cpf_III_Dn.chk
# apfd/def2TZVP scf=maxcycle=999 geom=check pop=nbo

1cpf_III_D

1 2


