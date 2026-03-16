%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=2gsm_III_Dn.chk
%chk=2gsm_III_D.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

2gsm_III_D_solv

1 2


