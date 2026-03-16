%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1aom_12_nbo.chk
%chk=1aom_12_solv.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1aom_12_solv

1 2


