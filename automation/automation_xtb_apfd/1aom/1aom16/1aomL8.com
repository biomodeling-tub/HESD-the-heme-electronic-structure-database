%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1aom_16_nbo.chk
%chk=1aom_16_solv.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1aom_16_solv

1 6


