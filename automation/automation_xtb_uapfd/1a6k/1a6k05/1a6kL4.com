%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1a6k_05_nbo.chk
%chk=1a6k_05_solv.chk
# uapfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1a6k_05_solv

0 5


