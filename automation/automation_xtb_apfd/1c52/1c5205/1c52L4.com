%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1c52_05_nbo.chk
%chk=1c52_05_solv.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1c52_05_solv

0 5


