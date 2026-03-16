%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1a6k_01_nbo.chk
%chk=1a6k_01_solv.chk
# apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read scrf=(smd,solvent=Chloroform)

1a6k_01_solv

0 1


