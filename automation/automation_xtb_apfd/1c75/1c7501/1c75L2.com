%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1c75_01_nbo.chk
%chk=1c75_01_solv.chk
# apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read scrf=(smd,solvent=Chloroform)

1c75_01_solv

0 1


