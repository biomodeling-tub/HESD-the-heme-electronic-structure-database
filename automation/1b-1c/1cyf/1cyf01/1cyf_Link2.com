%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cyf_II_Sn.chk
%chk=1cyf_II_S.chk
# apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read scrf=(smd,solvent=Chloroform) 

1cyf_II_S_solv 

0 1 


