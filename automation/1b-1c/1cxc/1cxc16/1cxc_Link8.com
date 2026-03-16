%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cxc_III_Hn.chk
%chk=1cxc_III_H.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1cxc_III_H_solv

1 6


