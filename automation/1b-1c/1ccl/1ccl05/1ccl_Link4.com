%mem=8GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1ccl_II_Qn.chk
%chk=1ccl_II_Q.chk
# apfd/def2TZVP geom=check guess=read scrf=(smd,solvent=Chloroform)

1ccl_II_Q_solv

0 5


