%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1ccl_II_S.chk
%chk=1ccl_II_Qn.chk
# apfd/def2TZVP scf=maxcycle=999 geom=check pop=nbo

1ccl_II_Q

0 5


