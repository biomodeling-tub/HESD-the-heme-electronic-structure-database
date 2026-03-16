%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1dp6_16_opt.chk
%chk=1dp6_16_nbo.chk
# opt uapfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

1dp6_16_nbo

1 6


