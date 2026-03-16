%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1aa4_05_opt.chk
%chk=1aa4_05_nbo.chk
# opt uapfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

1aa4_05_nbo

0 5


