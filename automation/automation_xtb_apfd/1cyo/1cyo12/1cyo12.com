%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cyo_12_opt.chk
%chk=1cyo_12_nbo.chk
# opt apfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

1cyo_12_nbo

1 2


