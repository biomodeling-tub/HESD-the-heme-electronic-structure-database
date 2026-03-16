%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1gjm_12_opt.chk
%chk=1gjm_12_nbo.chk
# opt uapfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

1gjm_12_nbo

1 2


