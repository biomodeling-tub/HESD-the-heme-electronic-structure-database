%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1cpt_05_opt.chk
%chk=1cpt_05_nbo.chk
# opt uapfd/def2TZVP guess=read scf=maxcycle=999 geom=check pop=nbo

1cpt_05_nbo

0 5


