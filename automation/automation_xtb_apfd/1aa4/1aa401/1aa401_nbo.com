%mem=16GB
%CPU=0-71
%gpucpu=0,1,2,3=0,1,2,3
%oldchk=1aa4_01_opt.chk
%chk=1aa4_01_nbo.chk
# apfd/def2TZVP int=(grid=ultrafine) geom=check guess=read pop=nbo scf=qc

1aa4_01_nbo

0 1


