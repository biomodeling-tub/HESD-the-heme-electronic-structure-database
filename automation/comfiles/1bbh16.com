--Link7--
%nprocs=12
%mem=16GB
%oldchk=1ebt_III_D.chk
%chk=1ebt_III_H.chk
# PBE1PBE/def2TZVP scf=maxcycle=999 empiricaldispersion=GD3BJ  geom=check pop=nbo

1ebt_III_H

1 6

--Link8--
%nprocs=12
%mem=16GB
%chk=1ebt_III_H.chk
# PBE1PBE/def2TZVP empiricaldispersion=GD3BJ int=(grid=ultrafine) geom=check guess=read scrf=(smd,solvent=Chloroform)

1ebt_III_H_solv

1 6

