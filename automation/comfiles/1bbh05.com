--Link3--
%nprocs=12
%mem=16GB
%oldchk=1bbh_II_S.chk
%chk=1bbh_II_Q.chk
# PBE1PBE/def2TZVP empiricaldispersion=GD3BJ scf=maxcycle=999 geom=check pop=nbo

1ebt_II_Q

0 5

--Link4--
%nprocs=12
%mem=16GB
%chk=1ebt_II_Q.chk
# PBE1PBE/def2TZVP empiricaldispersion=GD3BJ geom=check guess=read scrf=(smd,solvent=Chloroform)

1ebt_II_Q_solv

0 5

