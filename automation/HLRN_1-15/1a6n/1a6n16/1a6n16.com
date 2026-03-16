--Link7--
%nprocs=12
%mem=16GB
%oldchk=1a6n_III_D.chk
%chk=1a6n_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ scf=maxcycle=999 geom=check pop=nbo

1a6n_III_H

1 6

--Link8--
%nprocs=12
%mem=16GB
%chk=1a6n_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ geom=check guess=read scrf=(smd,solvent=Chloroform)

1a6n_III_H_solv

1 6
