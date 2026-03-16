--Link7--
%nprocs=12
%mem=16GB
%oldchk=1ac8_III_D.chk
%chk=1ac8_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ scf=maxcycle=999 geom=check pop=nbo

1ac8_III_H

1 6

--Link8--
%nprocs=12
%mem=16GB
%chk=1ac8_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ geom=check guess=read scrf=(smd,solvent=Chloroform)

1ac8_III_H_solv

1 6
