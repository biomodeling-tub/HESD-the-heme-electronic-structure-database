--Link5--
%nprocs=12
%mem=16GB
%oldchk=1a6m_II_S.chk
%chk=1a6m_III_D.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ scf=maxcycle=999 geom=check pop=nbo

1a6m_III_D

1 2

--Link6--
%nprocs=12
%mem=16GB
%chk=1a6m_III_D.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ geom=check guess=read scrf=(smd,solvent=Chloroform)

1a6m_III_D_solv

1 2
