--Link7--
%nprocs=12
%mem=16GB
%oldchk=102m_III_D.chk
%chk=102m_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ scf=maxcycle=999 geom=check pop=nbo

102m_III_H

1 6

--Link8--
%nprocs=12
%mem=16GB
%chk=102m_III_H.chk
# apfd/def2TZVP empiricaldispersion=GD3BJ geom=check guess=read scrf=(smd,solvent=Chloroform)

102m_III_H_solv

1 6
