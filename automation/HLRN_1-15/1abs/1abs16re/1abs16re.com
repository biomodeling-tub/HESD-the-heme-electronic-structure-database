--Link7--
%nprocs=12
%mem=16GB
%oldchk=1abs_III_D.chk
%chk=1abs_III_Dre.chk
# apfd/gen scf=maxcycle=999 geom=check

1abs_III_H

1 6

C N O H 0
6-31G*
****
Fe 0
TZVP
****
