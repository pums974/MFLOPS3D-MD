#! /bin/bash

points=16
blocs=2
precond=none
solv=gmres
iter=1200
pression=1
ordrev=2
ordrep=2

cd build ; make ; cd ..
cd bin
rm -f matrix_* *.nc
/usr/bin/time mpiexec -n 4 numactl -l ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds 10 -ts 0.01 -ksp_rtol 1.e-13 -pc_type $precond -ksp_type $solv -ntime $iter -nlt 2 -pt $pression -to $ordrev,$ordrep
cd ..
