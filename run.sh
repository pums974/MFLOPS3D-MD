#! /bin/bash

points=16
blocs=2
precond=none
solv=gmres
iter=120
pression=""
ordrev=2
ordrep=2

export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'

cd build ; make testnav || exit ; cd ..
echo
echo

cd bin
for pression in `echo "1"`; do
  rm -f matrix_* *.nc
  /usr/bin/time mpiexec -n 4 numactl -l \
    ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds 10 -ts 0.01 -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep \
              -ksp_rtol 1.e-13 -pc_type $precond -ksp_type $solv $* \
    | tail -n 14
done
gnuplot cav.gnu

cd ..

echo
echo

exit 0





points=16
blocs=2
precond=none
solv=gmres
iter=120
pression=4
ordrev=2
ordrep=2

export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'


cd build ; make testnav || exit ; cd ..
cd bin
rm -f matrix_* *.nc
/usr/bin/time mpiexec -n 4 numactl -l ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds 10 -ts 0.01 -ntime $iter -nlt 2 -pt $pression -to $ordrev,$ordrep -ksp_rtol 1.e-13 -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_type preonly $*
cd ..


exit 0

-pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_type preonly $*\
-pc_type lu -pc_factor_mat_solver_package mumps        -ksp_type preonly $*\
-pc_type $precond -ksp_type $solv $*\
gnuplot cav2.gnu
