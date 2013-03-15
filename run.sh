#! /bin/bash

points=16
blocs=2
#precond=none
#solv=lgmres pgmres
iter=80
#pression=""
ordrev=2
ordrep=2

export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'
nproc=`echo "$blocs^2" | bc -l `

echo
echo
cd build ; make testnav || exit ; cd ..
cd bin

for points in `echo "16"`; do
#for blocs in `echo "2 3 4 5"`; do
  rm -f matrix_* *.nc
echo
echo
echo $points
for pression in `echo "3"`; do
for nsub in `echo "1"`; do
echo
echo
  /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds 10 -ts 0.01 -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub \
              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2\
              $*   #  2>/dev/null   |  tail -n 11 #11
#              $*    > test_$pression_$blocs.out
#gnuplot cav2.gnu
done
done
done




cd ..

echo
echo

exit 0


              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2\
              -u_ksp_rtol 1.e-8 -u_pc_type bjacobi -u_ksp_type gmres -u_ksp_max_it 10\
              -p_ksp_rtol 1.e-8 -p_pc_type bjacobi -p_ksp_type gmres -p_ksp_max_it 60 -psm 1\
              -u_ksp_rtol 1.e-8 -u_pc_type pbjacobi -u_ksp_type gmres -u_ksp_max_it 20 \
              -p_ksp_rtol 1.e-8 -p_pc_type bjacobi -p_ksp_sub_pc_type ilu -p_sub_ksp_type gmres -p_sub_ksp_max_it 6 -p_sub_pc_type bjacobi -p_sub_sub_pc_type ilu -p_ksp_max_it 20 -psm 1 \





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

-pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_type preonly \
-pc_type lu -pc_factor_mat_solver_package mumps        -ksp_type preonly \
-pc_type $precond -ksp_type $solv \
gnuplot cav2.gnu
