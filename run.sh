#! /bin/zsh
trap 'exit 1' 2

points=20 # 30
blocs=1   # 4


zmodload zsh/mathfunc 
pression=3
ordrev=2
ordrep=2
ts=2e-3
nsub=1
iter=10000
re=3250 # re_tau = 180


export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'

echo
echo
cd build ; make -j8 testnav || exit ; cd ..
cd bin
echo 
echo

  export OMP_NUM_THREADS=1

  pointsx=20
  pointsy=20
  pointsz=20
   blocsx=1 #3
   blocsy=3 #3
   blocsz=1 #2

  nproc=`echo "$blocsx*$blocsy*$blocsz" | bc `

  if [ ! "`ls | grep -c 'matrix_'`" -eq "0" ]; then rm -f matrix_*;fi
  if [ ! "`ls | grep -c 'fort.'`" -eq "0" ];   then rm -f fort.*  ;fi
  if [ ! "`ls | grep -c '.nc'`" -eq "0" ];     then rm -f *.nc    ;fi
  
 /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $pointsx,$pointsy,$pointsz -dom $blocsx,$blocsy,$blocsz -period 1,0,1 -reynolds $re -ts $ts -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub -psm 1\
              -u_ksp_rtol 1.e-8 -u_pc_type bjacobi -u_ksp_type gmres -u_ksp_max_it 60 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-5 -p_ksp_max_it 50  -p_ksp_initial_guess_nonzero -p_ksp_type gmres  -p_pc_type ksp     -p_ksp_ksp_type chebyshev -p_ksp_ksp_chebychev_estimate_eigenvalues 0.1,1.1 -p_ksp_ksp_max_it 2 -p_ksp_ksp_norm_type none -p_ksp_pc_type bjacobi -p_ksp_sub_pc_type ilu \
              $* 

#gnuplot canal.gnu

cd ..

echo
echo

exit 0


              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2 \

              -u_ksp_rtol 1.e-8 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 60 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-8 -p_pc_type jacobi -p_ksp_type gmres -p_ksp_max_it 60 -p_ksp_initial_guess_nonzero -psm 1 \

              -u_ksp_rtol 1.e-8 -u_ksp_max_it 50  -u_ksp_initial_guess_nonzero -u_ksp_type lgmres -u_pc_type bjacobi -u_ksp_sub_pc_type ilu -u_sub_ksp_type gmres -u_sub_ksp_max_it 6 -u_sub_pc_type bjacobi -u_sub_sub_pc_type ilu \
              -p_ksp_rtol 1.e-5 -p_ksp_max_it 50  -p_ksp_initial_guess_nonzero -p_ksp_type gmres  -p_pc_type ksp     -p_ksp_ksp_type chebyshev -p_ksp_ksp_chebychev_estimate_eigenvalues 0.1,1.1 -p_ksp_ksp_max_it 2 -p_ksp_ksp_norm_type none -p_ksp_pc_type bjacobi -p_ksp_sub_pc_type ilu \

              -p_ksp_rtol 1.e-8 -p_ksp_max_it 50 -p_pc_type ml -p_ksp_type gmres -p_ksp_initial_guess_nonzero -p_pc_ml_CoarsenScheme MIS -p_pc_ml_Threshold 0.75 -p_mg_levels_0_pc_type bjacobi

