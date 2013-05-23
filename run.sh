#! /bin/zsh
trap 'exit 1' 2

points=20 # 30
blocs=1   # 4


zmodload zsh/mathfunc 
pression=1
ordrev=2
ordrep=2
ts=4e-3
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

if [   "$blocs" -eq "0" ] ; then 
  export OMP_NUM_THREADS=2

  pointsx=`echo "$points*9" | bc `
  pointsy=`echo "$points"   | bc `
  pointsz=`echo "$points*3" | bc `
   blocsx=1
   blocsy=1
   blocsz=1

  nproc=1

else
  export OMP_NUM_THREADS=1

  pointsx=`echo "$points"  | bc `
  pointsy=`echo "$points"  | bc `
  pointsz=`echo "$points"  | bc `
   blocsx=`echo "$blocs*5" | bc `
   blocsy=`echo "$blocs"   | bc `
   blocsz=`echo "$blocs*2" | bc `

  nproc=`echo "$blocsx*$blocsy*$blocsz" | bc `
fi

  if [ ! "`ls | grep -c 'matrix_'`" -eq "0" ]; then rm -f matrix_*;fi
  if [ ! "`ls | grep -c 'fort.'`" -eq "0" ];   then rm -f fort.*  ;fi
####  if [ ! "`ls | grep -c '\.nc'`" -eq "0" ];    then rm -f *.nc    ;fi
  
 /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $pointsx,$pointsy,$pointsz -dom $blocsx,$blocsy,$blocsz -period 1,0,1 -reynolds $re -ts $ts -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub \
              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2 \
              $* 

gnuplot canal.gnu

cd ..

echo
echo

exit 0


              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2 \

              -u_ksp_rtol 1.e-8 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 60 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-8 -p_pc_type jacobi -p_ksp_type gmres -p_ksp_max_it 60 -p_ksp_initial_guess_nonzero -psm 1 \
