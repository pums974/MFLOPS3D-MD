#! /bin/zsh
trap 'exit 1' 2
points=20
blocs=2
#precond=none
#solv=lgmres pgmres
zmodload zsh/mathfunc 
iter=80
#pression=""
ordrev=2
ordrep=2
ts=1e0
re=1e0
nsub=1
ligne=8

export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'
if [   "$blocs" -eq "1" ] ; then export OMP_NUM_THREADS=6;fi
if [ ! "$blocs" -eq "1" ] ; then export OMP_NUM_THREADS=1;fi

echo
echo
cd build ; make -j8 testnav || exit ; cd ..
cd bin
 

for points   in `echo "16"`  ; do
#for iter     in `echo "800"` ; do
#for ts       in `echo "1e0"` ; do
for re       in `echo "1e0"` ; do
for ordrep   in `echo "2"`   ; do
for ordrev   in `echo "2 "`   ; do
#if [ $ordrep -le $ordrev ];then
for pression in `echo "1"`  ; do
echo
echo
oldv=0
oldp=0
olds=0
oldd=0
echo "iter="$iter",re="$re", to="$ordrev","$ordrep", pt="$pression
#for ts in `echo "1 2 4 8 16 32 64 128 256 512 1024 2048 5096"`; do
#ts=`echo "$((1.0/($ts)))"`
#iter=`echo "$(($iters/(1*$ts)))"`
#for points   in `echo "16 24 36 54 81"`  ; do
for iter in `echo "20 40 80 160 320 640 1280 2560 5120 10240"`; do
ts=`echo "$((5.0/$iter))"`
#for blocs in `echo "1"`; do
  nproc=`echo "$blocs^3" | bc -l `
  if [ ! "`ls | grep -c 'matrix_'`" -eq "0" ]; then rm -f matrix_*;fi
  if [ ! "`ls | grep -c '\.nc'`" -eq "0" ];    then rm -f *.nc    ;fi
  if [ ! "`ls | grep -c 'fort.'`" -eq "0" ];   then rm -f fort.*  ;fi
  

#for pression in `echo " 1 3"`; do
#for nsub in `echo "1"`; do
lignes=$ligne
if [ ! "$blocs" -eq "1" ] ; then lignes=`echo "$ligne+4" | bc -l `;fi
if [ ! "$nsub"  -eq "1" ] ; then lignes=`echo "$ligne+2" | bc -l `;fi
#echo
#echo
#echo $pression $nsub
s=$iter
resultat=` /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $points,$points,$points -dom $blocs,$blocs,$blocs -period 0,0,0 -reynolds $re -ts $ts -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub \
              -u_ksp_rtol 1.e-12 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 60 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-12 -p_pc_type jacobi -p_ksp_sub_pc_type ilu -p_sub_ksp_type gmres -p_sub_ksp_max_it 6 -p_sub_pc_type bjacobi -p_sub_sub_pc_type ilu -p_ksp_max_it 60 -p_ksp_initial_guess_nonzero -psm 1 \
              $*    2>/dev/null   | tee run.out |  tail -n $lignes `

#resultat=`numactl -l ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds $re -ts $ts -ntime $iter \
#              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub \
#              $*    2>/dev/null   | tee run.out |  tail -n $lignes`
#echo $resultat
d=`echo $resultat | grep "error Div V" | awk '{print $5}'`
v=`echo $resultat | grep "error tot V" | awk '{print $5}'`
p=`echo $resultat | grep "error tot P" | awk '{print $5}'`

printf " %5i %E %E %E %E " $iter $ts $v $p $d
if [ ! "$olds" = "0" ]; then
ov=`echo $(((log($v)-log($oldv))/(log($olds)-log($s))))`
op=`echo $(((log($p)-log($oldp))/(log($olds)-log($s))))`
od=`echo $(((log($d)-log($oldd))/(log($olds)-log($s))))`
printf "%4.2f %4.2f %4.2f" $ov $op $od
fi
#echo
olds=$s
oldv=$v
oldp=$p
oldd=$d
#              $*    > test_$pression_$blocs.out
#gnuplot cav2.gnu
./visu.sh
done
done
#fi
done
done
done
done
done



cd ..

echo
echo

exit 0


              -u_pc_type lu -u_pc_factor_mat_solver_package mumps        -u_ksp_type preonly \
              -p_pc_type lu -p_pc_factor_mat_solver_package mumps        -p_ksp_type preonly -psm 2 \

              -u_ksp_rtol 1.e-12 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 10 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-12 -p_pc_type jacobi -p_ksp_type gmres -p_ksp_max_it 60 -p_ksp_initial_guess_nonzero -psm 1 \

              -u_ksp_rtol 1.e-12 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 10 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-12 -p_pc_type jacobi -p_ksp_sub_pc_type ilu -p_sub_ksp_type gmres -p_sub_ksp_max_it 6 -p_sub_pc_type bjacobi -p_sub_sub_pc_type ilu -p_ksp_max_it 10 -p_ksp_initial_guess_nonzero -psm 1 \





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
