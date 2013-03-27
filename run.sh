#! /bin/bash

points=100
popul=180
blocs=2
stretch_type=2
stretch_value1=1
stretch_value2=1
stretch_value3=1
export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'

echo
echo
cd build ; make -j8 teststretch || exit ; cd ..
cd build ; make -j8 testsolver_3d || exit ; cd ..
cd bin
#rm fort.*
echo 
echo 

#echo   "+------+-----+----------+----------+----------+--------+----------+----------+----------+--------+-----+"
#printf "|      |     |   Rmin   |   Rmax   |   Rmoy   | Rratio |   Smin   |   Smax   |   Smoy   | Sratio |  S  |\n"
#printf  "    nf    nm         alpha         err_s/err_r   err_r_b/err_r_m      err__sb/err_s_m\n"

for popul in `echo " 70 "`; do # population de l'echantillon
#for blocs in `echo "1 2 4 8"`; do   # points par domaines
#points=`echo 64/$blocs | bc `

for stretch_value1 in `echo "3 4 5 "`; do   # points par domaines
#echo 
#echo $stretch_value1
#echo
#for stretch_value2 in `echo "1.2 1.3"`; do   # points par domaines
#for stretch_value3 in `echo "1.5 2.5 3.5 4.5 5.5 6.5 "`; do   # points par domaines

err=0
solv=1

for points in `echo "16 32 64 128"`; do   # points par domaines

#printf  "  $points"
#printf "| %4d | %3d |" $points $popul
err1=`mpiexec -n $blocs numactl -l ./teststretch -dim $points,$popul,$popul -dom $blocs,1,1 -stretch_type $stretch_type -stretch_value $stretch_value1,$stretch_value2,$stretch_value3 $*    # >/dev/null `

  printf  " $stretch_value1  $points  $err1"

#if [ "$err1" = "              NaN" ]; then err1="0"; solv="0"; fi

#err2=`echo $err1 | sed 's/E/\\*10\\^/' | sed 's/+//' `
#err=`echo "$err + $err2" | bc -l`

solv1=`./testsolver_3d -dim $points,$points,$points -stretch_type $stretch_type -stretch_value  $stretch_value1,$stretch_value2,$stretch_value3 | head -n 1 | awk '{print $2}' `

  printf  "  $solv1\n"

#if [ "$solv1" = "cmplx" ]; then solv="0"; fi
#if [ "$solv1" = "ACML" ]; then solv="0"; fi

#echo "  $err1  $err  $solv1  $solv"
#  cat fort.1? | awk 'BEGIN{min=1e15;max=1e-15;moy=0}{if($4<min){min=$4};if($4>max){max=$4};moy=moy+$4}END{printf(" %4.2e | %4.2e | %4.2e | %6d |",min,max,moy/NR,max/min)}'
#  cat fort.2? | awk 'BEGIN{min=1e15;max=1e-15;moy=0}{if($4<min){min=$4};if($4>max){max=$4};moy=moy+$4}END{printf(" %4.2e | %4.2e | %4.2e | %6d |",min,max,moy/NR,max/min)}'
#  cat fort.3? | awk 'BEGIN{min=1e15;max=1e-15}{if($3<min){min=$3};if($3>max){max=$3}}END{printf(" %3d |\n",max/min)}'
#gnuplot stretch.gnu  -persist

done

#if [ "$solv1" = "cmplx" ]; then solv1="0"; fi
#if [ "$solv1" = "ACML" ]; then solv1="0"; fi

#solv1=`echo $solv1 | sed 's/E/\\*10\\^/' | sed 's/+//' `
#solv=`echo "$solv * $solv1" | bc -l `

##  printf  "  $stretch_value2  $stretch_value3  %4.2e  %4.2e\n" $err  $solv
#if [ ! "$solv" = "0" ]; then
#if [ $(echo "$solv < 0.001" | bc) -eq 1 ]; then
#solv=`echo $solv | sed 's/\./,/'`
#err=`echo $err | sed 's/\./,/'`
#  printf  "  $stretch_value2  $stretch_value3  %4.2e  %4.2e\n" $err  $solv
#fi
#fi
#done
done
done
#done
#echo   "+------+-----+----------+----------+----------+--------+----------+----------+----------+--------+-----+"
cd ..

echo
echo

exit 0

