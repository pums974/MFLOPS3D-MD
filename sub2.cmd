#!/bin/bash -l

# job name (default is name of pbs script file)
#PBS -N testcanal

# resource limits:  number of nodes and process per nodes 
#PBS -l nodes=1:ppn=48

# resource limits: amount of memory to be used
#PBS -l mem=100gb,vmem=100gb

# resource limits: max. wall clock time during which job can be running
#PBS -l walltime=10:00:00

# Go to submission directory
cd $PBS_O_WORKDIR

module load gnu/4.7.0
module load openmpi/1.4.5-gcc-4.7.0 

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=1

# exports perso ---------------------------------------------------------------

export points=20 # 30
export blocs=2
export pression=3
export ordrev=2
export ordrep=2
export ts=1e-3
export nsub=1
export iter=1000
export re=3250 # re_tau = 180


export     TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'
export  pointsx=`echo "$points"  | bc `
export  pointsy=`echo "$points"  | bc `
export  pointsz=`echo "$points"  | bc `
export   blocsx=`echo "$blocs*3" | bc `
export   blocsy=`echo "$blocs"   | bc `
export   blocsz=`echo "$blocs*2" | bc `
export    nproc=`echo "$blocsx*$blocsy*$blocsz" | bc `

mkdir $PBS_JOBID
cp ../testnav $PBS_JOBID
cd $PBS_JOBID
ln -s ../../gridi.nc gridi.nc  
ln -s ../../init.nc init.nc  

echo "***************************************************************************************************************"

 /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $pointsx,$pointsy,$pointsz -dom $blocsx,$blocsy,$blocsz -period 1,0,1 -reynolds $re -ts $ts -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep -nsub $nsub \
              -u_ksp_rtol 1.e-8 -u_pc_type jacobi -u_ksp_type gmres -u_ksp_max_it 60 -u_ksp_initial_guess_nonzero \
              -p_ksp_rtol 1.e-8 -p_pc_type jacobi -p_ksp_type gmres -p_ksp_max_it 100 -p_ksp_initial_guess_nonzero -psm 1 \
              > resultat.out 

echo "***************************************************************************************************************"
echo
