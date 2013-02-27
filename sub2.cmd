#!/bin/bash -l

# job name (default is name of pbs script file)
#PBS -N testconv

# resource limits:  number of nodes and process per nodes 
#PBS -l nodes=1:ppn=8

# resource limits: amount of memory to be used
#PBS -l mem=20gb,vmem=20gb

# resource limits: max. wall clock time during which job can be running
#PBS -l walltime=5:00:00

# Go to submission directory
cd $PBS_O_WORKDIR

module load gnu/4.7.0
module load openmpi/1.4.5-gcc-4.7.0 

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=1

# exports perso ---------------------------------------------------------------

#export points='16'
#export pression='1'
export blocs='2'
#export precond='none'
#export solv='gmres'
export iter='1000'
export ordrev='2'
export ordrep='2'

export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'
export nproc=`echo "$blocs^2" | bc -l `

mkdir $PBS_JOBID
cp bin/testnav $PBS_JOBID
cd $PBS_JOBID

for points in `echo "20 30 40 50"`; do
rm -f matrix_* *.nc
for pression in `echo "1 3"`; do
echo
echo "***************************************************************************************************************"
echo
  /usr/bin/time mpiexec -n $nproc numactl -l \
    ./testnav -dim $points,$points,$points -dom $blocs,$blocs,1 -period 0,0,0 -reynolds 100 -ts 0.1 -ntime $iter \
              -nlt 2 -pt $pression -to $ordrev,$ordrep \
              -u_ksp_rtol 1.e-13 -u_pc_type pbjacobi -u_ksp_type gmres -u_ksp_max_it 10\
              -p_ksp_rtol 1.e-13 -p_pc_type pbjacobi -p_ksp_type gmres -p_ksp_constant_null_space -p_ksp_max_it 60\
      > resultats_$pression_$points.out
echo
echo "***************************************************************************************************************"
echo
done
done

