#!/bin/bash -l

# job name (default is name of pbs script file)
#PBS -N test

# resource limits:  number of nodes and process per nodes 
#PBS -l nodes=1:ppn=4

# resource limits: amount of memory to be used
#PBS -l mem=5gb,vmem=5gb

# resource limits: max. wall clock time during which job can be running
#PBS -l walltime=1:00:00

# Go to submission directory
cd $PBS_O_WORKDIR

module load gnu/4.7.0
module load openmpi/1.4.5-gcc-4.7.0 

export OMP_SCHEDULE=dynamic
export OMP_NUM_THREADS=1
export TIME='Elapsed : %e ,User: %U ,CPU : %P ,Max Mem : %R , command : %C'

#export DIM='24,24,24'
#export DIM='32,32,32'
export DIM='16,16,16'
export DOM='2,2,1'
export REY='10'
export PROC='4'
export TS='0.01'
export NTIME='80'

mkdir $PBS_JOBID
cp testnav $PBS_JOBID
cd $PBS_JOBID


export DIM='10,10,10'
echo "*********************************************************************************************************************************"
rm -f matrix_* *.nc ; /usr/bin/time mpiexec -n $PROC numactl -l ./testnav -dim $DIM -dom $DOM -period 0,0,0 -reynolds $REY -ts $TS -ksp_rtol 1.e-8 -pc_type none -ksp_type gmres -ntime $NTIME -nlt 2 -pt 1 -to 2,2 > out.4
echo "*********************************************************************************************************************************"


#echo "*********************************************************************************************************************************"
#rm -f matrix_* *.nc ; /usr/bin/time mpiexec -n $PROC numactl -l ./testnav -dim $DIM -dom $DOM -period 0,0,0 -reynolds $REY -ts $TS -ksp_rtol 1.e-8 -pc_type none -ksp_type gmres -ntime $NTIME -nlt 2 -pt 1 -to 2,2 > out.32

#rm -f matrix_* *.nc ; /usr/bin/time mpiexec -n $PROC numactl -l ./testnav -dim $DIM -dom $DOM -period 0,0,0 -reynolds $REY -ts $TS -ksp_rtol 1.e-8 -pc_type none -ksp_type bcgs -ntime $NTIME -nlt 2 -pt 1 -to 2,2 > out.32

#rm -f matrix_* *.nc ; /usr/bin/time mpiexec -n $PROC numactl -l ./testnav -dim $DIM -dom $DOM -period 0,0,0 -reynolds $REY -ts $TS -ksp_rtol 1.e-8 -ksp_type lgmres -ntime $NTIME -nlt 2 -pt 1 -to 2,2 -pc_type bjacobi -ksp_sub_pc_type ilu -sub_ksp_type gmres -sub_ksp_max_it 6 -sub_pc_type bjacobi -sub_sub_pc_type ilu > out.32

#rm -f matrix_* *.nc ; /usr/bin/time mpiexec -n $PROC numactl -l ./testnav -dim $DIM -dom $DOM -period 0,0,0 -reynolds $REY -ts $TS -ksp_rtol 1.e-8 -ksp_type bcgs -ntime $NTIME -nlt 2 -pt 1 -to 2,2 -pc_type ksp -ksp_ksp_type chebyshev -ksp_ksp_chebychev_estimate_eigenvalues 0.1,1.1 -ksp_ksp_max_it 2 -ksp_ksp_norm_type none -ksp_pc_type bjacobi -ksp_sub_pc_type ilu > out.32
#echo "*********************************************************************************************************************************"

