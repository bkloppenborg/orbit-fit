#!/bin/bash
### Job name
#PBS -N fitrv
### Number of nodes and processors per node.
#PBS -l nodes=1:ppn=8
#PBS -q normal
#PBS -j oe

### Uncomment to send email when the job is completed:
##PBS -m abe
##PBS -M youremail@email.com

cd $PBS_O_WORKDIR
echo Host: $HOSTNAME
echo Date: $(date)
echo Dir: $PWD
echo This jobs runs on the following processors:
cat $PBS_NODEFILE
# openmpi test

echo fitrv_turb

mpirun -np 8 -machinefile $PBS_NODEFILE ./fitrv ../test_data/epsAur_rv_noeclipse.txt -gamma_min -10 -gamma_max 10 -K_min 5 -K_max 20 -tau_min 2445400 -tau_max 2455400 -T_min 8000 -T_max 12000 -rv_err 1 -turb -rv_s_min 0 -rv_s_max 20

