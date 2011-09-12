#!/bin/bash
### Job name
#PBS -N fitast
### Number of nodes and processors per node.
#PBS -l nodes=2:ppn=8
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
mpirun -np 16 -machinefile $PBS_NODEFILE ./fitast ../test_data/epsAur_ast_refit-pm.txt -tau_min 2428900 -tau_max 2442800 -T_min 9500 -T_max 10500 -alpha_min 1E-7 -alpha_max 1E-4 -e_min 0.252 -e_max 0.283 -inc_min 80 -inc_max 100 -omega_min 33 -omega_max 40 -motion -x0_min 75.4920 -x0_max 75.4923 -y0_min 43.8232 -y0_max 43.8235 -mu_x_min -1E-8 -mu_x_max 1E-8 -mu_y_min -1E-8 -mu_y_max 1E-8 -pi_min -1E-7 -pi_max 1E-7
