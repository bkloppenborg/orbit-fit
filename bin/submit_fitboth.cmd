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
mpirun -np 16 -machinefile $PBS_NODEFILE ./fitboth ../test_data/epsAur_rv_noeclipse.txt ../test_data/epsAur_ast_refit-pm.txt -gamma_min -10 -gamma_max 10 -K_min 5 -K_max 20 -tau_min 2445400 -tau_max 2455400 -T_min 8000 -T_max 12000 -rv_err 1 -turb -rv_s_min 0 -rv_s_max 20 -alpha_min 1E-7 -alpha_max 1E-4 -inc_min 80 -inc_max 100 -motion -x0_min 75.4920 -x0_max 75.4923 -y0_min 43.8232 -y0_max 43.8235 -mu_x_min -1E-5 -mu_x_max 1E-5 -mu_y_min -1E-5 -mu_y_max 1E-5 -pi_min -1E-5 -pi_max 1E-5