#!/bin/bash
### Job name
#PBS -N fitast
### Number of nodes and processors per node.
#PBS -l nodes=1:ppn=8
#PBS -q bigmem
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

echo fitast_ppmx-refit
AST_FILE=epsAur_ast_tycho2-normal.txt

mpirun -np 8 -machinefile $PBS_NODEFILE ./fitast ../test_data/${AST_FILE} -tau_min 2413932 -tau_max 2454878 -T_min 9896 -T_max 9898 -alpha_min 1E-3 -alpha_max 7E-2 -inc_min 89 -inc_max 89.5 -Omega_min 116 -Omega_max 118 -motion -x0_min -25 -x0_max -23 -y0_min 74 -y0_max 76 -mu_x_min -1E-4 -mu_x_max 1E-4 -mu_y_min -1E-4 -mu_y_max 1E-4 -pi_min -1E-1 -pi_max 1E-1 -omega_min 45 -omega_max 47 -e_min 0.240 -e_max 0.241 -noerr -ast_err 0.04


#mpirun -np 8 -machinefile $PBS_NODEFILE ./fitast ../test_data/${AST_FILE} -tau_min 2413932 -tau_max 2454878 -T_min 8000 -T_max 12000 -alpha_min 1E-3 -alpha_max 7E-2 -inc_min 85 -inc_max 95 -Omega_min 116 -Omega_max 118 -motion -x0_min -25 -x0_max -23 -y0_min 74 -y0_max 76 -mu_x_min -1E-4 -mu_x_max 1E-4 -mu_y_min -1E-4 -mu_y_max 1E-4 -pi_min -1E-1 -pi_max 1E-1 -omega_min 30 -omega_max 50 -e_min 0.260 -e_max 0.285 -noerr -ast_err 0.04
