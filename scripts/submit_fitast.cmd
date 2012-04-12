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

FITBOTH=../bin/fitboth

OUTPUT=fitast_ppmx-refit
AST_FILE=epsAur_ast_${OUTPUT}.txt
echo Running file: ${AST_FILE}

tau_min=2413932
tau_max=2454878

mpirun -np 8 -machinefile $PBS_NODEFILE ${FITBOTH} ../test_data/epsAur_rv_noeclipse.txt ../test_data/${AST_FILE} -o ${OUTPUT} -gamma_min -10 -gamma_max 10 -K_min 5 -K_max 20 -tau_min ${tau_min} -tau_max ${tau_max} -T_min 9000 -T_max 11000 -rv_err 1 -turb -rv_s_min 0 -rv_s_max 20 -alpha_min 1E-3 -alpha_max 7E-2 -inc_min 89 -inc_max 89.5 -Omega_min 116 -Omega_max 118 -motion -x0_min -25 -x0_max -23 -y0_min 74 -y0_max 76 -mu_x_min -1E-4 -mu_x_max 1E-4 -mu_y_min -1E-4 -mu_y_max 1E-4 -pi_min -1E-1 -pi_max 1E-1 -noerr -ast_err 0.04

