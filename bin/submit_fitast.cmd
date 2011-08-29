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
mpirun -np 16 -machinefile $PBS_NODEFILE ./fitast ./test_ast
