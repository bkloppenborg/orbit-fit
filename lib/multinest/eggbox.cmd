#!/bin/bash
### Job name
#PBS -N eggbox
### Number of nodes and processors per node.
#PBS -l nodes=1:ppn=8
#PBS -q small
#PBS -j oe

### Uncomment to send email when the job is completed:
##PBS -m abe
##PBS -M youremail@email.com

env

cd $PBS_O_WORKDIR
echo Host: $HOSTNAME
echo Date: $(date)
echo Dir: $PWD
echo This jobs runs on the following processors:
cat $PBS_NODEFILE
# openmpi test
mpirun -np 8 -machinefile $PBS_NODEFILE /home/bkloppen/lib/MultiNest_v2.12/eggboxC
