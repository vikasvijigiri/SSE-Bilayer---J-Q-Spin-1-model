#!/bin/bash

########################
### -- queue information
### medium - 3days
### long - 7days
#######################

#--------------------------------------------------------------------#
#SBATCH --partition bigq         ### parition
#SBATCH --mem 512M                   ### memory pool for all cores
#SBATCH --ntasks=1								 ### number of nodes
#SBATCH --cpus-per-task 1          ### number of cores
##SBATCH --nodelist=bose9
#SBATCH --time 3-00:00             ### time (D-HH:MM)
#SBATCH --output slurm.out         ##slurm.%N.%j.out   ##  STDOUT
#SBATCH --error  slurm.err         #%N.%j.err   ### STDERR
#SBATCH --job-name bilayer         ### name of the job
#--------------------------------------------------------------------#

module load gcc/gcc-10.1.0 

#make
#rm *.o

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

time ./main 
