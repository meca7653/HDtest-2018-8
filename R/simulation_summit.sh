#!/bin/bash
#SBATCH --export=NONE
#SBATCH -p shas                   # Send this job to the serial partition
#SBATCH -N 1                        # number of nodes
#SBATCH --ntasks-per-node=20
#SBATCH -t 0-6:00                  # wall time (D-HH:MM)
##SBATCH -A drzuckerman             # Account to pull cpu hours from (commented out)
#SBATCH -o simulation_summit.%j.out             # STDOUT (%j = JobId)
#SBATCH -e simulation_summit.%j.err             # STDERR (%j = JobId)
#SBATCH --array=1-6
##SBATCH --mail-type=END,FAIL        # notifications for job done & fail
R_PROFILE=/home/meca7653@colostate.edu/R/x86_64-pc-linux-gnu-library/3.3/snow/RMPISNOWprofile;
export R_PROFILE
module load R/3.4.3
ml gcc
ml intel
ml impi
mpirun Rscript --no-save $HOME/lustrefs/Chang/simulation_summit.R simulation_summit.${SLURM_ARRAY_TASK_ID}.out
