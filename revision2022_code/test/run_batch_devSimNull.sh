#!/bin/bash

#SBATCH --time=6:00:00       # walltime
#SBATCH --partition=cpu_short
#SBATCH --mem-per-cpu=4GB
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --array=1-100      # job array of input size


file=$1

module load matlab
# Slurm array must be max number of pairs per session
# If there are fewer pairs in session, these jobs will
# finish gracefully within MATLAB
srun -p cpu_medium matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "devAssembly_simVrev_null('$file',$SLURM_ARRAY_TASK_ID); exit;"
