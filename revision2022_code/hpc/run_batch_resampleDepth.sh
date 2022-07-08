#!/bin/bash

#SBATCH --time=1:00:00       # walltime
#SBATCH --partition=cpu_short,cpu_medium,cpu_long,fn_short,fn_medium,fn_long,gpu4_short,gpu4_medium,gpu4_long,cpu_dev,gpu4_dev
#SBATCH -J "assemblyTrack"       # job name
#SBATCH --mem-per-cpu=4GB    # 64GB for assembly detection
#SBATCH --nodes=1
#SBATCH --tasks=1


file=$1

module load matlab
# Slurm array must be max number of pairs per session
# If there are fewer pairs in session, these jobs will
# finish gracefully within MATLAB
matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "resample_depth('$file'); exit;"
