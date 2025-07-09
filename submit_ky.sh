#!/bin/bash
##NECESSARY JOB SPECIFICATIONS
#SBATCH --partition=shared
#SBATCH --job-name=Run          #Set the job name to "JobExample1"
#SBATCH --time=1:00:00          #Set the wall clock limit to 12 hours
### Array
#SBATCH --array=0-40
#SBATCH --ntasks-per-node=1

#SBATCH --mem=1000M              
#SBATCH --output=output/out.%j  

sbatch run_2d_abs.py ${SLURM_ARRAY_TASK_ID}