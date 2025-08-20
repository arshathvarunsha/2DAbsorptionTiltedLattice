#!/bin/bash

#SBATCH --export=NONE               
#SBATCH --get-user-env=L

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=demo-tamulauncher
#SBATCH --output=demo-tamulauncher.%j
#SBATCH --time=01:00:00            
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1000M

tamulauncher commands.in