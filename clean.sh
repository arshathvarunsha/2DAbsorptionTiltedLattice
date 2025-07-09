#!/bin/bash
##NECESSARY JOB SPECIFICATIONS
#SBATCH --partition=shared
#SBATCH --job-name=Run_Plot2D
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --output=output/clean.%j.out
#SBATCH --error=output/clean.%j.err


# Change to the directory where the script is located
cd $SLURM_SUBMIT_DIR

# Define omega value (default to 1 if not passed as argument)
OMEGA="3"

# Path to the Python script
SCRIPT="plotting_2d.py"

while true; do
    echo "[$(date)] Running the Python script with omega = $OMEGA"
    python3 "$SCRIPT" "$OMEGA"

    echo "Sleeping for 1 minute..."
    sleep 60
done