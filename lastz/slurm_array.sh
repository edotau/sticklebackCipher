#!/bin/sh
#SBATCH --mem=8G
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --job-name=lastz.sh
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=eric.au@duke.edu
set -e
SCRIPT=/data/lowelab/edotau/rabsTHREEspine/lastz/reciprocalBest/bepaTHREErabsSpine/lastz.sh
#uses a text file as the only input, must know the number of lines ahead of time
command=$(head -${SLURM_ARRAY_TASK_ID} $1 | tail -1)
srun $command

