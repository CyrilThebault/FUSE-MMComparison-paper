#!/bin/bash
#SBATCH --job-name=SimArray
#SBATCH --account=rrg-alpie
#SBATCH --output=/home/thebault/Out/SimArray-%A-%a.out
#SBATCH --error=/home/thebault/Err/SimArray-%A-%a.error
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00


module load r/4.5.0


Rscript 02_SimArrayDC.R
