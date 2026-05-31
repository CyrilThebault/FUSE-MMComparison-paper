#!/bin/bash
#SBATCH --job-name=MergeEval
#SBATCH --account=rrg-alpie
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cyril.thebault@ucalgary.ca
#SBATCH --output=/home/thebault/Out/MergeEval-%A-%a.out
#SBATCH --error=/home/thebault/Err/MergeEval-%A-%a.error
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00


module load r/4.5.0


Rscript 04_MergeEvalDC.R
