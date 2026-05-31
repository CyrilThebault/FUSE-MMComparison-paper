#!/bin/bash
#SBATCH --job-name=Eval
#SBATCH --account=rrg-alpie
#SBATCH --output=/home/thebault/Out/Eval-%A-%a.out
#SBATCH --error=/home/thebault/Err/Eval-%A-%a.error
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --array=1-559


module load r/4.5.0

## Instruction BASH permettant de récupérer la ligne correspondante au numéro d’ARRAY correspondant
input=${SLURM_ARRAY_TASK_ID}
## On lance le programme R avec en  argument le chemin vers le fichier correspondant au bassin SLURM_ARRAY_TASK_ID
Rscript 03_EvalDC.R ${input}
