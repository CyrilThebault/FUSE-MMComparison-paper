#!/bin/bash
#SBATCH --job-name=EvalFUSE
#SBATCH --partition=cpu2022,cpu2021,cpu2019
#SBATCH --output=/home/cyril.thebault/Out/Bassin-%A-%a.out
#SBATCH --error=/home/cyril.thebault/Err/Bassin-%A-%a.error
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --array=1-559


module load R/4.3.1
module load osgeo/gdal
module load lib/netcdf


input=${SLURM_ARRAY_TASK_ID}

Rscript 04_EvalFUSE.R ${input}
