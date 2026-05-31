#!/bin/bash
#SBATCH --job-name=PREPA_CALIB_FUSE
#SBATCH --partition=cpu2022,cpu2021,cpu2019
#SBATCH --output=/home/cyril.thebault/Out/Bassin-%A-%a.out
#SBATCH --error=/home/cyril.thebault/Err/Bassin-%A-%a.error
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00-02:00:00
#SBATCH --array=1-559

module restore FUSE

input=${SLURM_ARRAY_TASK_ID}

Rscript 01_FromCamelsToFUSE_Cluster.R ${input}