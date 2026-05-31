#!/bin/bash
#SBATCH --job-name=DC
#SBATCH --account=rrg-alpie
#SBATCH --output=/home/thebault/Out/DC-%A-%a.out
#SBATCH --error=/home/thebault/Err/DC-%A-%a.error
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03-00:00:00
#SBATCH --array=1-9000


module load r/4.5.0 

#! ----------------------------  Parallelization

# retrieve slurm array task id
input=${SLURM_ARRAY_TASK_ID}

#! --- Sub-parallelization (if the number of jobs in the array is greater than the limit allowed, i.e. 4000)

# Limit allowed
joblim=${SLURM_ARRAY_TASK_COUNT}

size=559000

# Number of subjoins per jobs
subjobs=$(echo "scale=0; ($size + $joblim - 1) / $joblim" | bc)

# Boundaries of the loop for subjobs
start=$(( (input - 1) * subjobs + 1 ))
end=$(( input * subjobs ))

# If the lower boundary is greater than the total number of jobs, exit the script
if [ $start -gt $size ]; then
  echo "Start value ($start) is greater than size ($size). Exiting."
  exit 1
fi

# If the upper boundary is greater than the total number of jobs, fix it to the total number of job value
if [ $end -ge $size ]; then
  end=$size
fi

# Loop for the subjobs
for ((id=start; id<=end; id++)); do

  echo "Running ID $id..."
  
  #timeout 60m Rscript 01_DC.R ${id}
  Rscript 01_DC.R ${id}

  # Capture the exit code
  exit_code=$?
  
  if [ $exit_code -eq 124 ]; then
    echo "ID $id timed out after 60 minutes."
  elif [ $exit_code -ne 0 ]; then
    echo "ID $id failed with exit code $exit_code."
  else
    echo "ID $id completed successfully."
  fi

done


echo "######### DONE! ########"