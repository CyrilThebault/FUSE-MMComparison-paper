#!/bin/bash
#SBATCH --job-name=FUSE_CAMELS
#SBATCH --partition=cpu2022,cpu2021,cpu2019
#SBATCH --error=/home/cyril.thebault/Err/ErrorFile_%A_%a.err
#SBATCH --output=/dev/null
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02-00:00:00
#SBATCH --array=1-4000

# See requirements file to see which modules are needed
module restore FUSE

#! ----------------------------  Param
dataset="CAMELS"
spatialisation="Lumped"
inputdata="daymet"
CritCal="KGECOMP"
TransfoCal="1"  # not used because Critical is set to KGECOMP which is fixed to (KGE(Q) + KGE(1/Q))/2

path_data="/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA"

#! ----------------------------  Basins
# Define the catchment list path
list_basins_path="$path_data/BDD/liste_BV_CAMELS_559.txt"

# Use cut to extract the first column and store it in an array
array_basins=($(cut -d ';' -f 1 "$list_basins_path"))

#! ----------------------------  Models

# Define the model list path
list_models_path="$path_data/FUSE/list_decision_78.txt"

# Use cut to extract the first column and store it in an array, skip the header
array_models=($(cut -d ';' -f 1 "$list_models_path" | tail -n +2))

#! ----------------------------  Basin-model pairs

# Initialize an empty array to store the combinations
combined=()

# Nested loops to create the combinations
for i in "${array_basins[@]}"; do
  for j in "${array_models[@]}"; do
    combined+=("$i,$j")
  done
done

# Check the size of the combined array
size=${#combined[@]}

#! ----------------------------  Parallelization

# retrieve slurm array task id
input=${SLURM_ARRAY_TASK_ID}

#! --- Sub-parallelization (if the number of jobs in the array is greater than the limit allowed, i.e. 4000)

# Limit allowed
joblim=${SLURM_ARRAY_TASK_COUNT}

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

  # Access a specific element (index equal to input-1 because arrays are zero-indexed)
  IFS=',' read -r basin model <<< "${combined[id-1]}"

  #! ----------------------------  Use RAM disk

  mkdir /dev/shm/$basin/output -p
  find "$path_data/FUSE/$dataset/$spatialisation/$inputdata/$CritCal/$TransfoCal/$basin" -mindepth 1 -maxdepth 1 ! -name "output" -exec cp -r {} /dev/shm/$basin \; 

  #! ----------------------------  Select the DecisionFile

  # define the decision file path
  decision_file_path="/dev/shm/$basin/settings/fuse_zDecisions/fuse_zDecisions_${model}.txt"

  # define the folder where the model is supposed to read the settings
  decision_file_output="/dev/shm/$basin/settings/"

  # copy the decision file to the right folder
  cp "$decision_file_path" "$decision_file_output"

  #! ----------------------------  Calibrate FUSE
 
  # define executable path
  exePATH="/home/cyril.thebault/Postdoc_Ucal/02_DATA/installs/fuse/bin/fuse.exe"

  # define file manager results path
  fmPATH="/dev/shm/$basin/${basin}_${model}.txt"

  # run FUSE for calibration
  "$exePATH" "$fmPATH" "$basin" calib_sce

  #! ----------------------------  Simulate FUSE

  # run FUSE for simulation
  "$exePATH" "$fmPATH" "$basin" run_best

  #! ----------------------------  Clean settings folder

  # remove the decision file
  rm "${decision_file_output}fuse_zDecisions_${model}.txt"

  #! ----------------------------  Send output to the project space

  cp -r /dev/shm/$basin/output/ $path_data/FUSE/$dataset/$spatialisation/$inputdata/$CritCal/$TransfoCal/$basin/
  rm -rf /dev/shm/$basin/output

done


echo "######### DONE! ########"