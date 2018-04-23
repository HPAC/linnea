#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-31]" # job name
#BSUB -o linnea_gen.%J.cout # job output 
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 10000 # memory in MB

module load python/3.6.0

# NEW_INDEX=$((JOB_INDEX-1))
# echo $NEW_INDEX
# DATE=$(date +'%Y_%m_%d_%H_%M_%S')

# SEED=${LSB_JOBINDEX}

cd /home/hb765588/linnea
source linnea_venv/bin/activate
cd /home/hb765588/linnea/results/generation
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}
# cp ../config.json .
python3 /home/hb765588/linnea/linnea/experiments/experiments_test.py ${LSB_JOBINDEX} 1

exit

