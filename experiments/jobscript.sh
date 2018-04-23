#!/usr/bin/env bash

#BSUB -J "linnea_gen"       # job name
#BSUB -o linnea_gen.%J.cout # job output 
#BSUB -W 0:15               # limits in hours:minutes
#BSUB -M 10000 # memory in MB

module load python/3.6.0

# NEW_INDEX=$((JOB_INDEX-1))
# echo $NEW_INDEX
DATE=$(date +'%Y_%m_%d_%H_%M_%S')

# SEED=${LSB_JOBINDEX}

cd /home/hb765588/linnea
source linnea_venv/bin/activate
cd /home/hb765588/linnea/results/generation
mkdir -p run_${DATE}
cd run_${DATE}
# cp ../config.json .
python3 /home/hb765588/linnea/linnea/experiments/experiments_test.py 1

exit

