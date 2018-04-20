#!/usr/bin/env bash

module load python/3.6.0

# NEW_INDEX=$((JOB_INDEX-1))
# echo $NEW_INDEX
DATE=$(date +'%d_%m_%Y_%H_%M_%S')

# SEED=${LSB_JOBINDEX}

cd /home/hb765588/linnea
source linnea_venv/bin/activate
cd /home/hb765588/linnea/results/generation
mkdir -p run_${DATE}
cd run_${DATE}
# cp ../config.json .
python3 /home/hb765588/linnea/linnea/experiments/experiments_test.py

exit

