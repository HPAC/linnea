#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-31]" # job name
#BSUB -o "linnea/results/generation/run_%J_c/out.out" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB

module load python/3.6.0

cd ${HOME}/linnea
source linnea_venv/bin/activate
cd ${HOME}/linnea/results/generation
mkdir -p run_${LSB_JOBID}_c
cd run_${LSB_JOBID}_c
python3 ${HOME}/linnea/linnea/experiments/time_generation.py ${LSB_JOBINDEX} 1 -c

exit