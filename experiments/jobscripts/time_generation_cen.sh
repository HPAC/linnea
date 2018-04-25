#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-31]" # job name
#BSUB -o "${HOME}/linnea/results/generation/run_%J_cen/out.out" # job output
#BSUB -W 24:00              # limits in hours:minutes
#BSUB -M 20000              # memory in MB

module load python/3.6.0

cd ${HOME}/linnea
source linnea_venv/bin/activate
cd ${HOME}/linnea/results/generation
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}
python3 ${HOME}/linnea/linnea/experiments/time_generation.py ${LSB_JOBINDEX} 1 -cen

exit

