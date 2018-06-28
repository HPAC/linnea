#!/usr/bin/env bash

#BSUB -J "linnea_time[1-31]" # job name
#BSUB -o "linnea/results/generation/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 4000               # memory in MB
#BSUB -P aices2
#BSUB -R model==Haswell_EP
##BSUB -x                   # exclusive access

module load python/3.6.0

cd ${HOME}/linnea
source linnea_venv/bin/activate
cd ${HOME}/linnea/results/generation
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}
python3 ${HOME}/linnea/linnea/experiments/experiments.py time_generation -j=${LSB_JOBINDEX} -r=10 -cen

exit