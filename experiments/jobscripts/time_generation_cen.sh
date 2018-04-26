#!/usr/bin/env bash

#BSUB -J "linnea_time[1-31]" # job name
#BSUB -o "linnea/results/generation/run_%J/cout.txt" # job output
#BSUB -W 48:00              # limits in hours:minutes
#BSUB -M 60000              # memory in MB
#BSUB -P aices
#BSUB -R model==SandyBridge_EP
#BSUB -N

module load python/3.6.0

cd ${HOME}/linnea
source linnea_venv/bin/activate
cd ${HOME}/linnea/results/generation
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}
python3 ${HOME}/linnea/linnea/experiments/experiments.py time_generation -j=${LSB_JOBINDEX} -r=1 -cen

exit

