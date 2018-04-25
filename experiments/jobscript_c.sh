#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-31]" # job name
#BSUB -o "/home/hb765588/linnea/results/generation/run_%J_c/out.out" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB

module load python/3.6.0

cd /home/hb765588/linnea
source linnea_venv/bin/activate
cd /home/hb765588/linnea/results/generation
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}
python3 /home/hb765588/linnea/linnea/experiments/run_experiments.py ${LSB_JOBINDEX} 1 -c

exit