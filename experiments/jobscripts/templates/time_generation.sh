#!/usr/bin/env bash

#BSUB -J "linnea_time[1-{jobs}]" # job name
#BSUB -o "linnea/results/generation/run_%J/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R model=={model}
{exclusive}

module load python/3.6.0

cd ${{HOME}}/linnea
source linnea_venv/bin/activate
cd ${{HOME}}/linnea/results/generation
mkdir -p run_${{LSB_JOBID}}
cd run_${{LSB_JOBID}}
python3 ${{HOME}}/linnea/linnea/experiments/experiments.py time_generation {name} -j=${{LSB_JOBINDEX}} -r={repetitions} -{strategy} -m={merging}

exit