#!/usr/bin/env bash

#BSUB -J "linnea_time[1-{jobs}]" # job name
#BSUB -o "linnea/results/{name}/generation/{strategy_name}/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R model=={model}
{exclusive}

echo "{name}${{LSB_JOBINDEX}}"

module load python/3.6.0

cd ${{HOME}}/linnea
source linnea_venv/bin/activate
cd ${{HOME}}/linnea/results/
mkdir -p {name}/generation/{strategy_name}
cd {name}/generation/{strategy_name}
python3 ${{HOME}}/linnea/linnea/experiments/experiments.py time_generation {name} -j=${{LSB_JOBINDEX}} -r={repetitions} -{strategy} -m={merging}

exit