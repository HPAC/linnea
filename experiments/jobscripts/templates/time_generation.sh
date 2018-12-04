#!/usr/bin/env bash

#BSUB -J "linnea_time[1-{jobs}]" # job name
#BSUB -o "{linnea_results_path}/{name}/generation/{strategy_name}/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R {model}
{exclusive}

echo $(printf "{name}%03d" $LSB_JOBINDEX)

module load python/3.6.0

source {linnea_virtualenv_path}/bin/activate
cd {linnea_results_path}
mkdir -p {name}/generation/{strategy_name}
cd {name}/generation/{strategy_name}
python3 {linnea_src_path}/experiments/experiments.py time_generation {name} -j=$LSB_JOBINDEX -r=1 -{strategy} -m={merging}

exit