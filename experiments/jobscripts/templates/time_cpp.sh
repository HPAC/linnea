#!/usr/bin/env bash

#BSUB -J "linnea_time_cpp[1-{jobs}]" # job name
#BSUB -o "{linnea_results_path}/{name}/execution/cpp/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R {model}
{exclusive}

module switch intel intel/18.0
module load gcc/7

cd {linnea_results_path}
mkdir -p {name}/execution/cpp
cd {name}/execution/cpp

expname=$(printf "{name}%03d" $LSB_JOBINDEX)
runner="{linnea_code_path}/${{expname}}/Cpp/build/${{expname}}"

if [ -f $runner ]; then
    $runner
else
    echo "File not found: $runner"
fi

exit
