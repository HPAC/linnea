#!/usr/bin/env bash

#BSUB -J "time_julia_{name}[1-{jobs}]" # job name
#BSUB -o "{linnea_results_path}/{name}/execution/julia/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R {model}
{exclusive}

module load gcc/7

cd {linnea_results_path}/
mkdir -p {name}/execution/julia
cd {name}/execution/julia

runner=$(printf "{output_code_path}/{name}%03d/Julia/runner.jl" $LSB_JOBINDEX)

if [ -f $runner ]; then
    {linnea_julia_path}/julia $runner
else
    echo "File not found: $runner"
fi

exit