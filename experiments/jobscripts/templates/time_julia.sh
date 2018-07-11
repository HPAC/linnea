#!/usr/bin/env bash

#BSUB -J "linnea_time_julia[1-{jobs}]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R model=={model}
{exclusive}

module load gcc/7

cd ${{HOME}}/linnea/results/execution
mkdir -p run_${{LSB_JOBID}}
cd run_${{LSB_JOBID}}

if [ -f ${{HOME}}/linnea/output/{name}${{LSB_JOBINDEX}}/Julia/runner.jl ]; then
    # echo "File found!: ${{HOME}}/linnea/output/{name}${{LSB_JOBINDEX}}/Julia/runner.jl"
    ${{HOME}}/julia/julia-0.6.3-haswell/julia ${{HOME}}/linnea/output/{name}${{LSB_JOBINDEX}}/Julia/runner.jl
else
    echo "File not found: ${{HOME}}/linnea/output/{name}${{LSB_JOBINDEX}}/Julia/runner.jl"
fi

exit