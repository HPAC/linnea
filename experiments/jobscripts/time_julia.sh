#!/usr/bin/env bash

#BSUB -J "linnea_time_julia[1-31]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB
#BSUB -P aices2
#BSUB -R model==Haswell_EP

# TODO commanline options for constructive/exhaustive, Matlab/C++/Julia?

module load gcc/7

cd ${HOME}/linnea/results/execution
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}

if [ -f ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Julia/runner.jl ]; then
    # echo "File found!: ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Julia/runner.jl"
    ${HOME}/julia/julia-0.6.3/julia ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Julia/runner.jl > julia_results_lamp_example${LSB_JOBINDEX}.txt 2>&1
else
    echo "File not found: ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Julia/runner.jl"
fi

exit