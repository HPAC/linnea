#!/usr/bin/env bash

#BSUB -J "linnea_time_julia[1-31]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB
#BSUB -P aices2


for i in {0..19}; do
    if [ -f ${HOME}/linnea/output/lamp_example${i}c/runner.jl ]; then
        echo "File found!"
        #JULIA_LOAD_PATH=/home/mc306023/projects/linnea/generators/MatrixGenerator.jl/src /home/mc306023/projects/linnea/install/julia/bin/julia ${DIR}/Julia/seed_${SEED}_matrix_chain_${i}/runner.jl > TEST_OUT_${SEED}_${i} 2>&1
    else
        echo "File not found: ${HOME}/linnea/output/lamp_example${i}c/runner.jl"
    fi
done

#exit