#!/usr/bin/env bash

#BSUB -J "linnea_time_julia[1-31]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB
#BSUB -P aices2

# TODO commanline options for constructive/exhaustive, Matlab/C++/Julia?

module load MISC
module load matlab/2018a



cd ${HOME}/linnea/results/execution
mkdir -p run_test
cd run_test
# mkdir -p run_${LSB_JOBID}
# cd run_${LSB_JOBID}

MATLAB_LOG_DIR=${HOME}/linnea/results/execution/run_${LSB_JOBID}

MATLABPATH=${HOME}/matrix_generators/MatrixGeneratorMatlab matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile matlab.log <"run ${HOME}/linnea/output/lamp_example1c/Matlab/runner;"

# for i in {1..31}; do
#     if [ -f ${HOME}/linnea/output/lamp_example${i}c/Julia/runner.jl ]; then
#         echo "File found!: ${HOME}/linnea/output/lamp_example${i}c/Julia/runner.jl"
#         # ${DIR}/Julia/seed_${SEED}_matrix_chain_${i}/runner.jl > TEST_OUT_${SEED}_${i} 2>&1
#     else
#         echo "File not found: ${HOME}/linnea/output/lamp_example${i}c/Julia/runner.jl"
#     fi
# done

#exit