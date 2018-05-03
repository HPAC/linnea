#!/usr/bin/env bash

#BSUB -J "linnea_time_matlab[1-31]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB
#BSUB -P aices2

# TODO commanline options for constructive/exhaustive, Matlab/C++/Julia?

module load MISC
module load matlab/2018a

cd ${HOME}/linnea/results/execution
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}

mkdir -p logs
export MATLAB_LOG_DIR=${HOME}/linnea/results/execution/run_${LSB_JOBID}/logs
export MATLABPATH=${HOME}/matrix_generators/MatrixGeneratorMatlab

if [ -f ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}c/Matlab/runner.m ]; then
    # it is not possible to run runner.m directly here because of how importing in Matlab works.
    matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile matlab.log <<EOF
    addpath('${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}c/Matlab/');
    runner;
    quit();
EOF
else
    echo "File not found: ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}c/Matlab/runner.m"
fi

exit