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

export MATLAB_LOG_DIR=${HOME}/linnea/results/execution/run_${LSB_JOBID}
export MATLABPATH=${HOME}/matrix_generators/MatrixGeneratorMatlab

# it is not possible to run runner.m directly here because of how importing in Matlab works.

for i in {1..31}; do
    if [ -f ${HOME}/linnea/output/lamp_example${i}c/Matlab/runner.m ]; then
        matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile matlab.log <<EOF
        run ${HOME}/linnea/output/lamp_example${i}c/Matlab/runner;quit();
EOF
    else
        echo "File not found: ${HOME}/linnea/output/lamp_example${i}c/Matlab/runner.m"
    fi
    if [ -f ${HOME}/linnea/output/lamp_example${i}e/Matlab/runner.m ]; then
        matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile matlab.log <<EOF
        run ${HOME}/linnea/output/lamp_example${i}e/Matlab/runner;quit();
EOF
    else
        echo "File not found: ${HOME}/linnea/output/lamp_example${i}e/Matlab/runner.m"
    fi
done

exit