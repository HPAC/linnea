#!/usr/bin/env bash

#BSUB -J "linnea_time_cpp[1-31]" # job name
#BSUB -o "linnea/results/execution/run_%J/cout.txt" # job output
#BSUB -W 2:00               # limits in hours:minutes
#BSUB -M 2000               # memory in MB
#BSUB -P aices2

# TODO commanline options for constructive/exhaustive, Matlab/C++/Julia?

module load gcc/5

cd ${HOME}/linnea/results/execution
mkdir -p run_${LSB_JOBID}
cd run_${LSB_JOBID}

if [ -f ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Cpp/build/lamp_example${LSB_JOBINDEX} ]; then
    ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Cpp/build/lamp_example${LSB_JOBINDEX}
else
    echo "File not found: ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}/Cpp/build/lamp_example${LSB_JOBINDEX}"
fi

exit

