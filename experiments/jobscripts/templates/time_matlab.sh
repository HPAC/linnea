#!/usr/bin/env bash

#BSUB -J "linnea_time_matlab[1-{jobs}]" # job name
#BSUB -o "linnea/results/{name}/execution/matlab/cout.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R {model}
{exclusive}

module load MISC
module load matlab/2018b

cd ${{HOME}}/linnea/results
mkdir -p {name}/execution/matlab
cd {name}/execution/matlab

mkdir -p logs
export MATLAB_LOG_DIR=${{HOME}}/linnea/results/execution/matlab/logs
export MATLABPATH=${{HOME}}/repositories/MatrixGeneratorMatlab

exppath=$(printf "${{HOME}}/linnea/output/{name}%03d/Matlab" $LSB_JOBINDEX)
runner="${{exppath}}/runner.m"

if [ -f $runner ]; then
    # it is not possible to run runner.m directly here because of how importing in Matlab works.
    matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile matlab.log <<EOF
    addpath('${{exppath}}');
    runner;
    quit();
EOF
else
    echo "File not found: $runner"
fi

exit