#!/usr/bin/env bash

#{directive} {flag_jobname} "time_matlab_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_matlab_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
{spec_model}
{spec_exclusive}

module load MISC
module load matlab/2018b

cd {linnea_results_path}
mkdir -p {name}/execution/matlab
cd {name}/execution/matlab

mkdir -p logs
export MATLAB_LOG_DIR={linnea_results_path}/execution/matlab/logs
export MATLABPATH={linnea_lib_path}/MatrixGeneratorMatlab

exppath=$(printf "{output_code_path}/{name}%03d/Matlab" ${var_array_idx})
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