#!/usr/bin/env bash

#{directive} {flag_jobname} "time_matlab_t{threads}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_matlab_t{threads}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{spec_exclusive}
{spec_node}

module purge
module load MISC
module load matlab/2019b

cd {linnea_results_path}
mkdir -p {name}/execution/matlab/t{threads}
cd {name}/execution/matlab/t{threads}

mkdir -p logs
export MATLAB_LOG_DIR={linnea_results_path}/execution/matlab/logs
export MATLABPATH={linnea_lib_path}/MatrixGeneratorMatlab

exppath=$(printf "{output_code_path}/{name}%03d/Matlab" ${var_array_idx})
runner="${{exppath}}/{runner_name}.m"

if [ -f $runner ]; then
    # it is not possible to run runner.m directly here because of how importing in Matlab works.
    matlab -nodisplay -nodesktop -nosplash -logfile matlab.log <<EOF
    addpath('${{exppath}}');
    {runner_name};
    quit();
EOF
else
    echo "File not found: $runner"
fi

exit