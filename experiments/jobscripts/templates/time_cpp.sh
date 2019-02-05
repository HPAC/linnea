#!/usr/bin/env bash

#{directive} {flag_jobname} "time_cpp_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_results_path}/{name}/execution/cpp/cout.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{exclusive}

module switch intel intel/18.0
module load gcc/7

cd {linnea_results_path}
mkdir -p {name}/execution/cpp
cd {name}/execution/cpp

expname=$(printf "{name}%03d" ${var_array_idx})
runner="{output_code_path}/${{expname}}/Cpp/build/${{expname}}"

if [ -f $runner ]; then
    $runner
else
    echo "File not found: $runner"
fi

exit
