#!/usr/bin/env bash

#{directive} {flag_jobname} "time_cpp_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_cpp_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
{spec_model}
{spec_exclusive}

module switch intel intel/{intel_version}
module load gcc/{gcc_version}

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
