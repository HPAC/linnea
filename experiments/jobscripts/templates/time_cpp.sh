#!/usr/bin/env bash

#{directive} {flag_jobname} "time_cpp_t{threads}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_cpp_t{threads}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_threads} {threads}
#{directive} {flag_tasks_per_core}{tasks_per_core}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{spec_exclusive}
{spec_node}

module switch intel/19.0 clang/{clang_version}
module load DEVELOP
module load LIBRARIES
module load intelmkl/{intel_version}
module load cmake/{cmake_version}

cd {linnea_results_path}
mkdir -p {name}/execution/cpp/t{threads}
cd {name}/execution/cpp/t{threads}

expname=$(printf "{name}%03d" ${var_array_idx})
runner="{output_code_path}/${{expname}}/Cpp/build/{runner_name}"

if [ -f $runner ]; then
    $runner
else
    echo "File not found: $runner"
fi

exit
