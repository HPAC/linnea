#!/usr/bin/env bash

#{directive} {flag_jobname} "time_cpp_t{threads}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_logs_path}/time_cpp_t{threads}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_nodes} {nodes}
#{directive} {flag_threads} {threads}
#{directive} {flag_tasks_per_core}{tasks_per_core}
{options}

module load DEVELOP
module load LIBRARIES
module unload intel/19.0
module load gcc/{gcc_version}
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
