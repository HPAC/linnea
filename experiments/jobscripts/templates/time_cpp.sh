#!/usr/bin/env bash

#{directive} {flag_jobname} "time_cpp_t{threads}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_cpp_t{threads}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{spec_exclusive}
{spec_node}

export GOMP_CPU_AFFINITY="1 3 5 7 9 11 13 15 17 19 21 23 0 2 4 6 8 10 12 14 16 18 20 22"

module switch intel intel/{intel_version}
module load gcc/{gcc_version}

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
