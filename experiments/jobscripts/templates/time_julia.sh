#!/usr/bin/env bash

#{directive} {flag_jobname} "time_{output_subdir}_t{threads}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_{output_subdir}_t{threads}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_execution}
#{directive} {flag_memory}{memory}
#{directive} {flag_threads} {threads}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{spec_exclusive}
{spec_node}

module switch intel/19.0 clang/{clang_version}
module load DEVELOP

cd {linnea_results_path}/
mkdir -p {name}/execution/{output_subdir}/t{threads}
cd {name}/execution/{output_subdir}/t{threads}

runner=$(printf "{output_code_path}/{name}%03d/Julia/{runner_name}.jl" ${var_array_idx})

if [ -f $runner ]; then
    {linnea_julia_path}/julia $runner
else
    echo "File not found: $runner"
fi

exit