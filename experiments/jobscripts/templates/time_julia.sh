#!/usr/bin/env bash

#{directive} {flag_jobname} "time_julia_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_results_path}/{name}/execution/julia/cout.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{exclusive}

module load gcc/7

cd {linnea_results_path}/
mkdir -p {name}/execution/julia
cd {name}/execution/julia

runner=$(printf "{output_code_path}/{name}%03d/Julia/runner.jl" ${var_array_idx})

if [ -f $runner ]; then
    {linnea_julia_path}/julia $runner
else
    echo "File not found: $runner"
fi

exit