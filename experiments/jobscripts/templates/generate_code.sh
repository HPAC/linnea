#!/usr/bin/env bash

#{directive} {flag_jobname} "gen_{strategy}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{output_code_path}/cout_{name}.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}


module load python/3.6.0

source {linnea_virtualenv_path}/bin/activate
python3 {linnea_src_path}/experiments/experiments.py generate_code {name} -j=${var_array_idx} -{strategy}

if {compile}; then
    module load cmake/3.10.1
    module switch intel intel/18.0
    module load gcc/7

    expname=$(printf "{name}%03d" ${var_array_idx})

    cd {output_code_path}/${{expname}}/Cpp
    mkdir -p build
    cd build
    rm -rf *
    cmake -DCMAKE_PREFIX_PATH={linnea_lib_path}/ ..
    make
fi

exit
