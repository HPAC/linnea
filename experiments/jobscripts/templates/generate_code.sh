#!/usr/bin/env bash

#{directive} {flag_jobname} "gen_{args}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/gen_{args}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}


module load python/{python_version}

source {linnea_virtualenv_path}/bin/activate
python3 {linnea_src_path}/experiments/experiments.py generate_code {name} -j=${var_array_idx} -{args}

if {compile}; then
    module load cmake/{cmake_version}
    module switch intel intel/{intel_version}
    module load gcc/{gcc_version}

    expname=$(printf "{name}%03d" ${var_array_idx})

    cd {output_code_path}/${{expname}}/Cpp
    mkdir -p build
    cd build
    rm -rf *
    cmake -DCMAKE_PREFIX_PATH={linnea_lib_path}/ ..
    make
fi

exit
