#!/usr/bin/env bash

#{directive} {flag_jobname} "gen_{name}{ref}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/gen_{name}{ref}{string_array_idx}.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
{options}

module load DEVELOP
module load LIBRARIES
module unload intel/19.0
module load gcc/{gcc_version}
module load intelmkl/{intel_version}
module load cmake/{cmake_version}
module load python/{python_version}

source {linnea_virtualenv_path}/bin/activate
python3 {linnea_src_path}/experiments/experiments.py generate_code {name} -j=${var_array_idx} {args}

expname=$(printf "{name}%03d" ${var_array_idx})

cd {output_code_path}/${{expname}}/Cpp
mkdir -p build
cd build
rm -rf *
cmake -DCMAKE_PREFIX_PATH={linnea_lib_path}/ ..
make

exit
