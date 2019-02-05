#!/usr/bin/env bash

#{directive} {flag_jobname} "time_gen_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_results_path}/{name}/generation/{strategy_name}/cout.txt"
#{directive} {flag_time} {time}
#{directive} {flag_memory}{memory}
#{directive} {flag_group} {group}
#{directive} {flag_model} {model}
{exclusive}

echo $(printf "{name}%03d" ${var_array_idx})

module load python/3.6.0

source {linnea_virtualenv_path}/bin/activate
cd {linnea_results_path}
mkdir -p {name}/generation/{strategy_name}
cd {name}/generation/{strategy_name}
python3 {linnea_src_path}/experiments/experiments.py time_generation {name} -j=${var_array_idx} -r=1 -{strategy} -m={merging}

exit