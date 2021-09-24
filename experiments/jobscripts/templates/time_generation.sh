#!/usr/bin/env bash

#{directive} {flag_jobname} "time_gen_{merging_label}_{name}{lsf_arrayjob}"
{slurm_arrayjob}
#{directive} {flag_output} "{linnea_output_path}/logs/time_gen_{merging_label}_{name}{string_array_idx}.txt"
#{directive} {flag_time} {time_generation}
#{directive} {flag_memory}{memory}
#{directive} {flag_tasks_per_core}{tasks_per_core}
{options}

echo $(printf "{name}%03d" ${var_array_idx})

module load DEVELOP
module load python/{python_version}

source {linnea_virtualenv_path}/bin/activate
cd {linnea_results_path}
mkdir -p {name}/generation/{dir_name}
cd {name}/generation/{dir_name}
python3 {linnea_src_path}/experiments/experiments.py time_generation {name} -j=${var_array_idx} -m={merging}

exit