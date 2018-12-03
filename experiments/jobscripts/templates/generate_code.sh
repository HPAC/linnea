#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-{jobs}]" # job name
#BSUB -o "{linnea_output_path}/cout_{name}.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R {model}


module load python/3.6.0

source {linnea_virtualenv_path}/bin/activate
python3 {linnea_src_path}/linnea/experiments/experiments.py generate_code {name} -j=$LSB_JOBINDEX -{strategy}

if {compile}; then
    module load cmake/3.10.1
    module switch intel intel/18.0
    module load gcc/7

    expname=$(printf "{name}%03d" $LSB_JOBINDEX)

    cd {linnea_code_path}/${{expname}}/Cpp
    mkdir -p build
    cd build
    rm -rf *
    cmake -DCMAKE_PREFIX_PATH={linnea_lib_path}/MatrixGeneratorCpp/ ..
    make
fi

exit
