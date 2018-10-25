#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-{jobs}]" # job name
#BSUB -o "linnea/output/cout_{name}.txt" # job output
#BSUB -W {time}:00            # limits in hours:minutes
#BSUB -M {memory}            # memory in MB
#BSUB -P {group}
#BSUB -R model=={model}


module load python/3.6.0

source ${{HOME}}/linnea/linnea_venv/bin/activate
python3 ${{HOME}}/linnea/linnea/experiments/experiments.py generate_code {name} -j=${{LSB_JOBINDEX}} -{strategy}

if {compile}; then
    module load cmake/3.10.1
    module load gcc/7

    cd ${{HOME}}/linnea/output/{name}${{LSB_JOBINDEX}}/Cpp
    mkdir -p build
    cd build
    rm -rf *
    cmake -DCMAKE_PREFIX_PATH=${{HOME}}/libraries/MatrixGeneratorCpp/ ..
    make
fi

exit
