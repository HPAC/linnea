#!/usr/bin/env bash

#BSUB -J "linnea_gen[1-31]" # job name
#BSUB -o "linnea/results/generation/run_%J/cout.txt" # job output
#BSUB -W 2:00              # limits in hours:minutes
#BSUB -M 2000              # memory in MB
#BSUB -P aices2

module load python/3.6.0

source ${HOME}/linnea/linnea_venv/bin/activate
cd ${HOME}/linnea/results/generation
mkdir -p run_${LSB_JOBID}
# cd run_${LSB_JOBID}
python3 ${HOME}/linnea/linnea/experiments/experiments.py generate_code -j=${LSB_JOBINDEX} -c

module load cmake/3.10.1
module load gcc/5

cd ${HOME}/linnea/output/lamp_example${LSB_JOBINDEX}c/Cpp
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=${HOME}/libraries/MatrixGeneratorCpp/ ..
make

exit

