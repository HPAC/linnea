#!/usr/bin/env bash

bsub < linnea/linnea/experiments/jobscript_c.sh

bsub -w "done(linnea_gen[1-31])" python3 /home/hb765588/linnea/linnea/experiments/collect_data.py

exit