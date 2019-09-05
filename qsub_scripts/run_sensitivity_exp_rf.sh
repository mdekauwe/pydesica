#!/bin/bash

#PBS -m ae
#PBS -P w35
#PBS -q normal
#PBS -M mdekauwe@gmail.com
#PBS -l mem=128GB
#PBS -l ncpus=128
#PBS -l walltime=02:00:00
#PBS -l wd
#PBS -j oe
#PBS -l other=gdata1

module load dot
source activate sci

# For example, if you request 128 CPUs in the normal queue, and 128GB of RAM,
# your job will be using 8 nodes (there are 16 CPU cores in each node in the
# normal queue on Raijin) and 16GB of RAM will be allocated on each of the
# 8 nodes.
python src/run_sensitivity_exp.py "rf" 128
#python src/run_sensitivity_exp.py "rf" 532 # multiples of 28, so 19 nodes
