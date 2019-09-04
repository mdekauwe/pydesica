#!/bin/bash

#PBS -m ae
#PBS -P w35
#PBS -q normalbw
#PBS -M mdekauwe@gmail.com
#PBS -l mem=16GB
#PBS -l ncpus=512
#PBS -l walltime=00:05:00
#PBS -l wd
#PBS -j oe
#PBS -l other=gdata1

module load dot
source activate sci

#python src/run_sensitivity_exp.py "rf" 16
python src/run_sensitivity_exp.py "rf" 512
