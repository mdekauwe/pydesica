#!/bin/bash
#PBS -m ae
#PBS -M mdekauwe\@gmail.com
#PBS -P w35
#PBS -q normal
#PBS -l walltime=00:30:00
#PBS -l ncpus=144
#PBS -l mem=288GB
#PBS -l wd
#PBS -j oe
#PBS -e logs/error_chunk_1_rows_0_4.txt
#PBS -o logs/log_chunk_1_rows_0_4.txt

ulimit -s unlimited
set -eu

row_start=0
row_end=4
col_start=0
col_end=27

cd $PBS_O_WORKDIR

# Process all of the valid points in parallel
export NP_DEBUG=1
./utils/node_parallel.sh python src/extract_forcing_timeseries_from_GSWP3.py {1} {2} :::: gswp3_land_sea/valid_points_"$row_start"_"$row_end"_"$col_start"_"$col_end"
