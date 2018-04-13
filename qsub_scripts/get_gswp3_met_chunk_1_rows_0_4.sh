#!/bin/bash
#PBS -m ae
#PBS -M mdekauwe\@gmail.com
#PBS -P w35
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -l ncpus=128
#PBS -l mem=16GB
#PBS -l wd
#PBS -j oe
#PBS -e logs/error.txt
#PBS -o logs/log.txt

ulimit -s unlimited

row_start=0
row_end=4
col_start=0
col_end=27

cd $PBS_O_WORKDIR

# Clear list of valid points
echo >> valid_points

# Work out which points are valid, writing them into the file 'valid_points'
for row in $(seq $row_start $row_end)
do
    for col in $(seq $col_start $col_end)
    do
        landsea=$(python gswp3_land_sea/check_nsw_gswp3_land_sea_mask.py $row $col gswp3_land_sea/nsw_gswp3_land_sea_mask.bin)
        if [ $landsea -eq 0 ]
        then
            echo $row $col >> valid_points
        endif
    done
done

# Process all of the valid points in parallel
export NP_PARALLEL_ARGS="--colsep ' '" # Required to split the input into columns

./utils/node_parallel.sh python src/extract_forcing_timeseries_from_GSWP3.py {1} {2} :::: valid_points
