#!/bin/bash
#PBS -m ae
#PBS -M mdekauwe\@gmail.com
#PBS -P w35
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=128
#PBS -l mem=64GB
#PBS -l wd
#PBS -j oe
#PBS -e logs/error.txt
#PBS -o logs/log.txt

ulimit -s unlimited

row_start=15
row_end=19
col_start=0
col_end=27

cd $PBS_O_WORKDIR

core=0
new_core=1
row=$row_start
while [ $row -le $row_end ]
do
    col=$col_start
    while [ $col -le $col_end ]
    do
        landsea=$(python gswp3_land_sea/check_nsw_gswp3_land_sea_mask.py $row $col gswp3_land_sea/nsw_gswp3_land_sea_mask.bin)
        if [ $landsea -eq 0 ]
        then
            pbsdsh -n $core python src/extract_forcing_timeseries_from_GSWP3.py $row $col
            let new_core=0
        else
            let new_core=1
        fi

        let col=col+1

        # only increment the core if we found a valid pixel to run
        if [ $new_core -eq 0 ]
        then
            let core=core+1
        fi

        if [ $core -ge $PBS_NCPUS ]
        then
            let core=0
            sleep 3m
        fi
    done
    let row=row+1
done

wait
