#!/bin/bash

row_start=0
row_end=4
col_start=0
col_end=27

core=0
new_core=0
row=$row_start
while [ $row -le $row_end ]
do
    col=$col_start
    while [ $col -le $col_end ]
    do
        landsea=$(python gswp3_land_sea/check_nsw_gswp3_land_sea_mask.py $row $col gswp3_land_sea/nsw_gswp3_land_sea_mask.bin)
        if [ $landsea -eq 0 ]
        then
            let new_core=1
        else
            let new_core=0
        fi

        # only increment the core if we found a valid pixel to run
        (( core += new_core ))

        echo $row, $col, $landsea $core

        let col=col+1
    done
    let row=row+1
done
