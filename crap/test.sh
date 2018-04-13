#!/bin/bash

row_start=15
row_end=19
col_start=0
col_end=27

row=$row_start
while [ $row -le $row_end ]
do
    col=$col_start
    while [ $col -le $col_end ]
    do
        landsea=$(python gswp3_land_sea/check_nsw_gswp3_land_sea_mask.py $row $col gswp3_land_sea/nsw_gswp3_land_sea_mask.bin)
        if [ $landsea -eq 0 ]
        then
            echo $landsea $row $col
        else
            echo $landsea $row $col
        fi

        let col=col+1

    done
    let row=row+1
done
