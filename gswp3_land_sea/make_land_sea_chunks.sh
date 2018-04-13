#!/bin/bash

# For each chunk, make a list of valid points in a text file.

col_start=0
col_end=27
st=0
en=19
step=5
for row_start in $(seq $st $step $en )
do
    row_end=$((row_start+$step-1))
    rm -f valid_points_"$row_start"_"$row_end"_"$col_start"_"$col_end"

    for row in $(seq $row_start $row_end)
    do
        for col in $(seq $col_start $col_end)
        do
            echo $row $col
            landsea=$(python check_nsw_gswp3_land_sea_mask.py $row $col nsw_gswp3_land_sea_mask.bin)
            if [ "$landsea" -eq 0 ]
            then
                echo $row $col >> valid_points_"$row_start"_"$row_end"_"$col_start"_"$col_end"
            fi
        done
    done
done
