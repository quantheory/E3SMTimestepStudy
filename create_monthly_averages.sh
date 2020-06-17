#!/bin/sh

CASE_NAMES="timestep_all_10s"
MONTHS="03 04 05 06 07 08 09 10 11 12"
OUTPUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    case_dir=/p/lscratchh/santos36/ACME/$case_name/run
    for month in $MONTHS; do
        date
        echo "On month $month."
        ncra $case_dir/$case_name.cam.h0.0002-$month-* $OUTPUT_DIR/$case_name.0002-$month.nc
    done
done
