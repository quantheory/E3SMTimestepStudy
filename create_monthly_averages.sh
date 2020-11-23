#!/bin/sh

CASE_NAMES="timestep_ctrl"
MONTHS="01"
OUTPUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    case_dir=/p/lscratchh/santos36/ACME/$case_name/run
    for month in $MONTHS; do
        date
        echo "On month $month."
        ncra $case_dir/$case_name.cam.h0.0007-$month-* $OUTPUT_DIR/$case_name.0007-$month.nc
    done
done
