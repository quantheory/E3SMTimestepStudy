#!/bin/sh

CASE_NAMES="timestep_ctrl"
MONTHS="03 04 05 06 07 08 09 10 11 12"
OUTPUT_DIR=/p/lscratchh/santos36/timestep_precip

times="00000 03600 07200 10800 14400 18000 21600 25200 28800 32400 36000 39600 43200 46800 50400 54000 57600 61200 64800 68400 72000 75600 79200 82800"

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    case_dir=/p/lscratchh/santos36/ACME/$case_name/run
    for month in $MONTHS; do
        date
        echo "On month $month."
        for time in $times; do
            echo "Time $time."
            outfile="$OUTPUT_DIR/$case_name.hourly.0003-$month-$time.nc"
            ncra -v PRECC,PRECL,PRECSC,PRECSL $case_dir/$case_name.cam.h0.0003-$month-*-$time.nc $outfile
            ncks -A -v lat,lon,area /p/lscratchh/santos36/timestep_monthly_avgs/timestep_ctrl.0001-01.nc $outfile
        done
    done
done

MONTHS="01 02"

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    case_dir=/p/lscratchh/santos36/ACME/$case_name/run
    for month in $MONTHS; do
        date
        echo "On month $month."
        for time in $times; do
            echo "Time $time."
            outfile="$OUTPUT_DIR/$case_name.hourly.0004-$month-$time.nc"
            ncra -v PRECC,PRECL,PRECSC,PRECSL $case_dir/$case_name.cam.h0.0004-$month-*-$time.nc $outfile
            ncks -A -v lat,lon,area /p/lscratchh/santos36/timestep_monthly_avgs/timestep_ctrl.0001-01.nc $outfile
        done
    done
done
