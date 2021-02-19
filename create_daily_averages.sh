#!/bin/sh

CASE_NAMES="timestep_presaer_cld_10s_lower_tau2"
#DAYS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"
DAYS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15"
OUTPUT_DIR=/p/lscratchh/santos36/timestep_daily_avgs

times="00000 03600 07200 10800 14400 18000 21600 25200 28800 32400 36000 39600 43200 46800 50400 54000 57600 61200 64800 68400 72000 75600 79200 82800"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    case_dir=/p/lscratchh/santos36/ACME/$case_name/run
    for day in $DAYS; do
        filenames=
        for time in $times; do
            filenames="$filenames $case_dir/$case_name.cam.h0.0001-01-$day-$time.nc"
        done
        echo "On day $day."
        ncra -O $filenames $OUTPUT_DIR/$case_name.0001-01-$day.nc
    done
done
