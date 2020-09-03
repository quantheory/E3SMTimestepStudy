#!/bin/sh

CASE_NAMES="timestep_ctrl timestep_all_10s"
MONTHS="01 02"
INPUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon
OUTPUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    for month in $MONTHS; do
        date
        echo "On month $month."
        ncra $INPUT_DIR/$case_name.00{02,03,04}-$month.nc $OUTPUT_DIR/${case_name}_${month}_climo.nc
    done
done

MONTHS="03 04 05 06 07 08 09 10 11 12"

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    for month in $MONTHS; do
        date
        echo "On month $month."
        ncra $INPUT_DIR/$case_name.00{01,02,03}-$month.nc $OUTPUT_DIR/${case_name}_${month}_climo.nc
    done
done

for case_name in $CASE_NAMES; do
    echo "On case $case_name"
    echo "Seasonal average DJF"
    ncra -w 31,31,28 $OUTPUT_DIR/${case_name}_{12,01,02}_climo.nc $OUTPUT_DIR/${case_name}_DJF_climo.nc
    echo "Seasonal average MAM"
    ncra -w 31,30,31 $OUTPUT_DIR/${case_name}_{03,04,05}_climo.nc $OUTPUT_DIR/${case_name}_MAM_climo.nc
    echo "Seasonal average JJA"
    ncra -w 30,31,31 $OUTPUT_DIR/${case_name}_{06,07,08}_climo.nc $OUTPUT_DIR/${case_name}_JJA_climo.nc
    echo "Seasonal average SON"
    ncra -w 30,31,30 $OUTPUT_DIR/${case_name}_{09,10,11}_climo.nc $OUTPUT_DIR/${case_name}_SON_climo.nc
    echo "Annual average"
    ncra -w 31,28,31,30,31,30,31,31,30,31,30,31 $OUTPUT_DIR/${case_name}_{01,02,03,04,05,06,07,08,09,10,11,12}_climo.nc $OUTPUT_DIR/${case_name}_ANN_climo.nc
done
