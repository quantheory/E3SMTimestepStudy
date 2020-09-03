#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_monthly_avgs
OUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon
CASE_NAMES="timestep_ctrl timestep_all_10s"
MONTHS="01 02"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    for month in $MONTHS; do
        date
        echo "On month $month"
        ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$case_name.0004-$month.nc
    done
done
