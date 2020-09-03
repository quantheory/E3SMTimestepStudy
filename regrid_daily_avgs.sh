#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_daily_avgs
OUT_DIR=/p/lscratchh/santos36/timestep_daily_avgs_lat_lon
CASE_NAMES="timestep_all_rad_10s"
DAYS="20 21 22 23 24 25 26 27 28 29 30"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    for day in $DAYS; do
        ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$case_name.0001-01-$day.nc
    done
done
