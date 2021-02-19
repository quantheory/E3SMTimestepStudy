#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_daily_avgs
OUT_DIR=/p/lscratchh/santos36/timestep_daily_avgs_lat_lon
CASE_NAMES="timestep_presaer_cld_10s_lower_tau2"
#DAYS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"
DAYS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    for day in $DAYS; do
        ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$case_name.0001-01-$day.nc
    done
done
