#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_daily_avgs
OUT_DIR=/p/lscratchh/santos36/timestep_daily_avgs_lat_lon
CASE_NAMES="timestep_presaer_ZM_10s timestep_MG2_10s timestep_CLUBB_10s timestep_dyn_10s"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$case_name.0001-01-*.nc
done
