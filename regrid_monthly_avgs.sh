#!/bin/sh

REMAP_FILE_NE30=/usr/gdata/climdat/maps/map_ne16np4_to_ne30np4_aave.20160601.nc
REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_monthly_avgs
OUT_DIR=/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon
CASE_NAMES="timestep_ctrl"
#MONTHS="01 02 03 04 05 06 07 08 09 10 11 12"
MONTHS="01"
YEAR="0007"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    filenames=
    for month in $MONTHS; do
        date
        echo "On month $month"
        filenames="$filenames $IN_DIR/$case_name.$YEAR-$month.nc"
#        ncremap -m $REMAP_FILE_NE30 -o $IN_DIR/$case_name.$YEAR-$month.nc $IN_DIR/$case_name.cam.h0.$YEAR-$month.nc
    done
    ls $filenames | ncremap -p bck -j 12 -m $REMAP_FILE -O $OUT_DIR
done
