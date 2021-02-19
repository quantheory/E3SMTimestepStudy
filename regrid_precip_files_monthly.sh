#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_precip
OUT_DIR=/p/lscratchh/santos36/timestep_precip_lat_lon
CASE_NAMES="timestep_ctrl timestep_all_10s"
#MONTHS="01 02 03 04 05 06 07 08 09 10 11 12"
MONTHS="01 02"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    for month in $MONTHS; do
        file_name=$case_name.freq.0004-$month.nc
        tmp_file_name=$case_name.freq.0004-$month-tmp.nc
        ncpdq -a nbins,ncol $IN_DIR/$file_name $IN_DIR/$tmp_file_name
        ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$tmp_file_name
        ncpdq -a lat,lon,nbins $OUT_DIR/$tmp_file_name $OUT_DIR/$file_name
        rm $IN_DIR/$tmp_file_name $OUT_DIR/$tmp_file_name
    done
done
