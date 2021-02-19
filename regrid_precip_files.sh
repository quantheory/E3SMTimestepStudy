#!/bin/sh

REMAP_FILE=/usr/gdata/climdat/maps/map_ne30np4_to_fv128x256_aave.20160301.nc
IN_DIR=/p/lscratchh/santos36/timestep_precip
OUT_DIR=/p/lscratchh/santos36/timestep_precip_lat_lon
CASE_NAMES="timestep_presaer_cld_10s_lower_tau2"

for case_name in $CASE_NAMES; do
    date
    echo "On case $case_name"
    ncpdq -a nbins,ncol $IN_DIR/$case_name.freq.short.d03-d15.nc $IN_DIR/$case_name.freq.short.d03-d15_tmp.nc
    ncremap -m $REMAP_FILE -O $OUT_DIR $IN_DIR/$case_name.freq.short.d03-d15_tmp.nc
    ncpdq -a lat,lon,nbins $OUT_DIR/$case_name.freq.short.d03-d15_tmp.nc $OUT_DIR/$case_name.freq.short.d03-d15.nc
    rm $IN_DIR/$case_name.freq.short.d03-d15_tmp.nc $OUT_DIR/$case_name.freq.short.d03-d15_tmp.nc
done
