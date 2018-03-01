#!/bin/bash

# This script overwrites the area_a and area_b fields in a SCRIP mapping
# file with the grid_area fields from the corresponding SCRIP grid files
# (using the grid_file_src and grid_file_dst global metadata).

if [ -z $1 ]; then
  echo "Usage: ./`basename $0` SCRIP_map.nc [-update_S]"
  echo "-update_S    Set this flag to also adjust S based on area"
  exit 0
fi

SRC_GRID=`ncks -M $1 | grep grid_file_src | sed  -r "s/^.+value = //"`
DST_GRID=`ncks -M $1 | grep grid_file_dst | sed  -r "s/^.+value = //"`

if [ -e $SRC_GRID ]; then
  echo "Source file is $SRC_GRID"
else
  echo "ERROR, can not find $SRC_GRID"
  stop 1
fi

if [ -e $DST_GRID ]; then
  echo Destination file is $DST_GRID
else
  echo "ERROR, can not find $DST_GRID"
  stop 1
fi

if [ "$2" == "-update_S" ]; then
  lupdate_S=True
else
  lupdate_S=False
fi

ncl overwrite_map_area.ncl map_file=\"$1\" src_file=\"$SRC_GRID\" dst_file=\"$DST_GRID\" lupdate_S=$lupdate_S
