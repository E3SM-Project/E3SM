#!/bin/bash
#
#
# Batch script to move and import mapping files to inputdata
# for several resolutions.
#

#----------------------------------------------------------------------

if [ -z "$CSMDATA" ]; then
  CSMDATA=/fis/cgd/cseg/csm/inputdata
fi

if [ ! -d "$CSMDATA" ]; then
   echo "Environment variable CSMDATA is not set to a valid directory!"
   exit 1
fi

mapdir="lnd/clm2/mappingdata/maps"
if [ ! -d "$CSMDATA/$mapdir" ]; then
   echo "Environment variable CSMDATA is not set to a valid inputdata directory!"
   exit 1
fi

if [ -z "$SVN_INP_DIR" ]; then
  SVN_INP_DIR=https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata
fi

if [ $# -gt 0 ]; then
   resols=""
   for arg in $@; do
      resols="$resols $arg"
   done
else
   echo "Run for all valid resolutions"
   resols=`../../bld/queryDefaultNamelist.pl -res list -silent`
fi
echo "Move and import mapping files for this list of resolutions: $resols"

#----------------------------------------------------------------------

for res in $resols; do
  echo "Move and import mapping files for: $res"
  dir=$mapdir/$res
  #----------------------------------------------------------------------
  files=(map_*${res}*_aave_da_c??????.nc)
  if [ ${#files[*]} -lt 2 ]; then
     echo "No mappingfiles found for $res"
     exit 2
  else
     if [ ! -d "$CSMDATA/$dir" ]; then
        echo "Create mapping directory: $CSMDATA/$dir"
        mkdir $CSMDATA/$dir
        svn mkdir $SVN_INP_URL/$dir -m "Create mapping directory for $res"
     fi
     for file in ${files[*]}; do
        echo "Copy and import file $file"
        cp -p $file $CSMDATA/$dir
        if [ $? -ne 0 ]; then
           echo "Problem copying file: $file"
           exit 3
        fi
        chmod 0444 $CSMDATA/$dir/$file
        if [ $? -ne 0 ]; then
           echo "Problem chmod on file: $file"
           exit 4
        fi
        svn import $CSMDATA/$dir/$file $SVN_INP_DIR/$dir/$file -m "Mapping file for $res"
        if [ $? -ne 0 ]; then
           echo "Problem doing svn import on file: $file"
           exit 4
        fi
     done
  fi
done
