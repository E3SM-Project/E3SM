#!/bin/bash

#KJE 1/2013 evanskj@ornl.gov
#This is the master script to set the parameters and paths to run the LIVV kit on NERSC euclid 
#Efforts funded by DOE BER PISCEES SciDAC project
#Currently it is designed specifically for the GLIDE dycore of the CISM model, because it is
#designed to read its output

# current environment to be loaded on NERSC's Carver machine 1/2/2013
#module unload netcdf/4.1.1
#module load netcdf/4.1.3
#module load ncar/6.0.0
#module load nco
#module load python

# and set env so that files are written with permissions to read the webpage
umask 002 

#user added comment of analysis to be performed
export COMMENT="evaluating code and test suite for CISM2.0 release"
export USERNAME='evanskj'

# flags to run the production analysis - you always run the test suite
RUN_GIS=0 
# ANT not activated yet, always set to 0. If 1, add -a to end of python command line at the bottom
RUN_ANT=0 

#specify location of executables 
#/reg_test and /livv needs to be placed in the subdirectory below
export TEST_FILEPATH="$GSCRATCH/higher-order"
export SCRIPT_PATH="$TEST_FILEPATH/livv"

#specify location where the html files will be sent 
export HTML_PATH="/project/projectdirs/piscees/www"
# providing a username creates a directory by that name in the location above in which all the web files will go
export HTML_LINK="portal.nersc.gov/project/~piscees"

if [ $RUN_GIS -eq 1 ] 
then
#Production run and production benchmark information
#directory of run
  export GIS_FILEPATH="$TEST_FILEPATH/gis_5km_long"
#configure file
  export GIS_CONFIG="$GIS_FILEPATH/gis_5km.config"
#sets the filepath to the benchmark netcdf file 
  export GIS_BENCH_NETCDF_FILEPATH="$GIS_FILEPATH/data"
  export GIS_BENCH_NETCDF_FILENAME="gis_5km.ice2sea.init.nc"
  export GIS_BENCH_NETCDF_FILE="$GIS_BENCH_NETCDF_FILEPATH/$GIS_BENCH_NETCDF_FILENAME"
#production run benchmark screen output for collecting convergence information
  export GIS_OUTPUT="$GIS_BENCH_NETCDF_FILEPATH/out.gnu"
# right now only plots the benchmark

#sets the filepath to the production output file to be analyzed
  export GIS_VAR_NETCDF_FILENAME="gis_5km.ice2sea.51-100.nc"
  export GIS_VAR_NETCDF_FILE="$GIS_BENCH_NETCDF_FILEPATH/$GIS_VAR_NETCDF_FILENAME"

#TODO once list of plots created, add feature to have user pick which plots to make, default provided
fi

if [ $RUN_ANT -eq 1 ] 
then
#  directory of run
  export ANT_FILEPATH="$TEST_FILEPATH/ant"
#  cofigure file
  export ANT_CONFIG="ant_5km.config"
#  production run screen output for collecting convergence information
  export ANT_OUTPUT="out.gnu"
fi

# date stamp of LIVV run to put with comments
NOW=$(date +"%m-%d-%Y-%r")
echo $NOW $COMMENT

# settings not generally altered, but leaving the option open for future extension
#location where the livv code is located
export PY_PATH="$SCRIPT_PATH/bin"
#location where the ncl directory of the ncl scripts and .nc files are located
export NCL_PATH="$SCRIPT_PATH/plots"
export DATA_PATH="$SCRIPT_PATH/data"

#command to run python script while inputting all of the files listed above
#NOTE: not all settings are required to run the python script
if [ $RUN_GIS -eq 1 ] 
then
python $PY_PATH/VV_main.py -p "$PY_PATH" -j "$HTML_PATH" -k "$NCL_PATH" -c "$GIS_CONFIG" -o "$GIS_OUTPUT" -s "$GIS_BENCH_NETCDF_FILE" -n "$GIS_VAR_NETCDF_FILE" -t "$TEST_FILEPATH" -i "$NOW" -m "$COMMENT" -u "$USERNAME" -d "$DATA_PATH" -l "$HTML_LINK" -g

else
python $PY_PATH/VV_main.py -p "$PY_PATH" -j "$HTML_PATH" -k "$NCL_PATH" -t "$TEST_FILEPATH" -i "$NOW" -m "$COMMENT" -u "$USERNAME" -d "$DATA_PATH" -l "$HTML_LINK" 
fi

#type "python VV_main -h" in the command line for a full list of options

chmod -R 2775 $HTML_PATH
chgrp -R piscees $HTML_PATH

