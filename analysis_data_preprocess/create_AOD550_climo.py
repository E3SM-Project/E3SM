# flake8: noqa
#!/bin/bash

#####################################################
# Purpose: Script for processing AOD_550 climatology
# for use in the ACME observations repository.
#
# Author: Susannah M. Burrows (susannah.burrows@pnnl.gov)
# Last modified: 21 Aug 2015
#####################################################

# All:
WHICHSEASONS="ANN DJF MAM JJA SON 01 02 03 04 05 06 07 08 09 10 11 12"
#
# Annual mean only:
# WHICHSEASONS="ANN"

# Original file name
#INFILE="/home/burr114/PreAndPostProcessingScripts/process_observational_datasets/AOD_550/input/sat_c_s_anet_aod_901.nc"
INFILE="/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/AOD_550/original_data/sat_c_s_anet_aod_901.nc"

# Add lat_bounds/lon_bounds attributes to lat/lon
#ncatted -O -h -a lat_bounds,lat,o,f,"-89.5,89.5" $INFILE $INFILE
#ncatted -O -h -a lon_bounds,lon,o,f,"0.5,359.5" $INFILE $INFILE

# Directory for output files
#OUTDIR="/home/burr114/PreAndPostProcessingScripts/process_observational_datasets/AOD_550/processed"
OUTDIR="/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/AOD_550/climatology"

FILLVAL=-999.0

# First, split the file into months.
# File name convention:
# climatology files:  <root>_<season>_<optional version>_climo.nc .
#                     Seasons can be three-letter seasons such as ANN, DFJ, MAR, or two digit months such as 01,02,...,12.
ncks -O -d time,0 $INFILE $OUTDIR/AOD_550_01_climo.nc
ncks -O -d time,1 $INFILE $OUTDIR/AOD_550_02_climo.nc
ncks -O -d time,2 $INFILE $OUTDIR/AOD_550_03_climo.nc
ncks -O -d time,3 $INFILE $OUTDIR/AOD_550_04_climo.nc
ncks -O -d time,4 $INFILE $OUTDIR/AOD_550_05_climo.nc
ncks -O -d time,5 $INFILE $OUTDIR/AOD_550_06_climo.nc
ncks -O -d time,6 $INFILE $OUTDIR/AOD_550_07_climo.nc
ncks -O -d time,7 $INFILE $OUTDIR/AOD_550_08_climo.nc
ncks -O -d time,8 $INFILE $OUTDIR/AOD_550_09_climo.nc
ncks -O -d time,9 $INFILE $OUTDIR/AOD_550_10_climo.nc
ncks -O -d time,10 $INFILE $OUTDIR/AOD_550_11_climo.nc
ncks -O -d time,11 $INFILE $OUTDIR/AOD_550_12_climo.nc

# Create averages of DJF, MAM, JJA, SON
# (possibly there is a slightly more elegant way, but this works - SMB)
ncra -O -w 31,31,28 $OUTDIR/AOD_550_12_climo.nc $OUTDIR/AOD_550_01_climo.nc $OUTDIR/AOD_550_02_climo.nc $OUTDIR/AOD_550_DJF_climo.nc
ncra -O -w 31,30,31 $OUTDIR/AOD_550_03_climo.nc $OUTDIR/AOD_550_04_climo.nc $OUTDIR/AOD_550_05_climo.nc $OUTDIR/AOD_550_MAM_climo.nc
ncra -O -w 30,31,31 $OUTDIR/AOD_550_06_climo.nc $OUTDIR/AOD_550_07_climo.nc $OUTDIR/AOD_550_08_climo.nc $OUTDIR/AOD_550_JJA_climo.nc
ncra -O -w 30,31,30 $OUTDIR/AOD_550_09_climo.nc $OUTDIR/AOD_550_10_climo.nc $OUTDIR/AOD_550_11_climo.nc $OUTDIR/AOD_550_SON_climo.nc
ncra -O -w 90,92,92,91 $OUTDIR/AOD_550_DJF_climo.nc $OUTDIR/AOD_550_MAM_climo.nc $OUTDIR/AOD_550_JJA_climo.nc $OUTDIR/AOD_550_SON_climo.nc  $OUTDIR/AOD_550_ANN_climo.nc

# Create additional backup of original file for safety
#cp processed/${ROOT}_${SEASON}_climo.nc processed/${ROOT}_${SEASON}_climo.nc.org

# Loop over seasons
for SEASON in $WHICHSEASONS
do

OBSFILE="${OUTDIR}/AOD_550_${SEASON}_climo.nc"

# Logic to insert a reasonable "climatological time axis", depending
# on the value of the SEASON variable.
#
# NOTE 1: Please only do this if your file is a climatology, and is missing
# any time information.  Otherwise, please comment out this code or modify
# as appropriate for your case.
#
# NOTE 2: The default values provided below assume that the December value
# is from the preceding year, so that DJF is composed of contiguous months
# This is the preferred approach for calculating DJF climatologies.
#
# NOTE 3: I am hard-coding some reasonable default values here.  There is
# likely a much better way to do this.  Perhaps it could be improved upon
# at a later date. (S.M. Burrows)


# File has climatology_bounds.  Use this instead of time_bnds.
NV="2"
if  [ "$SEASON" == "ANN" ]; then
    TIME="182.2125"           # Midpoint of a Gregorian year
#    TIMEBNDS="{0.0,365.2425}"  # One Gregorian year
    CLMBNDS="{0.0,5297.5149}" # (MD) Is this correct??
elif [ "$SEASON" == "DJF" ]; then
    TIME="321.7095"
    CLMBNDS="{275.1792,5113.3950}"     # Used DJF of 2000-2001 since no 1999 data for DJF 1999-2000
#    CLIMBNDS="2000-12-01","2014-03-01"
elif [ "$SEASON" == "MAM" ]; then
    TIME="46.0300"
    CLMBNDS="{0.0,5205.4549}"
#    CLIMBNDS="2000-03-01","2014-06-01"
elif [ "$SEASON" == "JJA" ]; then
    TIME="137.5895"
    CLMBNDS="{92.0599,5297.5149}"
#    CLMBNDS="2000-06-01","2014-09-01"
elif [ "$SEASON" == "SON" ]; then
    TIME="229.6495"
    CLMBNDS="{184.1199,5023.3365}"
#    CLIMBNDS="2000-09-01","2013-11-30"
elif [ "$SEASON" == "01" ]; then
    TIME="321.7142"
#    CLIMBNDS="2001-01-01","2014-02-01"
    CLMBNDS="{306.2041,5085.3768}"
elif [ "$SEASON" == "02" ]; then
    TIME="351.2334"
#    CLIMBNDS="2001-02-01","2014-03-01"
    CLMBNDS="{337.2243,5113.395}"
elif [ "$SEASON" == "03" ]; then
    TIME="15.5101"
#    CLIMBNDS="2000-03-01","2014-04-01"
    CLMBNDS="{0.0,5144.4152}"
elif [ "$SEASON" == "04" ]; then
    TIME="46.0300"
#    CLIMBNDS="2000-04-01","2014-05-01"
    CLMBNDS="{31.0202,5174.4347}"
elif [ "$SEASON" == "05" ]; then
    TIME="76.5498"
#    CLIMBNDS="2000-05-01","2014-06-01"
    CLMBNDS="{61.0397,5205.4549}"
elif [ "$SEASON" == "06" ]; then
    TIME="107.5700"
#    CLIMBNDS="2000-06-01","2014-07-01"
    CLMBNDS="{92.0599,5235.4744}"
elif [ "$SEASON" == "07" ]; then
    TIME="137.5895"
#    CLIMBNDS="2000-07-01","2014-08-01"
    CLMBNDS="{122.0794,5266.4946}"
elif [ "$SEASON" == "08" ]; then
    TIME="168.6097"
#    CLIMBNDS="2000-08-01","2014-09-01"
    CLMBNDS="{153.0996,5297.5149}"
elif [ "$SEASON" == "09" ]; then
    TIME="198.6292"
#    CLIMBNDS="2000-09-01","2014-10-01"
    CLMBNDS="{184.1198,4962.2919}"
elif [ "$SEASON" == "10" ]; then
    TIME="229.6494"
#    CLIMBNDS="2000-10-01","2013-11-01"
    CLMBNDS="{214.1393,4993.3121}"
elif [ "$SEASON" == "11" ]; then
    TIME="259.6689"
#    CLIMBNDS="2000-11-01","2013-12-01"
    CLMBNDS="{245.1595,5023.3316}"
elif [ "$SEASON" == "12" ]; then
    TIME="290.6891"
#    CLIMBNDS="2000-12-01","2014-01-01"
    CLMBNDS="{275.1792,5054.3518}"
else
    continue
fi

# Print informative message
echo "Setting time axis for season $SEASON"
echo "Time: $TIME"
#echo "Time bounds: $TIMEBNDS"
echo "Climatology bounds: $CLMBNDS"

# Set time axis values
ncap2 -O -h -s "time[time]=$TIME" $OBSFILE $OBSFILE

# Define bounds attribute
ncatted -O -a bounds,time,o,c,"time_bnds" -h $OBSFILE $OBSFILE

# Define a dimension with length 2 for use in defining the bounds attributes
ncap2 -h -O -s 'defdim("nv",2)' $OBSFILE $OBSFILE

# Set time_bnds values
# (MD) Set to -999.9 for missing value since time bounds are not in use
ncap2 -O -h -s "time_bnds[time,nv]={-999.9,-999.9}" $OBSFILE $OBSFILE # Swiched for climate bounds below

# Add lat_bnds/lon_bnds variables
#
# Note: This was particularly complicated to figure out in shell script,
# so it was decided to do this in Python and call it from here. If there
# is another (better) way to do this, please adjust this script accordingly. (MD)
#ncap2 -O -S latlon.nco $OBSFILE $OBSFILE

# Set _FillValue ("missing" value) for time and time_bnds
ncatted -O -h -a _FillValue,time,o,d,$FILLVAL $OBSFILE $OBSFILE

#ncatted -O -h -a _FillValue,time_bnds,o,f,$FILLVAL $OBSFILE $OBSFILE

# Example: add climatology_bounds attribute if file is a climatology.
#    Note that all climatology files should have this attribute, which
#    specifies the bounds over which the climatological means were calculated.
#
ncatted -O -h -a units,time,o,c,"days since 2000-03-01" $OBSFILE $OBSFILE

# Remove climatology_bounds
#
# Note: This was creating a problem with assigning new climatology bounds, so we are making our own.
ncks -O -x -v climatology_bounds $OBSFILE tmp.nc
mv tmp.nc $OBSFILE

# This data set has this attribute/variable
ncatted -O -a climatology,time,o,c,"climatology_bnds" -h $OBSFILE $OBSFILE
ncap2 -O -s "climatology_bnds[time,nv]=${CLMBNDS}" $OBSFILE $OBSFILE

# Set values to "missing" if not known
# (MD) Using climate bounds and adjusted time bounds accordingly
ncatted -O -h -a _FillValue,climatology_bnds,o,d,-999.9 $OBSFILE $OBSFILE

# Add the cell methods to each flux
# There is probably a lot better way to do this that I am not aware of (MD)
ncatted -O -h -a cell_methods,AOD_550,o,c,"time: mean within years time: mean over years" $OBSFILE $OBSFILE
ncatted -O -h -a cell_methods,AOD_550_ann,o,c,"time: mean over years" $OBSFILE $OBSFILE
# Climatology files should contain a "season" attribute, which should
# be either a three-letter season (e.g., ANN, DJF, MAR) or a two-digit
# month (e.g. 01, 02, ... 12)/mo
ncatted -O -h -a season,global,o,c,$SEASON $OBSFILE

# Variable cleanup
ncatted -O -a standard_name,lat,o,c,"latitude" $OBSFILE $OBSFILE
ncatted -O -a standard_name,lon,o,c,"longitude" $OBSFILE $OBSFILE
ncatted -O -a long_name,time,o,c,"time" -a standard_name,time,o,c,"time" $OBSFILE $OBSFILE
# Flip FLNS/FLNSC signs to match model output
#
# Note: Simply reversing valid_min/valid_max was
# throwing an error in ncview. By doing the below
# operation, observed output now matches model
# output more closely.
#ncap2 -O -s 'FLNS[\$time,\$lat,\$lon]*=-1' $OBSFILE $OBSFILE
#ncap2 -O -s 'FLNSC[\$time,\$lat,\$lon]*=-1' $OBSFILE $OBSFILE

#  [contact] : name and contact information (e.g., email, address,
#  phone number) of person(s) who should be contacted for more
#  information about the data.  This may/should include both the
#  original author of the dataset (if known) and an ACME POC who
#  processed or is familiar with the data.
ncatted -O -h -a contact,global,o,c,"Processed for ACME by:\nSusannah M. Burrows\n<susannah.burrows@pnnl.gov>\nAtmospheric Science and Global Change Division\nPacific Northwest National Laboratory\ntel: (509) 372-6183\nPo-Lun Ma\n<Po-Lun.Ma@pnnl.gov>\ntel:(509) 372-6936" $OBSFILE $OBSFILE

# [date] : or "date_created" but also all the changes should have date
# recorded).  Please use an unambiguous date format, e.g., "1 Apr
# 2015", not "2015-01-04" or "2015-04-01".
ncatted -O -h -a date,global,o,c,"21 Aug 2015" $OBSFILE $OBSFILE

# [version] : Model or dataset version number, if applicable.
ncatted -O -h -a version,global,o,c,"r1" $OBSFILE $OBSFILE

# [comment] : any additional comments or miscellaneous information
# about the data or the methods used to produce it.
ncatted -O -a comment,global,o,c,"Updated to improve compliance with CF and ACME conventions.\n   SMB, 21 Aug. 2015" -h $OBSFILE $OBSFILE

# [script_git-hash] : Hash that uniquely identifies the file
#$HOME/PreAndPostProcessingScripts/utils/add_git_hash_to_netcdf_metadata $OBSFILE


done
