#!/bin/csh -f

#*****************************************************************
# Check to see of 12 months of a year already exist, then
# get the monthly files from Mass Storage System if needed
#*****************************************************************

# This file reads in files from MSS 
# $DATE_FORMAT  form of date in history file name (eg. yyyy-mm), input
# $read_dir     case name of file to read   , input
# $BEG_READ     first year of data to read  , input
# $END_READ     last year of data to read   , input
# $FILE_HEADER  beginning of filename
# $PATH_MSS     directory on MSS where data resides
# $WKDIR      directory on dataproc where data will be put

if ($#argv != 5) then
  echo "usage: read_from_mss.csh $DATE_FORMAT $read_dir $BEG_READ $END_READ $month"
  exit
endif

set DATE_FORMAT = $1
set read_dir = $2
@ BEG_READ = $3
@ END_READ = $4
set month = $5

# Set msread password for b20.003 case
if ( $read_dir == b20.003 ) then
  set msspwd = '-rpwd ccsm1330'
else
  set msspwd = ' ' 
endif

echo GETTING MONTHLY FILES FROM THE MSS
echo THIS MIGHT TAKE SOME TIME
echo ' '

if ($BEG_READ < 1) then    # so we don't get a negative number
  echo ERROR: FIRST YEAR OF TEST DATA $BEG_READ MUST BE GT ZERO
  exit
endif

@ IYEAR = $BEG_READ
#-------------------------------------------------------
# Loop through years
#-------------------------------------------------------

while ($IYEAR <= $END_READ)

    set four_digit_year = `printf "%04d" {$IYEAR}`

#-------------------------------------------------------------
# Read a month of data from MSS
#-------------------------------------------------------------

    echo  'GETTING '{$PATH_MSS}{$FILE_HEADER}${four_digit_year}'*.nc'
    if (`which msrcp | wc -w` == 1) then
       msls {$PATH_MSS}{$FILE_HEADER}${four_digit_year}'-01.nc'
       if ( $status == 0 ) setenv TAR 0
    else
       hsi -q 'ls {$PATH_MSS}{$FILE_HEADER}${four_digit_year}-01.nc'
       if ( $status == 0 ) setenv TAR 0
    endif
    if (`which msrcp | wc -w` == 1) then
       msls {$PATH_MSS}{$FILE_HEADER}${four_digit_year}'.tar'
       if ( $status == 0 ) setenv TAR 1
    else
       hsi -q 'ls {$PATH_MSS}{$FILE_HEADER}${four_digit_year}.tar'
       if ( $status == 0 ) setenv TAR 1
    endif
    if ( $TAR == 0 ) then
       if (`which msrcp | wc -w` == 1) then
          msrcp $msspwd 'mss:'{$PATH_MSS}{$FILE_HEADER}${four_digit_year}'-'${month}'.nc' $WKDIR
       else
          pushd $WKDIR
          hsi -q 'cd {$PATH_MSS}; get {$FILE_HEADER}${four_digit_year}-${month}.nc'
          popd
       endif
    else
       pushd $WKDIR
       if (`which msrcp | wc -w` == 1) then
          msrcp $msspwd 'mss:'{$PATH_MSS}{$FILE_HEADER}${four_digit_year}'.tar' $WKDIR
       else
          hsi -q 'cd {$PATH_MSS}; get {$FILE_HEADER}${four_digit_year}.tar'
       endif
       tar -xvf {$FILE_HEADER}${four_digit_year}.tar {$FILE_HEADER}${four_digit_year}-${month}.nc
       /bin/rm {$FILE_HEADER}${four_digit_year}.tar
       popd
    endif

  @ IYEAR++                             # advance year
end        # End of IYEAR <= END_READ

echo MONTHLY FILES COPIED FROM THE MSS TO {$WKDIR}
echo ' '

end
