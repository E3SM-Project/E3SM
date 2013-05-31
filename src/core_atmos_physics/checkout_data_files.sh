#!/usr/bin/env sh

################################################################################
# File: checkout_data_files.sh
#
# The purpose of this script is to obtain lookup tables used by the WRF physics
#   packages. At present, the only method for acquiring these tables is through
#   the MPAS-Dev github repository using either git or svn.
#
# If none of the methods used in this script are successful in acquiring the 
#   tables, please contact the MPAS-A developers.
#
################################################################################

if [ -s physics_wrf/files/CAM_ABS_DATA.DBL ]; then
   echo "*** WRF physics tables appear to already exist; no need to obtain them again ***"
   exit 0
fi


#
# Try using 'git'
#
which git
if [ $? == 0 ]; then
   echo "*** trying git to obtain WRF physics tables ***"
   git clone git://github.com/MPAS-Dev/MPAS-Data.git
   if [ $? == 0 ]; then
      mv MPAS-Data/atmosphere/physics_wrf/files physics_wrf/
      rm -rf MPAS-Data
      exit 0
   else
      echo "*** failed to obtain WRF physics tables using git ***"
   fi
else
   echo "*** git not in path ***"
fi


#
# Try using 'svn'
#
which svn
if [ $? == 0 ]; then
   echo "*** trying svn to obtain WRF physics tables ***"
   svn checkout --non-interactive --trust-server-cert https://github.com/MPAS-Dev/MPAS-Data.git
   if [ $? == 0 ]; then
      mv MPAS-Data.git/trunk/atmosphere/physics_wrf/files physics_wrf/
      rm -rf MPAS-Data.git
      exit 0
   else
      echo "*** failed to obtain WRF physics tables using svn ***"
   fi
else
   echo "*** svn not in path ***"
fi


echo "***************************************************"
echo "Unable to obtain WRF physics tables by any means"
echo "***************************************************"

exit 1
