#!/bin/bash

# This script is for generating surface mass balance forcing in time.
# The time-variant perturbation is added on top of a background smb field
# input: ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc, smb_anomaly_04km.nc (from initMIP FTP site ftp://cryoftp1.gsfc.nasa.gov/), output.nc that contains the xtime field
# output: ais20km_smb_background_anomaly.nc

ncks -x -v sfcMassBal ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc ais20km_withoutsmb.nc
cp ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc ais20km_withsmb.nc

#ncrename -v x,x1 smb_anomaly_04km.nc
#ncrename -v y,y1 smb_anomaly_04km.nc
#ncrename -v asmb,acab smb_anomaly_04km.nc

# rename the coordinate and var name so that we can use the script interpolate_to_mpasli_grid.py for interpolation

~/Apps/MPAS-Tools/grid_gen/landice_grid_tools/interpolate_to_mpasli_grid.py -s smb_anomaly_04km.nc -d ais20km_withsmb.nc -m d
#change this script path as your own
mv ais20km_withsmb.nc ais20km_withsmbanomaly.nc


ncks -v sfcMassBal ais20km_withsmbanomaly.nc smb_anomaly.nc
ncks -v sfcMassBal ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc smb_Steve.nc

ncks -O -v xtime -d Time,0 output.nc xtime.nc

for i in {1..100..1}
do

    m=$(($i+99))
    m1=$(($i+199))

    if [ $i -le 40 ]
    then
        n=$(echo "scale=2; $i/40.0"|bc)
    else
        n=1
    fi


    ncap2 -O -s 'xtime="0000-01-01_00:00:00                                            "' xtime.nc ${m1}.nc

    ncap2 -O -s "sfcMassBal=sfcMassBal*$n" smb_anomaly.nc tmp.nc
    ncra -y ttl -v sfcMassBal tmp.nc smb_Steve.nc ${m}.nc

done

ncrcat -n 100,3,1 100.nc smb_all.nc
# combine the smb data for all years into a single file, smb_all.nc
ncrcat -n 100,3,1 200.nc xtime_all.nc
# combine the xtime for all years into a sinble file, xtime_all.nc
# note that here xtime is the same for each year. we may change its values using the script process_xtime.py
ncks -A -v sfcMassBal smb_all.nc ais20km_withoutsmb.nc

mv ais20km_withoutsmb.nc ais20km_smb_background_anomaly.nc
ncks -A -v xtime xtime_all.nc ais20km_smb_background_anomaly.nc

rm 1*.nc
rm 2*.nc
# remove temporary files

