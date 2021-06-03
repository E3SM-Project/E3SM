#!/bin/bash

# This script is for generating basal melt forcing in time.
# The time-variant perturbation is added on top of a background basal melt field
# input: ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc, output.nc that contains the xtime field
# output: ais20km_bmb_background_anomaly.nc

ncks -O -x -v floatingBasalMassBal ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc ais20km_withoutbmb.nc
ncks -O -v floatingBasalMassBal ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc bmb_anomaly.nc
ncks -O -v floatingBasalMassBal RignotBasalMelt.nc bmb_Rignot.nc

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

    ncap2 -O -s "floatingBasalMassBal=floatingBasalMassBal*$n" bmb_anomaly.nc tmp.nc
    # change the basal melt value for each year by b = b * yr/40

    #ncap2 -O -s "floatingBasalMassBal=floatingBasalMassBal*$n" bmb_anomaly.nc ${m}.nc 
    # if we don't want to add the perturbation on top of the background field, use this line
    ncra -O -y ttl -v floatingBasalMassBal tmp.nc bmb_Rignot.nc ${m}.nc
    # otherwise use this line

done

ncrcat -n 100,3,1 100.nc bmb_all.nc
# combine the bmb data for all years into a single file, bmb_all.nc
ncrcat -n 100,3,1 200.nc xtime_all.nc
# combine the xtime for all years into a sinble file, xtime_all.nc
# note that here xtime is the same for each year. we may change its values using the script process_xtime.py

ncap2 -s "floatingBasalMassBal=floatingBasalMassBal*(-1)" bmb_all.nc ais_bmb.nc
# let the value be negative
ncks -A -v floatingBasalMassBal ais_bmb.nc ais20km_withoutbmb.nc

mv ais20km_withoutbmb.nc ais20km_Rignot_anomaly.nc
mv ais20km_Rignot_anomaly.nc ais20km_bmb_background_anomaly.nc

ncks -A -v xtime xtime_all.nc ais20km_bmb_background_anomaly.nc

rm 1*.nc
rm 2*.nc
# remove temporary files

