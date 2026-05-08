#!/bin/bash

# this script starts with the 0.5x0.5 CMIP6 CEDS files for historical CO2 emissions
# the final 12 months of 2014 are extracted and processed
# the output is three text files that contain:
#    surface co2 emissions (excluding international shipping) at the surface
#    aircraft co2 emissions for the defined aggregate levels (currently two) 
#    international shipping co2 emissions at the surface
# the text files are csv tables with date,lon,lat,value
#    values are in co2_kg/m^2/s
#    lat,lon are in degrees, cell center
# the aircraft text file date includes a level label (yyyymmll)
#    this starts at zero for the lowest level

# the 25 input aircraft levels are in order starting from the surface
#    their defined heights are in the file
# define the aggregation of levels here so that changing it for the model simply means updating this file
#    the ehc code accepts up to two levels - so make that the default
#    the ehc can then sum the two levels if only one level is passed to eam

# the defined aircraft output levels are:
#    lo: 15 bottom levels including 8.845km
#    hi: 10 upper levels from 9.455 to 14.945
# note that GCAM outputs data up to ~11km, ?which is distributed to these 25 levels by ceds? 
#    11.285 is the highest large emission level
#    9.455 is where long-range emissions become noticable in the plots

# there are eight sectors in the other co2 file
# international shipping is id 7
# sector:ids = "0: Agriculture; 1: Energy; 2: Industrial; 3: Transportation; 4: Residential, Commercial, Other; 5: Solvents production and application; 6: Waste; 7: International Shipping"

# the final netcdf files that the text is extracted from are retained

date

# needed modules
#module load intel
#module load nco

# this gets what is needed also
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# eventually create options for the other three resolutions if needed and make arguments to the script

# some useful bash functions

# ncdmnsz $dmn_nm $fl_nm : What is dimension size?
function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }
# ncmax $var_nm $fl_nm : What is maximum of variable?
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
# ncmin $var_nm $fl_nm : What is minimum of variable?
function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }

proc_dir='/compyfs/inputdata/iac/giac/gcam/gcam_6_0/data/emission_processing/'

# set the desired output resolution here
# the original data are at 0.5x0.5 res with origin at -180,-90; and corners aligned with these limits 
# this is also how the nomask scrip file is defined for 0.5x0.5
map_file="/compyfs/inputdata/lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_nomask_to_0.9x1.25_nomask_aave_da_c121019.nc"

# output nc file names - make sure they match the map file out resolution
sfc_file_out=${proc_dir}'CO2-em-SFC-anthro_0.9x1.25_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'
air_file_out=${proc_dir}'CO2-em-AIR-2lvl-anthro_0.9x1.25_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'
ship_file_out=${proc_dir}'CO2-em-SHIP-anthro_0.9x1.25_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'

# final output text file names
sfc_text_out=${proc_dir}'CO2-em-SFC-anthro_0.9x1.25_input4MIPs_2014.csv'
air_text_out=${proc_dir}'CO2-em-AIR-2lvl-anthro_0.9x1.25_input4MIPs_2014.csv'
ship_text_out=${proc_dir}'CO2-em-SHIP-anthro_0.9x1.25_input4MIPs_2014.csv'

# source files
air_file="/compyfs/inputdata/atm/cam/ggas/CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-2017-08-30_gn_200001-201412.nc"
other_file="/compyfs/inputdata/atm/cam/ggas/CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_200001-201412.nc"

# intermediate files
air_file_1y=${proc_dir}'CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-2017-08-30_gn_2014.nc'
air_file_1y_agg=${proc_dir}'CO2-em-AIR-2lvl-anthro_input4MIPs_emissions_CMIP_CEDS-2017-08-30_gn_2014.nc'
other_file_1y=${proc_dir}'CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'
ship_file_1y=${proc_dir}'CO2-em-SHIP-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'
sfc_file_1y=${proc_dir}'CO2-em-SFC-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2014.nc'
air_file_lo=${proc_dir}'CO2-em-AIR-LO-anthro_input4MIPs_emissions_CMIP_CEDS-2017-08-30_gn_2014.nc'
air_file_lo_level=${proc_dir}'CO2-em-AIR-LO-level-temp.nc'
air_file_hi=${proc_dir}'CO2-em-AIR-HI-anthro_input4MIPs_emissions_CMIP_CEDS-2017-08-30_gn_2014.nc'
air_file_hi_level=${proc_dir}'CO2-em-AIR-HI-level-temp.nc'


# extract year 2014 from the original files, the last 12 months (nces?)
ncea -O -F -d time,169,180 ${air_file} ${air_file_1y}
ncea -O -F -d time,169,180 ${other_file} ${other_file_1y}

# aircraft level processing

# aggregate the aircraft data to two levels
ncwa -O -N -b -v CO2_em_AIR_anthro -a level -d level,0,14 ${air_file_1y} ${air_file_lo}
ncwa -O -N -b -v CO2_em_AIR_anthro -a level -d level,15,24 ${air_file_1y} ${air_file_hi}
# Make level record dimension - must put record dimension first for ncrcat to work properly
ncks -O --fix_rec_dmn time ${air_file_lo} ${air_file_lo_level}
ncks -O --mk_rec_dmn level ${air_file_lo_level} ${air_file_lo_level}
ncpdq -O -a level,time,lat,lon ${air_file_lo_level} ${air_file_lo_level}
# Make level record dimension 
ncks -O --fix_rec_dmn time ${air_file_hi} ${air_file_hi_level}
ncks -O --mk_rec_dmn level ${air_file_hi_level} ${air_file_hi_level}
ncpdq -O -a level,time,lat,lon ${air_file_hi_level} ${air_file_hi_level}
# concatenate along level
ncrcat -O ${air_file_lo_level} ${air_file_hi_level} ${air_file_1y_agg}
# revert time to record dimension
ncks -O --fix_rec_dmn level ${air_file_1y_agg} ${air_file_1y_agg}
ncks -O --mk_rec_dmn time ${air_file_1y_agg} ${air_file_1y_agg}
ncpdq -O -a time,level,lat,lon ${air_file_1y_agg} ${air_file_1y_agg}

rm ${air_file_lo}
rm ${air_file_hi}
rm ${air_file_lo_level}
rm ${air_file_hi_level}

# extract surface shipping data
ncea -O -d sector,7,7 ${other_file_1y} ${ship_file_1y}

# extract non-shipping surface data summed across remaining sectors
ncwa -O -N -v CO2_em_anthro -a sector -d sector,0,6 ${other_file_1y} ${sfc_file_1y}

# remap the data to the desired grid
ncremap -m ${map_file} ${sfc_file_1y} ${sfc_file_out}
ncremap -m ${map_file} ${air_file_1y_agg} ${air_file_out}
ncremap -m ${map_file} ${ship_file_1y} ${ship_file_out}

# convert time variable from days since 1750-01-01 to month in yyyymm format
# the time value is the day in the middle of the month, so the simple math below works
# only need to do this once
ncap2 -O -s "time=int(trunc(float(time/365.0)+1750.0)*100 + trunc(((float(time/365.0)+1750.0)-trunc(float(time/365.0)+1750.0))*365.0/30.0+1.0))" ${sfc_file_out} ${sfc_file_out}

# output text files
# need nested loop for this to work properly

num_time=$(ncdmnsz time ${sfc_file_out})
num_lat=$(ncdmnsz lat ${sfc_file_out})
num_lon=$(ncdmnsz lon ${sfc_file_out})

# write separate files then paste them together
# they need to be line by line and have the same length
# so still need the time-lat loops to create the full length time, lat, and lon files
# remove the blank lines from files
# the netcdf time, lat, and lon are the same for all files

tin=`ncks -C -H -v time -s "%i " ${sfc_file_out}`
ain=`ncks -C -H -v lat -s "%f " ${sfc_file_out}`
ncks -O -C -H -v lon -s "%f\n" ${sfc_file_out} > ${proc_dir}'temp1_lon.txt'
tr -s '\n' < ${proc_dir}'temp1_lon.txt' > ${proc_dir}'temp_lon.txt'
rm ${proc_dir}'temp1_lon.txt'

ncks -O -C -H -v CO2_em_anthro -s "%22.20f\n" ${sfc_file_out} > ${proc_dir}'temp_value.txt' 
tr -s '\n' < ${proc_dir}'temp_value.txt' > ${proc_dir}'sfc_value.txt'

ncks -O -C -H -v CO2_em_anthro -s "%22.20f\n" ${ship_file_out} > ${proc_dir}'temp_value.txt' 
tr -s '\n' < ${proc_dir}'temp_value.txt' > ${proc_dir}'ship_value.txt'

# lo air level
ncks -O -d level,0,0 -C -H -v CO2_em_AIR_anthro -s "%22.20f\n" ${air_file_out} > ${proc_dir}'temp_value.txt' 
tr -s '\n' < ${proc_dir}'temp_value.txt' > ${proc_dir}'air_lo_value.txt'

# hi air level
ncks -O -d level,1,1 -C -H -v CO2_em_AIR_anthro -s "%22.20f\n" ${air_file_out} > ${proc_dir}'temp_value.txt' 
tr -s '\n' < ${proc_dir}'temp_value.txt' > ${proc_dir}'air_hi_value.txt'

# this is the spacer for the time and lat variables
IFS=" "
tz=( $tin )
az=( $ain )

# make sure that these are new, clean files
> time.txt
> lat.txt
> lon.txt

for ((t=0 ; t<$num_time ; t++));
do

   for ((a=0 ; a<$num_lat ; a++));
   do

      # write repeated appropriate time and lat and lon
      yes ${tz[$t]} | head -n ${num_lon} >> ${proc_dir}'time.txt'
      # the -- tells that the next argument is an input, not an option, so that negative numbers can be used
      yes -- ${az[$a]} | head -n ${num_lon} >> ${proc_dir}'lat.txt'
      dd if=${proc_dir}'temp_lon.txt' bs=1M status=none >> ${proc_dir}'lon.txt'

   done
   
done

wc -l ${proc_dir}'sfc_value.txt'
wc -l ${proc_dir}'ship_value.txt'
wc -l ${proc_dir}'air_lo_value.txt'
wc -l ${proc_dir}'air_hi_value.txt'
wc -l ${proc_dir}'time.txt'
wc -l ${proc_dir}'lat.txt'
wc -l ${proc_dir}'lon.txt'

# surface csv file

# no quotes enters the text exactly
echo yyyymm,lon_deg,lat_deg,co2_kg/m2/s > ${sfc_text_out}
# - represents piped input
paste -d "," ${proc_dir}'time.txt' ${proc_dir}'lon.txt' | paste -d "," - ${proc_dir}'lat.txt' | paste -d "," - ${proc_dir}'sfc_value.txt' >> ${sfc_text_out}
$(echo wc -l ${sfc_text_out})


# ship csv file

echo yyyymm,lon_deg,lat_deg,co2_kg/m2/s > ${ship_text_out}
paste -d "," ${proc_dir}'time.txt' ${proc_dir}'lon.txt' | paste -d "," - ${proc_dir}'lat.txt' | paste -d "," - ${proc_dir}'ship_value.txt' >> ${ship_text_out}
$(echo wc -l ${ship_text_out})


# aircraft csv file

# write the low level, then the high level
# add a level tag to the first column for lo (00) and hi (01)
echo yyyymmll,lon_deg,lat_deg,co2_kg/m2/s > ${air_text_out}

sed "s/$/00/" ${proc_dir}'time.txt' > ${proc_dir}'temp_time.txt'
paste -d "," ${proc_dir}'temp_time.txt' ${proc_dir}'lon.txt' | paste -d "," - ${proc_dir}'lat.txt' | paste -d "," - ${proc_dir}'air_lo_value.txt' >> ${air_text_out}

sed "s/$/01/" ${proc_dir}'time.txt' > ${proc_dir}'temp_time.txt'
paste -d "," ${proc_dir}'temp_time.txt' ${proc_dir}'lon.txt' | paste -d "," - ${proc_dir}'lat.txt' | paste -d "," - ${proc_dir}'air_hi_value.txt' >> ${air_text_out}

$(echo wc -l ${air_text_out})

rm ${proc_dir}'time.txt'
rm ${proc_dir}'temp_time.txt'
rm ${proc_dir}'lat.txt'
rm ${proc_dir}'lon.txt'
rm ${proc_dir}'temp_lon.txt'
rm ${proc_dir}'sfc_value.txt'
rm ${proc_dir}'ship_value.txt'
rm ${proc_dir}'air_lo_value.txt'
rm ${proc_dir}'air_hi_value.txt'
rm ${proc_dir}'temp_value.txt'

date
