#!/bin/bash

# A Bash script to post-process E3SM 6 hourly (h2) instantaneous output to generate a text file storing Tropical Cyclone tracks
# TO-DO need stable tempestremap and tempestextremes MPI or non-MPI builds for executing this script. Current using builds from Paul Ullrich's home directory on Cori. Request an interactive session as follows: salloc --nodes=1 --partition=debug --time=30:00 -C haswell
# Note: The results directory needs to be on $SCRATCH on Cori for this script at this moment.

start="1880"
end="1881"
caseid="20180215.DECKv1b_H1.ne30_oEC.edison"

drc_in=/global/cfs/cdirs/e3smpub/E3SM_simulations/${caseid}/archive/atm/hist
result_dir=/global/cscratch1/sd/chengzhu/tc-analysis-test/
mkdir -p $result_dir
file_name=${caseid}_${start}_${end}

#res = 30/120 for ne30/ne120 grids
res=30

# Generate mesh files (.g).
/global/homes/p/paullric/tempestremap/bin/GenerateCSMesh --res $res --alt --file ${result_dir}outCSne$res.g
# Generate connectivity files (.dat)
/global/homes/p/paullric/tempestextremes/bin/GenerateConnectivityFile --in_mesh ${result_dir}outCSne$res.g --out_type CGLL --out_connect ${result_dir}connect_CSne${res}_v2.dat

# Get the list of files
cd ${drc_in};eval ls ${drc_in}/${caseid}.cam.h2.*{${start}..${end}}*.nc >${result_dir}inputfile_${file_name}.txt
cd ${drc_in};eval ls ${caseid}.cam.h2.*{${start}..${end}}*.nc >${result_dir}outputfile_${file_name}.txt

cd ${result_dir}
# Detection threshold including:
# The sea-level pressure (SLP) must be a local minimum; SLP must have a sufficient decrease (300 Pa) compared to surrounding nodes within 4 degree radius; The average of the 200 hPa and 500 hPa level temperature decreases by 0.6 K in all directions within a distance of 4Â° from the location o fSLP minima
/global/homes/p/paullric/tempestextremes/bin/DetectNodes --verbosity 0 --in_connect ${result_dir}connect_CSne${res}_v2.dat --closedcontourcmd "PSL,300.0,4.0,0;_AVG(T200,T500),-0.6,4,0.30" --mergedist 6.0 --searchbymin PSL --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2" --timestride 1 --in_data_list ${result_dir}inputfile_${file_name}.txt --out ${result_dir}out.dat

cat ${result_dir}out.dat0* > ${result_dir}cyclones_${file_name}.txt

# Stitch all candidate nodes in time to form tracks, with a maximum distance between candidates of 6.0, minimum time steps of 6, and with a maximum gap size of one (most consecutive time steps with no associated candidate). And there is threshold for wind speed, lat and lon.
/global/homes/p/paullric/tempestextremes/bin/StitchNodes --in_fmt "lon,lat,slp,wind" --in_connect ${result_dir}connect_CSne${res}_v2.dat --range 6.0 --mintime 6 --maxgap 1 --in ${result_dir}cyclones_${file_name}.txt --out ${result_dir}cyclones_stitch_${file_name}.dat --threshold "wind,>=,17.5,6;lat,<=,40.0,6;lat,>=,-40.0,6"
rm ${result_dir}cyclones_${file_name}.txt

# Generate histogram of detections
/global/homes/p/paullric/tempestextremes/bin/HistogramNodes --in ${result_dir}cyclones_stitch_${file_name}.txt --iloncol 2 --ilatcol 3 --out ${result_dir}cyclones_hist_${file_name}.nc

# Calculate relative vorticity
sed -i 's/.nc/_vorticity.nc/' ${result_dir}outputfile_${file_name}.txt
/global/homes/p/paullric/tempestextremes/bin/VariableProcessor --in_data_list ${result_dir}inputfile_${file_name}.txt --out_data_list ${result_dir}outputfile_${file_name}.txt --var "_CURL{4,0.5}(U850,V850)" --varout "VORT" --in_connect ${result_dir}connect_CSne${res}_v2.dat

/global/homes/p/paullric/tempestextremes/bin/DetectNodes --verbosity 0 --in_connect ${result_dir}connect_CSne${res}_v2.dat --closedcontourcmd "VORT,-5.e-6,4,0" --mergedist 2.0 --searchbymax VORT --outputcmd "VORT,max,0" --in_data_list ${result_dir}outputfile_${file_name}.txt --out ${result_dir}aew_out.dat --minlat -35.0 --maxlat 35.0
cat ${result_dir}aew_out.dat0* > ${result_dir}aew_${file_name}.txt

/global/homes/p/paullric/tempestextremes/bin/StitchNodes --in_fmt "lon,lat,VORT" --in_connect ${result_dir}connect_CSne${res}_v2.dat --range 3.0 --minlength 8 --maxgap 0 --min_endpoint_dist 10.0 --in ${result_dir}aew_${file_name}.txt --out ${result_dir}aew_stitch_5e-6_${file_name}.dat --threshold "lat,<=,25.0,8;lat,>=,0.0,8"
rm ${result_dir}aew_${file_name}.txt

/global/homes/p/paullric/tempestextremes/bin/HistogramNodes --in ${result_dir}aew_stitch_5e-6_${file_name}.dat --iloncol 2 --ilatcol 3 --nlat 256 --nlon 512 --out ${result_dir}aew_hist_${file_name}.nc
rm ${result_dir}*out.dat00*.dat
rm ${result_dir}${caseid}*.nc
rm ${result_dir}*.txt

exit
