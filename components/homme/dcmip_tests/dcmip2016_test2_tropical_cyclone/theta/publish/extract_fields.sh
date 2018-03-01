#!/bin/bash
#
#  copnvert data to format required for dcmip2016 intercomparison
#

#runtype="nonhydro"
runtype="hydro"

ncks -v time,ps     ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.PS.nc
ncks -v time,u      ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.U.nc
ncks -v time,v      ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.V.nc
ncks -v time,w      ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.W.nc
ncks -v time,omega  ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Omega.nc
ncks -v time,T      ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.T.nc
ncks -v time,Th     ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Theta.nc
ncks -v time,Q      ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Q.nc
ncks -v time,Q2     ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Qc.nc
ncks -v time,Q3     ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Qr.nc
ncks -v time,geo    ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.Geopotential.nc
ncks -v time,precl  ../movies/r50-prec0-pbl0-dcmip2016_test21.nc   acme.162-rjpbl.r50.L30.interp_latlon.${runtype}.PRECT.nc

ncks -v time,ps     ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.PS.nc
ncks -v time,u      ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.U.nc
ncks -v time,v      ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.V.nc
ncks -v time,w      ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.W.nc
ncks -v time,omega  ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Omega.nc
ncks -v time,T      ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.T.nc
ncks -v time,Th     ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Theta.nc
ncks -v time,Q      ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Q.nc
ncks -v time,Q2     ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Qc.nc
ncks -v time,Q3     ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Qr.nc
ncks -v time,geo    ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.Geopotential.nc
ncks -v time,precl  ../movies/r50-prec0-pbl1-dcmip2016_test21.nc   acme.162-bryanpbl.r50.L30.interp_latlon.${runtype}.PRECT.nc

