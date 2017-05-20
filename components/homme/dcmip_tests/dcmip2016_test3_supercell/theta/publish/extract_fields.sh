#!/bin/bash
#
#  copnvert data to format required for dcmip2016 intercomparison
#

ncks -v time,u   ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.U.nc
ncks -v time,v   ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.V.nc
ncks -v time,w   ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.W.nc
ncks -v time,T   ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.T.nc
ncks -v time,Th  ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.Theta.nc
ncks -v time,Q   ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.Q.nc
ncks -v time,Q2  ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.Qc.nc
ncks -v time,Q3  ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.Qr.nc
ncks -v time,geo ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.Geopotential.nc
ncks -v time,precl ../movies/dcmip2016_test3_r400.nc acme.163.r400.L40.interp_latlon.nonhydro.PRECT.nc

ncks -v time,u   ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.U.nc
ncks -v time,v   ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.V.nc
ncks -v time,w   ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.W.nc
ncks -v time,T   ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.T.nc
ncks -v time,Th  ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.Theta.nc
ncks -v time,Q   ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.Q.nc
ncks -v time,Q2  ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.Qc.nc
ncks -v time,Q3  ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.Qr.nc
ncks -v time,geo ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.Geopotential.nc
ncks -v time,precl ../movies/dcmip2016_test3_r200.nc acme.163.r200.L40.interp_latlon.nonhydro.PRECT.nc

ncks -v time,u   ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.U.nc
ncks -v time,v   ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.V.nc
ncks -v time,w   ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.W.nc
ncks -v time,T   ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.T.nc
ncks -v time,Th  ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.Theta.nc
ncks -v time,Q   ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.Q.nc
ncks -v time,Q2  ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.Qc.nc
ncks -v time,Q3  ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.Qr.nc
ncks -v time,geo ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.Geopotential.nc
ncks -v time,precl ../movies/dcmip2016_test3_r100.nc acme.163.r100.L40.interp_latlon.nonhydro.PRECT.nc

ncks -v time,u   ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.U.nc
ncks -v time,v   ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.V.nc
ncks -v time,w   ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.W.nc
ncks -v time,T   ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.T.nc
ncks -v time,Th  ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.Theta.nc
ncks -v time,Q   ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.Q.nc
ncks -v time,Q2  ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.Qc.nc
ncks -v time,Q3  ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.Qr.nc
ncks -v time,geo ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.Geopotential.nc
ncks -v time,precl ../movies/dcmip2016_test3_r50.nc acme.163.r50.L40.interp_latlon.nonhydro.PRECT.nc
