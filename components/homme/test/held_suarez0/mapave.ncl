load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  average HS output in time 
;  1deg, takes ~10min
;

begin


; NETCDF files to process
fnames =  systemfunc( \
  "ls movies/held_suarez01.nc;"  \
)
fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
print("files = "+fnames)
fcnt = dimsizes(fnames)
print ("number of uv input files = " + fcnt )


fuv=fuv_all[0];    file with u, v
nlev = dimsizes(fuv->lev)
nlat = dimsizes(fuv->lat)
nlon = dimsizes(fuv->lon)
ntimes = dimsizes(fuv->time)
print (" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon+" ntimes="+ntimes)


name = "hsmapave.nc"
system("/bin/rm "+name)           ; remove any pre-existing file
ncdf     = addfile(name ,"c")  ; open output netCDF file

ncdf->zeta = dim_avg_n_Wrap(fuv->zeta,0)
ncdf->div = dim_avg_n_Wrap(fuv->div,0)
;ncdf->omega = dim_avg_n_Wrap(fuv->omega,0)

ncdf->gw = fuv->gw
ncdf->ilev = fuv->ilev
ncdf->hyam = fuv->hyam
ncdf->hyai = fuv->hyai
ncdf->hybm = fuv->hybm
ncdf->hybi = fuv->hybi

ncdf@ne = fuv@ne

end

