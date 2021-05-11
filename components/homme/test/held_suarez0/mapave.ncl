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

if (isfilevar(fuv,"omega")) then
  ;ncdf->omega = dim_avg_n_Wrap(fuv->omega,0)
  omegap = new( (/nlat,nlon/), "double")
  copyatt( omegap, fuv->omega(0,0,:,:))

  print("INTERPOLATING omega to 500mb")
  omegap=0
  do n=0,ntimes-1
    pnew=500
    varlev = fuv->omega(n,:,:,:)
    interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    extrap = True  ; is extrapolation desired if data is outside the range of PS
    P0mb = 1000    ; ps in Pa, but this argument must be in mb
    var  = vinth2p(varlev,fuv->hyam,fuv->hybm,pnew,fuv->ps(n,:,:),interp,P0mb,1,extrap)          
    omegap(:,:) = omegap(:,:) + (/ var(0,:,:) /)/ntimes
  end do
  print("writing omega500")
  ncdf->omega500=omegap

  print("INTERPOLATING omega to 850mb")
  omegap=0
  do n=0,ntimes-1
    pnew=850
    varlev = fuv->omega(n,:,:,:)
    interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    extrap = True  ; is extrapolation desired if data is outside the range of PS
    P0mb = 1000    ; ps in Pa, but this argument must be in mb
    var = vinth2p(varlev,fuv->hyam,fuv->hybm,pnew,fuv->ps(n,:,:),interp,P0mb,1,extrap)          
    omegap(:,:) = omegap(:,:) + (/ var(0,:,:) /)/ntimes
  end do
  print("writing omega850")
  ncdf->omega850=omegap
  delete(omegap)
end if

print("writing dim_avg ps")
ncdf->ps = dim_avg_n_Wrap(fuv->ps,0)
print("writing dim_avg zeta")
ncdf->zeta = dim_avg_n_Wrap(fuv->zeta,0)
print("writing dim_avg div")
ncdf->div = dim_avg_n_Wrap(fuv->div,0)

ncdf->gw = fuv->gw
ncdf->ilev = fuv->ilev
ncdf->hyam = fuv->hyam
ncdf->hyai = fuv->hyai
ncdf->hybm = fuv->hybm
ncdf->hybi = fuv->hybi

ncdf@ne = fuv@ne

end

