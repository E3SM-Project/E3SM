load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin
;-- define file name
  diri  = "./"
;  fili  = "planar_nonhydro_gravity_wave1.nc"
  fili  = "planar_rising_bubble1.nc"

;-- open file and read variable 
  f     =  addfile(diri+fili, "r")  

;  var   = f->Th(0,:,{0},{-150000:147487})     ;-- first time step, latitude=40N, longitude=0-60E.
  lon_t = f->lon            ;-- longitude=0-60E
  lat_t = f->lat
  lev_t = f->lev                    ;-- currently 17 levels

nlev = dimsizes(lev_t)
nlat = dimsizes(lat_t)
nlon = dimsizes(lon_t)

print(lev_t)


;tt=0

aa=ind(lat_t.eq.0)

print(aa)

print(lat_t(aa))
print(lon_t(aa))

bb=dimsizes(aa)

print("bb="+bb)

vals=new((/nlev,bb/),double)


do ii=0,bb-1
;print("ii="+ii+", aa="+aa(ii))
;print("vals(0:nlev-1,ii)"+vals(0:nlev-1,ii))
;print("f->Th(tt,0:nlev-1,aa(ii))"+f->Th(tt,0:nlev-1,aa(ii)))

vals(0:nlev-1,ii)=f->Th(tt,0:nlev-1,aa(ii))
end do



;-- define workstation
  wks = gsn_open_wks("png","plot_slices-t"+tostring(tt))

;-- set resources
  res                          =  True
  res@tiMainString             = "DKRZ NCL Example: Slice plot at 40N"   ;-- title string

  cmap = read_colormap_file("MPL_GnBu")      ;-- read the color map

  res@cnFillOn                 =  True       ;-- turn on color fill
  res@cnFillPalette            = cmap(::-1,:)
  res@cnLineLabelsOn           =  False      ;-- turns off contour line labels
  res@cnInfoLabelOn            =  False      ;-- turns off contour info label
  res@lbOrientation            = "vertical"  ;-- vertical label bar
;  res@tiYAxisString            =  var@long_name+"  [hPa]"
 
  res@sfXArray                 =  lon_t(aa)      ;-- uses lon_t as plot x-axis
  res@sfYArray                 =  lev_t/100  ;-- uses lev_t in hPa as plot y-axis

;  res@trYReverse               =  True       ;-- reverses y-axis

;-- generate the plot
  plot = gsn_csm_contour(wks,vals,res)

end
