load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin

P0=100000.0


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

  hyam =   f->hyam         
  hybm =   f->hybm
  P0mb =   P0*0.01
  ps   =   f->ps  

  lev_p           = (/100, 200, 300, 400, 500, 600,700,800,900,1000/)
  lev_p@units     = "hPa"                        ; required for vinth2p

nlev = dimsizes(lev_t)
nlat = dimsizes(lat_t)
nlon = dimsizes(lon_t)

print(lev_t)

;tt=19
field=f->Th(tt,:,:)

;  vp = vinth2p (field(:,:), hyam,hybm, lev_p ,ps(tt,:), 2, P0mb, 1, False )


selectedcols=ind(lat_t.eq.0)

print(selectedcols)

print(lat_t(selectedcols))
print(lon_t(selectedcols))

nscols=dimsizes(selectedcols)
nplev=dimsizes(lev_p)

print("num select cols="+nscols)

vals=new((/nplev,nscols/),double)
psdummy=new((/nlat,1/),double)
psdummy(:,0)=ps(tt,:)
valsdummy=new((/nlev,nlat,1/),double)
valsdummy(:,:,0)=field(:,:)

vp = vinth2p (valsdummy, hyam,hybm, lev_p ,psdummy(:,:), 2, P0mb, 1, False )


do ii=0,nscols-1
;print("ii="+ii+", aa="+aa(ii))
;print("vals(0:nlev-1,ii)"+vals(0:nlev-1,ii))
;print("f->Th(tt,0:nlev-1,aa(ii))"+f->Th(tt,0:nlev-1,aa(ii)))

;vals(0:nlev-1,ii)=f->Th(tt,0:nlev-1,aa(ii))

;sample=field(:,aa(i))
;vals(:,ii)=interpolated(:)

vals(:,ii)=vp(:,selectedcols(ii),0)


end do



;-- define workstation
  wks = gsn_open_wks("png","plot-in-pressure-t"+tostring(tt))

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
 
  res@sfXArray                 =  lon_t(selectedcols)      ;-- uses lon_t as plot x-axis
;  res@sfYArray                 =  lev_t/100  ;-- uses lev_t in hPa as plot y-axis
  res@sfYArray                 =  lev_p  ;-- uses lev_t in hPa as plot y-axis

  res@trYReverse               =  True       ;-- reverses y-axis

;-- generate the plot
  plot = gsn_csm_contour(wks,vals,res)

end
