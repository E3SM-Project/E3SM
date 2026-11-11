;_____________________________________________________________________
;
; Plot potential temperature differences from background state
; works for unevenly spaced z-levels
;_____________________________________________________________________

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;_____________________________________________________________________
; dcmip3-1 model parameters (must match those in dcmip123_mod.F90)

g       = 9.80616               ; grav const
a       = 6371229.0             ; earth radius in meters
Rd      = 287.0                 ; dry gas const
cp      = 1004.5                ; heat capacity const pressure
kappa   = Rd/cp
pi      = 3.141592654
p_ref   = 100000.0              ; reference pressure

ztop    = 10000.0
X       = 125.0                 ; reduced Earth reduction factor
Om      = 0.0                   ; rotation rate of Earth
as      = a/X                   ; new radius of small Earth
u0      = 20.0                  ; max zonal wind speed
Teq     = 300.0                 ; temperature at equator
Peq     = 100000.0              ; reference ps at equator
lambdac = 2.0*pi/3.0            ; lon of pert center
d       = 5000.0                ; width for pert
phic    = 0.0                   ; lat of pert Center
delta_theta = 1.0               ; max amplitude of temp perturbation
Lz      = 20000.0               ; vertical wavelength of pert
N       = 0.010                 ; Brunt-Vaisala frequency
N2      = N*N                   ; Brunt-Vaisala frequency squared
bigG    = (g*g)/(N2*cp)         ; constant

;_____________________________________________________________________
; open file, read in data (time,lev,lat,long) in range [0 to n-1]

f       = addfile("./movies/dcmip2012_test31.nc","r")                  ; open ncl data file
lat1d   = f->lat*pi/180                                                 ; load 1d latitutde array
lon1d   = f->lon*pi/180                                                 ; load 1d longitude array
lev     = f->lev                                                        ; load 1d vertical-level array
time    = f->time
hyam    = f->hyam
hybm    = f->hybm
nlon    = getfilevardimsizes(f, "lon" )                                 ; get number of lon points
nlat    = getfilevardimsizes(f, "lat" )                                 ; get number of lat points
nlev    = getfilevardimsizes(f, "lev" )                                 ; get number of vertical levels
nt      = dimsizes(time)
equator = nlat/2                                                        ; lat index of equator for odd # latitudes

lat = new( (/nlat,nlon/), typeof(lat1d))
lon = new( (/nlat,nlon/), typeof(lon1d))

do j=0, nlon-1
  do i=0, nlat-1
    lat(i,j)=lat1d(i)
    lon(i,j)=lon1d(j)
  end do
end do

print("nlon = "+nlon+" nlat = "+nlat+" nlev = "+nlev+" nt = "+nt)                 ; display data dimensions
print("equator = "+lat(equator,0))

t_index = nt-1                                                          ; set time to max time
if(isvar("t")) then
  t_index = t                               ; use t from command line if possible
end if
print ("t_index="+t_index)

;_____________________________________________________________________
; define data fields

T     = f->T(0,:,:,:)                                                 ; read Temperature field
Ts0   = new( (/nlat,nlon/)      , typeof(T), -9999.0)                 ; define initial surf temp field
ps0   = new( (/nlat,nlon/)      , typeof(T), -9999.0)                 ; define initial surf pressure field
p0    = new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)                 ; define mean pressure at midpoints
T_mean= new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)                 ; define mean temperature field
theta_mean= new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)             ; define mean potential temperature field

;_____________________________________________________________________
; recompute initial conditions

Ts0 = bigG + (Teq-bigG)*exp( -(u0*N2/(4.d0*g*g))*(u0+2.d0*Om*as)*(cos(2.d0*lat)-1.d0))
ps0 = Peq * exp((u0/(4.0*bigG*Rd))*(u0+2.0*Om*as)*(cos(2.0*lat)-1.0))*(Ts0/Teq)^(cp/Rd)

; get height of eta levels from the pressure
h = new( (/nlev/), typeof(T))

do k=0, nlev-1
  p = hyam(k)*p_ref + hybm(k)*Peq
  h(k) = (-g/N2)*log( (Teq/bigG)*( (p/Peq)^(Rd/cp) - 1.d0  ) + 1.d0 )
end do

do k=0, nlev-1
  print( "height (" +k+ ")= " +h(k))
end do

do k=0, nlev-1
  p0(k,:,:)     = hyam(k)*p_ref + hybm(k)*ps0
  T_mean(k,:,:) = bigG*(1.d0 - exp(N2*h(k)/g))+ Ts0*exp(N2*h(k)/g)
end do

  ;_____________________________________________________________________
  ; plot contour field of lon vs height at the equator

  maxval= delta_theta/10
  minval=-delta_theta/10
  gap   = (maxval-minval)/10

  print("maxval = "+maxval)
  print("minval = "+minval)
  print("gap    = "+gap)

  res1                      = True
  res1@gsnMaximize          = True
  res1@gsnDraw              = False
  res1@gsnFrame             = False
  res1@gsnSpreadColors      = True
  res1@gsnCenterString      = ""
  res1@tiMainString         = ""
  res1@vpWidthF             = 2
  res1@vpHeightF            = 1
  res1@lbLabelBarOn         = True
  res1@lbOrientation        = "vertical"
  res1@lbLabelStride        = 1
  res1@cnFillOn             = True
  res1@cnLinesOn            = False
  res1@cnInfoLabelOn        = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = maxval
  res1@cnMinLevelValF       = minval
  res1@cnLevelSpacingF      = gap
  res1@sfYArray             = h
  res1@sfXArray             = lon1d*180.0/pi    ; use lon1d for x axis
  res1@trYReverse           = False ;True             ; reverses y-axis, false
  res1@tiYAxisString        = "height (m)"      ; y-axis title
  res1@tiXAxisString        = "Longitude"       ; x-axis title
  res1@gsnPaperOrientation  = "portrait"
  res1@gsnYAxisIrregular2Linear = True
  ;res1@trYTensionF          = 1.0

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 6
  wks_type@wkPaperWidthF    = 11
  wks_type@wkOrientation    = "portrait"

  if(isvar("t")) then
    filename =sprinti("./dcmip2012_test3_theta_diff_%0.4i",t_index)
  else
    filename ="./dcmip2012_test3_theta_diff"
  end if

  print("filename = "+filename)
  wks = gsn_open_wks(wks_type,filename)
  gsn_define_colormap(wks,"gui_default")  ; choose a colormap


;_____________________________________________________________________
; get potential-temperature differences at several times
do t_index = 0,nt-1  
  days2sec = 3600.0*24.0
  t_sec = time(t_index)*days2sec
  
  theta      = f->T(t_index,:,:,:)*(p_ref/p0)^kappa
  theta_mean = T_mean*(p_ref/p0)^kappa
  
  theta_diff = theta - theta_mean
  
  theta_diff_slice = theta_diff(0:nlev-1,equator,:)


  plot = new(1,graphic)
  res1@gsnStringFontHeightF     = 0.04
  res1@tmXBLabelFontHeightF     = 0.03
  res1@tmYLLabelFontHeightF     = 0.03
  res1@tiYAxisFontHeightF       = 0.03
  res1@tiXAxisFontHeightF       = 0.03

  pltTitle = sprintf("CAM-SE Test 3-1  Theta'  t=%5.2f s",t_sec)
  res1@gsnCenterString = pltTitle

  plot(0) = gsn_csm_contour(wks,theta_diff_slice,res1)

  draw(wks)
  frame(wks)
end do
end

