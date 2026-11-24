;_____________________________________________________________________
;
; Plot potential temperature differences from background state
;_____________________________________________________________________

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;_____________________________________________________________________
; dcmip3-1 model parameters (must match those in dcmip123_mod.F90)

ts      = 0.1                   ; timestep in seconds
interval= 100                   ; steps per write

g       = 9.80616               ; grav const
a       = 6371229.0             ; earth radius in meters
Rd      = 287.0                 ; dry gas const
cp      = 1004.5                ; heat capacity const pressure
kappa   = Rd/cp
pi      = 3.141592654
p_ref   = 100000.0              ; reference pressure

ztop    = 10000.0               ; top of the atmosphere in z coords
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

f     = addfile("./movies/dcmip2012_test31.nc","r")                       ; open ncl data file

lat1d   = f->lat*pi/180                                               ; load 1d latitutde array
lon1d   = f->lon*pi/180                                               ; load 1d longitude array
lev     = f->lev                                                      ; load 1d vertical-level array
time    = f->time

nlon  = getfilevardimsizes(f, "lon" )                                 ; get number of lon points
nlat  = getfilevardimsizes(f, "lat" )                                 ; get number of lat points
nlev  = getfilevardimsizes(f, "lev" )                                 ; get number of vertical levels
nt   = dimsizes(time)

lat = new( (/nlat,nlon/), typeof(lat1d))
lon = new( (/nlat,nlon/), typeof(lon1d))

do j=0, nlon-1
  do i=0, nlat-1
    lat(i,j)=lat1d(i)
    lon(i,j)=lon1d(j)
  end do
end do

print("nlon = "+nlon+" nlat = "+nlat+" nlev = "+nlev+" nt = "+nt)                 ; display data dimensions

t_index = nt-1                                                          ; set time to max time
if(isvar("t")) then
  t_index = t                               ; use t from command line if possible
end if
print ("t_index="+t_index)


;_____________________________________________________________________
; compute evenly spaced z-levels

dz = ztop / nlev
z  = new((/nlev/),double)

do k=0,nlev-1
  z(k) = ztop - dz*(k+0.5)
end do

;_____________________________________________________________________
; define data fields

T     = f->T(0,:,:,:)                                                 ; read Temperature field
Ts0   = new( (/nlat,nlon/)      , typeof(T), -9999.0)                 ; define initial surf temp field
ps0   = new( (/nlat,nlon/)      , typeof(T), -9999.0)                 ; define initial surf pressure field
p0    = new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)                 ; define mean temperature field
T_mean= new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)                 ; define mean temperature field
theta_mean= new( (/nlev,nlat,nlon/) , typeof(T), -9999.0)             ; define mean potential temperature field

;_____________________________________________________________________
; recompute initial conditions

Ts0 = bigG + (Teq-bigG)*exp( -(u0*N2/(4.d0*g*g))*(u0+2.d0*Om*as)*(cos(2.d0*lat)-1.d0))
ps0 = Peq * exp((u0/(4.0*bigG*Rd))*(u0+2.0*Om*as)*(cos(2.0*lat)-1.0))*(Ts0/Teq)^(cp/Rd)

do k=0, nlev-1
  p0(k,:,:)     = ps0*( (bigG/Ts0) * exp(-N2*z(k)/g) +1.d0 -(bigG/Ts0) )^(cp/Rd)
  T_mean(k,:,:) = bigG*(1.d0 - exp(N2*z(k)/g))+ Ts0*exp(N2*z(k)/g)
end do

theta_mean = T_mean*(p_ref/p0)^kappa

;_____________________________________________________________________
; get potential-temperature differences at several times

ts = 0.50                                                             ; timestep in seconds
interval = 100                                                        ; steps per write

t_sec =ts*interval*t_index                                            ; get time from time-index

theta = f->T(t_index,:,:,:)*(p_ref/p0)^kappa
theta_diff = theta - theta_mean

equator = nlat/2                                                      ; lat index of equator for odd # latitudes
theta_diff_slice = theta_diff(:,equator,:)

;w = f->w(t_index,:,equator,:)                                          ; read w field
omega = f->omega(t_index,:,equator,:)                                           ; read w field

;_____________________________________________________________________
; plot contour field of lon vs height at the equator

maxval= 2.0
minval=-2.0
gap   = (maxval-minval)/10

print("maxval = "+maxval)
print("minval = "+minval)
print("spacing = "+gap)

  res1                      = True
  res1@gsnMaximize          = True

  res1@gsnDraw              = True
  res1@gsnFrame             = False
  res1@cnFillOn             = True
  res1@cnLinesOn            = False

  res1@gsnSpreadColors      = True
  res1@lbLabelAutoStride    = True

  res1@gsnCenterString      = ""
  res1@tiMainString         = ""
  res1@vpWidthF             = 0.50
  res1@vpHeightF            = 0.25
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnInfoLabelOn        = False
  res1@cnLineLabelsOn       = False
  res1@lbLabelBarOn         = True
  res1@lbLabelStride        = 2

  res1@lbOrientation        = "horizontal"
  res1@sfYArray             = z             	; use z for y axis
  res1@sfXArray             = lon1d*180.0/pi  ; use lon1d for x axis
  res1@trYReverse           = False         	; reverses y-axis, false
  res1@tiYAxisString        = "height (m)"    ; y-axis title
  ;res1@tiXAxisString        = "Longitude"      ; x-axis title
  res1@cnMaxLevelValF       = maxval
  res1@cnMinLevelValF       = minval
  res1@cnLevelSpacingF      = gap
  res1@gsnPaperOrientation  = "portrait"

  pres                      = True
  pres@gsnMaximize          = True                     ; make ps, eps, pdf as large as possible
  pres@gsnPanelLabelBar     = True
  pres@pmLabelBarOrthogonalPosF = -0.03
  pres@amJust               = "TopLeft"
  pres@gsnFrame             = False
  pres@lbLabelStride        = 1
  pres@gsnPaperOrientation  = "portrait"

  ;wks_type                  = "png"
  ;wks_type@wkWidth          = 2500
  ;wks_type@wkHeight         = 2500

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  wks_type@wkOrientation    = "portrait"

if(isvar("t")) then
  name =sprinti("test31_omega_%0.4i",t_index)
else
  name ="./test31_omega"
end if

print("name = "+name)
wks = gsn_open_wks(wks_type,name)

gsn_define_colormap(wks,"gui_default")  ; choose a colormap
;gsn_define_colormap(wks,"BlWhRe")  ; choose a colormap

plot = new(1,graphic)

  res1@gsnStringFontHeightF     = 0.015
  res1@tmXBLabelFontHeightF     = 0.01
  res1@tmYLLabelFontHeightF     = 0.01
  res1@tiYAxisFontHeightF       = 0.01
  res1@tiXAxisFontHeightF       = 0.01
  ;res1@gsnCenterString = "omega  " +t_sec+ " sec "
  res1@tiMainString="omega " +t_sec+ " sec "
  res1@tiMainFontHeightF       = 0.01

plot(0) = gsn_csm_contour(wks,omega,res1)
; gsn_panel(wks,plot,(/1,1/),pres)

frame(wks)

end

