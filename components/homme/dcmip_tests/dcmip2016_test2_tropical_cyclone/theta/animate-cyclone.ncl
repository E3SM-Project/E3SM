;=======================================================================================
; This NCL code calculates radially-averaged tangential and radial wind components
; as well as T anomaly for DCMIP test case #2 (cyclone)
; this code requires the accompanying function set "radialAvg.ncl"
;
; Usage: User should modify "user options" for their particular data set. Currently,
; U, V, T, PS are required as variables.
; If variables are on constant Z surfaces, life is easy.
;
; Grepping for "MODELSPEC" will point to possible areas of the code in need of modification
; for model specific output
;
; Written by Colin Zarzycki (zarzycki@ucar.edu)
; Version 0.1 (6/5/2016) - DCMIP-2016 release
;=======================================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

;=======================================================================================
; User options
;=======================================================================================


name  = "animation_tc"
wks_type                  = "pdf" ; output format, popular options are x11, png, pdf, eps
wks_type@wkPaperHeightF   = 6.5
wks_type@wkPaperWidthF    = 11
wks = gsn_open_wks(wks_type,name)


filename="./movies/r50-prec0-pbl0-dcmip2016_test21.nc"
f     = addfile(filename,"r")
time  = f->time
nt    = dimsizes(time)-1

do ti=1,nt,1

data_on_constant_Z    =False  ; is data already on CONSTANT Z surfaces?
convert_hybridP_to_Z  =False  ; is data on hybrid pressure levels?
hasTimeIx             =True   ; does file have time index?
timeStep              =ti     ; If yes, what time index do you want to plot?
model="ACME-A"                ; used for mainStr, but also for model specific if statements
mainStr=model ;+" Day 8"        ; main string for plot titles

Uname ="u"                    ; Variable name of zonal wind (typically "U" or "ua")
Vname ="v"                    ; Variable name of meridional wind (typically "V" or "ua")
PSname="ps"                   ; Variable name of surface pressure (typically "PS","PSL","ps",or "slp")
Tname ="T"                    ; Variable name of air temperature (typically "T" or "ta")
Zname ="z"                    ; Variable name of height (typically "Z" or "za")
PRECTname="precl"             ; Variable name of height (typically "Z" or "za")

print("timeIndex = "+timeStep)

;=======================================================================================
; Other settings, required to remain constant for intercomparison
; Generally, don't touch unless experimenting/debugging
;=======================================================================================
stride=1                       ; subset in horizontal, leave @ 1 unless debugging
minLat=-40.                      ; max lat to read in (deg)
maxLat=80.                     ; min lat to read in (deg)
thetaWindMax=60.0              ; max theta wind to plot (m/s)
radMinMax=18.0                 ; min/max for radial plot (m/s)
tAnomMax=15.0                  ; max for tAnom plots (K)
offAnomDist=5.0                ; offset distance for anom calcs (deg)
zintmin=20.0                   ; min height for z interpolation (m)
zintmax=18000.0                ; max height for z interpolation (m)
nzint=100                      ; num points for z interpolation
Zmidpt=2500.0                  ; half of Z points BELOW this height, half ABOVE

;=======================================================================================
; Get file, coordinate variables, and U, V, PSL
;=======================================================================================
print("Loading data from file...")

lat  = f->lat({minLat:maxLat:stride})
lon  = f->lon(::stride)
lev  = f->lev(:)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
U =  f->$Uname$(timeStep,:,{minLat:maxLat:stride},::stride)
V =  f->$Vname$(timeStep,:,{minLat:maxLat:stride},::stride)
W =  f->w(timeStep,:,{minLat:maxLat:stride},::stride)

T =  f->$Tname$(timeStep,:,{minLat:maxLat:stride},::stride)
PS = f->$PSname$(timeStep,{minLat:maxLat:stride},::stride)
PRECT = f->$PRECTname$(timeStep,{minLat:maxLat:stride},::stride)
;Qv    = f->Q(timeStep,:,{minLat:maxLat:stride},::stride)
;Qv0   = f->Q(0,:,{minLat:maxLat:stride},::stride)
;Qv=Qv-Qv0
;Qv=Qv*1000.0
;Qv@units = "g/kg"

; If U and V are not m/s, please convert here
U@units="m/s"
V@units="m/s"
W@units="m/s"
T@units="K"
T@long_name="Temperature"

;=======================================================================================
; Find center of storm by minimizing PS
;=======================================================================================
print("Finding minimum PS location...")
a = new((/nlat,nlon/),typeof(PS))
a(:,:)    = PS(:,:)
a1D      = ndtooned(a)
dsizes_a = dimsizes(a)
a_indices  = ind_resolve(minind(a1D),dsizes_a) ; Resolve 1D indices to original array
psminlat = lat(a_indices(0,0))
psminlon = lon(a_indices(0,1))
delete([/a,a1D,dsizes_a,a_indices/])
print("... PS minimum found at lat: "+psminlat+" lon: "+psminlon)

;=======================================================================================
; Calculate temperature anomaly
;=======================================================================================
Tanom = T
Tref = T(:,{psminlat+offAnomDist},{psminlon+offAnomDist})
Tanom = T - conform(T,Tref,0)
Tanom@long_name="Temperature anomaly"

;=======================================================================================
; Do Z calculations/interpolations if necessary
;=======================================================================================

Z = f->geo(timeStep,:,{minLat:maxLat:stride},::stride)/9.81806

print("Generating constant Z levels to interpolate to")
ZlevPBL = fspan(zintmin,Zmidpt,nzint/2)
ZlevABL = fspan(Zmidpt,zintmax,nzint/2)

Zlev=(/100.,1000.,1500.,2500.,5000.,10000.,15000./)
Zlev@units = "m"
Zlev!0     = "lev"
Zlev&lev    = Zlev

;=======================================================================================
; Interpolate lat/lon variables to constant Z levels
;=======================================================================================
print("Interpolating to Z surfaces")
U_Z     = int2p_n_Wrap(Z,U,Zlev,2,0)
V_Z     = int2p_n_Wrap(Z,V,Zlev,2,0)
Tanom_Z = int2p_n_Wrap(Z,Tanom,Zlev,2,0)
;Qv_Z   = int2p_n_Wrap(Z,Qv,Zlev,2,0)
W_Z     = int2p_n_Wrap(Z,W,Zlev,2,0)
W_Z@long_name="vertical wind speed"

WIND_Z  = U_Z
WIND_Z  = sqrt(U_Z^2+V_Z^2)
WIND_Z@long_name="horizontal wind speed"

PRECT       = PRECT*8.64e7   ; convert m/s to mm/day
PRECT@long_name="Precipitation rate"
PRECT@units="mm/d"

print("Plotting...")

contour = new(3,"graphic")

;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

res                  = True
res@gsnDraw          = False
res@gsnFrame         = False
;res@gsnSpreadColors  = True        ; Span full color map
res@cnFillOn         = True        ; Turn on contour fill
res@cnLinesOn        = False
res@cnLineLabelsOn   = False
res@cnInfoLabelOn    = False
res@gsnAddCyclic     = True

latWidth1 = 40.0
lonWidth1 = 40.0

latWidth = 40.0
lonWidth = 40.0

res@mpOutlineOn = True
res@mpMaxLatF   = psminlat + latWidth
res@mpMinLatF   = psminlat - latWidth
res@mpMaxLonF   = psminlon + lonWidth
res@mpMinLonF   = psminlon - lonWidth
res@mpCenterLonF= psminlon

print(res@mpMaxLatF)
print(res@mpMinLatF)
print(res@mpMaxLonF)
print(res@mpMinLonF)

res0=res
res1=res
res2=res
res3=res

res0@cnFillPalette = "WhiteBlueGreenYellowRed"
res0@cnLevelSelectionMode = "ManualLevels"
res0@cnLevelSpacingF      =  1
res0@cnMinLevelValF       =  1
res0@cnMaxLevelValF       =  24
res0@mpMaxLatF   = psminlat + latWidth1
res0@mpMinLatF   = psminlat - latWidth1
res0@mpMaxLonF   = psminlon + lonWidth1
res0@mpMinLonF   = psminlon - lonWidth1

res0@tiMainString="Surface Wind Speed t="+sprintf("%3.1f", time(ti) )+" days"
contour(0) = gsn_csm_contour_map(wks,WIND_Z({100.0},:,:),res0)

res1@cnFillPalette        ="BlAqGrYeOrReVi200"
res1@cnLevelSelectionMode = "ManualLevels"
res1@cnLevelSpacingF      =  2e2
res1@cnMinLevelValF       =  989e2
res1@cnMaxLevelValF       =  1020e2
res1@tiMainString         ="Surface Pressure"
contour(1) = gsn_csm_contour_map(wks,PS,res1)

res2@cnFillPalette        ="temp_diff_18lev"
res2@cnLevelSelectionMode = "ManualLevels"
res2@cnLevelSpacingF      =  0.5
res2@cnMinLevelValF       = -5.0
res2@cnMaxLevelValF       =  5.0
res2@tiMainString         ="5km T anom."
contour(2) = gsn_csm_contour_map(wks,Tanom_Z({5000.0},:,:),res2)

;res3@cnFillPalette        = "BlueWhiteOrangeRed"; "gsdtol"
;res3@cnLevelSelectionMode = "ManualLevels"
;res3@cnLevelSpacingF      =  0.01
;res3@cnMinLevelValF       = -0.1
;res3@cnMaxLevelValF       =  0.1
;res3@tiMainString         ="W vel at 5km"
;contour(3) = gsn_csm_contour_map(wks,W_Z({5000.0},:,:),res3)


resP                     = True                ; modify the panel plot
resP@gsnMaximize         = True
resP@gsnPaperOrientation = "portrait"

gsn_panel(wks,contour,(/1,3/),resP)             ; now draw as one plot
end do

end
