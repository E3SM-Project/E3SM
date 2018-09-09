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

filename="./movies/dcmip2016_test21.nc"
f = addfile(filename,"r")

time = f->time
nt = dimsizes(time)-1

data_on_constant_Z=False     ; is data already on CONSTANT Z surfaces?
convert_hybridP_to_Z=False   ; is data on hybrid pressure levels?
hasTimeIx=True               ; does file have time index?
timeStep=nt                   ; If yes, what time index do you want to plot?
model="cam-se"              ; used for mainStr, but also for model specific if statements
mainStr=model+" Day 8"       ; main string for plot titles
out_type="pdf"                ; output format, popular options are x11, png, pdf, eps
;out_type="x11"                ; output format, popular options are x11, png, pdf, eps

Uname="u"                    ; Variable name of zonal wind (typically "U" or "ua")
Vname="v"                    ; Variable name of meridional wind (typically "V" or "ua")
PSname="ps"                  ; Variable name of surface pressure (typically "PS","PSL","ps",or "slp")
Tname="T"                    ; Variable name of air temperature (typically "T" or "ta")
Zname="z"                    ; Variable name of height (typically "Z" or "za")
PRECTname="precl"            ; Variable name of height (typically "Z" or "za")

;=======================================================================================

;=======================================================================================
; Other settings, required to remain constant for intercomparison
; Generally, don't touch unless experimenting/debugging
;=======================================================================================
stride=1                       ; subset in horizontal, leave @ 1 unless debugging
minLat=0.                      ; max lat to read in (deg)
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

;=======================================================================================
; Input data checks
;=======================================================================================
if (data_on_constant_Z .and. convert_hybridP_to_Z) then
  print("Both data_on_constant_Z and convert_hybridP_to_Z cannot be true, exiting...")
  exit
end if

;=======================================================================================
; Get file, coordinate variables, and U, V, PSL
;=======================================================================================
print("Loading data from file...")

lat = f->lat({minLat:maxLat:stride})
lon = f->lon(::stride)
lev = f->lev(:)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

if (hasTimeIx) then
  U =  f->$Uname$(timeStep,:,{minLat:maxLat:stride},::stride)
  V =  f->$Vname$(timeStep,:,{minLat:maxLat:stride},::stride)
  T =  f->$Tname$(timeStep,:,{minLat:maxLat:stride},::stride)
  PS = f->$PSname$(timeStep,{minLat:maxLat:stride},::stride)
  PRECT = f->$PRECTname$(timeStep,{minLat:maxLat:stride},::stride)
else
  U =  f->$Uname$(:,{minLat:maxLat:stride},::stride)
  V = f->$Vname$(:,{minLat:maxLat:stride},::stride)
  T = f->$Tname$(:,{minLat:maxLat:stride},::stride)
  PS = f->$PSname$({minLat:maxLat:stride},::stride)
  PRECT = f->$PRECTname$({minLat:maxLat:stride},::stride)
end if
; If U and V are not m/s, please convert here
U@units="m/s"
V@units="m/s"
T@units="K"
T@long_name="Temperature"
;=======================================================================================

;=======================================================================================
; Find center of storm by minimizing PS
;=======================================================================================
print("Finding minimum PS location...")
a = new((/nlat,nlon/),typeof(PS))
a(:,:) = PS(:,:)
a1D      = ndtooned(a)
dsizes_a = dimsizes(a)
a_indices  = ind_resolve(minind(a1D),dsizes_a) ; Resolve 1D indices to original array
psminlat = lat(a_indices(0,0))
psminlon = lon(a_indices(0,1))
delete([/a,a1D,dsizes_a,a_indices/])
print("... PS minimum found at lat: "+psminlat+" lon: "+psminlon)
;=======================================================================================


;=======================================================================================
; Calculate temperature anomaly
;=======================================================================================
Tanom = T
Tref = T(:,{psminlat+offAnomDist},{psminlon+offAnomDist})
Tanom = T - conform(T,Tref,0)
Tanom@long_name="Temperature anomaly"
;=======================================================================================


;=======================================================================================
; Do Z calculations/interpolations if necessary
;=======================================================================================
if (.not. data_on_constant_Z) then
  ;=======================================================================================
  ; Convert from hybrid levels to Z (CAM)
  ;=======================================================================================
  if (convert_hybridP_to_Z)
    print("Converting from hybrid P to Z...")
    ; If need_to_convert_P_to_Z is true, variables T,Q,hyam,hybm,hyai,hybm,P0 needed
    ; all hybrid coefs and TV need to be TOP TO BOTTOM!

    hyam=f->hyam
    hybm=f->hybm
    hyai=f->hyai
    hybi=f->hybi
; P0=f->P0
    P0=100000
    Tconv = f->T(timeStep,:,{minLat:maxLat:stride},::stride)
    Qconv = f->Q(timeStep,:,{minLat:maxLat:stride},::stride)

    ; Calculate TV from T (K) and Q (kg/kg)
    TV=Tconv*(1.+0.61*Qconv)

    ; PHIS is nlatxnlon array = 0.0 for DCMIP Test 2
    PHIS=PS
    PHIS=0.0

    Z = cz2ccm(PS,PHIS,TV,P0,hyam(::-1),hybm(::-1),hyai(::-1),hybi(::-1))
    Z@units="m"

    delete([/Tconv,Qconv,TV,PHIS/])

  else
    Z = f->geo(timeStep,:,{minLat:maxLat:stride},::stride)/9.81806
  end if


  ;=======================================================================================
  ; If all else fails, try to load Z directly from file
  ;=======================================================================================
  if(.not. isdefined("Z")) then
    if (isfilevar(f, Zname)) then
      print("Found Z on file...")
      Z = f->Z
    else
      print("FATAL: Z needs to be either loaded from a file or otherwise defined before continuing...")
      exit
    end if
  end if
  ;=======================================================================================

end if

;=======================================================================================
; Generate Z levels of interest
;=======================================================================================
print("Generating constant Z levels to interpolate to")
ZlevPBL = fspan(zintmin,Zmidpt,nzint/2)
ZlevABL = fspan(Zmidpt,zintmax,nzint/2)

Zlev=(/100.,1000.,1500.,2500.,5000.,10000.,15000./)
Zlev@units = "m"
Zlev!0     = "lev"
Zlev&lev = Zlev

;=======================================================================================
; Interpolate lat/lon variables to constant Z levels
;=======================================================================================
print("Interpolating to Z surfaces")
U_Z     = int2p_n_Wrap(Z,U,Zlev,2,0)
V_Z     = int2p_n_Wrap(Z,V,Zlev,2,0)
Tanom_Z = int2p_n_Wrap(Z,Tanom,Zlev,2,0)
;=======================================================================================

WIND_Z = U_Z
WIND_Z = sqrt(U_Z^2+V_Z^2)
WIND_Z@long_name="Horizontal wind"

PRECT = PRECT*8.64e7   ; convert m/s to mm/day
PRECT@long_name="Precipitation rate"
PRECT@units="mm/d"

print("Plotting...")

wks   = gsn_open_wks (out_type,"x-sections")               ; send graphics to PNG file
contour = new(4,"graphic")

gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

res                 = True
res@gsnDraw         = False
res@gsnFrame         = False
res@gsnSpreadColors  = True        ; Span full color map
res@cnFillOn         = True        ; Turn on contour fill
res@cnLinesOn        = False
res@cnLineLabelsOn   = False
res@cnInfoLabelOn    = False
res@gsnAddCyclic = True

latWidth = 10.0
lonWidth = 10.0
res@mpOutlineOn = False
res@mpMaxLatF = psminlat + latWidth
res@mpMinLatF = psminlat - latWidth
res@mpMinLonF = psminlon - lonWidth
res@mpMaxLonF = psminlon + lonWidth

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF      =  5.0
res@cnMinLevelValF       =  10.0
res@cnMaxLevelValF       =  80.0
res@tiMainString="1500m Horiz. Wind"
contour(0) = gsn_csm_contour_map(wks,WIND_Z({1500.0},:,:),res)  ; create the plot

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF      =  10.0
res@cnMinLevelValF       =  -60.0
res@cnMaxLevelValF       =  60.0
res@tiMainString="100m U-Wind"
contour(1) = gsn_csm_contour_map(wks,U_Z({100.0},:,:),res)  ; create the plot

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF      =  1.0
res@cnMinLevelValF       =  0.0
res@cnMaxLevelValF       =  10.0
res@tiMainString="5000m T anom."
contour(2) = gsn_csm_contour_map(wks,Tanom_Z({5000.0},:,:),res)  ; create the plot

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF      =  200.0
res@cnMinLevelValF       =  200.0
res@cnMaxLevelValF       =  4000.0
res@tiMainString="Precip. rate"
contour(3) = gsn_csm_contour_map(wks,PRECT(:,:),res)  ; create the plot

resP                     = True                ; modify the panel plot
resP@gsnMaximize         = True

gsn_panel(wks,contour,(/2,2/),resP)             ; now draw as one plot

end
