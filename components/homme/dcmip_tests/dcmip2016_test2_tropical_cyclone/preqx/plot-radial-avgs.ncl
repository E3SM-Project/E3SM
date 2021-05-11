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
load "radialAvg.ncl"

begin

;=======================================================================================
; User options
;=======================================================================================
filename="./movies/dcmip2016_test21.nc"
f = addfile(filename,"r")

time = f->time
nt = dimsizes(time)-1

data_on_constant_Z=False     ; is data already on CONSTANT Z surfaces?
convert_hybridP_to_Z=False    ; is data on hybrid pressure levels?
hasTimeIx=True               ; does file have time index?
timeStep=nt                  ; If yes, what time index do you want to plot?
model="PREQX"               ; used for mainStr, but also for model specific if statements
mainStr=model+" Day 8"       ; main string for plot titles
;outType="x11"                ; output format, popular options are x11, png, pdf, eps
outType="pdf"                ; output format, popular options are x11, png, pdf, eps

Uname="u"                    ; Variable name of zonal wind (typically "U" or "ua")
Vname="v"                    ; Variable name of meridional wind (typically "V" or "ua")
PSname="ps"                  ; Variable name of surface pressure (typically "PS","PSL","ps",or "slp")
Tname="T"                    ; Variable name of air temperature (typically "T" or "ta")
Zname="z"                    ; Variable name of height (typically "Z" or "za")
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
else
  U =  f->$Uname$(:,{minLat:maxLat:stride},::stride)
  V = f->$Vname$(:,{minLat:maxLat:stride},::stride)
  T = f->$Tname$(:,{minLat:maxLat:stride},::stride)
  PS = f->$PSname$({minLat:maxLat:stride},::stride)
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
  end if

  Z = f->geo(timeStep,:,{minLat:maxLat:stride},::stride)/9.81806

  ;=======================================================================================
  ; Generate uniform Z levs to interpolate to
  ; Extra density (half of pts) is specified in PBL, controlled by Zmidpt
  ;=======================================================================================
  print("Generating constant Z levels to interpolate to")
  ZlevPBL = fspan(zintmin,Zmidpt,nzint/2)
  ZlevABL = fspan(Zmidpt,zintmax,nzint/2)

  Zlev=new(nzint-1,"float")
  Zlev(0:(nzint/2)-1)=ZlevPBL
  Zlev((nzint/2)-1:nzint-2)=ZlevABL
  Zlev@units = "m"
  Zlev!0     = "lev"
  Zlev&lev = Zlev
  delete([/ZlevPBL,ZlevABL/])
  print(Zlev)

  ;=======================================================================================

  ;=======================================================================================
  ; Interpolate lat/lon variables to constant Z levels
  ;=======================================================================================
  print("Interpolating from Z model levels to constant Z surfaces")
  U_Z     = int2p_n_Wrap(Z,U,Zlev,2,0)
  V_Z     = int2p_n_Wrap(Z,V,Zlev,2,0)
  Tanom_Z = int2p_n_Wrap(Z,Tanom,Zlev,2,0)
  ;=======================================================================================
else
  print("Data already on Z surfaces, copying to new vars")
  U_Z     = U
  V_Z     = V
  Tanom_Z = Tanom
end if

;=======================================================================================
; Calculate radial and tangential wind components on native grid 
;=======================================================================================
print("Separating winds into radial and tangential components given PS min loc...")
vComps = calcWindComponents(U_Z,V_Z,lat,lon,psminlat,psminlon)
v_rad = vComps[0]
v_theta = vComps[1]
delete(vComps)
;=======================================================================================


;=======================================================================================
; Perform radial integrations
;=======================================================================================
print("Do radial integrations...")
rad_v_theta  = radialAvg3D(v_theta,lat,lon,v_theta&lev,psminlat,psminlon,600.,True)
rad_v_rad    = radialAvg3D(v_rad  ,lat,lon,v_rad&lev  ,psminlat,psminlon,600.,True)
rad_t_anom   = radialAvg3D(Tanom_Z,lat,lon,Tanom_Z&lev  ,psminlat,psminlon,600.,True)
;=======================================================================================



print("Plotting...")

res  = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnSpreadColors  = True        ; Span full color map
res@cnFillOn         = True        ; Turn on contour fill
res@cnLinesOn        = False
res@cnLineLabelsOn   = False
res@cnInfoLabelOn    = False
res@tiYAxisString    = "Height (m)"
res@gsnYAxisIrregular2Linear = True
res@gsnXAxisIrregular2Linear = True

wks   = gsn_open_wks (outType,"t_anom_"+model)
gsn_define_colormap(wks,"matlab_jet")
res_rad = res
res_rad@cnLevelSelectionMode = "ExplicitLevels"
res_rad@cnLevels = ispan(0,toint(tAnomMax),1)
plot = gsn_csm_contour(wks,rad_t_anom,res_rad)
draw(plot)
frame(wks)

delete(plot)
delete(wks)
delete(res_rad)

wks   = gsn_open_wks (outType,"v_theta_"+model)
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
res_rad = res
res_rad@tiMainString     = mainStr+" tangential wind"
res_rad@cnLevelSelectionMode = "ExplicitLevels"
res_rad@cnLevels = fspan(0.0,thetaWindMax,21)
plot = gsn_csm_contour(wks,rad_v_theta,res_rad)
draw(plot)
frame(wks)

delete(plot)
delete(wks)
delete(res_rad)

wks   = gsn_open_wks (outType,"v_rad_"+model)
gsn_define_colormap(wks,"hotcolr_19lev")
res_rad = res
res_rad@tiMainString     = mainStr+" radial wind"
res_rad@cnLevelSelectionMode = "ExplicitLevels"
res_rad@cnLevels = ispan(toint(-radMinMax),toint(radMinMax),1)
plot = gsn_csm_contour(wks,rad_v_rad,res_rad)
draw(plot)
frame(wks)

delete(plot)
delete(wks)
delete(res_rad)



end
