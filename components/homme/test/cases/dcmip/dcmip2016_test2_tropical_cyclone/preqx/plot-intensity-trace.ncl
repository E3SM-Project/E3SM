;=======================================================================================
; This NCL 
;
; Written by Colin Zarzycki (zarzycki@ucar.edu)
; Version 0.1 (6/5/2016) - DCMIP-2016 release
;=======================================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;=======================================================================================

top_to_bottom = True       ; is data top to bottom? (false for bottom to top)
model         = "PREQX"

; needs to be a catted file ordered time, lev, lat, lon (or time, lat, lon for PS)
filename="./movies/dcmip2016_test21.nc"

;=======================================================================================

maxLat=90.0
minLat=0.0

;=======================================================================================

a = addfile(filename,"r")
sfcLev=dimsizes(a->lev)-1                  ; dimsizes(a->lev)-1 subsets bottom level
PS = a->ps(:,{minLat:maxLat},:)              ; just take out NH since we know storm is there
PS = PS/100.                                 ; convert from Pa to hPa/mb
PS@long_name="Surface pressure"
PS@units="hPa"
var_vtime = dim_min_n_Wrap(PS,(/1,2/))     ; for now, just find min over lat/lon domain at each time

wks = gsn_open_wks ("pdf","ps")               ; send graphics to PNG file

res                  = True                     ; plot mods desired
res@gsnDraw          = False
res@gsnFrame         = False
res@tiMainString     = "DCMIP2 "+model+" PS v. time"          ; add title
res@tiXAxisString    = "Time (days)"
res@tiYAxisString    = var_vtime@long_name+" ("+var_vtime@units+")"
res@xyLineThicknessF = 5.0
res@xyLineColor      ="Red"

plot  = gsn_csm_xy (wks,var_vtime&time,var_vtime,res) ; create plot
draw(plot)
frame(wks)

UBOT=a->u(:,sfcLev,{minLat:maxLat},:)
VBOT=a->v(:,sfcLev,{minLat:maxLat},:)
WIND=UBOT
WIND=sqrt(UBOT^2+VBOT^2)
WIND@long_name="Lowest level wind"
WIND@units="m/s"
var_vtime = dim_max_n_Wrap(WIND,(/1,2/))

wks = gsn_open_wks ("pdf","wind")               ; send graphics to PNG file

res@tiMainString     = "DCMIP2 "+model+" WIND vs time"          ; add title
res@tiXAxisString    = "Time (days)"
res@tiYAxisString    = var_vtime@long_name+" ("+var_vtime@units+")"
res@xyLineThicknessF = 5.0
res@xyLineColor      ="Blue"

plot  = gsn_csm_xy (wks,var_vtime&time,var_vtime,res) ; create plot
draw(plot)
frame(wks)

end
