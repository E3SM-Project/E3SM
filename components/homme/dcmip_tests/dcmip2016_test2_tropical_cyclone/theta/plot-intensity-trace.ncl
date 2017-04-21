;=======================================================================================
; This NCL 
;
; Written by Colin Zarzycki (zarzycki@ucar.edu)
; Version 0.1 (6/5/2016) - DCMIP-2016 release
;=======================================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
begin

;=======================================================================================

whichVar="WIND"          ; PS or WIND
top_to_bottom=True     ; is data top to bottom? (false for bottom to top)
model="CAM-SE"

; needs to be a catted file ordered time, lev, lat, lon (or time, lat, lon for PS)
filename="/glade/u/home/zarzycki/scratch/cam5_work/dcmip2/ne120_kessler_bryan/dcmip2.cam.h1.2000-01-01-00000.nc.nc"

;=======================================================================================

maxLat=90.0
minLat=0.0

;=======================================================================================

a = addfile(filename,"r")

; figure out what is surface level for wind speed...
if (top_to_bottom) then
  sfcLev=dimsizes(a->lev)-1                  ; dimsizes(a->lev)-1 subsets bottom level
else
  sfcLev=0
end if

if (whichVar .eq. "PS") then
  PS=a->PS(:,{minLat:maxLat},:)              ; just take out NH since we know storm is there
  PS=PS/100.                                 ; convert from Pa to hPa/mb
  PS@long_name="Surface pressure"
  PS@units="hPa"
  var_vtime = dim_min_n_Wrap(PS,(/1,2/))     ; for now, just find min over lat/lon domain at each time
else if (whichVar .eq. "WIND") then
  UBOT=a->U(:,sfcLev,{minLat:maxLat},:)
  VBOT=a->V(:,sfcLev,{minLat:maxLat},:)
  WIND=UBOT
  WIND=sqrt(UBOT^2+VBOT^2)
  WIND@long_name="Lowest level wind"
  WIND@units="m/s"
  var_vtime = dim_max_n_Wrap(WIND,(/1,2/))
else
  print("unsupported var choice... exiting...")
  exit
end if
end if

wks   = gsn_open_wks ("x11","xy")               ; send graphics to PNG file

res                  = True                     ; plot mods desired
res@gsnDraw=False
res@gsnFrame=False
res@tiMainString     = "DCMIP2 "+model+" "+whichVar+" v. time"          ; add title

res@tiXAxisString        = "Time (days)"  
res@tiYAxisString        = var_vtime@long_name+" ("+var_vtime@units+")"

res@xyLineThicknessF=5.0
res@xyLineColor="Blue"

plot  = gsn_csm_xy (wks,var_vtime&time,var_vtime,res) ; create plot
draw(plot)
frame(wks)

end
