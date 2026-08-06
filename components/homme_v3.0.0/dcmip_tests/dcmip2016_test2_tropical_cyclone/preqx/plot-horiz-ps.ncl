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



nlat = dimsizes(f->lat)
nlon = dimsizes(f->lon)
nlev = dimsizes(f->lev)

t10  = closest_val( 10.0, f->time)
t9  = closest_val( 9.0, f->time)
t8  = closest_val( 8.0, f->time)
t0  = closest_val( 0.0, f->time)




print("Plotting...")

wks   = gsn_open_wks ("pdf","psmap")               ; send graphics to PNG file
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

res@mpOutlineOn = False
res@mpMaxLatF = 90.
res@mpMinLatF = -45.
res@mpCenterLonF = 160.
res@mpMinLonF = 160.-135.
res@mpMaxLonF = 160+135.


res@cnLevelSelectionMode = "ManualLevels"
dx=20.0
clev=20
res@cnLevelSpacingF      =  dx
res@cnMinLevelValF       =  101500-(clev+.5)*dx
res@cnMaxLevelValF       =  101500+(clev+.5)*dx

res@tiMainString="PS t="+f->time(t0)
ps = f->ps(t0,:,:)
print("ps min/max = "+min(ps)+","+max(ps))
contour(0) = gsn_csm_contour_map(wks,ps(:,:),res)  ; create the plot

res@tiMainString="PS t="+f->time(t8)
ps = f->ps(t8,:,:)
contour(1) = gsn_csm_contour_map(wks,ps(:,:),res)  ; create the plot

res@tiMainString="PS t="+f->time(t9)
ps = f->ps(t9,:,:)
contour(2) = gsn_csm_contour_map(wks,ps(:,:),res)  ; create the plot

res@tiMainString="PS t="+f->time(t10)
ps = f->ps(t10,:,:)
contour(3) = gsn_csm_contour_map(wks,ps(:,:),res)  ; create the plot

resP                     = True                ; modify the panel plot
resP@gsnMaximize         = True

gsn_panel(wks,contour,(/2,2/),resP)             ; now draw as one plot

end
