;-------------------------;
; Author: Christoph Erath ;
; Date: 27th October 2011 ;
; make a movie:           ;
; convert test.ps test.gif
; convert -delay 100 -size 1024x1024 -page +1+2 test.gif -loop 0 animation.gif
;________________________;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ;use this order only 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   ;of loading

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  fnames =  "cslam_boomerang1.nc"  
  output_format = "ps"                ; change for different outputs
  wks = gsn_open_wks(output_format,"cslam_boomerang1")  ; Open an X11 window


  f = addfile(fnames,"r") 
  time=f->time
  printVarSummary(time)

  lat1d=doubletofloat(f->phys_lat);*90./acos(0.)
  lon1d=doubletofloat(f->phys_lon);*90./acos(0.)
  
  ntimes = dimsizes(f->time);

  nxy = dimsizes(lon1d);
  print ("ntimes="+ntimes+" nxy="+nxy)

  firsttime=1
  lasttime=ntimes-1
  timeinc=1
  
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")  ; Change color map
  wks@gsnPaperOrientation= "landscape" 
    ;************************************************
    ; add time stamp to plot
    ;************************************************
 ;     infoTimeStamp(wks,0.025,fnames)
;frame (wks)  


  res                      = True            ; Set some plot options

 ; res@gsnMaximize          = True            ; Maximize plot in frame
  res@tiMainString = "What the hell are you doing here?"
  res@gsnAddCyclic          = True         
  res@cnFillOn              = True              ; color Fill 
  ;res@cnFillMode           = "RasterFill"       ; Raster Mode
  res@cnLinesOn             =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"  ; Set contour levels manually
  ;res@cnMinLevelValF       = 195             ; Minimum contour level
  ;res@cnMaxLevelValF       = 328             ; Maximum contour level
  ;res@cnLevelSpacingF      = 2.25            ; Contour level spacing


  res@gsnSpreadColors      = True            ; Span full color map
  res@gsnSpreadColorStart  = 14              ; Start at color index 14
  res@gsnSpreadColorEnd    = -3              ; Stop at 3rd color from end

  res@lbBoxLinesOn         = False           ; Turn off labelbar box lines
  res@lbLabelAutoStride    = True            ; Control labelbar labels.
  
  ;res@gsnDraw = False
  ;res@gsnFrame = False
  ; plotsu = new(dimsizes(time),graphic)
  ; plotsv = new(dimsizes(time),graphic)
  res@txString="Haha"
  
  res@sfXArray            = lon1d        ; Required to tell NCL where to
  res@sfYArray            = lat1d        ; overlay data on globe.
  ; plotval = new(dimsizes(time),graphic)
  
  k=1  ; first level only
    do t = firsttime,lasttime,timeinc
    val=f->c2(t,k,:)
    res@cnMinLevelValF = 0.09
    res@cnMaxLevelValF = 1
    res@cnLevelSpacingF = 0.05
    
    res@tiMainString        = "c at level = "+k+", time="+time(t)+" days"    ;
    print("plotting var. min/max="+min(val)+" "+max(val))
    plotval = gsn_csm_contour_map(wks,val,res)  ; Call the gsn_csm function for  
    print("Contour complete: "+res@tiMainString)  
  end do
  ;pan = True
  ;pan@txString = "Haha"
  ;gsn_panel(wks,plotval,(/3,1/),pan)  
end 