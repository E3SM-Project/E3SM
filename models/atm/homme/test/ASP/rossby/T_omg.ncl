load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;; Temp and Omega (vertical velocity)

begin
; fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-4-0-0-low-L26/movies/asp_RH1.nc")
;fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium-L26/movies/asp_rossby1.nc")
 fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium_high-L26/movies/asp_rossby1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  PS = f->ps
  t850 = f->T
  omg850 = f->omega 
  
; output_format = "X11"
; output_format = "ncgm"
   output_format = "pdf"
;  output_format@wkOrientation = "landscape"

  ntimes = dimsizes(f->time)
  nlev = dimsizes(f->lev)
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)


; for HOMME dycore.  
  if (isatt(f,"ne")) then
    NE=f@ne
    NPTS = f@nv
    deg=360.0 / (4*NE*(NPTS-1))  ; average grid spacing at equator
    hdeg=sprintf("%4.1f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    deg=360.0/nlon
  end if


; change the levels and or times to plot here
  firsttime=15
  lasttime=15 ; ntimes-1
  timeinc=1


  time = f->time
  wks = gsn_open_wks(output_format,"T_omg_rossby_1deg")
  print ("Writing file T_omg"+"."+output_format)

  do n=firsttime,lasttime,timeinc

    res                     = True
;   res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnSpreadColors     = True         ; Use full colormap, but start
;   res@gsnSpreadColorStart =  35           ; at color index 24.
;   res@gsnSpreadColorEnd   = -35           ; at color index 24.
    
    res@cnFillOn            = True         ; Turn on contour fill
    res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLinesOn           = True        ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@cnLevelSelectionMode = "ManualLevels"
    
    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpShapeMode          = "FreeAspect"
    
    res@tmYLLabelFont = 21    ;  21 = helvetica
    res@tmXBLabelFont = 21    ;  22 = helvetica-bold 
    res@lbLabelFont = 21
    res@lbLabelBarOn = True
    res@lbLabelAutoStride = True

; res@vpKeepAspect = True   ;;; multi-panel
; res@vpWidthF = 0.38
; res@vpHeightF = 0.19

    
    res@lbLabelBarOn = True
    res@mpCenterLonF         = 180


    ps = PS(n,:,:)
    if ( "Pa" .ne. ps@units) then
      print("converting ps to Pascals")
      ps=ps*100
      ps@units="Pa"
    end if
    
    hyam = f->hyam
    hybm = f->hybm
    pnew = (/ 850 /)
    
    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = False
    P0mb = 1000
    
  gsn_define_colormap(wks,"gui_default")

;************************************************
; Multi Panel plot
;************************************************
; plot = new(2,graphic)
; res@gsnDraw = False
; res@gsnFrame = False
; pres = False


   res@tiMainString    = "HOMME/ASP_RH, Res:"+hdeg+"deg, Day="+time(n)+"."

    var = vinth2p(t850(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("t850 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 281.46
    res@cnMaxLevelValF = 281.95
    res@cnLevelSpacingF = 0.035
    var@long_name = "temperature at 850 hPa"
    var@units = "K"

    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
;   plot(0) = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    
    var = vinth2p(omg850(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("omg850 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
;   res@cnMinLevelValF = 8000.0
;   res@cnMaxLevelValF = 9900.0
;   res@cnLevelSpacingF = 100.
    res@cnMinLevelValF = -0.02  
    res@cnMaxLevelValF = 0.02   
    res@cnLevelSpacingF = 0.0025
    var@long_name = "vertcal pressure velociy at 850 hPa"
    var@units = "Pa/s"

    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
;   plot(1) = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
;   gsn_panel(wks,plot,(/2,1/),pres)


end do 
end

