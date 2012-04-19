load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;  Surface Pressure and GP height at 500mb

begin
; fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-4-0-0-low-L26/movies/asp_RH1.nc")
;fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium-L26/movies/asp_rossby1.nc")
 fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium_high-L26/movies/asp_rossby1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  PS = f->ps
  t850 = f->T
  z500 = f->geo
; omg850 = f->omega 
  z500 = z500/9.80616
  
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
  wks = gsn_open_wks(output_format,"Ps_Z500_rossby_1deg")
  print ("Writing file UV"+"."+output_format)

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

    
;   thick=4
;   res@tmXBMajorThicknessF = thick
;   res@tmYLMajorThicknessF = thick
;   res@tmXBMinorThicknessF = thick 
;   res@tmYLMinorThicknessF = thick 
;   res@tmBorderThicknessF = thick
;  res@tmXBLabelFontHeightF = .015
;
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


   res@tiMainString    = "HOMME/ASP_RH, Res:"+hdeg+"deg, Day="+time(n)+"."

;   var = vinth2p(t850(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
;   print ("t850 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
;   res@cnMinLevelValF = 281.4
;   res@cnMaxLevelValF = 282.0
;   res@cnLevelSpacingF = 0.04
;   var@long_name = "temperature at 850 hPa"
;   var@units = "K"
;   contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

   PS = PS / 100.
    print ("Ps min/max = "+min(PS(n,:,:))+", "+max(PS(n,:,:)) )
    res@cnMinLevelValF = 955.
    res@cnMaxLevelValF = 1025.0
    res@cnLevelSpacingF = 5.0
    PS@long_name = "surface pressure"
    PS@units = "hPa"
    contour = gsn_csm_contour_map(wks,PS(n,:,:),res)     ; Draw a contour plot.  ;



    
    pnew = (/ 500 /)
    var = vinth2p(z500(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("Z min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 5200.0
    res@cnMaxLevelValF = 5740.0
    res@cnLevelSpacingF = 40.
    var@long_name = "gepotential height at 500 hPa"
    var@units = "m"

    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

end do 
end

