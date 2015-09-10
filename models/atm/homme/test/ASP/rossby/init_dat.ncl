load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"



begin
; fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-4-0-0-low-L26/movies/asp_RH1.nc")
 fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium-L26/movies/asp_rossby1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  PS = f->ps
  U = f->u
  V = f->v

  
  output_format = "X11"
; output_format = "ncgm"
;  output_format = "pdf"
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
  firsttime=0
  lasttime=0 ; ntimes-1
  timeinc=1


  time = f->time
  wks = gsn_open_wks(output_format,"UV")
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
    
; res@lbLabelAutoStride = True
  gsn_define_colormap(wks,"gui_default")

   res@tiMainString    = "HOMME/ASP_RH, Res:"+hdeg+"deg, Day="+time(n)+"."
    PS = PS / 100.
     res@cnMinLevelValF = 955
     res@cnMaxLevelValF = 1035
     res@cnLevelSpacingF = 5
;   print ("V min/max = "+min(PS(n,:,:))+", "+max(PS(n,:,:)) )
;   PS@long_name = "Initial Surface Pressure"
;   PS@units = "hPs"
   contour = gsn_csm_contour_map(wks,PS(n,:,:),res)     ; Draw a contour plot.  ;

;  res@tiMainString    = "HOMME/ASP_RH, Res:"+hdeg+"deg, Day="+time(n)+"."

    var = vinth2p(U(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("U min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = .0
    res@cnMaxLevelValF = 26.
    res@cnLevelSpacingF = 2.
    var@long_name = "zonal wind at 850 hPa"
    var@units = "m/s"

    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    
    var = vinth2p(V(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("V min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = -16.0
    res@cnMaxLevelValF = 16.
    res@cnLevelSpacingF = 2.
    var@long_name = "meridional wind at 850 hPa"
    var@units = "m/s"

   contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;



end do 
end

