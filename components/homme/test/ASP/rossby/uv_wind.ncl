load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Two plots per page (U & V)

begin
;fnames =  systemfunc("ls /ptmp/rnair/preqx/HOMME-5-0-0-medium-L26/movies/asp_mountain1.nc")
 fnames =  systemfunc("ls /ptmp/rnair/rossby/HOMME-4-0-0-medium_high-L26/movies/asp_rossby1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  PS = f->ps
  U = f->u
  V = f->v

  
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
  firsttime=0
  lasttime=30 ; ntimes-1
  timeinc=15


  time = f->time
  wks = gsn_open_wks(output_format,"UV_Rossby_1deg")
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
  res@vpKeepAspect = True
  res@vpWidthF = 0.38
  res@vpHeightF = 0.19

  plot = new(2,graphic)
  res@gsnDraw = False
  res@gsnFrame = False
  pres = False



   res@tiMainString    = "HOMME/ASP_Rossby, Res:"+hdeg+"deg, Day="+time(n)+"."

    var = vinth2p(U(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("U min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 26.
    res@cnLevelSpacingF = 2.
    var@long_name = "zonal wind at 850 hPa"
    var@units = "m/s"

;   contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    plot(0) = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

    
    var = vinth2p(V(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("V min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = -16.0
    res@cnMaxLevelValF = 16.
    res@cnLevelSpacingF = 2.
    var@long_name = "meridional wind at 850 hPa"
    var@units = "m/s"

;   contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    plot(1) = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

    gsn_panel(wks,plot,(/2,1/),pres)  ;2in1 plots


end do 
end

