load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


;; Geopotential heights for the Rossby/Mountain test case

begin
; fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-4-0-0-low-L26/movies/asp_RH1.nc")
;fnames =  systemfunc("ls /ptmp/rnair/mountain/HOMME-5-0-0-medium-L26/movies/asp_mountain1.nc")
 fnames =  systemfunc("ls /ptmp/rnair/mountain/HOMME-5-0-0-medium_high-L26/movies/asp_mountain1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  PS = f->ps
  z700 = f->geo
; omg850 = f->omega 
  z700 = z700/9.80616
  
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
  firsttime=5
  lasttime=25 ; ntimes-1
  timeinc=10


  time = f->time
  wks = gsn_open_wks(output_format,"Zgp_mountain_1deg")
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
    pnew = (/ 700 /)
    
    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = False
    P0mb = 1000
    
  gsn_define_colormap(wks,"gui_default")
   res@tiMainString    = "HOMME/ASP_mountain, Res:"+hdeg+"deg, Day="+time(n)+"."

    var = vinth2p(z700(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("Z min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 2400.0
    res@cnMaxLevelValF = 3300.0
    res@cnLevelSpacingF = 100.0
    var@long_name = "geopotential height at 700 hPa"
    var@units = "m"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

end do 
end

