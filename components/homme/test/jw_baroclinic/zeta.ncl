load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
  fnames =  systemfunc("ls movies/jw_baroclinic?.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  fps=fall[0];   file with ps
  f=fall[0];     f has zeta
  fT=fall[1];    file with T, Q, Q2
  fQ=fall[1];    file with Q
  

; for HOMME dycore.  
  if (isatt(f,"ne")) then
    NE=f@ne
    NPTS = f@np
    deg=360.0 / (4*NE*(NPTS-1))
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    hdeg=""
  end if 
  
;  output_format = "ncgm"
  output_format = "pdf"
  output_format@wkOrientation = "landscape"


  ntimes = dimsizes(fps->time)
  nlev = dimsizes(fps->lev)
  nlat = dimsizes(fps->lat)
  nlon = dimsizes(fps->lon)
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)



; change the levels and or times to plot here
  n9 = closest_val(9,fps->time)

  firsttime=n9
  lasttime=n9
  timeinc=1


  time = f->time
  wks = gsn_open_wks(output_format,"zeta")
  print ("Writing file zeta"+"."+output_format)

  do n=firsttime,lasttime,timeinc

    res                     = True
    res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnSpreadColors     = True         ; Use full colormap, but start
    res@gsnSpreadColorStart =  35           ; at color index 24.
    res@gsnSpreadColorEnd   = -35           ; at color index 24.
    
    res@cnFillOn            = True         ; Turn on contour fill
    res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLinesOn           = True        ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@cnLevelSelectionMode = "ManualLevels"
    
    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = 0
    res@mpShapeMode          = "FreeAspect"
    
    res@tmYLLabelFont = 21    ;  21 = helvetica
    res@tmXBLabelFont = 21    ;  22 = helvetica-bold 
    res@lbLabelFont = 21
    res@lbLabelBarOn = True
    
    
    thick=4
    res@tmXBMajorThicknessF = thick
    res@tmYLMajorThicknessF = thick
    res@tmXBMinorThicknessF = thick 
    res@tmYLMinorThicknessF = thick 
    res@tmBorderThicknessF = thick
;  res@tmXBLabelFontHeightF = .015
;
;   ZOOM into region
    res@mpCenterLonF         = 200
    res@mpMinLonF = 120
    res@mpMaxLonF = 270
    res@mpMinLatF = 25
    res@mpMaxLatF = 75
;    res@mpMinLatF = -15   ; when rotated 45 degrees
;    res@mpMaxLatF = 35
    
    res@vpHeightF           = .3
    

    ps = fps->ps(n,:,:)
    copyatt( ps, fps->ps(n,:,:))
    if ( "Pa" .ne. ps@units) then
      print("converting ps to Pascals")
      ps=ps*100
      ps@units="Pa"
    end if
    print ("ps min/max = "+min(ps(:,:))+", "+max(ps(:,:)))
    hyam = f->hyam
    hybm = f->hybm
    pnew = (/ 850 /)
    
    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = False
    P0mb = 1000
    
    var = vinth2p(f->zeta(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    var=var*1e5
    print ("zeta 1e5 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    
    
    if ( time(n) .le. 7  ) then
       res@cnMinLevelValF = -2.
       res@cnMaxLevelValF = 5.
       res@cnLevelSpacingF = 1.
    else
       res@cnMinLevelValF = -5.
       res@cnMaxLevelValF = 35.
       res@cnLevelSpacingF = 5.
    end if 

;   gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;   gsn_reverse_colormap(wks)
    gsn_define_colormap(wks,"BlWhRe")

    res@tiMainString        = ":F21:HOMME "+hdeg+":F34:0"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
    
    
    if ( True ) then
      res@lbLabelBarOn = True

      ;var = vinth2p(fT->T(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      ;print ("T min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
      ;res@cnLevelSelectionMode = "EqualSpacedLevels"
      ;res@tiMainString        = "T"
      ;contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

      var = vinth2p(fT->Q(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      print ("Q min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
      res@cnMinLevelValF = .5
      res@cnMaxLevelValF = .8
      res@cnLevelSpacingF = .025
      res@tiMainString        = "Q1"
      contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

      var = vinth2p(fT->Q2(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      print ("Q min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
      res@cnMinLevelValF = 0.0
      res@cnMaxLevelValF = .6
      res@cnLevelSpacingF = .05
      res@tiMainString        = "Q2"
      contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    end if
end do 
end

