load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
  fnames =  systemfunc("ls movies/dcmip2012_test41.nc")
  ;fnames =  systemfunc("ls -r x*/movies/dcmip2012_test41.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];     f has zeta
  
  NE=f@ne
  NPTS = f@np
  deg=360.0 / (4*NE*(NPTS-1))
  hdeg=sprintf("%4.2f", deg)
  print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  
;  output_format = "ncgm"
  output_format = "pdf"


  ntimes = dimsizes(f->time)
  nlev = dimsizes(f->lev)
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)


; assume data is output 1 per eath day, and day 15 is last day in file
; day 9, 10, 11, 12, 13, 14, 15 
; output day 9, day 15
  timeinc=6
  lasttime = ntimes-1
  firsttime=lasttime-timeinc
;  lasttime=ntimes-1 
  lasttime=firsttime


  wks = gsn_open_wks(output_format,"zeta")
  print ("Writing file zeta"+"."+output_format)

  do n=firsttime,lasttime,timeinc
  do kf = 0,fcnt-1
    f=fall[kf];     f has zeta

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
    res@mpShapeMode          = "FreeAspect"
    
    res@gsnPaperOrientation = "Portrait"
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
  res@tmXBLabelFontHeightF = .013
  res@tmYLLabelFontHeightF = .013

;    res@gsnMajorLonSpacing = 90.0
;    res@gsnMajorLatSpacing = 90.0
;    res@gsnMinorLonSpacing = 30.0
;    res@gsnMinorLatSpacing = 30.0
    

    res@mpCenterLonF         = 200
    res@mpMinLonF = 120
    res@mpMaxLonF = 270
    res@mpMinLatF = 25
    res@mpMaxLatF = 75

    res@vpHeightF           = .3


    time = f->time
;  compute scaling.  ASSUME first output in file is day 7
    xscale = 7d0/time(0)
    print("time = "+time(n)+" Earth days="+time(n)*xscale)
    time=time*xscale
    timestr=sprintf("%4.2f", time(n))
    xscalestr = sprintf("%.0f",xscale)
    

    ps = f->ps(n,:,:)
    copyatt( ps, f->ps(n,:,:))
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
    var=var*1e5/xscale
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
;  gsn_define_colormap(wks,"cosam")

    res@tiMainString        = "X="+xscalestr+"  "+timestr+" days"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
    
    
  end do 
  end do 
end

