load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"



begin
  fnames =  systemfunc("ls movies/asp_baroclinic1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	
  print ("files = "+fnames)
  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  
; output_format = "ncgm"
   output_format = "pdf"
;  output_format@wkOrientation = "portrait"

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
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    deg=360.0/nlon
  end if

  psname = "ps"
  if (isfilevar(f,"PS")) then
     psname = "PS"
  end if


; change the levels and or times to plot here
  firsttime=9
  lasttime=ntimes-1
  timeinc=4


  time = f->time
  wks = gsn_open_wks(output_format,"zeta")
  print ("Writing file zeta"+"."+output_format)

  do n=firsttime,lasttime,timeinc
    print ("time = "+n)
    zeta = f->zeta(n,:,:,:)
    ps = f->$psname$(n,:,:)

    psunits=stringtochar(ps@units)
    if ( "Pa" .ne. psunits(0:1) ) then
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
    extrap = True
    P0mb = 1000
    print ("ps min/max = "+min(ps(:,:))+", "+max(ps(:,:)) )
    var = vinth2p(zeta,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("zeta 850 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    var=var*1e5

    
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
    res@mpCenterLonF         = 180

;   ZOOM into region
    res@mpCenterLonF         = 200
    res@mpMinLonF = 120
    res@mpMaxLonF = 270
    res@mpMinLatF = 25
    res@mpMaxLatF = 75

    if ( time(n) .le. 7  ) then
       res@cnMinLevelValF = -2.
       res@cnMaxLevelValF = 5.
       res@cnLevelSpacingF = 1.
    else
       res@cnMinLevelValF = -5.
       res@cnMaxLevelValF = 35.
       res@cnLevelSpacingF = 5.
    end if 

    gsn_define_colormap(wks,"BlWhRe")
    res@tiMainString        = ":F21:HOMME "+hdeg+":F34:0"
    print("drawing contour: "+res@tiMainString)
    contour = gsn_csm_contour_map(wks,var(0,:,:),res) 

    

end do 
end

