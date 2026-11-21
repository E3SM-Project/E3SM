load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"



begin
;  fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-3-0-56-medium_high-L60/movies/*.nc")
;  fnames =  systemfunc("ls ~/scratch1/preqx/HOMME-3-0-56-medium-L60/movies/*.nc")
  fnames =  systemfunc("ls movies/*.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	
  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics

  
  output_format = "ncgm"
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
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    deg=360.0/nlon
  end if


; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1
  timeinc=1

  time = f->time
  wks = gsn_open_wks(output_format,"q")
  print ("Writing file q"+"."+output_format)
  gsn_define_colormap(wks,"BlWhRe")

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
    res@lbLabelBarOn = True
    res@trYReverse = True

    deg_plot = 270   + 360*time(n)/12
    if (deg_plot .gt. 360) then
      deg_plot = deg_plot - 360
    end if
    klon=doubletoint( (deg_plot*nlon)/360 )
    if (klon .ge. nlon ) then
        klon=nlon-1
    end if
    
    q1 = f->Q(n,:,:,klon)
    q2 = f->Q2(n,:,:,klon)

    print ("time="+time(n)+" q1 min/max = "+min(q1)+", "+max(q1) )
    res@tiMainString        = "Time="+time(n)
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = 0.10

    ;res@sfYArray = (/f->lev/)
    ;res@trYAxisType = "LinearAxis"
    ;res@trYLog = True
    contour = gsn_csm_contour(wks,q1,res)     ; Draw a contour plot.  ;

;    q1b = f->Q(n,:,nlat/2+1,:)
;    contour = gsn_csm_contour(wks,q1b,res)     ; Draw a contour plot.  ;

;    q1c = f->Q(n,37,:,:)
;    contour = gsn_csm_contour(wks,q1c,res)     ; Draw a contour plot.  ;


    res@tiMainString        = "Time="+time(n)
    res@cnMinLevelValF =0.0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = 0.05
    print ("time="+time(n)+" q2 min/max = "+min(q2)+", "+max(q2) )
    contour = gsn_contour(wks,q2,res)     ; Draw a contour plot.  ;


end do 
end

