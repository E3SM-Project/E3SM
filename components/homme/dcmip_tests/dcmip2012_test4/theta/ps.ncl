load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


;***********************************************************************;
; NCL script to plot PS
;***********************************************************************;
begin

fnames =  systemfunc("ls movies/dcmip2012_test41.nc")
;fnames =  systemfunc("ls -r x*/movies/dcmip2012_test41.nc")
 f      = addfiles(fnames,"r")       ; Open netCDF files.	


 NE=f[0]@ne
 NPTS = f[0]@np
 deg=360.0 / (4*NE*(NPTS-1))
 hdeg=sprintf("%4.2f", deg)
 print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)

 ntimes = dimsizes(f[0]->time)
 nlev = dimsizes(f[0]->lev)
 nlat = dimsizes(f[0]->lat)
 nlon = dimsizes(f[0]->lon)
 print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

; assume data is output 1 per eath day, and day 15 is last day in file
; day 9, 10, 11, 12, 13, 14, 15 
; output day 9, day 15
  timeinc=6
  lasttime=ntimes-1
  firsttime=lasttime-timeinc




  field = "ps"
  if ( isfilevar(f[0],"PS")) then
    field = "PS"
  end if 

  lat = f[0]->lat
  lon = f[0]->lon



  output_format = "pdf"

  fcnt = dimsizes(fnames)
  print ("number of input files = " + fcnt )

  var = new((/nlat,nlon/),"double")
  diffl2 = new((/fcnt,lasttime+1/),"double")
  diffl2=0

  print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

  time = f[0]->time
  wks = gsn_open_wks(output_format,field)
  print ("Writing file "+field+"."+output_format)


do n=firsttime,lasttime,timeinc
do kf = 0,fcnt-1
   time = f[kf]->time

;  compute scaling.  ASSUME first output in file is day 7
   xscale = 7d0/time(0)
   print("time = "+time(n)+" Earth days="+time(n)*xscale)
   time=time*xscale
   timestr=sprintf("%4.2f", time(n))
   xscalestr = sprintf("%.0f",xscale)


;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
   gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)

  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 20           ; at color index 24.


  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  res@mpGridAndLimbOn  = False
  res@mpGridSpacingF  = 30.0
  res@cnLevelSelectionMode = "ManualLevels"

  res@mpOutlineOn          = False
  res@mpFillOn             = False
  res@mpCenterLonF         = 0
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
  
  res@gsnMajorLonSpacing = 90.0
  res@gsnMajorLatSpacing = 90.0
  res@gsnMinorLonSpacing = 30.0
  res@gsnMinorLatSpacing = 30.0

;  res@mpProjection = "CylindricalEquidistant"
;  res@mpLimitMode = "LatLon"
 res@mpCenterLonF         = 180

  var(:,:)= f[kf]->$field$(n,:,:)
  copyatt( var, f[kf]->$field$(n,:,:) )

  if ( "Pa" .eq. var@units) then
     var=var/100
     var@units="mb"
  end if

  res@cnMinLevelValF = 940.
  res@cnMaxLevelValF = 1020.
  res@cnLevelSpacingF = 10.0

  res@tiMainString        = "X="+xscalestr+"  "+timestr+" days"

  ;gsn_define_colormap(wks,"rainbow")                ; choose colormap
  ;gsn_define_colormap(wks,"wh-bl-gr-ye-re")
  ;gsn_define_colormap(wks,"WhViBlGrYeOrRe")
  gsn_define_colormap(wks,"cosam")
  print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) + "  time="+timestr+ " Xscale="+xscalestr  )
  contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
  print("Contour complete: "+res@tiMainString)

end do  ; loop over file names
end do  ; loop over time level


end

