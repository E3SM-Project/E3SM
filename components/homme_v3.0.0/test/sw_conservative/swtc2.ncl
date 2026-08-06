load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"




;***********************************************************************;
; NCL script to plot "geop"
;***********************************************************************;
begin

;fnames =  systemfunc("\ls -t movies/swtc21.nc")
fnames =  systemfunc("\ls ~/scratch1/swtc2-NE10-t240-nu34e15/movies/swtc21.nc; " +\
"\ls ~/scratch1/swtc2-NE10-t240-nu0e15/movies/swtc21.nc; " +\
"ls ~/scratch1/swtc2-serendipity-ne10-nu32e15/movies/swtc21.nc")
 f      = addfiles(fnames,"r")       ; Open netCDF files.	
print("fnames = "+fnames)


  field = "geop"
;  output_format = "ncgm"
  output_format = "pdf"
;  output_format = "X11"

  fcnt = dimsizes(fnames)
  print ("number of input files = " + fcnt )


  wks = gsn_open_wks(output_format,"h")
  print ("Writing file h"+output_format)


do kf = 0,fcnt-1

  if (isatt(f[0],"ne")) then
    NE=f[0]@ne
    NPTS = f[0]@np
    deg=360.0 / (4*NE*(NPTS-1))
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  end if 

lat = f[kf]->lat
lon = f[kf]->lon

ntimes = dimsizes(f[kf]->time)
nlat = dimsizes(f[kf]->lat)
nlon = dimsizes(f[kf]->lon)
temp = dimsizes(f[kf]->geop)
nlev=temp(1);

  print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

; change the levels and or times to plot here
  firsttime=ntimes-1
  lasttime=ntimes-1
  timeinc=1

  var= f[kf]->$field$(0,0,:,:)


; exact solution
  alpha=0  ; rotation angle radians
  csalpha = cos(alpha)
  snalpha = sin(alpha)
  rearth = 6.376D6  
  omega = 7.292D-5
  PI = 4*atan(1d0)
  u0 = 2.0D0*PI*rearth/(12.0*86400.0)
  coef = rearth*omega*u0 + (u0^2)/2.0D0
  g            = 9.80616D0 
  hexact=var
  do j=0,nlat-1
    do i=0,nlon-1
      snlat = sin(lat(j)*PI/180)
      cslat = cos(lat(j)*PI/180)
      cslon = cos(lon(i)*PI/180)
      hexact(j,i)= 2.94D4 -coef*( -cslon*cslat*snalpha + snlat*csalpha )^2
      hexact(j,i)=hexact(j,i)/g
    end do
  end do




do k=0,nlev-1
do n=firsttime,lasttime,timeinc

;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  res                     = True
 res@gsnPaperOrientation = "portrait"
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
  res@mpGridSpacingF  = 10.0
  res@cnLevelSelectionMode = "ManualLevels"

  res@mpOutlineOn          = False
  res@mpFillOn             = False
  res@mpCenterLonF         = 0
  res@mpShapeMode          = "FreeAspect"

    res@mpProjection = "CylindricalEquidistant"
;    res@mpProjection = "Orthographic"
;    res@mpCenterLatF         = 90

;    res@mpLimitMode = "LatLon"
;    res@mpCenterLonF         = 0
;  res@mpMinLatF = 0
;  res@mpMaxLatF = 90
  ;res@vpWidthF            = 1           ; change aspect ratio of plot
  ;res@vpHeightF           = .3


  print("reading time: "+f[kf]->time(n))
  varmin = 1.e9
  varmax = -1.
  var(:,:)= f[kf]->$field$(n,k,:,:)
  copyatt( var, f[kf]->$field$(n,k,:,:) )
  copyatt( hexact, f[kf]->$field$(n,k,:,:) )

  res@cnMinLevelValF = 0000.
  res@cnMaxLevelValF = 4000.
  res@cnLevelSpacingF = 100.

  ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
  ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 


;    print("plotting min/max: "+min(var(:,:))+" "+max(var(:,:)))
     res@tiMainString        = "HOMME "+hdeg+"~F34~0"+"~F21~  min="+sprintf("%10.2e",min(var(:,:)))+" max="+sprintf("%5.1f",max(var(:,:))) 
     contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;


  res@cnMinLevelValF =  -2.0
  res@cnMaxLevelValF = 2.0
  res@cnLevelSpacingF =  .05

  res@cnLinesOn           = False        ; Turn off contour lines

  var=var-hexact
  res@tiMainString        = "DAY"+f[kf]->time(n)+" ERROR min="+sprintf("%10.2e",min(var))+" max="+sprintf("%10.2e",max(var)) 
  contour = gsn_csm_contour_map(wks,var,res)     ; Draw a contour plot.  ;
  print("k="+k+" n="+n+" Contoured: "+res@tiMainString)
  

end do  ; loop over time level
end do  ; loop over levels (different I.C.s)
end do  ; loop over file names



end

