load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"




;***********************************************************************;
; NCL script to plot "geop"
;***********************************************************************;
begin

fnames =  systemfunc("\ls -t movies/swtc61.nc")
 f      = addfiles(fnames,"r")       ; Open netCDF files.	
print("fnames = "+fnames)


plotsol=1
plotref=0

if (isatt(f[0],"ne")) then
   NE=f[0]@ne
   NPTS = f[0]@nv
   deg=360.0 / (4*NE*(NPTS-1))
   hdeg=sprintf("%4.2f", deg)
   print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
end if 




lat = f[0]->lat
lon = f[0]->lon

  ntimes = dimsizes(f[0]->time)
  nlat = dimsizes(f[0]->lat)
  nlon = dimsizes(f[0]->lon)
  temp = dimsizes(f[0]->geop)
  nlev=temp(1);

; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1
  timeinc=1


  field = "geop"
;  output_format = "ncgm"
  output_format = "pdf"
;  output_format = "X11"

  fcnt = dimsizes(fnames)
  print ("number of input files = " + fcnt )
   

  var = new((/nlat,nlon/),"double")
  print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

  time = f[0]->time
  wks = gsn_open_wks(output_format,"rossby")
  print ("Writing file "+field+"."+output_format)


do kf = 0,fcnt-1
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
;  gsn_reverse_colormap(wks)
; Any time you don't want a plot drawn, just set the resources:
;  res@gsnDraw  = False
;  res@gsnFrame = False
;
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


;  res@gsnAddCyclic = False                  ; This is not global data.
    res@mpProjection = "CylindricalEquidistant"
;    res@mpLimitMode = "LatLon"
    res@mpCenterLonF         = 0
;  res@mpMinLatF = 0
;  res@mpMaxLatF = 90

  ;res@vpWidthF            = 1           ; change aspect ratio of plot
  ;res@vpHeightF           = .3


  print("reading time: "+f[kf]->time(n))

  varmin = 1.e9
  varmax = -1.
  var(:,:)= f[kf]->$field$(n,k,:,:)
  copyatt( var, f[kf]->$field$(n,k,:,:) )

  res@cnMinLevelValF = 8000. 
  res@cnMaxLevelValF = 11000
  res@cnLevelSpacingF = 100.

  ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
  ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 


  if (plotsol .eq. 1) then
    res@tiMainString        = "HOMME "+hdeg+":F34:0"+":F21:  min="+sprintf("%10.2e",min(var(:,:)))+" max="+sprintf("%5.1f",max(var(:,:))) 
    
     contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
     print("k="+k+" n="+n+" Contoured: "+res@tiMainString)
   end if



end do  ; loop over time level
end do  ; loop over levels (different I.C.s)
end do  ; loop over file names



end

