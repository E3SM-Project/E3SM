
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  fname =  "~/scratch1/cslam/results/cslam_solidbody1.nc"
  f     = addfile(fname,"r")       ; Open netCDF file.	

;  output_format = "pdf"
;  output_format = "x11"
  output_format = "ncgm"

  lat1d=doubletofloat(f->phys_lat);*90./acos(0.)
  lon1d=doubletofloat(f->phys_lon);*90./acos(0.)

  ntimes = dimsizes(f->time);
  nxy = dimsizes(lon1d);

; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1
  timeinc=10


 print ("ntimes="+ntimes+" nxy="+nxy)


  time = f->time
 
  wks = gsn_open_wks(output_format,"c")

  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                     = True
  ;res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.

    res@sfXArray            = lon1d        ; Required to tell NCL where to
    res@sfYArray            = lat1d        ; overlay data on globe.

    res@cnLevelSelectionMode = "ManualLevels"

    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = -90
    res@mpShapeMode          = "FreeAspect"

;    res@mpMinLatF = 20
;    res@mpMaxLatF = 80
;    res@vpHeightF           = .3

    k=0  ; first level only
    do t = firsttime,lasttime,timeinc
      var=f->c(t,k,:)

      res@cnMinLevelValF = .9
      res@cnMaxLevelValF = 2.1
      res@cnLevelSpacingF = .02
      
      res@tiMainString        = "c at level = "+k+", time="+time(t)+" days"    ;
      print("plotting var. min/max="+min(var)+" "+max(var))
      contour = gsn_csm_contour_map(wks,var,res)     ; Draw a contour plot.
      print("Contour complete: "+res@tiMainString)
    end do
    

end

