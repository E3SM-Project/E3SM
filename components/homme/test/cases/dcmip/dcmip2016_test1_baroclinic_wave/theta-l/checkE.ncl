load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
setvalues NhlGetWorkspaceObjectId() 
   "wsMaximumSize" : 50000000
end setvalues 


;  output_format = "ncgm"
  output_format = "pdf"
;  output_format = "X11"

  f = addfile("movies/dcmip2016_test11.nc","r")


  
  wks = gsn_open_wks(output_format,"prec-vor.pdf")
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)
  nlev = dimsizes(f->lev)
  ntimes = dimsizes(f->time)
  
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)
  
  
  
  res=True
  res@cnLevelSelectionMode = "ManualLevels"
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "portrait"
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 50           ; at color index 24.
    
    
  res@cnFillOn            = True         ; Turn on contour fill
   res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
;  res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
;    res@cnFillMode          = "CellFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False     ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
;  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  
  res@mpOutlineOn          = False
  res@mpFillOn             = False
;    res@mpShapeMode          = "FreeAspect"


res@gsnMajorLatSpacing = 45
res@gsnMajorLonSpacing = 90
res@gsnMinorLatSpacing = 15
res@gsnMinorLonSpacing = 30
; res@tmXBLabelFontHeightF = .015
;  res@tmYLLabelFontHeightF = .015



    k=21
    eta = (f->hyam(k)+f->hybm(k))
    eta_mb = 1000*eta
    print("model level:  k="+k+" eta_mb="+eta_mb)

;  snapshots
    ; precip: 4e-7 meters/s
    ; convert to mm/day:  
    xfac = 1000*24*3600d0

    gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    res@lbLabelStride = 2
    res@cnMinLevelValF = 1d0
    res@cnMaxLevelValF = 37d0
    res@cnLevelSpacingF = 2d0

    field = f->precl(ntimes-1,:,:)
    field = field*xfac
    field@units = "mm/day"
    print("precl snapshot min/max ="+min(field)+" "+max(field))
    print("levels: [" + res@cnMinLevelValF +","+res@cnMaxLevelValF +"] spacing="+res@cnLevelSpacingF)
    contour = gsn_csm_contour_map(wks,field,res)   

;    gsn_define_colormap(wks,"MPL_viridis")
    gsn_define_colormap(wks,"MPL_RdBu")
    gsn_reverse_colormap(wks)         
    res@lbLabelStride = 6
    res@cnMinLevelValF =  -8d-5
    res@cnMaxLevelValF =  8d-5
    res@cnLevelSpacingF = 2*res@cnMaxLevelValF/24
    field = f->zeta(ntimes-1,k,:,:)
;    field = field*1e5
;    field@units = "1e-5 1/s"
    print("zeta snapshot min/max ="+min(field)+" "+max(field))
    print("levels: [" + res@cnMinLevelValF +","+res@cnMaxLevelValF +"] spacing="+res@cnLevelSpacingF)
    res@tiMainString=""
    contour = gsn_csm_contour_map(wks,field,res)   

    


delete(wks)
end


