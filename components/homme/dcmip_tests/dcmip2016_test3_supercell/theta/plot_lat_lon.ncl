load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin


  output_format                  = "pdf"
  output_format@wkPaperHeightF   = 8
  output_format@wkPaperWidthF    = 11
  output_format@wkOrientation    = "portrait"
  wks  = gsn_open_wks(output_format,"output_latlon.pdf")   ; send graphics to PNG file
  gsn_merge_colormaps(wks,"temp_diff_18lev","sunshine_9lev")

  filename= "./movies/dcmip2016_test31.nc"
  f       = addfile(filename,"r")
  ntimes  = dimsizes(f->time)
  lev     = f->lev
  lon     = f->lon

do time=0,ntimes-1,1

  print("time="+time+" of "+ (ntimes-1))

  geo  = f->geo(time,::-1,:,:)
  var  = f->w (time,::-1,:,:)
  var0 = f->w (0   ,::-1,:,:)

  if(lon(0) .eq. 0) then

    geo  = lonFlip(geo )
    var  = lonFlip(var)
    var0 = lonFlip(var0)

  end if

  ;_________________________________________________
  ; Interpolate fields from eta coord to z coord

  z_eta = geo/9.81/1000.0
  z     = lev
  z     = ispan(250,19750,500)/1000.0

  var  = linint1_n_Wrap(z_eta,var,False,z,0,0)
  var0 = linint1_n_Wrap(z_eta,var0,False,z,0,0)

  slice = var(9,:,:)

  res                     = True          ; plot mods desired
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "Portrait"

  res@gsnSpreadColorStart = 8
  res@gsnSpreadColorEnd = 19

;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF = 1.0
;  res@cnMinLevelValF = -3
;  res@cnMaxLevelValF = 3

;  res@trYMinF = -40.0
;  res@trYMaxF =  40.0
;  res@trXMinF = -20.0
;  res@trXMaxF =  20.0

; res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF = 2.0
;  res@cnMinLevelValF = -8.0
;  res@cnMaxLevelValF = 24.0

  res@cnFillOn            = True          ; turn on color
  res@cnFillPalette       = "gui_default" ; set color map
  res@vpWidthF            = 0.5           ; change aspect ratio of plot
  res@vpHeightF           = 0.3
  res@gsnMaximize         = True          ; maximize plot size

  plot = gsn_csm_contour(wks,slice,res)       ; contour the variable

end do

end
