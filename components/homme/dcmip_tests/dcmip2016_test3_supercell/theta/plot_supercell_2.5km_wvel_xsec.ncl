;***********************************************
; plot_supercell_2.5km_wvel_xsec.ncl
;
; Version 1.0
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

;****************************************************************************
; DO NOT EDIT ABOVE THIS LINE
;****************************************************************************
; Read in data and store as follows:
;
;   Output:
;     w_2500m(lon, lat)
;     2D vertical velocity field at 2.5km altitude at 7200s (units m/s)
;
;   Dimensions:
;     lon:   arbitrary
;     lat:   arbitrary
;----------------------------------------------------------------------------

  data_file = "./movies/dcmip2016_test31.nc"
  data      = addfile(data_file,"r")

  time = data->time
  nt = dimsizes(time)
  print("nt="+nt)
  min_per_day = 60.0*24
  t1 = closest_val( 120.0/min_per_day, time)
  t1 = where(t1 .gt. (nt-1), nt-1, t1)
  print("t1="+t1)

  zm = data->geo(t1,::-1,:,:)/9.80616
  w  = data->w  (t1,::-1,:,:)

  w_interp = linint1_n_Wrap(zm,w, False,2500.0,0,0)
  w_2500m  = w_interp(0,:,:)

  lon = data->lon
  if(lon(0) .eq. 0) then
    w_2500m  = lonFlip(w_2500m )
  end if

;****************************************************************************
; DO NOT EDIT BELOW THIS LINE
;****************************************************************************

;------------------------------------------------
; Workspace
;------------------------------------------------
  wks = gsn_open_wks("pdf","2.5km_wvel_xsec")

  gsn_merge_colormaps(wks,"temp_diff_18lev","sunshine_9lev")

;------------------------------------------------
; Plot resources [options]
;------------------------------------------------
  plot = new(1,graphic)

; Resources:
  res                     = True

  res@gsnDraw             = False
  res@gsnFrame            = False

  res@vpKeepAspect        = True
  res@vpWidthF            = 0.5
  res@vpHeightF           = 0.6
  res@gsnMaximize         = True
  res@gsnPaperOrientation = "portrait"

  res@cnFillOn          = True
  res@cnLinesOn         = True
  res@cnLineLabelsOn    = False
  res@lbLabelAutoStride = True
  res@cnInfoLabelOn     = False

  res@trXMinF           = -45.0
  res@trXMaxF           =  30.0

  res@trYMinF           = -45.0
  res@trYMaxF           =  45.0

  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@tmXBMode          = "Explicit"
  res@tmXBValues        = (/-45,-30,-15,0,15,30/)
  res@tmXBLabels        = (/"45W","30W","15W","0","15E","30E"/)

  res@tmXBMinorOn       = True
  res@tmXBMinorValues   = (/-40,-35,-25,-20,-10,-5,5,10,20,25/)

  res@tmXTLabelsOn      = True
  res@tmXTOn            = True

  res@tmYLLabelsOn      = True
  res@tmYLOn            = True

  res@tmYLMode          = "Explicit"
  res@tmYLValues        = (/-45,-30,-15,0,15,30,45/)
  res@tmYLLabels        = (/"45S","30S","15S","0","15N","30N","45N"/)

  res@tmYLMinorOn       = True
  res@tmYLMinorValues   = (/-40,-35,-25,-20,-10,-5,5,10,20,25,35,40/)

  res@gsnSpreadColors = True

  res@gsnContourZeroLineThicknessF = 0.0
  res@gsnContourNegLineDashPattern = 1

  res@gsnSpreadColorStart = 7
  res@gsnSpreadColorEnd = 19
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 1.0
  res@cnMinLevelValF = -10.0
  res@cnMaxLevelValF = 20.0

  plot(0) = gsn_csm_contour(wks, w_2500m(:,:), res)

;------------------------------------------------
; Panel plot
;------------------------------------------------

  pres = True
  pres@gsnPanelBottom = 0.1
  pres@gsnPanelTop = 0.9
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  gsn_panel(wks, plot, (/1,1/), pres)

end

