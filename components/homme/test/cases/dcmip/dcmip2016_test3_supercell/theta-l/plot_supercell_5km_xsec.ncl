;***********************************************
; plot_supercell_5km_xsec.ncl
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
;     w_5km(time, lon, lat)
;     2D vertical velocity field at 5km altitude (units m/s)
;
;   Output:
;     qr_5km(time, lon, lat)
;     2D rain water mixing ratio at 5km altitude (units g/kg)
;
;   Dimensions:
;     time:  4 elements at 30, 60, 90 and 120 minutes
;     lon:   arbitrary
;     lat:   arbitrary
;----------------------------------------------------------------------------

  data_file = "./movies/dcmip2016_test31.nc"
  data      = addfile(data_file,"r")

  time = data->time
  nt = dimsizes(time)
  print("nt="+nt)

  min_per_day = 60*24

  t1  = closest_val( 30.0/min_per_day, time)
  t2  = closest_val( 60.0/min_per_day, time)
  t3  = closest_val( 90.0/min_per_day, time)
  t4  = closest_val(120.0/min_per_day, time)

  times  = (/t1,t2,t3,t4/)
  times = where(times .gt. (nt-1), 0*times+(nt-1), times)
  print("times="+times)
  print(time(times)*min_per_day)

  zm = data->geo(times,::-1,:,:)/9.80616
  w  = data->w  (times,::-1,:,:)
  qr = data->Q3 (times,::-1,:,:)

  w_interp = linint1_n_Wrap(zm,w, False,5000.0,0,1)
  w_5km    = w_interp(:,0,:,:)

  qr_interp= linint1_n_Wrap(zm,qr,False,5000.0,0,1)
  qr_5km   = qr_interp(:,0,:,:)
  qr_5km = qr_5km*1000

  delete(qr_5km@long_name)
  delete(qr_5km@units)

  lon = data->lon
  if(lon(0) .eq. 0) then
    w_5km  = lonFlip(w_5km )
    qr_5km = lonFlip(qr_5km)
  end if

  ;------------------------------------------------
  ; Workspace
  ;------------------------------------------------
  wks = gsn_open_wks("pdf","5km_xsec")
  gsn_merge_colormaps(wks,"temp_diff_18lev","sunshine_9lev")

  ;------------------------------------------------
  ; Plot resources [options]
  ;------------------------------------------------
  plot   = new(2,graphic)
  wplot  = new(3,graphic)
  qrplot = new(3,graphic)

  ; Resources:
  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False

  res@vpKeepAspect        = True
  res@vpWidthF            = 0.4
  res@vpHeightF           = 0.8
; res@gsnSpreadColors     = True
  res@gsnMaximize         = True
  res@gsnPaperOrientation = "portrait"

  res@cnFillOn          = True
  res@cnLinesOn         = True
  res@cnLineLabelsOn    = False
  res@lbLabelAutoStride = True
  res@cnInfoLabelOn     = False
  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@tmXTLabelsOn      = False
  res@tmXTOn            = False

  res@gsnSpreadColors = True

  res@trYMinF = -40.0
  res@trYMaxF =  40.0

  res@gsnContourZeroLineThicknessF = 0.0
  res@gsnContourNegLineDashPattern = 1

; res@lbLabelBarOn = False
  res@tmYLLabelsOn = False

  res@lbOrientation         = "Vertical"      ; ncl default is vertical
  res@pmLabelBarSide        = "Right"          ; default is right

  resw = res
  resw@gsnSpreadColorStart = 8
  resw@gsnSpreadColorEnd = 19
  resw@cnLevelSelectionMode = "ManualLevels"
  resw@cnLevelSpacingF = 2.0
  resw@cnMinLevelValF = -8.0
  resw@cnMaxLevelValF = 24.0

  resqr = res
  resqr@gsnSpreadColorStart = 21
  resqr@gsnSpreadColorEnd = 30
  resqr@cnLevelSelectionMode = "ManualLevels"
  resqr@cnLevelSpacingF = 1.0
  resqr@cnMinLevelValF = 1.0
  resqr@cnMaxLevelValF = 12.0

  ;------------------------------------------------
  ; Vertical velocity plots
  ;------------------------------------------------

  resw@trXMinF = -20.0
  resw@trXMaxF =  20.0
  plot(0) = gsn_csm_contour(wks, w_5km(0,:,:), resw)

  resw@trXMinF = -19.99
  resw@trXMaxF =  20.0
  wplot(0) = gsn_csm_contour(wks, w_5km(1,:,:), resw)

  resw@trXMinF = -15.0
  resw@trXMaxF =  25.0
  wplot(1) = gsn_csm_contour(wks, w_5km(2,:,:), resw)

  resw@trXMinF = -10.0
  resw@trXMaxF =  30.0
  wplot(2) = gsn_csm_contour(wks, w_5km(3,:,:), resw)

  mergedw = gsn_attach_plots(plot(0), (/wplot(0),wplot(1),wplot(2)/), False, False)

  ;------------------------------------------------
  ; Rain water plots
  ;------------------------------------------------

  resqr@trXMinF = -20.0
  resqr@trXMaxF =  20.0
  plot(1) = gsn_csm_contour(wks, qr_5km(0,:,:), resqr)

  resqr@trXMinF = -19.99
  resqr@trXMaxF =  20.0
  qrplot(0) = gsn_csm_contour(wks, qr_5km(1,:,:), resqr)

  resqr@trXMinF = -15.0
  resqr@trXMaxF =  25.0
  qrplot(1) = gsn_csm_contour(wks, qr_5km(2,:,:), resqr)

  resqr@trXMinF = -10.0
  resqr@trXMaxF =  30.0
  qrplot(2) = gsn_csm_contour(wks, qr_5km(3,:,:), resqr)

  mergedqr = gsn_attach_plots(plot(1), (/qrplot(0),qrplot(1),qrplot(2)/), False, False)

  ;------------------------------------------------
  ; Panel plot
  ;------------------------------------------------

  pres = True
  pres@gsnPanelBottom = 0.1
  pres@gsnPanelTop = 0.9
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  gsn_panel(wks, plot, (/2,1/), pres)

end

