;***********************************************
; plot_supercell_prect.ncl
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
;     prect_rate(time)
;     time series precipitation rate (units m/s)
;
;     time_s(time)
;     time of each sample (units s)
;
;   Dimensions:
;     time:  arbitrary
;----------------------------------------------------------------------------

precip  = asciiread("./measurement_prect_rate.txt", -1, "double")
time_s  = asciiread("./measurement_time.txt", -1, "double")

; convert from m/s to mm/day
precip = precip*(1000.0)*(24.0*3600)

;****************************************************************************
; DO NOT EDIT BELOW THIS LINE
;****************************************************************************

;------------------------------------------------
; Workspace
;------------------------------------------------
  wks = gsn_open_wks("pdf","max_precip")

;------------------------------------------------
; Plot resources [options]
;------------------------------------------------
; Resources:
  res                     = True

  res@gsnDraw             = False
  res@gsnFrame            = False

  res@vpKeepAspect        = True
  res@vpWidthF            = 0.6
  res@vpHeightF           = 0.4
  res@gsnMaximize         = True
  res@gsnPaperOrientation = "portrait"

  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@tmXTLabelsOn      = False
  res@tmXTOn            = False

  res@tmXBMinorOn       = True
  res@tmYLMinorOn       = True

  res@xyLineThicknessF  = 3.0
;res@xyMarkLineMode    = "MarkLines"

  resw = res

;------------------------------------------------
; Vertical velocity plot
;------------------------------------------------

  resw@tiXAxisOn = True
  resw@tiXAxisString = "time (s)"

  resw@tiYAxisOn = True
  resw@tiYAxisString = "max precipitation rate (mm/day)"

  resw@trXMinF = 0.0
  resw@trXMaxF = 7200.0
  
  resw@trYMinF = 0.0
  resw@trYMaxF = 7000.0

  resw@tmXBMode = "Manual"
  resw@tmXBTickStartF = 0.0
  resw@tmXBTickEndF = 7200.0
  resw@tmXBTickSpacingF = 1800.0

; resw@tmYLMode = "Explicit"
;  resw@tmYLTickStartF = 0.0
;  resw@tmYLTickEndF = 50.0

  plot = gsn_csm_xy(wks, time_s, precip, resw)

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

