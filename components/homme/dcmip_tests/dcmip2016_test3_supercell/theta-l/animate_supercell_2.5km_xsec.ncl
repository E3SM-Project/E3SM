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
;     w_slice(time, lon, lat)
;     2D vertical velocity field at 5km altitude (units m/s)
;
;   Output:
;     qr_slice(time, lon, lat)
;     2D rain water mixing ratio at 5km altitude (units g/kg)
;
;   Dimensions:
;     time:  4 elements at 30, 60, 90 and 120 minutes
;     lon:   arbitrary
;     lat:   arbitrary
;----------------------------------------------------------------------------

;data_file = "./movies/dcmip2016_test31.nc"
data_file="./movies/dcmip2016_test3_r50.nc"
;data_file="./movies/dcmip2016_test3_r100.nc"
;data_file="./movies/dcmip2016_test3_r400.nc"
  data      = addfile(data_file,"r")

  time = data->time
  nt = dimsizes(time)
  print("nt="+nt)

  min_per_day = 60*24

;------------------------------------------------
; Workspace
;------------------------------------------------
name  = "animation_2.5km_xsec"
wks_type                  = "pdf"
wks_type@wkPaperHeightF   = 6.5
wks_type@wkPaperWidthF    = 11
;wks_type@wkAntiAlias      = True

wks = gsn_open_wks(wks_type,name)

;gsn_merge_colormaps(wks,"temp_diff_18lev","sunshine_9lev")

altitude = 2500.0

t1  = 0 ;closest_val( 15.0/min_per_day, time)

do ti=t1,nt-1,1
  minutes =time(ti)*min_per_day

  print("ti="+ti+" time in minutes="+minutes)

  zm = data->geo(ti,::-1,:,:)/9.80616
  w  = data->w  (ti,::-1,:,:)
  qr = data->Q3 (ti,::-1,:,:)
  qc = data->Q2 (ti,::-1,:,:)

  w_interp = linint1_n_Wrap(zm,w, False,altitude,0,0)
  w_slice    = w_interp(0,:,:)

  qr_interp = linint1_n_Wrap(zm,qr,False,altitude,0,0)
  qr_slice  = qr_interp(0,:,:)
  qr_slice  = qr_slice*1000

  qc_interp = linint1_n_Wrap(zm,qc,False,altitude,0,0)
  qc_slice  = qc_interp(0,:,:)
  qc_slice  = qc_slice*1000

  delete(qr_slice@long_name)
  delete(qr_slice@units)
delete(qc_slice@long_name)
delete(qc_slice@units)

  lon = data->lon
  if(lon(0) .eq. 0) then
    w_slice  = lonFlip(w_slice )
    qr_slice = lonFlip(qr_slice)
  end if


  ;------------------------------------------------
  ; Plot resources [options]
  ;------------------------------------------------
  plot   = new(3,graphic)
;wplot  = new(3,graphic)
;  qrplot = new(3,graphic)

  ; Resources:
  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@gsnMaximize         = False
  res@gsnPaperOrientation = "auto"

  res@cnConstFLabelOn   = False
  res@cnFillOn          = True

  res@cnLinesOn         = False
  res@cnLineLabelsOn    = False
  res@lbLabelAutoStride = True
  res@cnInfoLabelOn     = False
  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@trYMinF         = -35.0
  res@trYMaxF         =  35.0
  res@vpHeightF       = .35
  res@vpWidthF        = .25
  res@trXMinF         = -25.0
  res@trXMaxF         =  25.0
  res@txFont   = "helvetica-bold"
  res@gsnContourZeroLineThicknessF = 0.0
  res@gsnContourNegLineDashPattern = 1

  res@lbLabelBarOn          = True
  res@lbOrientation         = "Horizontal"      ; ncl default is vertical
  res@pmLabelBarSide        = "Bottom"          ; default is right

  resw = res
  resw@cnLevelSelectionMode   = "ManualLevels"
  resw@cnLevelSpacingF        = 2
  resw@cnMinLevelValF         = -24 ;8.0
  resw@cnMaxLevelValF         = 24.0
  resw@tiYAxisString          = "latitude"
  resw@gsnCenterString        = "vertical velocity (m/s) "+sprintf("%3.0f", minutes)+" min"
  resw@cnFillPalette          = "temp_diff_18lev"; "WhBlGrYeRe"; "BlWhRe";  "BlueWhiteOrangeRed"; "BlueRed"; "cmp_b2r";

  resqr = res
  resqr@cnLevelSelectionMode  = "ManualLevels"
  resqr@cnLevelSpacingF       = 1.0
  resqr@cnMinLevelValF        = 1.0
  resqr@cnMaxLevelValF        = 10; 12.0
  resqr@tiYAxisString         = "latitude"
  resqr@gsnCenterString       = "rain water (g/kg) "+sprintf("%3.0f", minutes)+" min"
  resqr@cnFillPalette          = "sunshine_9lev"; "WhBlGrYeRe"; "BlWhRe";  "BlueWhiteOrangeRed"; "BlueRed"; "cmp_b2r";

  resqc = res
  resqc@cnLevelSelectionMode  = "ManualLevels"
  resqc@cnLevelSpacingF       = 0.1
  resqc@cnMinLevelValF        = 0.0
  resqc@tiYAxisString         = "latitude"
  resqc@gsnCenterString       = "cloud water (g/kg) "+sprintf("%3.0f", minutes)+" min"

;resqc@cnFillPalette         ="WhiteBlue"
;  resqc@cnMaxLevelValF        = 2

  resqc@cnFillPalette         ="radar"
  resqc@cnMaxLevelValF        = 5

  ;------------------------------------------------
  ; Vertical velocity plot
  ;------------------------------------------------
  plot(0) = gsn_csm_contour(wks, qc_slice(:,:), resqc)

  plot(1) = gsn_csm_contour(wks, w_slice(:,:), resw)

  ;------------------------------------------------
  ; Rain water plot
  ;------------------------------------------------
  plot(2) = gsn_csm_contour(wks, qr_slice(:,:), resqr)


  ;------------------------------------------------
  ; Panel plot
  ;------------------------------------------------

  pres                       = True
  pres@gsnMaximize           = True
  pres@gsnPaperOrientation = "portrait"

  gsn_panel(wks, plot, (/1,3/), pres)

end do
delete(data)

;system("convert -density 100 -strip animation_2.5km_xsec.pdf tempimage-%02d.png")
;system("convert -layers OptimizePlus -delay 6 tempimage-*.png -loop 1 anim.gif")
;system("rm tempimage-*.png")

end

