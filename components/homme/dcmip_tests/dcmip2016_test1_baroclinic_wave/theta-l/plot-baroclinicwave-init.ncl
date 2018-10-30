;***********************************************
; plot_baroclinicwave_init.ncl
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

;************************************************
; Input parameters
;************************************************

  datafilename = "./movies/dcmip2016_test11.nc"

;************************************************
; Initialize workspace
;************************************************

  wks = gsn_open_wks("pdf","plot_baroclinicwave_init")

;************************************************
; Plot resources [options]
;************************************************

  res                     = True
 ;res@mpCenterLonF        = 180
  res@cnFillOn            = True       ; color              [default=False]
  res@cnLinesOn           = True       ; draw contour lines [default=True]
  res@lbLabelAutoStride   = True
  res@vpKeepAspect        = True
  res@vpWidthF            = 0.38       ; user specified shape
  res@vpHeightF           = 0.19
  res@gsnSpreadColors     = True       ; use all colors
  res@gsnSpreadColorStart = 2          ; default=2
  res@gsnSpreadColorEnd   = 23        ; final color of "gui-default"
;  res@gsnSpreadColorStart = 2          ; default=2
;  res@gsnSpreadColorEnd   = 11        ; final color of "cosam"

  res@gsnDraw             = False      ; default=True
  res@gsnFrame            = False      ; default=True

  res_p                   = res        ; pressure plot resources
  res_ps                  = res        ; surface pressure plot resources
  res_u                   = res        ; zonal velocity plot resources
  res_up                  = res        ; zonal velocity pert plot resources
  res_t                   = res        ; temperature plot resources
  res_qv                  = res        ; water vapor plot resources
  res_theta               = res        ; potential temperature plot resources

;************************************************
; Panel plot
;************************************************
  plot = new(6,graphic)
  gsn_merge_colormaps(wks,"gui_default", "BlRe")
;  gsn_merge_colormaps(wks,"cosam", "BlRe")

  datafile = addfile(datafilename, "r")

  lev   = datafile->lev
  lon   = datafile->lon

  nlon = dimsizes(lon)

  t     = datafile->T(0,::-1,:,:)
  theta = datafile->Th(0,::-1,:,:)
  u     = datafile->u(0,::-1,:,:)
  qv    = datafile->Q(0,::-1,:,:)
  p     = datafile->pnh(0,::-1,:,:)
  geo   = datafile->geo(0,::-1,:,:)
  qv    = qv*1000.0
  p     = p/100.0
  ps    = datafile->ps(0,:,:)

delete(u@long_name)
delete(p@long_name)
delete(qv@long_name)
delete(theta@long_name)
delete(t@long_name)

  print("lon(0)="+lon(0))
  if(lon(0) .eq. 0) then
    print("lonFlip")
    t     = lonFlip(t)
    theta = lonFlip(theta)
    u     = lonFlip(u)
    qv    = lonFlip(qv)
    p     = lonFlip(p)
    ps    = lonFlip(ps)
    geo   = lonFlip(geo)
  end if


  up = u
  do i=0,nlon-1
    up(0,:,i) = up(0,:,i) - u(0,:,0)
  end do

  ;_________________________________________________
  ; Interpolate fields from eta coord to z coord

  z_eta = geo/9.81
  z     = lev
  z     = ispan(500,29500,1000)

  print("interpolating from geopotential levels to evenly spaced z")
  theta = linint1_n_Wrap(z_eta,theta,False,z,0,0)
  t     = linint1_n_Wrap(z_eta,t ,False,z,0,0)
  u     = linint1_n_Wrap(z_eta,u ,False,z,0,0)
  qv    = linint1_n_Wrap(z_eta,qv,False,z,0,0)
  p     = linint1_n_Wrap(z_eta,p ,False,z,0,0)

; ---
;  z         = lev * 30.0
  p&lev     = z/1000
  t&lev     = z/1000
  theta&lev = z/1000
  qv&lev    = z/1000
  u&lev     = z/1000

; ---
  res_t@cnLevelSelectionMode= "ManualLevels"
  res_t@cnLevelSpacingF   = 10.
  res_t@cnMinLevelValF    = 180.
  res_t@cnMaxLevelValF    = 310.

  res_t@tiYAxisString     = "Altitude (km)"
  res_t@gsnCenterString   = "Temperature"
  res_t@gsnRightString    = "(K)"
  plot(0) = gsn_csm_contour(wks,t(:,:,0),res_t)

  res_theta@cnLevelSelectionMode= "ExplicitLevels"
  res_theta@cnLevels = (/250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 380., 400., 450., 500., 550., 600., 650., 700./)

  res_theta@tiYAxisString     = "Altitude (km)"
  res_theta@gsnCenterString   = "Potential Temperature"
  res_theta@gsnRightString    = "(K)"
  plot(1) = gsn_csm_contour(wks,theta(:,:,0),res_theta)

; ---

  res_u@cnLevelSelectionMode= "ManualLevels"
  res_u@cnLevelSpacingF   = 2.0
  res_u@cnMinLevelValF    = 2.0
  res_u@cnMaxLevelValF    = 26.0

  res_u@tiYAxisString     = "Altitude (km)"
  res_u@gsnCenterString   = "Zonal velocity"
  res_u@gsnRightString    = "(m/s)"
  plot(2) = gsn_csm_contour(wks,u(:,:,0),res_u)

  res_qv@cnLevelSelectionMode= "ManualLevels"
  res_qv@cnLevelSpacingF   = 2.0
  res_qv@cnMinLevelValF    = 2.0
  res_qv@cnMaxLevelValF    = 18.0

  res_qv@tiYAxisString     = "Altitude (km)"
  res_qv@gsnCenterString   = "Specific Humidity"
  res_qv@gsnRightString    = "(g/kg)"
  plot(3) = gsn_csm_contour(wks,qv(:,:,0),res_qv)

; ---

  res_up@cnLevelSelectionMode= "ManualLevels"
  res_up@cnLevelSpacingF   = 0.1
  res_up@cnMinLevelValF    = 0.1
  res_up@cnMaxLevelValF    = 1.0

  res_up@trXMinF = 0.
  res_up@trXMaxF = 40.
  res_up@trYMinF = 30.
  res_up@trYMaxF = 50.

  res_up@tiYAxisString     = "Latitude"
  res_up@gsnCenterString   = "Zonal wind perturbation"
  res_up@gsnRightString    = "(m/s)"
  plot(4) = gsn_csm_contour(wks,up(0,:,:),res_up)

  res_p@cnLevelSelectionMode= "ExplicitLevels"
  res_p@cnLevels = (/1.e1, 2.e1, 5.e1, 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 7.e2, 8.e2, 9.e2, 1.e3/)

  res_p@tiYAxisString     = "Altitude (km)"
  res_p@gsnCenterString   = "Pressure"
  res_p@gsnRightString    = "(hPa)"
  plot(5) = gsn_csm_contour(wks,p(:,:,0),res_p) 

 ; res_ps@cnLevelSelectionMode= "ManualLevels"
 ; res_ps@cnLevelSpacingF   = 10.
 ; res_ps@cnMinLevelValF    = 950.
 ; res_ps@cnMaxLevelValF    = 1050.

 ; res_ps@gsnCenterString   = "Surfae Pressure"
 ; res_ps@gsnRightString    = "(hPa)"
 ; plot(5) = gsn_csm_contour(wks,ps(:,:),res_ps) 

; ---

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.1

  gsn_panel(wks,plot,(/3,2/),pres)

  frame(wks)


end
