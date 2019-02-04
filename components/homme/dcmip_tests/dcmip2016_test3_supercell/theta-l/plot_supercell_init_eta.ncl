;***********************************************
; plot_supercell_init.ncl
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

 datafilename = "./movies/dcmip2016_test31.nc"


;************************************************
; Initialize workspace
;************************************************

  wks = gsn_open_wks("pdf","init_eta")

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

  res@trYMaxF = 20.0

  res_p                   = res        ; pressure plot resources
  res_pd                  = res        ; pressure plot resources
  res_u                   = res        ; zonal velocity plot resources
  res_t                   = res        ; temperature plot resources
  res_qv                  = res        ; water vapor plot resources
  res_theta               = res        ; potential temperature plot resources
  res_thetad              = res        ; potential temperature diff plot resources

;************************************************
; Panel plot
;************************************************
  plot = new(6,graphic)
  gsn_merge_colormaps(wks,"gui_default", "BlWhRe")
;  gsn_merge_colormaps(wks,"cosam", "BlRe")

  datafile = addfile(datafilename, "r")

  lev   = datafile->lev
  qv    = datafile->Q  (0,::-1,:,:)
  geo   = datafile->geo(0,::-1,:,:)
  lat   = datafile->lat

  theta = datafile->Th (0,::-1,:,:)
  t     = datafile->T  (0,::-1,:,:)
  u     = datafile->u  (0,::-1,:,:)
  p     = datafile->pnh(0,::-1,:,:)

  p  = p/100.0  ; convert to hectopascals
  qv = qv*1000.0; convert to g/kg

  delete(u@long_name)
  delete(p@long_name)
  delete(qv@long_name)
  delete(theta@long_name)

  nlat = dimsizes(theta(0,:,0))
  ilateq = nlat/2

  print("lat at equator = "+lat(ilateq)+" (should be zero)")

  ;_________________________________________________
  ; Interpolate fields from eta coord to z coord

  z_eta = geo/9.81
  z     = lev
  z     = ispan(250,19750,500)

  ;print("z="+z)
  ;print("z_eta="+z_eta(:,0,0))
  print("interpolating from geopotential levels to evenly spaced z")
  theta = linint1_n_Wrap(z_eta,theta,False,z,0,0)
  t     = linint1_n_Wrap(z_eta,t ,False,z,0,0)
  u     = linint1_n_Wrap(z_eta,u ,False,z,0,0)
  p     = linint1_n_Wrap(z_eta,p ,False,z,0,0)
  qv    = linint1_n_Wrap(z_eta,qv,False,z,0,0)


  thetad = theta
  do i=0,nlat-1
    thetad(:,i,:) = theta(:,i,:) - theta(:,ilateq,:)
  end do

  pd = p
  do i=0,nlat-1
    pd(:,i,:) = p(:,i,:) - p(:,ilateq,:)
  end do

  z_eta = geo/9.81

; ---
;  z           = lev * 20.0
;  z           = geo(:,0,0) / 9.81 / 1000.0
  p&lev       = z
  pd&lev      = z
  t&lev       = z
  theta&lev   = z
  thetad&lev  = z
  qv&lev      = z
  u&lev       = z

; ---

  res_u@cnLevelSelectionMode= "ManualLevels"
  res_u@cnLevelSpacingF   = 2.5
  res_u@cnMinLevelValF    = -15.0
  res_u@cnMaxLevelValF    = 15.0

  res_u@tiYAxisString     = "Altitude (km)"
  res_u@gsnCenterString   = "Zonal velocity"
  res_u@gsnRightString    = "(m/s)"
  plot(0) = gsn_csm_contour(wks,u(:,:,0),res_u)

  res_qv@cnLevelSelectionMode= "ManualLevels"
  res_qv@cnLevelSpacingF   = 2.0
  res_qv@cnMinLevelValF    = 2.0
  res_qv@cnMaxLevelValF    = 14.0

  res_qv@tiYAxisString     = "Altitude (km)"
  res_qv@gsnCenterString   = "Specific Humidity"
  res_qv@gsnRightString    = "(g/kg)"
  plot(1) = gsn_csm_contour(wks,qv(:,:,0),res_qv)

; ---

  res_theta@cnLevelSelectionMode= "ExplicitLevels"
  res_theta@cnLevels = (/300., 310., 320., 330., 340., 350., 360., 380., 400., 420., 440., 460., 480., 500./)

  res_theta@tiYAxisString     = "Altitude (km)"
  res_theta@gsnCenterString   = "Virtual Potential Temperature (~F5~q~F~~B~v~E~  )"
  res_theta@gsnRightString    = "(K)"
  plot(2) = gsn_csm_contour(wks,theta(:,:,0),res_theta)

; ---

  res_thetad@cnLevelSelectionMode= "ManualLevels"
  res_thetad@cnLevelSpacingF   = 0.2
  res_thetad@cnMinLevelValF    = -2.0
  res_thetad@cnMaxLevelValF    = 2.0

  res_thetad@trXMinF = 0.
  res_thetad@trXMaxF = 90.

  res_thetad@gsnSpreadColorStart = 24
  res_thetad@gsnSpreadColorEnd = 150

  res_thetad@tiYAxisString     = "Altitude (km)"
  res_thetad@gsnCenterString   = "~F5~q~F~~B~v~E~  - ~F5~q~F~~B~v,eq~E~"
  res_thetad@gsnRightString    = "(K)"
  plot(3) = gsn_csm_contour(wks,thetad(:,:,180),res_thetad)

; ---

 res_p@cnLevelSelectionMode= "ExplicitLevels"
 res_p@cnLevels = (/1.e1, 2.e1, 5.e1, 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 7.e2, 8.e2, 9.e2, 1.e3/)

  res_p@tiYAxisString     = "Altitude (km)"
  res_p@gsnCenterString   = "Pressure (p)"
  res_p@gsnRightString    = "(hPa)"
  plot(4) = gsn_csm_contour(wks,p(:,:,0),res_p)

; ---

  res_pd@cnLevelSelectionMode= "ManualLevels"
  res_pd@cnLevelSpacingF   = 0.1
  res_pd@cnMinLevelValF    = -1.0
  res_pd@cnMaxLevelValF    = 1.0

  res_pd@trXMinF = 0.
  res_pd@trXMaxF = 90.

  res_pd@gsnSpreadColorStart = 24
  res_pd@gsnSpreadColorEnd = 150

  res_pd@tiYAxisString     = "Altitude (km)"
  res_pd@gsnCenterString   = "p - p~B~eq~E~"
  res_pd@gsnRightString    = "(hPa)"
  plot(5) = gsn_csm_contour(wks,pd(:,:,180),res_pd)

;  res_t@cnLevelSelectionMode= "ManualLevels"
;  res_t@cnLevelSpacingF   = 0.1
;  res_t@cnMinLevelValF    = -1.0
;  res_t@cnMaxLevelValF    =  1.0

;  res_t@tiYAxisString     = "Altitude (km)"
;  res_t@gsnCenterString   = "Temperature"
;  res_t@gsnRightString    = "(K)"
;  plot(0) = gsn_csm_contour(wks,t(:,:,0),res_t)

; ---

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.1

  gsn_panel(wks,plot,(/3,2/),pres)

  frame(wks)


end
