;***********************************************
; plot_tropicalcyclone_init.ncl
;***********************************************
;  https://github.com/ClimateGlobalChange/DCMIP2016/

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

;************************************************
; Input parameters
;************************************************

datafilename = "./movies/dcmip2016_test21.nc"

;************************************************
; Initialize workspace
;************************************************

  wks = gsn_open_wks("pdf","init.pdf")

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
res@gsnYAxisIrregular2Linear = True

  res_p                   = res        ; pressure plot resources
  res_ps                  = res        ; surface pressure plot resources
  res_u                   = res        ; zonal velocity plot resources
  res_pp                  = res        ; zonal velocity pert plot resources
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
  u     = datafile->u  (0,::-1,:,:)
  t     = datafile->T  (0,::-1,:,:)
  theta = datafile->Th (0,::-1,:,:)
  qv    = datafile->Q  (0,::-1,:,:)
  p     = datafile->pnh(0,::-1,:,:)
  geo   = datafile->geo(0,::-1,:,:)
  ps    = datafile->ps (0,::-1,:  )

  nlon = dimsizes(datafile->lon)
  ilon0   = 0          ; location of perturbation
  ilon1   = nlon/2     ; +180 degrees

; original dcmip 2016 NCL script plots:
;  T, theta,qv  at ilon1     
;  u,p          at ilon0
;  pp  perturbation from ilon0 and ilon1
;  but then we dont see the purturbation in T and theta, so
;  we plot at ilon0
;

  print("max u="+max(u))
  print("max u xsection="+max(u(:,:,ilon0))+" lon="+lon(ilon0))

  ;_________________________________________________
  ; Interpolate fields from eta coord to z coord

  z_eta = geo/9.81/1000
  z     = lev
  z     = 1.0d-3*ispan(500,29550,1000)

  ;print("z="+z)
  print("z_eta="+z_eta(:,0,0))
  print("interpolating from geopotential levels to evenly spaced z")
  theta = linint1_n_Wrap(z_eta,theta,False,z,0,0)
  t     = linint1_n_Wrap(z_eta,t ,False,z,0,0)
  u     = linint1_n_Wrap(z_eta,u ,False,z,0,0)
  p     = linint1_n_Wrap(z_eta,p ,False,z,0,0)
  qv    = linint1_n_Wrap(z_eta,qv,False,z,0,0)

  delete(u@long_name)
  delete(p@long_name)
  delete(qv@long_name)
  delete(theta@long_name)
  delete(t@long_name)

  p  = p/100
  qv = qv*1000
  pp = p(:,:,ilon0)-p(:,:,ilon1)

  copy_VarCoords(p(:,:,0), pp)
;print("pp="+pp)

  nlon = dimsizes(u(0,0,:))
  up = u
  do i=0,nlon-1
    up(0,:,i) = up(0,:,i) - u(0,:,0)
  end do

; lev   = datafile->lev
;  rho   = datafile->Rho
;  t     = datafile->T
;  theta = datafile->Theta
;  u     = datafile->U
;  qv    = datafile->RhoQv / rho * 1000.0
;  copy_VarCoords(rho, qv)
;  p     = rho * 287.0 * t / 100.0
;  copy_VarCoords(rho, p)

;  pp    = p(0,:,:,0) - p(0,:,:,ilon0)
;  copy_VarCoords(p(0,:,:,0), pp)
; ps    = (p(0,0,:,:) * lev(1) - p(0,1,:,:) * lev(0)) / (lev(1) - lev(0))
;  copy_VarCoords(rho(0,0,:,:), ps)
;  nlon = dimsizes(u(0,0,0,:))
;  up = u
;  do i=0,nlon-1
;    up(0,0,:,i) = up(0,0,:,i) - u(0,0,:,0)
;  end do

; ---
;  z         = lev; * 30.0
;z         = geo(:,0,0) / 9.81 / 1000.0
;print("z="+z)
;  p&lev     = z
;  t&lev     = z
;  theta&lev = z
;  qv&lev    = z
;  u&lev     = z
;  pp&lev    = z



; ---
  res_t@cnLevelSelectionMode= "ManualLevels"
  res_t@cnLevelSpacingF   = 10.
  res_t@cnMinLevelValF    = 200.
  res_t@cnMaxLevelValF    = 310.

  res_t@tiYAxisString     = "Altitude (km)"
  res_t@gsnCenterString   = "Temperature"
  res_t@gsnRightString    = "(K)"
  plot(0) = gsn_csm_contour(wks,t(:,:,ilon0),res_t)

  res_theta@cnLevelSelectionMode= "ExplicitLevels"
  res_theta@cnLevels = (/250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 380., 400., 450., 500., 550., 600., 650., 700./)

  res_theta@tiYAxisString     = "Altitude (km)"
  res_theta@gsnCenterString   = "Potential Temperature"
  res_theta@gsnRightString    = "(K)"
  plot(1) = gsn_csm_contour(wks,theta(:,:,ilon0),res_theta)

; ---

  res_u@cnLevelSelectionMode= "ManualLevels"
  res_u@cnLevelSpacingF   = 2.0
  res_u@cnMinLevelValF    = -26.0
  res_u@cnMaxLevelValF    = 26.0

  res_u@trXMinF = 0.
  res_u@trXMaxF = 20.
  res_u@tiYAxisString     = "Altitude (km)"
  res_u@gsnCenterString   = "Zonal velocity"
  res_u@gsnRightString    = "(m/s)"
  plot(2) = gsn_csm_contour(wks,u(:,:,ilon0),res_u)

  res_qv@cnLevelSelectionMode= "ManualLevels"
  res_qv@cnLevelSpacingF   = 2.0
  res_qv@cnMinLevelValF    = 2.0
  res_qv@cnMaxLevelValF    = 18.0

  res_qv@tiYAxisString     = "Altitude (km)"
  res_qv@gsnCenterString   = "Specific Humidity"
  res_qv@gsnRightString    = "(g/kg)"
  plot(3) = gsn_csm_contour(wks,qv(:,:,ilon0),res_qv)

; ---

  res_pp@cnLevelSelectionMode= "ManualLevels"
  res_pp@cnLevelSpacingF   = 1.0
  res_pp@cnMinLevelValF    = -10.0
  res_pp@cnMaxLevelValF    = -1.0

  res_pp@trXMinF = 0.
  res_pp@trXMaxF = 20.

  res_pp@tiYAxisString     = "Altitude (km)"
  res_pp@gsnCenterString   = "Pressure perturbation"
  res_pp@gsnRightString    = "(hPa)"
  plot(4) = gsn_csm_contour(wks,pp(:,:),res_pp) 

 res_p@cnLevelSelectionMode= "ExplicitLevels"
 res_p@cnLevels = (/1.e1, 2.e1, 5.e1, 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 7.e2, 8.e2, 9.e2, 1.e3/)

  res_p@tiYAxisString     = "Altitude (km)"
  res_p@gsnCenterString   = "Background pressure"
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
