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

  wks = gsn_open_wks("pdf","init_eta2")

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

  res_u                   = res        ; zonal velocity plot resources
  res_t                   = res        ; temperature plot resources
  res_td                  = res        ; temperature diff plot resources
  res_thetap              = res        ; theta perturbation plot resources

;************************************************
; Panel plot
;************************************************
  plot = new(4,graphic)
  gsn_merge_colormaps(wks,"gui_default", "BlWhRe")
;  gsn_merge_colormaps(wks,"cosam", "BlRe")

  datafile = addfile(datafilename, "r")

  lev   = datafile->lev
  qv    = datafile->Q  (0,::-1,:,:)
  geo   = datafile->geo(0,::-1,:,:)
  lat   = datafile->lat
  lon   = datafile->lon

  t     = datafile->T  (0,::-1,:,:)
  theta = datafile->Th (0,::-1,:,:)
  u     = datafile->u  (0,::-1,:,:)
  p     = datafile->pnh(0,::-1,:,:)

;  theta = datafile->Q2  (0,::-1,:,:)
;  t     = datafile->Q3  (0,::-1,:,:)
;  p     = datafile->Q4  (0,::-1,:,:)
;  u     = datafile->Q5  (0,::-1,:,:)


  delete(t@long_name)

  p     = p/100.0  ; convert to hectopascals
  qv    = qv*1000.0; convert to g/kg
  z_eta = geo/9.81


  ;_________________________________________________
  ; Interpolate var field from eta coord to z coord

  nlev = dimsizes(lev)
  z  = new((/nlev/),double)
  z = ispan(250,19750,500)

  print("z="+z)
  print("z_eta="+z_eta(:,0,0))

  t_i     = linint1_n_Wrap(z_eta,t,False,z,0,0)
  theta_i = linint1_n_Wrap(z_eta,theta,False,z,0,0)
  nlat = dimsizes(theta(0,:,0))
  ilateq = nlat/2

  td = t
  do i=0,nlat-1
    td(:,i,:) = t_i(:,i,:) - t_i(:,ilateq,:)
  end do
  copy_VarCoords(t, td)

  thetap = theta_i(:,:,0) - theta_i(:,:,180)
  copy_VarCoords(theta(:,:,0), thetap)
  print("max thetap ="+max(thetap))
  print("min thetap ="+min(thetap))

; ---
  z           = z / 1000.0
  p&lev       = z
  t&lev       = z
  td&lev      = z
  theta&lev   = z
  qv&lev      = z
  u&lev       = z
  thetap&lev  = z

; ---

  res_t@cnLevelSelectionMode= "ManualLevels"
  res_t@cnLevelSpacingF   = 10.0
  res_t@cnMinLevelValF    = 200.0
  res_t@cnMaxLevelValF    = 310.0

  res_t@tiYAxisString     = "Altitude (km)"
  res_t@gsnCenterString   = "Temperature (T)"
  res_t@gsnRightString    = "(K)"
;  plot(0) = gsn_csm_contour(wks,t(:,:,0),res_t)
;plot(0) = gsn_csm_contour(wks,t(:,:,0),res_t)
plot(0) = gsn_csm_contour(wks,t_i(:,:,0),res_t)

  res_td@cnLevelSelectionMode= "ManualLevels"
  res_td@cnLevelSpacingF   = 0.2
  res_td@cnMinLevelValF    = -2.0
  res_td@cnMaxLevelValF    = 2.0

  res_td@trXMinF = 0.
  res_td@trXMaxF = 90.

  res_td@gsnSpreadColorStart = 24
  res_td@gsnSpreadColorEnd = 150

  res_td@tiYAxisString     = "Altitude (km)"
  res_td@gsnCenterString   = "T - T~B~eq~E~"
  res_td@gsnRightString    = "(K)"
;  plot(1) = gsn_csm_contour(wks,td(:,:,0),res_td)
plot(1) = gsn_csm_contour(wks,td(:,:,180),res_td)


; ---

  res_thetap@cnLevelSelectionMode= "ManualLevels"
  res_thetap@cnLevelSpacingF   = 0.2
  res_thetap@cnMinLevelValF    = 0.2
  res_thetap@cnMaxLevelValF    = 3.0

  res_thetap@gsnSpreadColorStart = 75
  res_thetap@gsnSpreadColorEnd = 150

  res_thetap@tiYAxisString     = "Altitude (km)"
  res_thetap@gsnCenterString   = "Theta perturbation"
  res_thetap@gsnRightString    = "(K)"
  plot(2) = gsn_csm_contour(wks,thetap(:,:),res_thetap)

; ---

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.1

  gsn_panel(wks,plot,(/2,2/),pres)

  frame(wks)


end
