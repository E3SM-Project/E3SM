;================================================;
; plot-slice.ncl
; David Hall (University of Colorado, Boulder)
;================================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin

  f1    = addfile("./movies/dcmip2012_test1_21.nc","r")
  lat   = f1->lat
  lon   = f1->lon
  lev   = f1->lev
  time  = f1->time

  nlat  = getfilevardimsizes(f1, "lat" )
  nlon  = getfilevardimsizes(f1, "lon" )
  nlev  = getfilevardimsizes(f1, "lev" )
  nt    = dimsizes(time)

  t_index = nt-1        ; set time to max time
  if(isvar("t")) then
    t_index = t         ; use t from command line if possible
  end if
  print ("t_index="+t_index)

  midlon  = nlon/2;-1
  q2      = f1->Q2 (t_index,:,:,midlon)
  print ("longitude = "+lon(midlon))

  ; Get height from geopotential

  g = 9.80616d0; Gravity (m s^2)
  z = f1->geo(0,:,0,midlon)/g

  ; Remove tracer units

  q2@units     = ""
  q2@long_name = ""

  res1                      = True            ; Plot resource
  res1@gsnDraw              = False          	; panel plot
  res1@gsnFrame             = False           ; don't draw yet
  res1@lbLabelBarOn         = False           ; label bar
  res1@vpWidthF             = 2               ; set aspect ratio 2:1
  res1@vpHeightF            = 1
  res1@sfXArray             = lat           	; uses lat as plot x-axis
  res1@sfYArray             = z/1000.0        ; vertical axis = z in km

  res1@cnInfoLabelOn        = False
  res1@cnFillOn             = True            ; fill contours
  res1@cnLinesOn            = True            ; draw contour lines
  res1@cnLineLabelsOn       = False           ; don't draw labels
  res1@cnMonoLineColor      = True            ; single color lines
  res1@cnLineColor          = "black"         ; black lines
  res1@cnLineThicknessF     = 0.01             ; thin lines
  res1@cnLineDrawOrder      = "PreDraw"
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = 1.0             ; max contour level
  res1@cnMinLevelValF       = 0.0             ; min contour level
  res1@cnLevelSpacingF      = 0.1             ; contour level spacing

  res1@tiYAxisString        = "height (km)"   ; y-axis title
  res1@tiXAxisString        = "latitude"      ; x-axis title
  res1@txFont               = 22              ; helvetica-bold
  res1@tiMainFont           = 22              ; helvetica-bold
  res1@tiMainFontHeightF    = 0.06            ; set font sizes
  res1@tiXAxisFontHeightF   = 0.05
  res1@tiYAxisFontHeightF   = 0.05
  res1@tmXBLabelFontHeightF = 0.04
  res1@tiXAxisOffsetYF      = 0.015           ; shift x label
  res1@tiYAxisOffsetXF      = 0.015           ; shift y label
  res1@gsnYAxisIrregular2Linear = True

  res1@gsnCenterString      = "DCMIP2012   test 1-2   1dg   L"+ (nlev) +"   t = "+time(t_index)*24.0+" hr"

  ; panel properties

  pres                      = True            ; Panel resource
  pres@gsnMaximize          = True            ; maximize panel size
  pres@gsnFrame             = False           ; don't advance yet
  pres@gsnPaperOrientation  = "portrait"      ; paper orientation
  pres@gsnPanelLabelBar     = True            ; add common label bar
  pres@lbLabelFontHeightF   =0.020            ; set labelbar font height

  ; Create workstation for file output

  name  = "dcmip2012_test1-2_L"+(nlev)+"_t"+t_index
  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  wks = gsn_open_wks(wks_type,name)
  wks@wkAntiAlias      = True

  ; Use rainbow colormap with a few of color substitutions

  ;gsn_define_colormap(wks,"gui_default")      ; choose default colormap
  ;colors      = gsn_retrieve_colormap(wks)    ; retrieve color map for editing
  ;colors(2 ,:)= (/ 0.4, 0.4, 0.4 /)           ; dark gray for undershoots
  ;colors(4 ,:)= (/ 1.0, 1.0, 1.0 /)           ; white     for tracer-free region
  ;colors(23,:)= (/ 0.6, 0.0, 0.0 /)           ; brick-red for overshoots

  ; WBGYR
  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")      ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)   ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                 ; white     for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)          ; gray for undershoots
  colors(255,:)= (/ 0.4, 0.0, 0.0 /)          ; brick-red for overshoots
  gsn_define_colormap(wks,colors)
	  
  ; Add colored contour plot

  plot    = new (1, graphic)                  ; define plot - need 4 panels
  plot(0) = gsn_csm_contour(wks,q2(:,:),res1)	; plot q2 at 12 hours

  ; Set properties for horizontal lines

  gsres                     = True            ; graphics resource
  gsres@tfPolyDrawOrder     = "PostDraw"
  gsres@gsLineThicknessF    = 0.1
  gsres@gsLineColor         = "gray50"
  z_line                    = new(nlev,graphic)
  plx                       = (/-180,180/)
  
  ; Add line every 5th vertical level

  do k = 0,nlev-1,5
    ply = (/z(k)/1000.0 ,z(k)/1000.0/)
    z_line(k)   = gsn_add_polyline (wks,plot(0),plx, ply, gsres)
  end do

  ; draw panel with contour and lines

  print("writing "+name)

	gsn_panel(wks,plot,(/1,1/),pres)
	frame(wks)

end
