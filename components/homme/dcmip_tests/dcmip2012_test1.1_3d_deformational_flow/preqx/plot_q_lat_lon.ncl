;================================================;
;  plot-q-lat-lon
;  plot lattitude-longitude cross section of tracer q
;  at time_index t
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

begin

  f    = addfile("./movies/dcmip2012_test1_11.nc","r")
  lat  = f->lat
  lon  = f->lon
  lev  = f->lev
  time = f->time
  nlat = getfilevardimsizes(f, "lat" )
  nlon = getfilevardimsizes(f, "lon" )
  nlev = getfilevardimsizes(f, "lev" )
  nt   = dimsizes(time)

  ; get time_index t from command line

  t_index = nt-1
  if(isvar("t")) then
    t_index = t
  end if

  ; get tracer_id q from command line

  tracer_id=1
  if(isvar("q")) then
    tracer_id = q         ; use q from command line if possible
  end if

  ; get zoom factor from command line

  zoom_factor = 1.0
  if(isvar("zoom"))
    zoom_factor = zoom
  end if

  ; find level closest to 4900 m

  g = 9.80616d0; Gravity (m s^2)
  heights = f->geo(0,:,0,0)/g
  
;  level = (nlev-1) - closest_val(4900,heights(::-1))
  level = (nlev-1) - closest_val(5100,heights(::-1))

  ; use level index from command line if specified
  if(isvar("l"))
    level = l
  end if

  print("height = "+heights(level))
  z = round( heights(level), 3)

  print ("tracer_id:    q    = " + tracer_id)
  print ("t_index:      t    = " + t_index)
  print ("zoom_factor:  zoom = " + zoom_factor)
  print ("level:        l    = " + level)
  print ("eta:               = " + lev(level) )
  print ("z:                 = " + z )

  ; load chosen field

  if(tracer_id .eq. 1) then
    q_i = f->Q (t_index,level,:,:)

  else if (tracer_id .eq. 2) then
    q_i = f->Q2(t_index,level,:,:)

  else if (tracer_id .eq. 3) then
    q_i = f->Q3(t_index,level,:,:)

  else if (tracer_id .eq. 4) then
    q_i = f->Q4(t_index,level,:,:)

  end if
  end if
  end if
  end if

  ; delete units and long names

  q_i@units = " "
  q_i@long_name = " "

  ; display min and max values
  minval = min(q_i)
  maxval = max(q_i)
  print("minval="+minval+" maxval="+maxval)

  ; round to two decimal places?
  ;q_i = decimalPlaces(q_i,2,True)
  ;printMinMax(q_i,0)

  ; Setup the plot

  res1                      = True
  res1@gsnDraw              = False        		; panel plot
  res1@gsnFrame             = False           ; don't draw yet
  res1@gsnCenterString      = ""
  res1@vpWidthF             = 2
  res1@vpHeightF            = 1
  res1@sfXArray             = lon           	; uses lon as plot x-axis
  res1@sfYArray             = lat         		; uses lat for y axis
  res1@tiYAxisString        = "latitiude"     ; y-axis title
  res1@tiXAxisString        = "longitude"      ; x-axis title
  res1@lbLabelBarOn         = False           ; label bar

  res1@cnLinesOn            = False ;True
  res1@cnFillOn             = True
  res1@cnInfoLabelOn        = False           ; don't give each
  res1@cnLineLabelsOn       = False           ; panel its own
  res1@cnMonoLineColor      = True           ; single color lines
  res1@cnLineColor          = "black"         ; black lines
  res1@cnLineDrawOrder      = "PreDraw"
  res1@cnLineThicknessF     = 0.1
  res1@cnFillOpacityF       = 1.0
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = 1.00            ; max contour color label
  res1@cnMinLevelValF       = 0.00            ; min contour color label
  res1@cnLevelSpacingF      = 0.10            ; contour color spacing

  ; text properties

  res1@tmXBMajorLengthF     = 0.03
  res1@txFont               = 22              ; helvetica-bold
  res1@tiMainFont           = 22              ; helvetica-bold
  res1@tiMainFontHeightF    = 0.045            ; set font sizes
  res1@tiXAxisFontHeightF   = 0.05
  res1@tiYAxisFontHeightF   = 0.05
  res1@tmXBLabelFontHeightF = 0.04
  res1@tiXAxisOffsetYF      = 0.010           ; shift x label
  res1@tiYAxisOffsetXF      = 0.010           ; shift y label
	res1@tiMainString         = "DCMIP 1-1  q"+tracer_id+"  t="+sprintf("%3.1f", time(t_index))+" days"+ "  z="+z+"m  range=("+sprintf("%5.4f", minval)+","+sprintf("%5.4f", maxval)+")"

  ; zoom in on the center of the plot

  res1@trXMinF =  180.0 - 180.0/zoom_factor
  res1@trXMaxF =  180.0 + 180.0/zoom_factor
  res1@trYMinF =    0.0 -  90.0/zoom_factor
  res1@trYMaxF =    0.0 +  90.0/zoom_factor

  ; panel properties

  pres                      = True            ; Panel resource
  pres@gsnMaximize          = True            ; maximize panel size
  pres@gsnFrame             = False           ; don't advance yet
  pres@gsnPaperOrientation  = "portrait"      ; paper orientation
  pres@gsnPanelLabelBar     = True            ; add common label bar
  pres@lbLabelFontHeightF   = 0.020           ; set labelbar font height

  ; set output type and name

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  filename ="dcmip2012_test1-1_q"+tracer_id+"_t"+t_index
	wks  = gsn_open_wks(wks_type,filename)
  print("writing "+filename)
  
  ; use rainbow colormap with a few of color substitutions

  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")      ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)   ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                 ; white     for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)          ; gray for undershoots
  colors(255,:)= (/ 0.4, 0.0, 0.0 /)          ; brick-red for overshoots
  gsn_define_colormap(wks,colors)

  ; plot countours and label bar

  plot    = new (1, graphic)               			; define plot - need 4 panels
  plot(0) = gsn_csm_contour(wks,q_i(:,:),res1)		; plot Q1
  gsn_panel(wks,plot,(/1,1/),pres)			; 2x2 plot
  frame(wks)

end




