;================================================;
;  test11-height-lon.ncl
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  ;=================================================;
  ; open file and read in data (time,lev,lat,long) from 0 to n-1
  ;=================================================;

  ; replace your models file name below

  f    = addfile("./movies/dcmip2012_test1_11.nc","r")

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
  time = f->time

	nlat = getfilevardimsizes(f, "lat" )
	nlon = getfilevardimsizes(f, "lon" )
	nlev = getfilevardimsizes(f, "lev" )
  nt   = dimsizes(time)

  ; get time index from command line

  t_index = nt-1
  if(isvar("t")) then
    t_index = t         ; use t from command line if possible
  end if
  print ("t_index="+t_index)

  ; get zoom factor from command line
  
  zoom_factor = 1.0
  if(isvar("zoom"))
    zoom_factor = zoom
  end if

	; We require plots of q1, q2, q3 and q4 at the equator.
  
  equator = nlat/2

  ; Input tracers

	q1 = f->Q (t_index,:,equator,:)
	q2 = f->Q2(t_index,:,equator,:)
	q3 = f->Q3(t_index,:,equator,:)
	q4 = f->Q4(t_index,:,equator,:)

  ; Delete tracer units and long names

	q1@units = " "
	q2@units = " "
	q3@units = " "
	q4@units = " "

	q1@long_name = " "
	q2@long_name = " "
	q3@long_name = " "
	q4@long_name = " "

  ; Set up height coordinate
  ; Note that for CAM-FV pressure coordinates are used,
  ; therefore the height levels are upside-down
  ; i.e. for 60 vertical levels index 0 is at 11900 m
  ; and index 59 is at 100 m.

	z = new((/nlev/),double)
	dz = 12000.0/nlev

	do kk=0,nlev-1
		z(kk) = (12000.0-0.5*dz)-(kk)*dz
  end do

  ; We now produce the plot

	plot  = new (4, graphic)  ; define plot - need 4 panels

	res1                      = True
	res1@gsnDraw              = False     ; panel plot
	res1@gsnFrame             = False			; don't draw yet
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = ""
	res1@vpWidthF             = 2
	res1@vpHeightF            = 1
  res1@lbLabelBarOn         = False			; label bar

  res1@cnFillOn             = True
	res1@cnLinesOn            = False
  res1@cnInfoLabelOn        = False			; don't give each
	res1@cnLineLabelsOn       = False			; panel its own
  res1@cnLevelSelectionMode = "ManualLevels"
	res1@cnMaxLevelValF       = 0.91            ; max contour color label
	res1@cnMinLevelValF       = 0.01            ; min contour color label
	res1@cnLevelSpacingF      = 0.05            ; contour color spacing
  res1@cnFillPalette        = "WhBlGrYeRe" ;"WhiteBlueGreenYellowRed";
	;pres@lbLabelStride            = 1

  ; zoom in on the center of the plot

  res1@trXMinF =  180.0 - 180.0/zoom_factor
  res1@trXMaxF =  180.0 + 180.0/zoom_factor
  res1@trYMinF =    5.0 -  5.0/zoom_factor
  res1@trYMaxF =    5.0 +  7.0/zoom_factor

	pres                          = True
	pres@gsnMaximize              = True  
	pres@gsnPanelLabelBar         = True			; Communal label bar
	pres@gsnFrame                 = False
  pres@gsnPaperOrientation      = "portrait"
  pres@gsnPanelXWhiteSpacePercent = 0

	pres@pmLabelBarOrthogonalPosF = 0.010
  pres@lbLabelFontHeightF       = 0.012

	res1@sfXArray                 = lon           		; uses lon as plot x-axis
	res1@sfYArray                 = z/1000.0             		; uses z  for y axis
	res1@trYReverse               = False          		; reverses y-axis, false
	res1@tiYAxisString            = ""   	       		; y-axis title
	res1@tiXAxisString            = ""             		; x-axis title

 	pltTitle="dcmip 1-1     vertical equatorial slice     t="+sprintf("%3.1f", time(t_index))+" days"
	pres@txString                 = pltTitle
  pres@txFontHeightF            = 0.015
  pres@txFont                   = 22             ; helvetica-bold
;  pres@txPerimOn =True

  pres@tiMainOffsetYF=1
  ;res1@txFont                   = 22             ; helvetica-bold
  res1@tmXBLabelFontHeightF     = 0.04
  res1@tiMainFontHeightF        = 0.04            ; change font heights
  ;res1@tiMainFont               = 22
  ;res1@gsnStringFontHeightF     = 0.04

  ; Change the output type and name

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  filename ="test11_heightlon_t"+t_index
	wks  = gsn_open_wks(wks_type,filename)
  print("writing "+filename)

  ;gsn_define_colormap(wks,"gui_default")            ; Colormap
	;res1@gsnSpreadColors      = True

	res1@gsnCenterString = "q1"
	plot(0) = gsn_csm_contour(wks,q1(:,:),res1)       ; plot Q1

	res1@gsnCenterString = "q2"
  plot(1) = gsn_csm_contour(wks,q2(:,:),res1)       ; plot Q2

	res1@gsnCenterString = "q3"
	plot(2) = gsn_csm_contour(wks,q3(:,:),res1)       ; plot Q3

	res1@gsnCenterString = "q4"
  plot(3) = gsn_csm_contour(wks,q4(:,:),res1)       ; plot Q4

  pres@gsnPanelXF = (/-0.50, .55, -0.50, .55/)    ; Adjust rightmost plots
  pres@gsnPanelYF = (/-1,-1,.49,.49/)             ; Adjust bottommost plots.

	gsn_panel(wks,plot,(/2,2/),pres)			; 2x2 plot

  ; Add height and longitude labels

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.015
	gsn_text_ndc(wks,"Height (km)",0.010,0.52,txres3)

	;txres3@txAngleF       = 0.
	;gsn_text_ndc(wks,"Longitude",0.52,0.23,txres3)

	frame(wks)

end  
