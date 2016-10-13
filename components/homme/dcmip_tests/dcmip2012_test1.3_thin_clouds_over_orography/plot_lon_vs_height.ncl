;================================================;
; DMCIP Test 1-3
; plot_tracers.ncl
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  f    = addfile("movies/dcmip2012_test1_31.nc","r")
	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
  time = f->time

	nlat = getfilevardimsizes(f, "lat" )
	nlon = getfilevardimsizes(f, "lon" )
	nlev = getfilevardimsizes(f, "lev" )
  nt   = dimsizes(time)

  ; get time_index t from command line

  interval = 24.0;hours
  ind2days = 1.0*interval/24
  t_index = nt-1
  if(isvar("t")) then
    t_index = t
  end if
  print ("t_index: t = " + t_index)

  ; Load tracer data and surface pressure

  equator = nlat/2
	q1 = f->Q4(t_index,:,equator,:)
	ps = f->ps(t_index,  equator,:)

  ; round to 4th decimal place, to prevent plotting errors like q=-1.0e-14
  q1=round(q1*1e4,2)/1e4

  print("min q = "+min(q1))
  print("max q = "+max(q1))

  ; Set up reference height coordinate

  numz = 4*nlev
	refz = new((/numz/),double)
	dz   = 12000.0/numz

	do kk=0,numz-1
		refz(kk) = (12000.0-0.5*dz)-(kk)*dz
  end do

	; We want to interpolate to the constant height levels, 
	; so we need to calculate P and Z at the models levels

  ; Calculate height at each point based on
  ; z = Rd T0/g * ln(p0/p)
  ; and p = a p0 * b ps

	g      = 9.80616
  radius = 6371229.
  T0     = 300.
  Rd     = 287.04
  cp     = 1004.64
  H      = Rd * T0 / g

  hyam = f->hyam
  hybm = f->hybm

	P0 = 100000.0

  z = new((/nlev,nlon/),double)
  P = new((/nlev,nlon/),double)

	do i =0, nlon-1
		do k =0, nlev-1

		P(k,i) = hyam(k)*P0 + hybm(k)*ps(i)
		z(k,i) = H * log(P0/P(k,i))

		end do
	end do
  zs = H * log(P0/ps)

	; Interpolate to height levels between 2000 and 10000m.


	zz  = new((/numz/),double)

	do kk = 0, numz-1
		zz(kk) = 000.0 + kk*dz;200.0
	end do

  ; Simple Interpolation

  qq      = new((/numz,nlon/),double)
  qqpart  = new((/numz/),double)
  q1part  = new((/nlev/),double)
  zpart   = new((/nlev/),double)
  q1part2 = new((/nlev/),double)
  zpart2  = new((/nlev/),double)

  q1n     = new((/nlev,nlon/),double)
  q2n     = new((/nlev,nlon/),double)

	do i = 0, nlon-1

		q1part(:) = q1(:,i)			; put in 1D column
		zpart(:)  = z(:,i)
		q1part2   = q1part(::-1)		; Reverse Q
		zpart2    = zpart(::-1)			; Reverse z
		;qqpart    = ftcurv(zpart2,q1part2,zz)	; Interpolate to zz
    qqpart = linint1(zpart2,q1part2,False,zz,0)
		qq(:,i)   = qqpart(:)       ; Make 2D again

	end do

  ; Remove tracer units and long-names from plot

	q1@units  = " "
	qq@units  = " "
	q1@long_name  = " "
	qq@long_name  = " "

  ; Setup the plots

	plot  = new (2, graphic)               			; define plot - need 4 panels

	res1                      = True
	res1@gsnDraw              = False         	; panel plot
	res1@gsnFrame             = False           ; don't draw yet
	;res1@gsnSpreadColors      = True
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = ""
	res1@vpWidthF             = 2
	res1@vpHeightF            = 1
	res1@lbLabelBarOn         = False
  res1@sfXArray             = lon             ; use lon as plot x-axis
	res1@sfYArray             = lev             ; use z for y axis
	res1@trYReverse           = True            ; flip y-axis

  res1@cnFillOn             = True
	res1@cnLinesOn            = False
  res1@cnInfoLabelOn        = False
	res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = 1.0
	res1@cnMinLevelValF       = 0.0
	res1@cnLevelSpacingF      = 0.1
  res1@cnLineThicknessF     = 0.2
  ;res1@cnLineDrawOrder      = "PreDraw"

  res1@tiYAxisString        = "~F33~h"+"~F21~ level"  ; y-axis title
	res1@tiXAxisString        = "longitude"     ; x-axis title
  res1@txFont               = 22              ; helvetica-bold
  res1@tiMainFont           = 22              ; helvetica-bold
  res1@tiMainFontHeightF    = 0.06            ; set font sizes
  res1@tiXAxisFontHeightF   = 0.05
  res1@tiYAxisFontHeightF   = 0.05
  res1@tmXBLabelFontHeightF = 0.04
  ;res1@tiXAxisOffsetYF      = 0.015           ; shift x label
  ;res1@tiYAxisOffsetXF      = 0.015           ; shift y label
  res1@gsnYAxisIrregular2Linear = True

	res2                      = True
	res2@gsnDraw              = False           ; panel plot
	res2@gsnFrame             = False           ; don't draw yet
	res2@lbLabelAutoStride    = True
  res2@lbLabelBarOn         = False
  res2@gsnCenterString      = ""
	res2@tiMainString         = ""
	res2@vpWidthF             = 2
	res2@vpHeightF            = 1
  res2@sfXArray             = lon             ; uses lon as plot x-axis
	res2@sfYArray             = zz              ; uses zz  for y axis
	res2@trYReverse           = False           ; reverses y-axis, false

  res2@cnFillOn             = True
	res2@cnLinesOn            = False
  res2@cnInfoLabelOn        = False           ; don't give each
	res2@cnLineLabelsOn       = False           ; panel its own
  res2@cnLevelSelectionMode = "ManualLevels"
  res2@cnMaxLevelValF       = 1.0             ; max contour color label
	res2@cnMinLevelValF       = 0.0             ; min contour color label
	res2@cnLevelSpacingF      = 0.1             ; contour color spacing
  ;res2@cnLineDrawOrder      = "PreDraw"
  res2@cnLineThicknessF     = 0.01
  res2@cnMonoLineColor      = True

  res2@tiYAxisString        = "z level"        ; y-axis title
	res2@tiXAxisString        = "longitude"     ; x-axis title
  res2@txFont               = 22              ; helvetica-bold
  res2@tiMainFont           = 22              ; helvetica-bold
  res2@tiMainFontHeightF    = 0.06            ; set font sizes
  res2@tiXAxisFontHeightF   = 0.05
  res2@tiYAxisFontHeightF   = 0.05
  res2@tmXBLabelFontHeightF = 0.04
  ;res2@tiXAxisOffsetYF      = 0.015           ; shift x label
  ;res2@tiYAxisOffsetXF      = 0.015           ; shift y label

	pres                      = True
	pres@gsnMaximize          = True
	pres@gsnFrame             = False
  pres@gsnPaperOrientation  = "portrait"
  pres@gsnPanelTop          = 1
	;pres@gsnPanelBottom       = 0.10
  pres@gsnPanelLabelBar     = True            ; add common label bar
  pres@lbLabelFontHeightF   = 0.020           ; set labelbar font height
  pres@gsnPanelYWhiteSpacePercent = 0.0

  ; set output type and name

  wks_type="pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 6

  filename ="dcmip2012_test1-3_L60_t"+t_index
	wks  = gsn_open_wks(wks_type,filename)
  print ("writing file: " + filename)

  ; WBGYR
  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")      ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)   ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                 ; white     for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)          ; gray for undershoots
  colors(255,:)= (/ 0.4, 0.0, 0.0 /)          ; brick-red for overshoots
  gsn_define_colormap(wks,colors)

  days = t_index*ind2days

	res1@gsnCenterString = "DCMIP 1-3   1 dg L60 "+sprintf("%6.2f", days)+" days"
	plot(0) = gsn_csm_contour(wks,q1(:,:),res1)		; plot Q on model levels

	res2@gsnCenterString = "DCMIP 1-3   1 dg L60   "+sprintf("%6.2f", days)+" days"
	plot(1) = gsn_csm_contour(wks,qq(:,:),res2)		; plot Q on height levels

  ; Add mountain range to plot

  gsres                     = True                    ; poly res
  gsres@gsLineColor         = "black"
  gsres@tfPolyDrawOrder     = "postdraw"
  gsres@gsLineOpacityF      = 0.5
  
  zs(0)=0
  zs(nlon-1)=0
  dummy = gsn_add_polygon (wks,plot(1),lon, zs, gsres)

  ; Add vertical-level contours to plot

  gsres@gsLineColor = "black"
  eta_line  = new(nlev,graphic)
  z_line    = new(nlev,graphic)
  plx       = (/lon(0),lon(nlon-1)/)

  do k = 0,nlev-1,10
    ply = (/lev(k),lev(k)/)
    eta_line(k) = gsn_add_polyline (wks,plot(0),plx, ply, gsres)
    z_line(k)   = gsn_add_polyline (wks,plot(1),lon, z(k,:), gsres)
  end do

	gsn_panel(wks,plot,(/2,1/),pres)			; 2x2 plot

	frame(wks)

end  
