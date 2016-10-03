; plot_z_lon
;
; Plots vertical cross section from interpolated data
;
; usage plot_z_lon t=(time_index) n=(variable_number)
; n = (0,1,2) => (omega,u,v)


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  T0        = 300.d0          ; temperature (K)
  lapse     = 0.0065d0        ; temperature lapse rate (K/m)
  ztop      = 12000.d0        ; model top (m)
  Rd        = 287.0d0         ; Ideal gas const dry air (J kg^-1 K^1)
  g         = 9.80616d0       ; Gravity (m s^2)
  Cp        = 1004.5d0        ; Specific heat capacity (J kg^-1 K^1)
  pi        = 4.d0*atan(1.d0) ; pi
  p0        = 100000.0        ; reference pressure

  exponent  = g/(Rd*lapse)
  exponent_r= 1.0/exponent

  ;_________________________________________________
  ; open file and read in data(time,lev,lat,long)

  f     = addfile("./movies/dcmip2012_test2_01.nc","r")
	lat   = f->lat
	lon   = f->lon
	lev   = f->lev
  time  = f->time
  hyam  = f->hyam
  hybm  = f->hybm

	nlat  = getfilevardimsizes(f, "lat" )
	nlon  = getfilevardimsizes(f, "lon" )
	nlev  = getfilevardimsizes(f, "lev" )
  nt    = dimsizes(time)

  ;_________________________________________________
  ; get time index from command line

  t_index = nt-1
  if(isvar("t")) then
    t_index = t
  end if
  print ("t_index: t = " + t_index)

	equator   = nlat/2

  ;_________________________________________________
  ; load data

  var_num = 0
  if(isvar("n")) then
    var_num = n
  end if

  omega = f->omega(t_index,:,equator,:)
  u     = f->u(t_index,:,equator,:)
  v     = f->v(t_index,:,equator,:)

  print("max u    ="+max(u)+" min u="+min(u))
  print("max v    ="+max(v)+" min v="+min(v))
  print("max omega="+max(omega)+" min omega="+min(omega))

  if(var_num .eq. 0) then
    var       = f->omega (t_index,:,equator,:)
    var_name  = "omega"
  end if

  if(var_num .eq. 1) then
    var       = f->u(t_index,:,equator,:)
    var_name  = "u"
  end if

  if(var_num .eq. 2) then
    var       = f->v(t_index,:,equator,:)
    var_name  = "v"
  end if

  ps  = f->ps(t_index,equator,:)

  norm = max(abs(var))

  amax = 1.2
  step = amax/6

  var = var/norm
	;var@long_name = " "

  ;_________________________________________________
  ; Set up height coordinate with constant z spacing

	z  = new((/nlev/),double)
	dz = 12000.0/nlev

	do kk=0,nlev-1
		z(kk) = (12000.0-0.5*dz)-(kk)*dz
  end do

  ; Get surface elevation curve

  zs = (T0/lapse)*(1.0-(ps/p0)^exponent_r)

  ;_________________________________________________
  ; Get pressure and z position of eta-level midpoints

  z_m = new((/nlev,nlon/),double)
  p_m = new((/nlev,nlon/),double)

  do k =0, nlev-1
		p_m(k,:) = hyam(k)*p0 + hybm(k)*ps
		z_m(k,:) = (T0/lapse)*(1.0-(p_m(k,:)/p0)^exponent_r)
  end do

  ;_________________________________________________
  ; Interpolate var field from eta coord to z coord

  var_z = new( (/nlev,nlon/), double)

	do i = 0, nlon-1
    z_col       = z_m(::-1,i) ; single column of z   field, inverted
    var_col     = var(::-1,i) ; single column of var field, inverted
		var_z(:,i)  = ftcurv(z_col, var_col, z) ; interpolate to evenly spaced z using splines
	end do

  ;_________________________________________________
  ; Setup the plot

	res1                      = True
	res1@gsnDraw              = False
	res1@gsnFrame             = False       ; don't draw yet
  res1@gsnMaximize          = True
	res1@cnFillOn             = True
	res1@cnLinesOn            = False       ; turn off contour lines
  res1@vpWidthF             = 1
	res1@vpHeightF            = 0.5
	res1@cnInfoLabelOn        = False
	res1@cnLineLabelsOn       = False
	res1@lbLabelBarOn         = True
  res1@sfXArray             = lon         ; use lon as plot x-axis
	res1@sfYArray             = z           ; use z for y axis
	res1@trYReverse           = False       ; reverses y-axis, false
  res1@tiXAxisString        = "longitude" ; x-axis title
  res1@tiYAxisString        = "height"    ; y-axis title
  res1@gsnPaperOrientation  = "portrait"
  res1@lbBoxLinesOn         = False     ; Turn off labelbar box lines
  res1@lbLabelStride        = 1
  res1@lbAutoManage         = False
  res1@lbTitlePosition      = "Bottom"
  res1@lbOrientation        = "Vertical"
  res1@lbTitleString        = sprintf("*%7.3f", norm)
  res1@cnFillPalette        = "temp_diff_18lev"; "WhBlGrYeRe"; "BlWhRe";  "BlueWhiteOrangeRed"; "BlueRed"; "cmp_b2r";
  res1@cnSpanFillPalette    = True
  res1@tiXAxisFontHeightF   = 0.020
  res1@tmXBLabelFontHeightF = 0.015
  res1@lbTitleFontHeightF   = 0.020
  res1@lbLabelFontHeightF   = 0.020
  res1@tiMainFont           = 22
  res1@tiMainFontHeightF    = 0.020
  res1@gsnCenterStringFontHeightF = 0.025			; to set the gsnLeft/Center/RightString font heights,
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = amax-0.5*step;  0.9; ; max contour color label
	res1@cnMinLevelValF       =-amax+0.5*step; -0.9; ; min contour color label
  res1@cnLevelSpacingF      = step;           0.1; ; contour color spacing
  res1@gsnCenterString      =""
	res1@tiMainString         = "DCMIP 2012 Test 2-0  "+var_name+"   t = "+sprintf("%5.2f", time(t_index)) + " days"

  wks_type="pdf"
  wks_type@wkPaperWidthF    = 11
  wks_type@wkPaperHeightF   = 7.5
  filename ="dcmip2012_test2_0_"+var_name+"_t"+t_index
	wks  = gsn_open_wks(wks_type,filename)
  print("outputfilename = "+filename+".pdf")

  gsres                     = True
  gsres@gsLineColor         = "black"
  gsres@gsFillColor         = "black"
  gsres@tfPolyDrawOrder     = "PostDraw"
  gsres@gsLineThicknessF    = 5.0

  zs(0)=0
  zs(nlon-1)=0

  plot1     = gsn_csm_contour (wks,var(:,:),res1)
  mountains = gsn_add_polygon (wks, plot1,lon, zs, gsres)
  mt        = gsn_add_polyline(wks, plot1, lon, zs, gsres)

  draw(plot1)

	frame(wks)

end

