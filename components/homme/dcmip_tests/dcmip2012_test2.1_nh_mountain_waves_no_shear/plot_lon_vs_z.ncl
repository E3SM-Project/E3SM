load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  ;_____________________________________________________________________
  ; open file and read in data(time,lev,lat,long) from 0 to n-1

  data_dir ="./movies/dcmip2012_test2_11.nc"
  print ("data_dir="+data_dir)
  f    = addfile(data_dir,"r")

  lat  = f->lat
  lon  = f->lon
  lev  = f->lev
  time = f->time

  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  nlev = dimsizes(lev)
  nt   = dimsizes(time)

  ti = nt-1                              ; set t to max time
  if(isvar("t")) then
    ti = t                               ; use t from command line if possible
  end if

  X               = 500   ; small-planet scale factor
  dt              = 0.1   ; seconds per step for hydrostatic
  output_interval = 100 ;1000  ; steps per output
  sec_per_day     = 3600.0*24


  t_sec = time(ti)*sec_per_day
  earth_days = time(ti)*X
  print ("ti="+ti+" t= "+t_sec+"s earth_days="+earth_days)

  equator = nlat/2                          ; The equator is at nlat/2

  ; Select var_choice
  ; 1 - Temperature Pert
  ; 2 - Zonal Velocity Pert
  ; 3 - Vertical Pressure Velocity

  var_choice = 1
  if(isvar("n")) then
    var_choice = n
  end if

  if (var_choice .eq. 1) then               ; Select T
    var1 = f->T(ti,:,equator,:) - 300
  end if

  if (var_choice .eq. 2) then               ; Select u
    var1 = f->u(ti,:,equator,:) - 20.0
  end if

  if (var_choice .eq. 3) then               ; Select omega
    var1 = f->omega(ti,:,equator,:)-0.0
  end if

  if (var_choice .eq. 4) then               ; Select v
    var1 = f->v(ti,:,equator,:)-0.0
  end if

  ;var1@units = ""
  ;var1@long_name = " "

  ;_____________________________________________________________________
  ; Get z positions of eta levels

  Teq     = 300.d0      ; Temperature at Equator
  Rd      = 287.0d0     ; Ideal gas const dry air (J kg^-1 K^1)
  g       = 9.80616d0   ; Gravity (m s^2)
  H       = Rd*Teq/g
  z_eta   = -H * log(lev)


  epsilon = 0.001 ;meter
  do i=1,nlev-1
    if( z_eta(i) .eq. z_eta(i-1) ) then
      z_eta(i) = z_eta(i)-epsilon
    end if
  end do

;print("z_eta="+z_eta)

  ;_____________________________________________________________________
  ; Get evenly spaced z levels

  z   = new((/nlev/),double)
  dz  = 30000.0/nlev
  do kk=0,nlev-1
    z(kk) = (30000.0-0.5*dz)-(kk)*dz
  end do
  
  ;_________________________________________________
  ; Interpolate from z_eta to evenly spaced z


  var_z = new( (/nlev,nlon/), double)

	do i = 0, nlon-1

    z_col       = z_eta(::-1)                 ; single column of z field, inverted
    var_col     = var1 (::-1,i)               ; single column of var field, inverted
		var_z(:,i)  = ftcurv(z_col, var_col, z) ; interpolate to evenly spaced z using splines

	end do


  ;_____________________________________________________________________
  plot  = new (1, graphic)                          ; define plot - need 3 panels
  res1                          = True
  res1@gsnDraw                  = True              ; panel plot
  res1@gsnFrame                 = False             ; don't draw yet
  res1@gsnMaximize              = True
  res1@cnFillOn                 = True
  res1@cnLinesOn                = True
  res1@gsnContourLineThicknessesScale = 0.1
  ;res1@gsnSpreadColors          = True
  res1@lbLabelAutoStride        = True
  res1@gsnCenterString          = ""
  res1@tiMainString             = ""
  res1@vpWidthF                 = 2.0
  res1@vpHeightF                = 1.0
  res1@cnInfoLabelOn            = False
  res1@cnLineLabelsOn           = False
  res1@lbLabelBarOn             = True
  res1@lbOrientation            = "horizontal"
  res1@lbLabelStride            = 1
  res1@gsnPaperOrientation      = "portrait"
  res1@sfXArray                 = lon           		; uses lon as plot x-axis
  res1@sfYArray                 = z/1000.0          ; uses z for y axis
  res1@trYReverse               = False          		; reverses y-axis, false
  res1@tiYAxisString            = "height (km)"     ; y-axis title
  res1@tiXAxisString            = ""                ; x-axis title
  res1@gsnStringFontHeightF     = 0.04
  res1@tmXBLabelFontHeightF     = 0.03
  res1@tmYLLabelFontHeightF     = 0.03
  res1@tiYAxisFontHeightF       = 0.03
  res1@cnLevelSelectionMode     = "AutomaticLevels"
  ;res1@cnFillPalette           ="temp_diff_18lev"; "BlWhRe"; "BlueWhiteOrangeRed";  "BlueRed"; "cmp_b2r";
  res1@lbOrientation        = "Vertical"

;  res1@gsnContourZeroLineThicknessF = 0.0


  print("max = "+max(var1)+" min="+min(var1))
  prefix = "dcmip2012_test2_1_"

  if (var_choice .eq. 1) then         ; Setup for T'

    pltTitle = sprintf("CAM-SE Test 2-1  T'  t=%.1f s",t_sec)
    fileName = sprinti(prefix + "T_t%i",ti)

    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnMaxLevelValF  = 1.2				; max contour color label
    res1@cnMinLevelValF  =-1.2				; min contour color label
    res1@cnLevelSpacingF = 0.2				; contour color spacing

  end if

  if (var_choice .eq. 2) then		; Setup for U'

    pltTitle = sprintf("CAM-SE Test 2-1  U'  t=%.1f s",t_sec)
    fileName = sprinti(prefix + "u_t%i",ti)
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnMaxLevelValF  = 2				; max contour color label
    res1@cnMinLevelValF  =-2				; min contour color label
    res1@cnLevelSpacingF = 0.4				; contour color spacing

  end if

  if (var_choice .eq. 3) then		; Setup for OMEGA

    pltTitle = sprintf("CAM-SE Test 2-1  Omega  t=%.1f s",t_sec)
    fileName = sprinti(prefix + "omega_t%i",ti)
     res1@cnLevelSelectionMode = "ManualLevels"
     res1@cnMaxLevelValF  = 30.0				; max contour color label
     res1@cnMinLevelValF  =-30.0				; min contour color label
     res1@cnLevelSpacingF = 5.0				; contour color spacing

  end if

  if (var_choice .eq. 4) then		; Setup for V

    pltTitle = sprintf("CAM-SE Test 2-1  V  t=%.1f s",t_sec)
    fileName = sprinti(prefix + "v_t%i",ti)
     res1@cnLevelSelectionMode = "ManualLevels"
     res1@cnMaxLevelValF  = 1.8				; max contour color label
     res1@cnMinLevelValF  =-1.8				; min contour color label
     res1@cnLevelSpacingF = 0.4				; contour color spacing

  end if
  print("filename = "+fileName)

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  wks_type@wkOrientation    = "portrait"
  wks = gsn_open_wks(wks_type,fileName)
  gsn_define_colormap(wks,"gui_default")  		; Colormap

  res1@gsnCenterString = pltTitle

;  plot(0) = gsn_csm_contour(wks,var1(:,:),res1)		; plot var1
  plot(0) = gsn_csm_contour(wks,var_z(:,:),res1)		; plot var1

  frame(wks)

end

