;================================================;
;	Example ncl script to produce the set of
;	lat-lon plots for DCMIP-2016 test case 1
; Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan) for DCMIP-2012 
;
; Modified by Colin Zarzycki for DCMIP-2016
;
; v1.01 - 6/7/16 
; v1.02 - 6/8/16 - Typo in var_choice options
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

; NOTE, if your model does not output T850 OR P at level
; midpoints OR you cannot calculate P from hybrid levels, you
; will need to make model specific additions denoted by
; CALCPHERE below

; PRECL needs to be in m/s for contours to plot correctly

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,lon) from 0 to n-1
;=================================================;

  ; NOTE, that this file needs to be a catted file including all times at daily interval
  ; if you are outputting 1 file per time, you can run "ncrcat dcmip1.files.*.nc cat.nc"
  ; for one file

	f = addfile("./movies/dcmip2016_test11.nc","r")

  ; Input useful parameters

	lat  = f->lat
        lon  = f->lon
	lev  = f->lev
	nlat = getfilevardimsizes(f, "lat" )
	nlon = getfilevardimsizes(f, "lon" )
	nlev = getfilevardimsizes(f, "lev" )

  ; We want plots of temperature, vertical velocity and relative vorticity 
  ; at the 850 hPa level. We also want surface pressure.

  ; Select var_choice
  ; 1 - Temperature at 850 hPa
  ; 2 - Large-scale precipitation
  ; 3 - Surface Pressure

        if(isvar("var_choice")) then
        else
           var_choice = 1
        end if
	
	if (var_choice .eq. 1) then ; Select T
		if (isfilevar(f,"T850")) then
			varload = f->T850

		else	

			; we interp T to the 850 hPa level
			varload0 = f->T				; T is a 4D var
			varload  = f->T(:,0,:,:)
			; Calculate P using hybrid coefficients
			; If you have pressure as a variable just do
			; P = f->P
      ; CALCPHERE
;     hyam = f->hyam
;  		hybm = f->hybm
;  		P0   = f->P0
;			ps   = f->ps
;			P    = 0.0*varload0
;			P = pres_hybrid_ccm(PS,P0,hyam,hybm)
;			P!0 = "time"
;			P!1 = "lev"
;			P!2 = "lat"
;			P!3 = "lon"

      P = f->pnh
			; interp to 850 hPa level
			plevel = 85000.0
			varload = int2p_n(P,varload0,plevel,2,1)

		end if

	else if (var_choice .eq. 2) then			; Select Large-scale precipitation

		varload = f->precl
    varload=varload*8.64e7		

	else if (var_choice .eq. 3) then			; Select PS

		varload = f->ps

	end if
	end if	
	end if

  ; We want the output at days 6, 9, 12 and 15. If your model
  ; output is in 6 hourly intervals this becomes indices 24, 36, 48 and 60

time = f->time
nt = dimsizes(time)
print("nt="+nt)

day6   = closest_val( 6, time)
day9   = closest_val( 9, time)
day12  = closest_val( 12, time)
day15  = closest_val( 15, time)

	var1 = varload(day6,:,:)
	var2 = varload(day9,:,:)
	var3 = varload(day12,:,:)
	var4 = varload(day15,:,:)

  ; We don't want long-name in the plots

	var1@long_name = " "
	var2@long_name = " "
	var3@long_name = " "
	var4@long_name = " "

  ; Delete loaded data

	delete(varload)

		; We now produce the plot

	plot  = new (4, graphic)               			; define plot - need 4 panels

	res1                      = True
	res1@gsnDraw              = False          		; panel plot 
	res1@gsnFrame             = False			; don't draw yet
	res1@cnFillOn             = True
	res1@cnLinesOn            = True
	res1@gsnSpreadColors      = True
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = ""
	res1@vpWidthF             = 0.38
	res1@vpHeightF            = 0.19
	res1@cnLevelSelectionMode = "ManualLevels"
	res1@cnInfoLabelOn        = False			; don't give each
	res1@cnLineLabelsOn       = False			; panel its own
	res1@lbLabelBarOn         = False			; label bar

	pres                          = True
	pres@gsnMaximize              = True  
	pres@gsnPanelLabelBar         = True			; Communal label bar
	pres@gsnPanelLeft             = 0.1
	pres@gsnPanelRight            = 0.9
	pres@pmLabelBarOrthogonalPosF = -0.03
	pres@gsnFrame                 = False
	pres@lbLabelStride            = 1

	res1@sfXArray        = lon            			; uses lon as plot x-axis
	res1@sfYArray        = lat            			; uses lat for y axis
	res1@trYReverse      = False          			; reverses y-axis, false
	res1@tiYAxisString   = ""   	       			; y-axis title
	res1@tiXAxisString   = ""             			; x-axis title

	if (var_choice .eq. 1) then 	;================ Temperature plot ====================
 
	pltTitle="PREQX Test 16-1, T850 hPa" 			; Plot title if required
	pres@txString = pltTitle

		; Change the output type and name

	;wks  = gsn_open_wks("X11","acme-test16-1-latlonT850")   	; output using X11
	wks  = gsn_open_wks("pdf","preqx-test16-1latlonT850")   	; output using eps
	gsn_define_colormap(wks,"gui_default")  

	res1@cnMaxLevelValF  = 300.0				; max contour color label
	res1@cnMinLevelValF  = 230.0				; min contour color label
	res1@cnLevelSpacingF = 10.0				; contour color spacing			; choose a colormap
 
	else if (var_choice .eq. 2) then 	;================ Large-scale precipitation plot ====================
 
	pltTitle="PREQX Test 16-1, Large-scale precipitation" 	; Plot title if required
	pres@txString = pltTitle

	var1@units = " "
	var2@units = " "
	var3@units = " "
	var4@units = " "

	res1@cnLinesOn            = False				; Turn lines off for clearer plot

		; Change the output type and name

	;wks  = gsn_open_wks("X11","acme-test16-1latlonPRECL")   	; output using X11
	wks  = gsn_open_wks("pdf","preqx-test16-1latlonPRECL")   	; output using eps
	gsn_define_colormap(wks,"gui_default")  

		; Note that the maximum might be larger than these contour spacing values

	res1@cnMaxLevelValF  = 50.0				; max contour color label
	res1@cnMinLevelValF  = 0.0				; min contour color label
	res1@cnLevelSpacingF = 4.0				; contour color spacing			; choose a colormap
 
	else if (var_choice .eq. 3) then 	;================ Surface Pressure plot ====================
 
	pltTitle="PREQX Test 16-1, PS" 			; Plot title if required
	pres@txString = pltTitle

		; Convert to hPa 

	var1=var1/100.0
	var2=var2/100.0
	var3=var3/100.0
	var4=var4/100.0

	var1@units = "hPa "
	var2@units = "hPa "
	var3@units = "hPa "
	var4@units = "hPa "

		; Change the output type and name

	;wks  = gsn_open_wks("X11","acme-test16-1latlonPS")   	; output using X11
	wks  = gsn_open_wks("pdf","preqx-test16-1latlonPS")   	; output using eps
	gsn_define_colormap(wks,"gui_default")  

	res1@cnMaxLevelValF  = 1020.0				; max contour color label
	res1@cnMinLevelValF  = 920.0				; min contour color label
	res1@cnLevelSpacingF = 10.0				; contour color spacing			; choose a colormap

	end if				;============== end plot choice ===============
	end if				;============== end plot choice ===============
	end if				;============== end plot choice ===============

	res1@gsnCenterString = "t = 6"
	plot(0) = gsn_csm_contour(wks,var1(:,:),res1)		; plot var1

	res1@gsnCenterString = "t = 9"
     	plot(1) = gsn_csm_contour(wks,var2(:,:),res1)		; plot var2

	res1@gsnCenterString = "t = 12"
	plot(2) = gsn_csm_contour(wks,var3(:,:),res1)		; plot var3

	res1@gsnCenterString = "t = 15"
     	plot(3) = gsn_csm_contour(wks,var4(:,:),res1)		; plot var4

	gsn_panel(wks,plot,(/2,2/),pres)			; 2x2 plot

		; Add latitude and longitude labels

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.02
	gsn_text_ndc(wks,"Latitude",0.08,0.49,txres3)

	txres3@txAngleF       = 0.
	gsn_text_ndc(wks,"Longitude",0.5,0.22,txres3)

	frame(wks)

end






