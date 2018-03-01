;================================================;
;  test22-height-lon.ncl
;
;	Example ncl script to produce the height-lon 
;	plots for test 22.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

	; This produces the height-lon plots for test 22

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

		; replace your models file name below
name = "dcmip2012_test2_21"
title = "NonH, rsplit=0 "

var_choice = 2

f=addfile(name+".nc","r")	

	; Input useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

		; We want to plot temperature perturbation 
		; Also, for reference: zonal and vertical velocity

		; Select var_choice 	1 - Temperature Pert
		;			2 - Zonal Velocity Pert
		; 			3 - Vertical Pressure Velocity


	if (var_choice .eq. 1) then				; Select T

		varload = f->T(:,:,:,:)

	else if (var_choice .eq. 2) then			; Select U

		varload = f->u(:,:,:,:)

	else if (var_choice .eq. 3) then			; Select OMEGA

		varload = f->omega(:,:,:,:)

	end if
	end if	
	end if	

		; We want the plots at the equator at times 2400, 3600 and 7200 s
		; This corresponds to indices 24, 36 and 72
		; The equator is at nlat/2 

	equator = nlat/2

		; Input the variables at this time

	;time2400 = 24
	;time3600 = 36
	;time7200 = 72

       finaltime = 1
	var1 = varload(finaltime,:,equator,:)
	;var2 = varload(time3600,:,equator,:)
	;var3 = varload(time7200,:,equator,:)

	if (var_choice .eq. 1) then				; Select T
       mint = min(var1)
       maxt = max(var1) 
	else if (var_choice .eq. 2) then			; Select U
       minu = min(var1)
       maxu = max(var1) 
	else if (var_choice .eq. 3) then			; Select OMEGA
       mino = min(var1)
       maxo = max(var1) 
	end if
	end if	
	end if	


		; Set up height coordinate (assuming constant spacing)
		; Note that for CAM-FV pressure coordinates are used, 
		; therefore the height levels are upside-down 
		; i.e. for 60 vertical levels index 0 is at 29750 m 
		; and index 59 is at 250 m.

	z = new((/nlev/),double)
	dz = 30000.0/nlev

	do kk=0,nlev-1

		z(kk) = (30000.0-0.5*dz)-(kk)*dz  

  	end do	

		; Reference temp/velocity at equator

	vareq = 0.0*var1

	if (var_choice .eq. 1) then

		; Temperature at equator Teq

		vareq(:,:) = 300.0

	else if (var_choice .eq. 2) then

		; Zonal Velocity at Equator Ueq

		do k=0, nlev-1

			vareq(k,:) = 20.0*sqrt(2.0*0.00025*doubletofloat(z(k)) + 1.0 )

		end do

	else 

		; Set to zero

		vareq(:,:) = 0.0

	end if
	end if

	var1 = var1 - vareq
	;var2 = var2 - vareq
	;var3 = var3 - vareq

		; Delete loaded data

	delete(varload)

		; We don't want long-name in the plots

	var1@long_name = " "
	;var2@long_name = " "
	;var3@long_name = " "

		; We now produce the plot

	plot  = new (1, graphic)               			; define plot - need 3 panels

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
	pres@gsnPanelLeft             = 0.05
	pres@gsnPanelRight            = 0.95
	pres@lbOrientation            = "horizontal"
	pres@pmLabelBarOrthogonalPosF = -0.0
	pres@gsnFrame                 = False
	pres@lbLabelStride            = 1

	res1@sfXArray                 = lon           		; uses lon as plot x-axis
	res1@sfYArray                 = z             		; uses z for y axis
	res1@trYReverse               = False          		; reverses y-axis, false
	res1@tiYAxisString            = " "   	       		; y-axis title
	res1@tiXAxisString            = " "             	; x-axis title

	if (var_choice .eq. 1) then				; Setup for T'
 
		pltTitle             =title+" T' min/max="+sprintf("%5.2f", mint)+"/"+sprintf("%5.2f", maxt)		; Plot title if required
		pres@txString        = pltTitle

		res1@cnMaxLevelValF  = 1.5				; max contour color label
		res1@cnMinLevelValF  = -1.5				; min contour color label
		res1@cnLevelSpacingF = 0.25				; contour color spacing

			; Change the output type and name

		wks  = gsn_open_wks("pdf",name+"-T")   	; output using X11
		;wks  = gsn_open_wks("eps","cam-fv-test22heightlon1")   ; output using eps
		;wks  = gsn_open_wks("png","cam-fv-test22heightlon1")   ; output using png
		gsn_define_colormap(wks,"gui_default")  		; Colormap

	else if (var_choice .eq. 2) then			; Setup for U'
 
		pltTitle             =title+" U' min/max="+sprintf("%5.2f", minu)+"/"+sprintf("%5.2f", maxu)		; Plot title if required
		pres@txString        = pltTitle

		res1@cnMaxLevelValF  = 3.0				; max contour color label
		res1@cnMinLevelValF  = -3.0				; min contour color label
		res1@cnLevelSpacingF = 0.5				; contour color spacing

			; Change the output type and name

		wks  = gsn_open_wks("pdf",name+"-U")   	; output using X11
		;wks  = gsn_open_wks("eps","cam-fv-test22heightlon2")   ; output using eps
		;wks  = gsn_open_wks("png","cam-fv-test22heightlon2")   ; output using png
		gsn_define_colormap(wks,"gui_default")  		; Colormap
		
	else if (var_choice .eq. 3) then			; Setup for OMEGA
 
		pltTitle             =title+"OMEGA" 		; Plot title if required
		pres@txString        = pltTitle

		res1@cnMaxLevelValF  = 48.0				; max contour color label
		res1@cnMinLevelValF  = -48.0				; min contour color label
		res1@cnLevelSpacingF = 8.0				; contour color spacing

			; Change the output type and name

		wks  = gsn_open_wks("pdf",name)   	; output using X11
		;wks  = gsn_open_wks("eps","cam-fv-test22heightlon3")   ; output using eps
		;wks  = gsn_open_wks("png","cam-fv-test22heightlon3")   ; output using png
		gsn_define_colormap(wks,"gui_default")  		; Colormap

	end if
	end if
	end if

	res1@gsnCenterString = title
	plot(0) = gsn_csm_contour(wks,var1(:,:),res1)		; plot var1

	;res1@gsnCenterString = "b) t = 3600 s"
     	;plot(1) = gsn_csm_contour(wks,var2(:,:),res1)		; plot var2

	;res1@gsnCenterString = "c) t = 7200 s"
	;plot(2) = gsn_csm_contour(wks,var3(:,:),res1)		; plot var3

	gsn_panel(wks,plot,(/1,1/),pres)			; 1x1 plot

		; Add height label

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.02
	gsn_text_ndc(wks,"Height (m)",0.16,0.48,txres3)

	frame(wks)




end
