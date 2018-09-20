;================================================;
;  test200-height-lon.ncl
;
;	Example ncl script to produce the height-lon 
;	plots for test 200.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

	; This produces the height-lon plots for test 200

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

		; replace your models file name below

;	f    = addfile("../cam-fv.200.medium.L30.latlon.hydro.4th_order_div_damping.nc","r")

;  f = addfile("d2-hydro.nc", "r")
;  f= addfile("d2preqx-ne30.nc","r")


;  name = "d2preqx-ne30"
;  name = "d2thetah-ne30"
;   name = "d2thetah-ne30-upwind"
;  title = "Hydrostatic eqns, original formulation"
;  title = "Hydrostatic eqns, Exner p., upwinding" 
;   name = "nonh-ne30"
;   title = "Nonhydrostatic eqns"
;  name = "d2non-ne30-discrete-phi"
;  title = "Nonhydro, discrete phi"
  name = "movies/dcmip2012_test2_01.nc"
  title = ""


  f = addfile(name+".nc","r")
        minuplot = -.25
        maxuplot = .25
		; Input useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

		; We want to plot zonal and vertical velocity at the equator
		; Note that CAM-FV uses vertical pressure velocity

	varload  = f->u(:,:,:,:)
        maxu = max(varload)
        minu = min(varload)
        print("max u    ="+maxu+" min u="+minu)                                   
	varload2 = f->omega(:,:,:,:)

		; We want the plots at the equator at time 6 days
		; For 6hr intervals this corresponds to index 24
		; The equator is at nlat/2 

	equator = nlat/2

		; Input the variables at this time

       day6 = closest_val( 6.0, f->time) 
       if (day6 .gt. dimsizes(f->time)-1 ) then
         day6 = dimsizes(f->time) -1 
       end if
       ;print("times = "+f->time)
       print("plotting n="+day6+" time(n) = "+f->time(day6)+"d")


	var1 = varload(day6,:,equator,:)
	var2 = varload2(day6,:,equator,:)

		; Delete loaded data

	delete(varload)
	delete(varload2)

		; Set up height coordinate (assuming constant spacing)
		; Note that for CAM-FV pressure coordinates are used, 
		; therefore the height levels are upside-down 
		; i.e. for 30 vertical levels index 0 is at 11800 m 
		; and index 29 is at 200 m.

        maxu = max(var1)
        minu = min(var1)
        print("crosssection: max u    ="+maxu+" min u="+minu)                                   


	z = new((/nlev/),double)
	dz = 12000.0/nlev

	do kk=0,nlev-1

		z(kk) = (12000.0-0.5*dz)-(kk)*dz  

  	end do	

		; We don't want long-name in the plots

	var1@long_name = " "
	var2@long_name = " "

		; We now produce the plot

	plot  = new (1, graphic)               			; define plot - need 2 panels

	res1                      = True
	res1@gsnDraw              = False          		; panel plot 
	res1@gsnFrame             = False			; don't draw yet
	res1@cnFillOn             = True
	res1@cnLinesOn            = False			; turn off contour lines
	res1@gsnSpreadColors      = True
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = ""
	res1@vpWidthF             = 0.38
	res1@vpHeightF            = 0.19
	res1@cnLevelSelectionMode = "ManualLevels"
	res1@cnInfoLabelOn        = False			; give each
	res1@cnLineLabelsOn       = False			; panel its own
	res1@lbLabelBarOn         = True			; label bar

	pres                          = True
	pres@gsnMaximize              = True  

	res1@sfXArray                 = lon           		; uses lon as plot x-axis
	res1@sfYArray                 = z             		; uses z for y axis
	res1@trYReverse               = False          		; reverses y-axis, false
	res1@tiYAxisString            = " "   	       		; y-axis title
	res1@tiXAxisString            = " "             	; x-axis title

	pltTitle             ="test 200, t=6 days, min/max="+sprintf("%5.2f", minu)+"/"+sprintf("%5.2f", maxu) 		; Plot title if required
	pres@txString        = pltTitle

		; Change the output type and name

	;wks  = gsn_open_wks("X11","cam-fv-test200heightlon")   	; output using X11
	wks  = gsn_open_wks("pdf",name)   ; output using eps
	;wks  = gsn_open_wks("png","cam-fv-test200heightlon")   ; output using png
	gsn_define_colormap(wks,"gui_default")  		; Colormap

		; Add height label

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.02
	gsn_text_ndc(wks,"Height (m)",0.03,0.48,txres3)

	txres3@txAngleF       = 0.
	gsn_text_ndc(wks,"Longitude",0.5,0.02,txres3)

		; Contour spacing for u

	res1@cnMaxLevelValF  = maxuplot;0.1				; max contour color label
	res1@cnMinLevelValF  = minuplot; -0.1				; min contour color label
	res1@cnLevelSpacingF = (maxuplot-minuplot)/20.	; contour color spacing
		
	res1@gsnCenterString = "a) U"
	plot(0) = gsn_csm_contour(wks,var1(:,:),res1)		; plot var1

		; Contour spacing for w

;	res1@cnMaxLevelValF  = 0.0003				; max contour color label
;	res1@cnMinLevelValF  = -0.0003				; min contour color label
;	res1@cnLevelSpacingF = 0.0001				; contour color spacing

	;res1@gsnCenterString = "b) OMEGA"
     	;plot(1) = gsn_csm_contour(wks,var2(:,:),res1)		; plot var2
       ; this is with omega
	;gsn_panel(wks,plot,(/2,1/),pres)			; 2x1 plot

	gsn_panel(wks,plot,(/1,1/),pres)			; 2x1 plot
	frame(wks)




end
