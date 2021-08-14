
  casename = "U10_compy_F2010_ne30pg2_r05_oECv3_e3smv2_1d66ec33"
  FilePath = "/compyfs/zhan391/cond_diag_scratch/"+casename + "/run/"
  FileName = casename+".eam.h0.0001-01.nc"

  FillValueInDataFile = 0.
  FillValue4Postpro   = -999.

  cnd = "cnd"+(/"01","02","03","04"/)

  species = "ncl"
  metric_name = (/"U10","U10","U10","U10"/)
  qoi_name    = (/"SFncl_a2_CHEM","SFncl_a1_CHEM","SFncl_a3_CHEM"/) 
  qoi_longname = (/"Aitken mode emis.","accumulation mode emis.","coarse mode emis."/) 

  cnMin = (/0., 0., 0e-10/)
  cnMax = (/1.,15., 1e-8/)
  cnInt = (/.1, 1., 1e-9/)

  ncnd   = dimsizes(cnd)
  nqoi   = dimsizes(qoi_name)


  ;-----------------------------
  ; for plotting
  ;-----------------------------
  wks  = gsn_open_wks("pdf","U10_and_SS_emission_maps")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False


  res@mpFillOn    = False                       ; turn map fill off
  res@mpOutlineOn = True                        ; turn the map outline on
  res@cnFillOn    = True                        ; turn on color fill
  res@cnLinesOn   = False                       ; turn off the contour lines
  res@cnLineLabelsOn = False                    ; turn the line labels off

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF =  res@gsnStringFontHeightF * 1.1

  res@lbOrientation = "vertical"
  res@lbLabelFontHeightF = res@gsnStringFontHeightF * 0.8

  res@mpProjection = "Robinson"
  res@mpDataBaseVersion    = "LowRes"    ; better map outlines
  res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
  res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
  res@mpPerimOn            = False          ; turn off map perimeter
  res@mpFillDrawOrder      = "PostDraw"     ; draw map fill last

  ;-----------------------------
  ; open data file
  ;-----------------------------
  File = addfile( FilePath+FileName,"r")
  area = tofloat(File->area)
  lat  = File->lat
  lon  = File->lon

  res@sfXArray = lon
  res@sfYArray = lat

  ncol = dimsizes(lon)

  ;-----------------------------
  ; start data processing and plotting
  ;-----------------------------
  cndplots = new( (/ncnd,2/),"graphic" )        ; frq and composite mean of the metric
  qoiplots = new( (/ncnd,nqoi,2/), "graphic" )  ; composite mean and contribution to total 
                                                ; (i.e., sum over all conditions)

  do icnd = 0,ncnd-1
 
     print("") 
     print("condition "+icnd)

     ; read and plot frequency (the "flag" array)

     varname = cnd(icnd)+"_"+metric_name(icnd)+"_flag"
     var2d = File->$varname$(0,:)
     delete(var2d@long_name)

     res@cnLevelSelectionMode = "ManualLevels"
     res@cnMinLevelValF = 0.1 
     res@cnMaxLevelValF = 1. 
     res@cnLevelSpacingF = 0.1 

     ; compute global mean of the frequency of occurrence. 0 is considered 0 so 
     ; we use here the original flag array.

     glbavg = sum(var2d*area)/sum(area)
     res@gsnCenterString = varname
     res@tiMainString = "Frequency of occurrence of "+cnd(icnd)+" (glbavg = "+sprintf("%5.2g",glbavg)+")"

     ; Masked version of the flag variable.
     ;
     ; Fillvalue of EAM output variables from the conditional sampling tool was set to 0 
     ; to allow for the calculation of frequency of occurrence and composite mean. 
     ; Here we need a copy of the "flag" variable with 0's marked as missing values 
     ; to avoid division by zero when "flag" is used in the denominator. 
     ;  
     ; It looks like NCL does not support 0 as fillvalue in plotting, so we use
     ; a non-zero FillValue4Postpro for the masked variable

     flag_masked = where(var2d.eq.FillValueInDataFile, FillValue4Postpro, var2d)
     flag_masked@_FillValue = FillValue4Postpro 

     ; when plotting, we treat 0 as missing values for the flag array to be consistant with 
     ; other plots created later in this script

     cndplots(icnd,0) =  gsn_csm_contour_map(wks,flag_masked,res)

     ; also created a masked version of the cell area array, to be used for the calculation of
     ; global averages 

     area_masked = where( ismissing(flag_masked), FillValue4Postpro, area )
     area_masked@_FillValue = FillValue4Postpro

     ;---------------------------------------------
     ; read and plot composite mean of the metric
     ;---------------------------------------------
     varname = cnd(icnd)+"_"+metric_name(icnd)

     delete(var2d)
     var2d = File->$varname$(0,:)
     var2d_compmean_masked = where( ismissing(flag_masked), FillValue4Postpro, var2d/flag_masked )
     var2d_compmean_masked@_FillValue = FillValue4Postpro

     ;average of the composite mean

     res@tiMainString = "Composite avg. of "+metric_name(icnd) +" (m s~S~-1~N~)" +", "+cnd(icnd) 
     res@gsnCenterString = varname+" / "+cnd(icnd)+"_"+metric_name(icnd)+"_flag" 

     res@cnLevelSelectionMode = "ManualLevels"
     res@cnMinLevelValF  = 0.
     res@cnMaxLevelValF  = 15.
     res@cnLevelSpacingF = 1. 

     cndplots(icnd,1) =  gsn_csm_contour_map(wks,var2d_compmean_masked,res)

     ;---------------------
     ; read and plot QoIs
     ;---------------------
     do iqoi = 0,nqoi-1

        varname = cnd(icnd)+"_"+qoi_name(iqoi)

        delete(var2d)
        var2d = File->$varname$(0,:)

        ; here we mask out not only the cells with frq = 0 but also cells with zero emission
        ; as the latter seem to cause plotting issues.

       ;var2d_compmean_masked = where( ismissing(flag_masked).or.(var2d.eq.FillValueInDataFile), \
        var2d_compmean_masked = where( ismissing(flag_masked), \
                                       FillValue4Postpro, var2d/flag_masked )
        var2d_compmean_masked@_FillValue = FillValue4Postpro

        res@tiMainString =  "Composite avg. of "+qoi_longname(iqoi)+" (kg m~S~-2~N~ s~S~-1~N~)" +", "+cnd(icnd)
        res@gsnCenterString = varname+" / "+ cnd(icnd)+"_"+metric_name(icnd)+"_flag"

        res@cnLevelSelectionMode = "ExplicitLevels"
         if (species.eq."ncl") then

           if (iqoi.eq.2) then
           res@cnLevels = (/1e-12,2e-12,5e-12, \
                            1e-11,2e-11,5e-11,1e-10,2e-10,5e-10,1e-9/)
           else
           res@cnLevels = (/1e-15,2e-15,5e-15,1e-14,2e-14,5e-14,1e-13,2e-13,5e-13,1e-12,2e-12,5e-12, \
                            1e-11,2e-11,5e-11,1e-10,2e-10,5e-10,1e-9,2e-9,5e-9/)
           end if

         end if

         if (species.eq."dst") then
           res@cnLevels = (/1e-12,2e-12,5e-12, \
                            1e-11,2e-11,5e-11,1e-10,2e-10,5e-10,1e-9,2e-9,5e-9, \
                            1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7/)
         end if

        qoiplots(icnd,iqoi,0) =  gsn_csm_contour_map(wks,var2d_compmean_masked,res)

        delete(res@cnLevels)

        ; contribution to total

        totname = "cnd04_"+qoi_name(iqoi) 
        total = File->$totname$(0,:)

        ; need a masked version of total to avoid division by zero

        total_masked = where( total.eq.FillValueInDataFile, FillValue4Postpro, total )
        total_masked@_FillValue = FillValue4Postpro 

        ; calculate the contibution of each condition to total
        contrib = where( ismissing(total_masked), FillValue4Postpro, var2d/total_masked)
        contrib@_FillValue = FillValue4Postpro


        glbavg = sum(var2d*area)/sum(total*area)

        res@tiMainString = "Contribution to total "+qoi_longname(iqoi)+" (glbavg = "+sprintf("%5.2g",glbavg)+")"
        res@gsnCenterString = varname+" / "+totname

        res@cnLevelSelectionMode = "ManualLevels"
        res@cnMinLevelValF = 0.1
        res@cnMaxLevelValF = 1. 
        res@cnLevelSpacingF = 0.1 
        qoiplots(icnd,iqoi,1) =  gsn_csm_contour_map(wks,contrib,res)

     end do
  end do 

     resP = True
   ; resP@gsnMaximize = True
     resP@gsnPanelYWhiteSpacePercent = 10
     resP@gsnPanelXWhiteSpacePercent = 8
     resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.3
     gsn_panel(wks,ndtooned(cndplots),(/ncnd,2/),resP)

 ;do iqoi=0,nqoi-1
  iqoi=nqoi-1
     gsn_panel(wks,ndtooned(qoiplots(:,iqoi,:)),(/ncnd,2/),resP)
 ;end do


