
  ldebug = True
  ldebug = False

  RHIcrit = (/125,-1/)
  casename = "RHI"+RHIcrit(0)+"_compy_F2010_ne30pg2_r05_oECv3_e3smv2_1d66ec33"
  FilePath = "/compyfs/zhan391/cond_diag_scratch/"+casename + "/run/"
  FileName = casename+".eam.h0.0001-01.nc"

  FillValue = 0.

   cnd_name = (/"cnd01","cnd02"/)
  flag_name = cnd_name+"_RHI_flag"

  QoI_name = (/"RHI","Q","QSATI"/)
  QoI_longname = (/"RH","Q","QSATI"/)
  QoI_unit     = (/"percent","kg kg~S~-1~N~","kg kg~S~-1~N~"/)
  QoI_format   = (/"0@*+^sg","0@*+^se","0@*+^se"/)

  chkpt    = (/"PBCDIAG","RAD", "PACEND","DYNEND","DEEPCU",\
              "CLDMAC01","CLDMIC01",\
              "CLDMAC02","CLDMIC02",\
              "CLDMAC03","CLDMIC03",\
              "CLDMAC04","CLDMIC04",\
              "CLDMAC05","CLDMIC05",\
              "CLDMAC06","CLDMIC06" \
              /)

  label = chkpt

  color    = (/194,109, 96, 81, 185,\
               162,  16, \
               154,  22, \
               146,  28, \
               138,  34, \
               130,  40, \
               132,  48  \
             /)
  marker   = (/ 16,16, 16, 16, 16, \
                4, 4, \
                5, 5, \ 
                6, 6, \ 
                7, 7, \ 
                9, 9, \ 
                3, 3  \ 
             /)
  dash     = (/0,0,0,0,0, \
               0,0,\
               5,5,\
               0,0,\
               5,5,\
               0,0,\
               5,5 \
             /)

  ncnd = dimsizes(cnd_name)
  nqoi   = dimsizes(QoI_name)
  nchkpt = dimsizes(chkpt)

 ;--------
  val_inc_suff     = (/"","_inc"/)
  val_inc_longname = (/""," increment"/)
  refvalue         = (/RHIcrit(0),0./)
  trX_maxmin_scale = (/1.1, 1.25/)

  nsuff = dimsizes(val_inc_suff)

  legendonoff = (/"Always","Never"/)
  legendonoff = (/"Always","Always"/)
 ;--------

  File = addfile( FilePath+FileName,"r")

  area = tofloat(File->area)

   levS = 100.
   levE = 1000.
   lev = File->lev({levS:levE})
  nlev = dimsizes(lev)

  composite_mean = new( (/nsuff,nqoi,nchkpt,nlev/),"float" )
  composite_mean@_FillValue = FillValue

;---------------------------------------------------------------------------------
  wks = gsn_open_wks("pdf","RHI_profiles")
  gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")

  resXY = True

 ;resXY@trYLog = True
  resXY@trYReverse = True
  resXY@trYMaxF = 1000.
  resXY@trYMinF = 100.
  resXY@tiMainFont = 22 ;helvetica-bold

  resXY@tmYLMode = "Explicit"
  resXY@tmYLValues = (/1000.,850.,700.,600.,500.,400.,300.,200.,100./)
  resXY@tmYLLabels = (/1000.,850.,700.,600.,500.,400.,300.,200.,100./)
  resXY@tiYAxisString = "Norminal pressure (hPa)" 

  resXY@tmLabelAutoStride  = True

  resXY@tmXTOn = False
  resXY@tmYROn = True

  cidmin = 3
  cidmax = 197
  cidint = (cidmax - cidmin)/nchkpt

     plot = new( (/ncnd,nsuff,nqoi/), "graphic" )
  frqplot = new( ncnd, "graphic" )


do icnd = 0,ncnd-1

;---------------------------------------------------------------------------------
; Read in the "flag" variable and cell area; create tmp arrays
;---------------------------------------------------------------------------------

  flag3d = File->$flag_name(icnd)$(0,{levS:levE},:)
  area3d = conform(flag3d,area,1)

  flag_glbavg = dim_sum(flag3d*area3d)/dim_sum(area3d)

  flag_glbavg_masked = flag_glbavg
  flag_glbavg_masked@_FillValue = FillValue

; plot the global mean of the "flag" variable

  resXY@gsnDraw = True
  resXY@gsnFrame = True

  resXY@tiMainString = "Frequency of occurrence of "+cnd_name(icnd)
  resXY@gsnCenterString = "glbavg( "+flag_name(icnd)+" )"
  resXY@pmLegendDisplayMode = "Never" 
  resXY@xyLineColor = 17
  resXY@xyLineThicknessF = 1.5 
  resXY@tiXAxisString = "Frequency (fraction)" 

  frqplot(icnd) = gsn_csm_xy(wks,flag_glbavg,lev,resXY)

;---------------------------------------------------------------------------------
; prepare for reading and plotting QoIs
;---------------------------------------------------------------------------------

  resXY@gsnDraw = False
  resXY@gsnFrame = False

  resXY@xyLineThicknessF       = 1.5 
 ;resXY@xyMarkLineMode         = "MarkLines"
 ;resXY@xyMarkers              = marker(0:nchkpt-1)
  resXY@xyLineColors           = color(0:nchkpt-1)
  resXY@xyDashPatterns         = dash(0:nchkpt-1) 

  resXY@xyExplicitLegendLabels = label
  resXY@lgPerimOn              = False      ; No legend perimeter.

  resXY@lgItemOrder = ispan( nchkpt-1,0,1 )

 resXY@lgJustification        = "TopRight"
 resXY@pmLegendParallelPosF   = 0.5       ; Move to right
 resXY@pmLegendWidthF         = 0.4       ; Decrease width
 resXY@pmLegendHeightF        = 0.35      ; Decrease height
 resXY@lgBoxMinorExtentF      = 0.1       ; Shorten the legend lines
 resXY@lgLabelFontHeightF     = 0.012     ; Change the font size

  do isuff = 0,nsuff-1

     resXY@pmLegendDisplayMode = legendonoff(isuff) 
     resXY@gsnXRefLine = refvalue(isuff)

     do iqoi = 0,nqoi-1

        resXY@tiMainString = "Composite avg. of "+QoI_name(iqoi)+val_inc_suff(isuff)+", "+cnd_name(icnd)
        resXY@gsnCenterString = "glbavg( "+cnd_name(icnd)+"_"+QoI_name(iqoi)+"_*"+val_inc_suff(isuff)+" )"+\
                                " / glbavg( "+cnd_name(icnd)+"_RHI_flag )"
        resXY@tiXAxisString = QoI_longname(iqoi)+val_inc_longname(isuff)+" ("+QoI_unit(iqoi)+")"

        resXY@pmLegendOrthogonalPosF = -0.78

        resXY@pmLegendParallelPosF = 0.5 
        if (iqoi.ne.2 .and. isuff.eq.nsuff-1 )then
           resXY@pmLegendParallelPosF = 1.1 
        end if


        resXY@tmXBFormat = QoI_format(iqoi) 

        do ichkpt = 0,nchkpt-1

           ;---------------------------------------------------------------------------------
           ; read and process data
           ;---------------------------------------------------------------------------------
           varname = cnd_name(icnd)+"_"+QoI_name(iqoi) +"_"+ chkpt(ichkpt)+val_inc_suff(isuff)

           qoi = File->$varname$(0,{levS:levE},:) ; read data
           print(min(qoi) + max(qoi))
           ; for each vertical layer, calculate the global horizonal mean

           qoi_glbavg = dim_sum(qoi*area3d)/dim_sum(area3d)

           ; the composite (event) mean of a vertical layer is the global (horizonal) mean of the QoI
           ; divided by the global (horizonal) mean frequency of occurrence of the event in that layer.
           ; (Here the use of flag_glbavg_masked instead of flag_glbavg helps avoid division by zero
           ; if the event does not occurr in any grid cell in a layer.) 

           composite_mean(isuff,iqoi,ichkpt,:) = where( ismissing(flag_glbavg_masked), FillValue, qoi_glbavg/flag_glbavg_masked )

           ; debug
           if (ldebug) then
              print("")
              print("flag_glbavg_masked = "+flag_glbavg_masked+", ismissing = "+ismissing(flag_glbavg_masked))
              print(qoi_glbavg/flag_glbavg_masked+" -> "+glbavg(isuff,iqoi,ichkpt,:))
              exit
           end if
           print("done with varname = "+varname)
        end do

        ; single QoI: plot the profiles at all checkpoints

        zmin = min(composite_mean(isuff,iqoi,:,:))
        zmax = max(composite_mean(isuff,iqoi,:,:))
        print(zmax)
        print(zmin)
        if (zmin.lt.0) then
           resXY@trXMinF = zmin * trX_maxmin_scale(isuff)
        else
           resXY@trXMinF = zmin / trX_maxmin_scale(isuff)
        end if
        resXY@trXMaxF = zmax * trX_maxmin_scale(isuff)

        plot(icnd,isuff,iqoi) = gsn_csm_xy(wks,composite_mean(isuff,iqoi,:,:),lev,resXY)

     end do 
     end do 

  ; plot all QoIs of the same condition on the same page

  resP = True
  resP@gsnMaximize = False
  resP@gsnPanelYWhiteSpacePercent = 10
  resP@gsnPanelXWhiteSpacePercent = 3 
  resP@gsnPanelRight= 0.9 
  gsn_panel(wks,ndtooned(plot(icnd,:,:)),(/nsuff,nqoi/),resP)

  delete(resXY@lgItemOrder)
  delete(resXY@trXMaxF)
  delete(resXY@trXMinF)
end do

  ; combine a subset of the panels and print on a separate page 

  figure4paper     = new(6,"graphic")
  figure4paper(0)  = frqplot(0)
  figure4paper(1)  = plot(1,0,0)  ; cnd02_RHI
  figure4paper(2)  = plot(0,0,0)  ; cnd01_RHI
  figure4paper(3:) = plot(0,1,:)  ; cnd01_xxx_inc

  gsn_panel(wks,figure4paper,(/2,3/),resP)
