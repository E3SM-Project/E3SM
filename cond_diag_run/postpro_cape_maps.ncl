
  CaseName = "cape_compy_F2010_ne30pg2_r05_oECv3_e3smv2_1d66ec33"
  FilePath = "/compyfs/zhan391/cond_diag_scratch/"+CaseName+"/run/"
  FileName = CaseName+ ".eam.h0.0001-01.nc"
  QoI_name = "cnd01_"+(/"CAPE"/)

  chkpt    = (/"RAD",      "PACEND", "DYNEND",  "DEEPCU", "STCLD"/)
  label    = (/"radiation","AC physics","dynamics","deep Cu","CLUBB+MG2"/) 

  dtime    = 1800.
  v_suff   = ""             ; "" for original QoI;  "_v" for vertical integral of 3D field

  inc_suff = (/"_inc",""/)  ; "_inc" is increment, "" gives state
  inc_fac  = (/1.,    1./)  ; set to 1./dtime if converting to tendency

  nqoi   = dimsizes(QoI_name)
  nchkpt = dimsizes(chkpt)
  ninc   = dimsizes(inc_suff)

  glbavg = new( (/nchkpt,nqoi/), "float" )

  ;-----------------------------
  ; for plotting
  ;-----------------------------
  plot = new( (/nqoi,nchkpt/),"graphic" )
  wks  = gsn_open_wks("pdf","cape_maps")
  gsn_define_colormap(wks,"testcmap")

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@mpProjection         = "Robinson"
  res@mpDataBaseVersion    = "LowRes"    ; better map outlines
  res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
  res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
  res@mpPerimOn            = False          ; turn off map perimeter
  res@mpFillDrawOrder      = "PostDraw"     ; draw map fill last
  res@mpGeophysicalLineColor = "gray"

  res@mpFillOn    = False                       ; turn map fill off
  res@mpOutlineOn = True                        ; turn the map outline on
  res@cnFillOn    = True                        ; turn on color fill
  res@cnLinesOn   = False                       ; turn off the contour lines
  res@cnLineLabelsOn = False                    ; turn the line labels off



  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF    = 0.025

  res@lbOrientation = "vertical"
  res@lbLabelFontHeightF = res@gsnStringFontHeightF * 0.8
 ;res@pmLabelBarWidthF    = 0.1      ; default is shorter
  res@pmLabelBarHeightF   = 0.41       ; default is taller
 ;res@pmLabelBarOrthogonalPosF = 0.15
 ;res@pmLabelBarParallelPosF = 0.6

  ;-----------------------------
  ; open data file
  ;-----------------------------
  File = addfile( FilePath+FileName,"r")
  area = tofloat(File->area)
  lat  = File->lat
  lon  = File->lon

  res@sfXArray = lon
  res@sfYArray = lat

  do iqoi = 0,nqoi-1     ; QoI loop

     do iinc = 0,ninc-1  ; increment and/or state

        ;-----------------------------------------------------------
        ; Set contour levels (different for state and increment)
        ;-----------------------------------------------------------
        if (inc_suff(iinc).eq."_inc") then

           add_string = QoI_name(iqoi)+" increment due to "
           res@cnLevelSelectionMode = "ExplicitLevels"
           res@cnLevels = (/-100,-50,-20,-10,-5,-2,2,5,10,20,50,100/)

        else

           add_string = QoI_name(iqoi)+" after "
           res@cnLevelSelectionMode = "ManualLevels"
           res@cnMinLevelValF =  30.
           res@cnMaxLevelValF =  300.
           res@cnLevelSpacingF = 30.

        end if

        ;----------------------------------------------------
        ; Read data, plot global map, calculate global mean
        ;----------------------------------------------------
        do ichkpt = 0,nchkpt-1 ; different atm. processes

           varname = QoI_name(iqoi)+v_suff+"_"+ chkpt(ichkpt)+inc_suff(iinc)
           var2d   = File->$varname$(0,:)*inc_fac(iinc)

           glbavg(ichkpt,iqoi) = sum(var2d*area)/sum(area)

           res@tiMainString = add_string + label(ichkpt)
           res@gsnCenterString = varname+", glbavg = "+glbavg(ichkpt,iqoi)
           plot(iqoi,ichkpt) = gsn_csm_contour_map(wks,var2d,res)    

        end do

        ;-----------------------------------------------------------------
        ; Create paneled plot: all atm. processes on the same page;
        ; Increment and state or different QoIs each has their own page.
        ;-----------------------------------------------------------------
        resP = True
        resP@gsnMaximize = True
        resP@gsnPanelYWhiteSpacePercent = 8 
       ;resP@gsnPanelMainString = QoI_name(iqoi)
        resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.3

        gsn_panel(wks,plot(iqoi,:),(/nchkpt,1/),resP)

     end do 

     print("")
     print("global mean = "+dim_sum(glbavg)+" after "+chkpt)

  end do


