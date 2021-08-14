
  species  = "dst"
  casename = "dst_compy_F2010_ne30pg2_r05_oECv3_e3smv2_1d66ec33"
  FilePath = "/compyfs/zhan391/cond_diag_scratch/"+casename + "/run/"
  FileName = casename+".eam.h0.0001-01.nc"
  QoI_name = "cnd01_"+(/"dst_a3"/)
  tendfac  = (/1./)
  chkpt    = (/"CFLXAPP","AERDRYRM","PBCINI", "STCLD","AERWETRM"/) 
  label    = (/"Emission","Dry removal","Resolved transport", "Turb. mixing and activation","Wet removal"/) 

  dtime    = 1800.
  v_suff   = "_v"    ; or "" for non-integrals
  inc_suff = "_inc"  ; or "" for state
  inc_fac  = 1./dtime
  unit     = "kg m~S~-2~N~ s~S~-1~N~"

  if (species.eq."dst") then
     cnLevels = (/ -1e-8,-5e-9, -1e-9, -5e-10,-1e-10,-5e-11,\
                    5e-11, 1e-10,5e-10, 1e-9,  5e-9,  1e-8 /) 
  end if

  nqoi   = dimsizes(QoI_name)
  nchkpt = dimsizes(chkpt)

  glbavg_of_modes = new( (/nchkpt,nqoi/), "float" )

  ;-----------------------------
  ; constants
  ;-----------------------------
  pi             = 3.14159265358979323846  ; from E3SM
  earth_radius   = 6.37122e6               ; from E3SM
  earth_sfc_area = 4.*pi*earth_radius^2

  kg2Tg  = 1e-9
  sec2yr = 86400.*365

  ;-----------------------------
  ; for plotting
  ;-----------------------------
  plot = new( (/nqoi,nchkpt/),"graphic" )
  wks  = gsn_open_wks("pdf",species+"_maps")
  gsn_define_colormap(wks,"testcmap")

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@mpProjection         = "Robinson"
  res@mpDataBaseVersion    = "LowRes"       ; better map outlines
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
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF    = 0.025

  res@lbOrientation = "vertical"
  res@lbLabelFontHeightF = res@gsnStringFontHeightF * 0.8
 ;res@pmLabelBarWidthF    = 0.1        ; default is shorter
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

  do iqoi = 0,nqoi-1

     res@cnLevels = cnLevels * tendfac(iqoi)

     print("")
     do ichkpt = 0,nchkpt-1

        varname = QoI_name(iqoi)+v_suff+"_"+ chkpt(ichkpt)+inc_suff
        var2d   = File->$varname$(0,:)*inc_fac
       ;delete(var2d@long_name)

        ; The quantities written out was q*dp/g/dt, unit = kg/kg * (kg/m^3) * m / s = kg/m^2/s

        glbavg_of_modes(ichkpt,iqoi) = sum(var2d*area)/sum(area) * earth_sfc_area * kg2Tg * sec2yr 
        print(QoI_name(iqoi)+": global mean "+sprintf("%7.2f",glbavg_of_modes(ichkpt,iqoi)) \
                                             +" Tg/yr from "+chkpt(ichkpt))

        res@tiMainString = label(ichkpt)+ " ("+unit+")"
        res@gsnCenterString = varname+" / dtime" 
        plot(iqoi,ichkpt) = gsn_csm_contour_map(wks,var2d,res)    

     end do

     resP = True
     resP@gsnMaximize = True
     resP@gsnPanelYWhiteSpacePercent = 8 
     resP@gsnPanelMainString = QoI_name(iqoi)
     resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.6
     gsn_panel(wks,plot(iqoi,:),(/nchkpt,1/),resP)

  end do 

 ;print("")
 ;print("global mean = "+dim_sum(glbavg_of_modes)+" after "+chkpt)


