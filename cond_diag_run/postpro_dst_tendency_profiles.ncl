
  species  = "dst"
  casename = "dst_compy_F2010_ne30pg2_r05_oECv3_e3smv2_1d66ec33"
  FilePath = "/compyfs/zhan391/cond_diag_scratch/"+casename + "/run/"
  FileName = casename+".eam.h1.0001-01-01-00000.nc"
  QoI_name = "cnd01_"+(/"dst_a3"/)
  chkpt    = (/"AERDRYRM","PBCINI", "STCLD"/) 
  label    = (/"Dry removal","Resolved transport","Turb. mixing and activation"/) 

  dtime    = 1800.
  v_suff   = ""    ; or "" for non-integrals
  inc_suff = "_inc"  ; or "" for state
  inc_fac  = 1./dtime
  unit     = "kg kg~S~-1~N~ s~S~-1~N~"

  nqoi   = dimsizes(QoI_name)
  nchkpt = dimsizes(chkpt)

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
  wks  = gsn_open_wks("pdf",species+"_profiles")
  gsn_define_colormap(wks,"StepSeq25")
  polycolor = (/15,20,26,10/)

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@vpHeightF = 0.4

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF    = 0.025
  res@tiXAxisFontHeightF   = res@gsnStringFontHeightF
  res@tiYAxisFontHeightF   = res@gsnStringFontHeightF

  res@tiMainFont = 22 ;helvetica-bold

  res@trYReverse = True
 ;res@trYLog = True
 ;res@trYMaxF = 1050.
  res@trYMaxF = 1000.
  res@trYMinF = 400.

  res@trXMaxF = 6E-16
  res@trXMinF =-2E-16 

  if (species.eq."dst") then
  res@trXMaxF =  20E-14
  res@trXMinF = -20E-14 
  end if

  res@gsnXRefLine = 0. 

  res@tmXTOn = False
  res@tmYROn = False

  res@xyMarkLineMode = "MarkLines"
  res@xyMarker = 1 
 ;res@xyMarkerSizeF = 0.007
  res@xyLineColor = "black"

  ;-----------------------------
  ; open data file
  ;-----------------------------
  File = addfile( FilePath+FileName,"r")
  area = tofloat(File->area)
  lev  = File->lev
  lev@long_name = "Norminal pressure (hPa)"

  ntime = dimsizes(File->time)
  nlev  = dimsizes(File->lev)

  glbavg_of_modes = new( (/nchkpt,nqoi,nlev,ntime/), "float" )

  dum = new( (/nqoi,nchkpt,nlev/), "graphic" )
  xp = new( 2*nlev, "float" )
  yp = new( 2*nlev, "double" )


  do iqoi = 0,nqoi-1

     print("")
     do ichkpt = 0,nchkpt-1

        varname = QoI_name(iqoi)+"_"+ chkpt(ichkpt)+inc_suff

        do itime = 2,ntime-1

           var3d   = File->$varname$(itime,:,:)*inc_fac

           do ilev = 0,nlev-1
              glbavg_of_modes(ichkpt,iqoi,ilev,itime) = sum(var3d(ilev,:)*area)/sum(area)
           end do
           print(varname+" glbavg, itime:"+itime)
          ;print(glbavg_of_modes(ichkpt,iqoi,:,itime))
        end do

        res@tiMainString = label(ichkpt) +" ("+unit+")"
        res@gsnCenterString = "glbavg( "+varname+" / dtime )"
        res@tiXAxisString = "Mass mixing ratio tendency ("+unit+")"

        zavg = dim_avg(glbavg_of_modes(ichkpt,iqoi,:,:))
        zstd = dim_stddev(glbavg_of_modes(ichkpt,iqoi,:,:))

        plot(iqoi,ichkpt) = gsn_csm_xy(wks,zavg,lev,res)    

      ;-- show +- sigma range as horizontal lines ----
        polyres = True
        polyres@gsLineColor = polycolor(ichkpt)-1
        polyres@gsLineThicknessF = 3 
        polyres@tfPolyDrawOrder   = "Predraw"
        do ilev = 0,nlev-1
           xx = (/zavg(ilev)-zstd(ilev),zavg(ilev)+zstd(ilev)/)
           yy = (/lev(ilev),lev(ilev)/)
           dum(iqoi,ichkpt,ilev) = gsn_add_polyline(wks,plot(iqoi,ichkpt),xx,yy, polyres)
        end do

      ;-- show +- sigma range as shading ----
      ; do k = 0,nlev-1
      ;    xp(k)          = zavg(k) + zstd(k)
      ;    yp(k)          = lev(k) 
      ;    yp(2*nlev-1-k) = lev(k) 
      ;    xp(2*nlev-1-k) = zavg(k) - zstd(k) 
      ; end do
      ; polyres = True
      ; polyres@gsFillColor = polycolor(ichkpt)
      ; polyres@tfPolyDrawOrder   = "Predraw"
      ; dum(iqoi,ichkpt,0) = gsn_add_polygon(wks,plot(iqoi,ichkpt),xp,yp, polyres)
      ;-----------

     end do

     resP = True
     resP@gsnMaximize = False ;True
     resP@gsnPanelYWhiteSpacePercent = 5 
     resP@gsnPanelMainString = QoI_name(iqoi)
     resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.8
     gsn_panel(wks,plot(iqoi,:),(/nchkpt,1/),resP)

  end do 
  print("")


