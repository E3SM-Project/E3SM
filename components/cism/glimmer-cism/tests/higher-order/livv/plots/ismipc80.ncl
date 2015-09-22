
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/higher_order/reg_test/ismip-hom-c/80km/bench/ishom.c.80km.JFNK.out.nc","r")
;  VAR1    = addfile("/global/scratch/sd/evanskj/higher_order/reg_test/ismip-hom-c/80km/data/ishom.c.80km.JFNK.out.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/higher_order/livv/plots/ismipa"


;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  vstock   = STOCK->vvel(0,1,:,:)
  ustock   = STOCK->uvel(0,1,:,:)
  velstock   = STOCK->velnorm(0,1,:,:)

  vvarbl  = VAR1->vvel(0,1,:,:)
  vdiff1 = new((/40,40/),float);
  vdiff1 = vvarbl(:,:)-vstock(:,:);

  uvarbl  = VAR1->uvel(0,1,:,:)
  udiff1 = new((/40,40/),float);
  udiff1 = uvarbl(:,:)-ustock(:,:);

  velvarbl  = VAR1->velnorm(0,1,:,:)
  veldiff1 = new((/40,40/),float);
  veldiff1 = velvarbl(:,:)-velstock(:,:);

;  maxv1=max(thkstock(:,:));
;  print(maxv1);
;  printVarSummary(thkstock)
;  maxv50=max(thkvarbl(:,:));
;  print(maxv50);
;  printVarSummary(thkvarbl);

;  thkstock@long_name   = "Benchmark"
;  thkdiff1@long_name   = "Difference with 1 proc"

  vstock@long_name   = "Benchmark"
  vdiff1@long_name   = "Difference with 1 proc"

  ustock@long_name   = "Benchmark"
  udiff1@long_name   = "Difference with 1 proc"

  velstock@long_name   = "Benchmark"
  veldiff1@long_name   = "Difference with 1 proc"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plotv = new(2,graphic)
  plotu = new(2,graphic)
  plotvel = new(2,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@cnInfoLabelOn       = False         ; turn off cn info label
;  res@gsnMaximize         = True          ; make as large as possible

  plotv(0) = gsn_contour(wks,vstock,res) 
  plotv(1) = gsn_contour(wks,vdiff1,res) 
  plotu(0) = gsn_contour(wks,ustock,res) 
  plotu(1) = gsn_contour(wks,udiff1,res) 
  plotvel(0) = gsn_contour(wks,velstock,res) 
  plotvel(1) = gsn_contour(wks,veldiff1,res) 

;************************************************
; create panel
;************************************************

  resv = True                             ; panel resources
  resv@gsnFrame = False                 
  resv@gsnPanelTop = 1.0                 
  resv@gsnPanelBottom = 0.66                 
;  resv@gsnMaximize = True                 ; make as large as possible
  resv@txString    = "V velocity"          ; title string (make blank to remove)
  resv@gsnPanelFigureStrings = (/"bench)","diff)"/)   ; add strings to panel

  resvel = True                             ; panel resources
  resu = True                             ; panel resources
  resu@gsnFrame = False                 
  resu@gsnPanelTop = 0.62                 
  resu@gsnPanelBottom = 0.33                
;  resu@gsnMaximize = True                 ; make as large as possible
  resu@txString    = "U velocity"          ; title string (make blank to remove)
  resu@gsnPanelFigureStrings = (/"bench)","diff)"/)   ; add strings to panel

  resvel = True                             ; panel resources
  resvel@gsnFrame = False                 
  resvel@gsnPanelTop = 0.29                 
  resvel@gsnPanelBottom = 0.0                 
  resvel@gsnPanelLabelBar = False		; make a common contour bar
;  resvel@gsnMaximize = True                 ; make as large as possible
  resvel@txString    = "Norm of Velocity"   ; title string (make blank to remove)
  resvel@gsnPanelFigureStrings = (/"bench)","diff)"/)   ; add strings to panel

  gsn_panel(wks,plotv,(/1,2/),resv)                 
  gsn_panel(wks,plotu,(/1,2/),resu)                 
  gsn_panel(wks,plotvel,(/1,2/),resvel)                  

  frame(wks)

end

