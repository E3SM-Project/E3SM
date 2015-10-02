
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/cism_tests/dome30/bench/dome.out.nc","r")
;  VAR1    = addfile("/global/scratch/sd/evanskj/cism_tests/dome30/data/dome.1.nc","r")
;  VAR4    = addfile("/global/scratch/sd/evanskj/cism_tests/dome30/data/dome.4.nc","r")
;  VAR9    = addfile("/global/scratch/sd/evanskj/cism_tests/dome30/data/dome.9.nc","r")
;  VAR15    = addfile("/global/scratch/sd/evanskj/cism_tests/dome30/data/dome.15.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/cism_tests/livv/plots/dome30"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  thkstock   = STOCK->thk(10,:,:)
  velstock   = STOCK->velnorm(10,1,:,:)

  thkvarbl  = VAR1->thk(10,:,:)
  thkdiff1 = new((/30,30/),float);
  thkdiff1 = thkvarbl(:,:)-thkstock(:,:);

  thkvarbl  = VAR4->thk(10,:,:)
  thkdiff4 = new((/30,30/),float);
  thkdiff4 = thkvarbl(:,:)-thkstock(:,:);

  thkvarbl  = VAR9->thk(10,:,:)
  thkdiff9 = new((/30,30/),float);
  thkdiff9 = thkvarbl(:,:)-thkstock(:,:);

  thkvarbl  = VAR15->thk(10,:,:)
  thkdiff15 = new((/30,30/),float);
  thkdiff15 = thkvarbl(:,:)-thkstock(:,:);

  velvarbl  = VAR1->velnorm(10,1,:,:)
  veldiff1 = new((/29,29/),float);
  veldiff1 = velvarbl(:,:)-velstock(:,:);

  velvarbl  = VAR4->velnorm(10,1,:,:)
  veldiff4 = new((/29,29/),float);
  veldiff4 = velvarbl(:,:)-velstock(:,:);

  velvarbl  = VAR9->velnorm(10,1,:,:)
  veldiff9 = new((/29,29/),float);
  veldiff9 = velvarbl(:,:)-velstock(:,:);

  velvarbl  = VAR15->velnorm(10,1,:,:)
  veldiff15 = new((/29,29/),float);
  veldiff15 = velvarbl(:,:)-velstock(:,:);

;  maxv1=max(thkstock(:,:));
;  print(maxv1);
;  printVarSummary(thkstock)
;  maxv50=max(thkvarbl(:,:));
;  print(maxv50);
;  printVarSummary(thkvarbl);

  thkstock@long_name   = "Benchmark"
  thkdiff1@long_name   = "Difference 1 proc"
  thkdiff4@long_name   = "Difference 4 proc"
  thkdiff9@long_name   = "Difference 9 proc"
  thkdiff15@long_name  = "Difference 15 proc"

  velstock@long_name   = "Benchmark"
  veldiff1@long_name   = "Difference 1 proc"
  veldiff4@long_name   = "Difference 4 proc"
  veldiff9@long_name   = "Difference 9 proc"
  veldiff15@long_name  = "Difference 15 proc"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plotthk = new(5,graphic)
  plotvel = new(5,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@cnInfoLabelOn       = False         ; turn off cn info label
;  res@gsnMaximize         = True          ; make as large as possible

  plotthk(0) = gsn_contour(wks,thkstock,res) 
  plotthk(1) = gsn_contour(wks,thkdiff1,res) 
  plotthk(2) = gsn_contour(wks,thkdiff4,res) 
  plotthk(3) = gsn_contour(wks,thkdiff9,res) 
  plotthk(4) = gsn_contour(wks,thkdiff15,res) 

  plotvel(0) = gsn_contour(wks,velstock,res) 
  plotvel(1) = gsn_contour(wks,veldiff1,res) 
  plotvel(2) = gsn_contour(wks,veldiff4,res) 
  plotvel(3) = gsn_contour(wks,veldiff9,res) 
  plotvel(4) = gsn_contour(wks,veldiff15,res) 

;************************************************
; create panel
;************************************************

  resthk = True                             ; panel resources
  resthk@gsnFrame = False                 
  resthk@gsnPanelTop = 0.9                 
  resthk@gsnPanelBottom = 0.45                 
;  resthk@gsnMaximize = True                 ; make as large as possible
  resthk@txString    = "Thickness"          ; title string (make blank to remove)
  resthk@gsnPanelFigureStrings = (/"bench)","1)","4)","9)","15)"/)   ; add strings to panel

  resvel = True                             ; panel resources
  resvel@gsnFrame = False                 
  resvel@gsnPanelTop = 0.45                 
  resvel@gsnPanelLabelBar = True		; make a common contour bar
;  resvel@gsnMaximize = True                 ; make as large as possible
  resvel@txString    = "Norm of Velocity"   ; title string (make blank to remove)
  resvel@gsnPanelFigureStrings = (/"bench)","1)","4)","9)","15)"/)   ; add strings to panel

  gsn_panel(wks,plotthk,(/1,5/),resthk)                 
  gsn_panel(wks,plotvel,(/1,5/),resvel)                  

  frame(wks)

end

