
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/diagnostic/bench/dome.out.nc","r")
;  VAR1    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/diagnostic/data/dome.1.nc","r")
;  VAR4    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/diagnostic/data/dome.4.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/higher-order/livv/plots/dome30d"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  velstock   = STOCK->velnorm(0,1,:,:)
  thkstock   = STOCK->thk(0,:,:)
  tmpstock   = STOCK->tempstag(0,2,:,:)

  velvarbl  = VAR1->velnorm(0,1,:,:)
  veldiff1 = new((/30,30/),float);
  veldiff1 = velvarbl(:,:)-velstock(:,:);

  velvarbl  = VAR4->velnorm(0,1,:,:)
  veldiff4 = new((/30,30/),float);
  veldiff4 = velvarbl(:,:)-velstock(:,:);

;  thkvarbl  = VAR1->thk(0,:,:)
;  thkdiff1 = new((/31,31/),float);
;  thkdiff1 = thkvarbl(:,:)-thkstock(:,:);

;  thkvarbl  = VAR4->thk(0,:,:)
;  thkdiff4 = new((/31,31/),float);
;  thkdiff4 = thkvarbl(:,:)-thkstock(:,:);

;  maxv1=max(velstock(:,:));
;  maxv1=max(velstock(:,:));
;  maxv1=max(velstock(:,:));
;  print(maxv1);
;  printVarSummary(velstock)
;  maxv50=max(velvarbl(:,:));
;  print(maxv50);
;  printVarSummary(velvarbl);

  vs= velstock@scale_factor * velstock
  vd1= velstock@scale_factor * veldiff1
  vd4= velstock@scale_factor * veldiff4
  vs@long_name   = "Benchmark"
  vd1@long_name   = "Diff 1 proc"
  vd4@long_name   = "Diff 4 proc"

;  ts= thkstock@scale_factor * thkstock
;  td1= thkstock@scale_factor * thkdiff1
;  td4= thkstock@scale_factor * thkdiff4
;  ts@long_name   = "Scaled Benchmark"
;  td1@long_name   = "Scaled Difference 1 proc"
;  td4@long_name   = "Scaled Difference 4 proc"

;  tmpstock@long_name   = "Benchmark"
;  tmpdiff1@long_name   = "Difference 1 proc"
;  tmpdiff4@long_name   = "Difference 4 proc"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plotvel = new(3,graphic)
  plotthk = new(3,graphic)
;  plottmp = new(3,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@cnInfoLabelOn       = False         ; turn off cn info label
;  res@gsnMaximize         = True          ; make as large as possible

  plotvel(0) = gsn_contour(wks,vs,res) 
  plotvel(1) = gsn_contour(wks,vd1,res) 
  plotvel(2) = gsn_contour(wks,vd4,res) 

;  plotthk(0) = gsn_contour(wks,ts,res) 
;  plotthk(1) = gsn_contour(wks,td1,res) 
;  plotthk(2) = gsn_contour(wks,td4,res) 
;  plottmp(0) = gsn_contour(wks,tmpstock,res) 
;  plottmp(1) = gsn_contour(wks,tmpdiff1,res) 
;  plottmp(2) = gsn_contour(wks,tmpdiff4,res) 

;************************************************
; create panel
;************************************************

  resvel = True                             ; panel resources
  resvel@gsnFrame = False                 
  resvel@gsnPanelTop = 1.0                 
  resvel@gsnPanelBottom = 0.60                 
  resvel@gsnPanelLabelBar = True		; make a common contour bar
;  resvel@gsnMaximize = True                 ; make as large as possible
  resvel@txString    = "Norm of Velocity"   ; title string (make blank to remove)
  resvel@gsnPanelFigureStrings = (/"bench)","1)","4)","9)","15)"/)   ; add strings to panel

  resthk = True                             ; panel resources
  resthk@gsnFrame = False                 
  resthk@gsnPanelTop = 0.6                 
  resthk@gsnPanelBottom = 0.0                 
  resthk@gsnPanelLabelBar = True		; make a common contour bar
;  resthk@gsnMaximize = True                 ; make as large as possible
  resthk@txString    = "Thickness"   ; title string (make blank to remove)
  resthk@gsnPanelFigureStrings = (/"bench)","1)","4)","9)","15)"/)   ; add strings to panel

;  restmp = True                             ; panel resources
;  restmp@gsnFrame = False                 
;  restmp@gsnPanelTop = 0.3                 
;  restmp@gsnPanelBottom = 0.0                 
;  restmp@gsnPanelLabelBar = True		; make a common contour bar
;  restmp@txString    = "Temperature"   ; title string (make blank to remove)
;  restmp@gsnPanelFigureStrings = (/"bench)","1)","4)","9)","15)"/)   ; add strings to panel

  gsn_panel(wks,plotvel,(/1,3/),resvel)                  
;  gsn_panel(wks,plotthk,(/2,3/),resthk)                  
;  gsn_panel(wks,plottmp,(/3,3/),restmp)                  

  frame(wks)

end

