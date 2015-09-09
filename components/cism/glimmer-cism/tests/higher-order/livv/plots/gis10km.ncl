
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/cism_tests/gis_10km/bench/gis_10km.seacism.nc","r")
;  VAR    = addfile("/global/scratch/sd/evanskj/cism_tests/gis_10km/data/gis_10km.seacism.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/cism_tests/livv/plots/gis10km"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  thkstock  = STOCK->thk(1,:,:)
  velstock  = STOCK->velnorm(1,1,:,:)

  thkvarbl  = VAR->thk(1,:,:)
  velvarbl  = VAR->velnorm(1,1,:,:)

;  maxv1=max(thkstock(:,:));
;  print(maxv1);
;  printVarSummary(thkstock)
;  maxv50=max(thkvarbl(:,:));
;  print(maxv50);
;  printVarSummary(thkvarbl);

  thkdiff = new((/280,150/),float);
  thkdiff = thkvarbl(:,:)-thkstock(:,:);

  veldiff = new((/279,149/),float);
  veldiff = velvarbl(:,:)-velstock(:,:);

  thkstock@long_name = "Initial, year 1"
  thkvarbl@long_name = "Simulation, year 1"
  thkdiff@long_name  = "Difference"
  velstock@long_name = "Initial, year 1"
  velvarbl@long_name = "Simulation, year 1"
  veldiff@long_name  = "Difference"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plotthk = new(3,graphic)
  plotvel = new(3,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnInfoLabelOn       = False         ; turn off cn info label
  res@vpWidthF            = 0.25           ; change aspect ratio of plot
  res@vpHeightF           = 0.6

  thkstock@_FillValue     = 9.96921e+36 ; set missing values for ncl v5.0
  thkvarbl@_FillValue     = 9.96921e+36 ; set missing values for ncl v5.0

  plotthk(0) = gsn_contour(wks,thkstock,res) 
  plotthk(1) = gsn_contour(wks,thkvarbl,res)
  plotthk(2) = gsn_contour(wks,thkdiff,res) 

  plotvel(0) = gsn_contour(wks,velstock,res) 
  plotvel(1) = gsn_contour(wks,velvarbl,res)
  plotvel(2) = gsn_contour(wks,veldiff,res) 

;************************************************
; create panel
;************************************************

  resthk = True                             	      ; panel resources
  resthk@gsnMaximize = True                 	      ; make as large as possible
  resthk@txString    = "Thickness"          	      ; title string (make blank to remove)

  resthk@gsnPanelLabelBar = True                      ; add common colorbar
  resthk@gsnFrame = False 
  resthk@gsnPanelTop = 1.0
  resthk@gsnPanelBottom = 0.5
  resthk@gsnPanelFigureStrings = (/"Bench)","Sim)","Diff)"/)   ; add strings to panel
  resthk@gsnMaximize = True                           ; make as large as possible

  resvel = True                                       ; panel resources
  resvel@gsnFrame = False 
  resvel@gsnPanelTop = 0.45
  resvel@gsnPanelBottom = 0.0
  resvel@gsnPanelLabelBar = True                ; make a common contour bar
  resvel@gsnMaximize = True                 	      ; make as large as possible
  resvel@txString    = "Norm of Surface Velocity"     ; title string (make blank to remove)

  resvel@gsnPanelLabelBar = True                      ; add common colorbar
  resvel@gsnPanelFigureStrings = (/"Bench)","Sim)","Diff)"/)   ; add strings to panel
  resvel@gsnMaximize = True                           ; make as large as possible

  gsn_panel(wks,plotvel,(/1,3/),resvel)                  ; now draw as one plot
  gsn_panel(wks,plotthk,(/1,3/),resthk)                  ; now draw as one plot

  frame(wks)

end

