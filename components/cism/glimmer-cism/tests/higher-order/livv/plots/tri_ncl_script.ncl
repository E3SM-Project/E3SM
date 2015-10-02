
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/cism_tests/gis_5km_long/data/gis_5km.ice2sea.init.nc","r")
;  VAR    = addfile("/global/scratch/sd/evanskj/cism_tests/gis_5km_long/data/gis_5km.ice2sea.51-100.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/cism_tests/livv/plots/gis5km_thk"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  thkstock  = STOCK->thk(0,:,100:380)
  thkvarbl  = VAR->thk(0,:,100:380)

;  maxv1=max(thkstock(:,:));
;  print(maxv1);
;  printVarSummary(thkstock)
;  maxv50=max(thkvarbl(:,:));
;  print(maxv50);
;  printVarSummary(thkvarbl);

  thkdiff = new((/601,281/),float);
  thkdiff = thkvarbl(:,:)-thkstock(:,:);

  thkstock@long_name = "Initial"
  thkvarbl@long_name = "Simulation"
  thkdiff@long_name  = "Difference"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plot = new(3,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnInfoLabelOn       = False         ; turn off cn info label
  res@vpWidthF            = 0.25           ; change aspect ratio of plot
  res@vpHeightF           = 0.7

  thkstock@_FillValue     = 9.96921e+36 ; set missing values for ncl v5.0
  thkvarbl@_FillValue     = 9.96921e+36 ; set missing values for ncl v5.0

  resP = True                             ; panel resources
  resP@gsnMaximize = True                 ; make as large as possible
  resP@txString    = "Thickness"          ; title string (make blank to remove)

  plot(0) = gsn_contour(wks,thkstock,res) 
  plot(1) = gsn_contour(wks,thkvarbl,res)
  plot(2) = gsn_contour(wks,thkdiff,res) 

;************************************************
; create panel
;************************************************

  resP@gsnPanelLabelBar = True                      ; add common colorbar
  resP@gsnPanelFigureStrings = (/"a)","b)","c)"/)   ; add strings to panel
  resP@gsnMaximize = True                           ; make as large as possible
  gsn_panel(wks,plot,(/1,3/),resP)                  ; now draw as one plot

end

