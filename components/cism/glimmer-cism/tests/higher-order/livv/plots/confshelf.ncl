
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/confined-shelf/bench/confined-shelf.gnu.JFNK.nc","r")
;  VAR1     = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/confined-shelf/data/confined-shelf.gnu.JFNK.nc","r")
;  PNG      = "/global/scratch/sd/evanskj/higher-order/livv/plots/conf"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  velstock   = STOCK->velnorm(0,1,:,:)
  velvarbl  = VAR1->velnorm(0,1,:,:)
  vs= velstock@scale_factor * velstock
  vv= velvarbl@scale_factor * velvarbl

  vs@long_name   = "Scaled Benchmark"
  vv@long_name   = "Scaled Test Run on 2 Processors"

  veldiff1 = new((/43,42/),float);
  veldiff1 = velvarbl(:,:)-velstock(:,:);

  veldiff1@long_name   = "Difference on 2 procs"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plotvel = new(3,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@cnInfoLabelOn       = False         ; turn off cn info label

  plotvel(0) = gsn_contour(wks,vs,res) 
  plotvel(1) = gsn_contour(wks,vv,res) 
  plotvel(2) = gsn_contour(wks,veldiff1,res) 

;************************************************
; create panel
;************************************************

  resvel = True                             ; panel resources
  resvel@gsnFrame = False                 
  resvel@gsnPanelTop = 1.0                 
  resvel@gsnPanelLabelBar = True		; make a common contour bar
  resvel@txString    = "Norm of Velocity"   ; title string (make blank to remove)
  resvel@gsnPanelFigureStrings = (/"bench)","test)","diff)"/)   ; add strings to panel

  gsn_panel(wks,plotvel,(/1,3/),resvel)                  

  frame(wks)

end
