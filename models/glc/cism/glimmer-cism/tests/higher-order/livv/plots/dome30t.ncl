
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; stand alone read in for testing
;  STOCK    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/evolving/bench/dome.out.nc","r")
;  VAR9    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/evolving/data/dome.9.nc","r")
;  VAR15    = addfile("/global/scratch/sd/evanskj/higher-order/reg_test/dome30/evolving/data/dome.15.nc","r")
;  PNG    = "/global/scratch/sd/evanskj/higher-order/livv/plots/dome30e"

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("png",PNG)                 ; open a png file

  tmpstock   = STOCK->tempstag(10,1,:,:)

  tmpvarbl  = VAR9->tempstag(10,1,:,:)
  tmpdiff9 = new((/31,31/),float);
  tmpdiff9 = tmpvarbl(:,:)-tmpstock(:,:);

  tmpvarbl  = VAR15->tempstag(10,1,:,:)
  tmpdiff15 = new((/31,31/),float);
  tmpdiff15 = tmpvarbl(:,:)-tmpstock(:,:);

  tmpstock@long_name   = "Benchmark"
  tmpdiff9@long_name   = "Difference 9 proc"
  tmpdiff15@long_name  = "Difference 15 proc"

;************************************************
; create plots
;************************************************

  wks = gsn_open_wks("png",PNG)            ; open a png file
  gsn_define_colormap(wks,"wgne15")        ; choose a colormap
  plottmp = new(3,graphic)
  
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnDraw             = False         ; don't draw
  res@gsnFrame            = False         ; don't advance frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@cnInfoLabelOn       = False         ; turn off cn info label
;  res@gsnMaximize         = True          ; make as large as possible

  plottmp(0) = gsn_contour(wks,tmpstock,res) 
  plottmp(1) = gsn_contour(wks,tmpdiff9,res) 
  plottmp(2) = gsn_contour(wks,tmpdiff15,res) 

;************************************************
; create panel
;************************************************

  restmp = True                             ; panel resources
  restmp@gsnFrame = False                 
  restmp@gsnPanelTop = 1.0                 
  restmp@gsnPanelBottom = 0.62                 
  restmp@lbLabelAngleF = 45			; angle of contour bar labels
  restmp@gsnPanelLabelBar = True		; make a common contour bar
  restmp@txString    = "Temperature"   ; title string (make blank to remove)
  restmp@gsnPanelFigureStrings = (/"bench)","9)","15)"/)   ; add strings to panel

  gsn_panel(wks,plottmp,(/1,3/),restmp)                  

  frame(wks)

end

