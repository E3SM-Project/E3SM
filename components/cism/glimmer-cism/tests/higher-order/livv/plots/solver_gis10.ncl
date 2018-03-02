load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; Load the NCL file that contains the
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  ; Load the NCL file that contains the
begin

;  nfile = "/global/scratch/sd/evanskj/higher-order/livv/data/newton_gis10km2.asc"
;  lfile = "/global/scratch/sd/evanskj/higher-order/livv/data/fgmres_gis10km2.asc"
;  nbfile = "/global/scratch/sd/evanskj/higher-order/livv/data/newton_gis10km2b.asc"
;  lbfile = "/global/scratch/sd/evanskj/higher-order/livv/data/fgmres_gis10km2b.asc"
;  PNG   = "/global/scratch/sd/evanskj/higher-order/livv/plots/gis10km_iter"

  ndata = asciiread(nfile,(/10/),"integer"); replace 10 with -1 to have a variable end
  ldata = asciiread(lfile,(/10/),"integer"); replace 10 with -1 to have a variable end
  nbdata = asciiread(nbfile,(/10/),"integer"); replace 10 with -1 to have a variable end
  lbdata = asciiread(lbfile,(/10/),"integer"); replace 10 with -1 to have a variable end

  tstart = 1
  tend = 10
  tstep = 1
  solve_file  = new((/4,10/),float)
  solve_file(0,:)  = ndata(:);
  solve_file(1,:)  = ldata(:);
  solve_file(2,:)  = nbdata(:);
  solve_file(3,:)  = lbdata(:);
  time  = new((/10/),float)
  do i = 1, 10;
    time(i-1) = i;
  end do 
  wks = gsn_open_wks("png",PNG); 

  res                        = True

  res@tiMainString           = "Iterations to converge GIS 10km" 
  res@tiXAxisString          = "Time Steps"
  res@tiYAxisString          = "Iterations"
  res@tiXAxisFontHeightF     = 0.02        ; Change the font size.
  res@tiYAxisFontHeightF     = 0.02

  res@xyLineColors           = (/2,3,4/)     ; Set the line colors.
  res@xyLineThicknessF       = 2.0         ; Double the width.

;  res@xyYStyle = "Log"   ; Convert Y axis to logarithmic

;  res@xyLabelMode            = "Custom"    ; Label XY curves.
;  res@xyExplicitLabels       = (/"U","V"/) ; Labels for curves
;  res@xyLineLabelFontHeightF = 0.02        ; Font size and color
;  res@xyLineLabelFontColor   = 2           ; for line labels

  plot = gsn_csm_xy(wks,time,solve_file,res) ; Draw an XY plot with 2 curves.

end
