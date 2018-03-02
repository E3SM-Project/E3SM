;================================================;
;  conwomap_2.ncl
;
; Concepts illustrated:
;   - Drawing a simple filled contour plot
;   - Selecting a different color map
;   - Changing the size/shape of a contour plot using viewport resources
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; ================================================;
begin
;=================================================;
; open file and read in data
;=================================================;
  fp    = addfile("/tmp/work/4ue/SEACISM/cism_tests/gis_5km/data/gis_5km.ice2sea.1-50.nc","r")
  velnorm1  =  fp->velnorm(1,0,:,100:380)                
  velnorm50  = fp->velnorm(50,0,:,100:380)                

;=================================================;

  wks  = gsn_open_wks("ps","gis")   	  	  ; open a ncgm file
  gsn_define_colormap(wks,"wgne15")          ; choose a colormap
  plot = new(3,graphic)

  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
 
  res@vpWidthF            = 0.3           ; change aspect ratio of plot
  res@vpHeightF           = 0.6

  velnorm1@_FillValue            = 9.96921e+36 ; set missing values for ncl v5.0
  velnorm50@_FillValue           = 9.96921e+36 ; set missing values for ncl v5.0

  maxv1=max(velnorm1(:,:));
  print(maxv1);
  printVarSummary(velnorm1)
  maxv50=max(velnorm50(:,:));
  print(maxv50);
  printVarSummary(velnorm50)

   veldiff = new((/600,281/),float);
   veldiff = velnorm50(:,:)-velnorm1(:,:);

  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels 
;  res@cnLevels = (/ 0.0025,0.005,0.0075,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.125,0.15,0.2,0.25,0.3,0.35,0.4,0.465/) ; set levels 
  res@cnLevels = (/ -1e-3,-7.5e-4,-5e-4,-2.5e-4,-1e-4,1e-4,2.5e-4,5e-4,7.5e-4,1e-3,2e-3,3e-3,4e-3,5e-3,7.5e-3,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.125,0.15,0.2,0.25,0.3,0.35,0.4,0.465/) ; set levels 

  res@tiMainString              = "" ; main title

  plot(0) = gsn_contour(wks,velnorm1,res)       ; contour the variable
  plot(1) = gsn_contour(wks,velnorm50,res)       ; contour the variable
  plot(2) = gsn_contour(wks,veldiff,res)       ; contour the variable

  resP = True                          ; panel resources
;  resP@gsnPanelLaberBar  = True        ; common contour label bar
  resP@gsnMaximize = True              ; make as large as possible
  resP@txString    = "Velocity Norm"  ; (optional) common title 
  gsn_panel(wks,plot,(/1,3/),resP)

  end
