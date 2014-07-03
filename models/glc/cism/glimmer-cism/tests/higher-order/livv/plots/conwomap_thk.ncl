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
;  varbl    = addfile("/tmp/work/4ue/SEACISM/cism_tests/gis_5km/data/gis_5km.ice2sea.51-100.nc","r")
;  stock    = addfile("/tmp/work/4ue/SEACISM/cism_tests/gis_5km/data/gis_5km.ice2sea.init.nc","r")
  thk1  =  stock->thk(1,:,100:380)                
  thk50  = varble->thk(0,:,100:380)                

;=================================================;

;  wks  = gsn_open_wks("png","conwomap_thk.ncl.png")   	  	  ; open a ncgm file
  wks  = gsn_open_wks("png",PNG)   	  	  ; open a ncgm file
  gsn_define_colormap(wks,"wgne15")          ; choose a colormap
  plot = new(3,graphic)

  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn           = False         ; turn off contour lines
  res@gsnFrame            = False           ; don't advance frame
 
  res@vpWidthF            = 0.3           ; change aspect ratio of plot
  res@vpHeightF           = 0.6

  thk1@_FillValue            = 9.96921e+36 ; set missing values for ncl v5.0
  thk50@_FillValue           = 9.96921e+36 ; set missing values for ncl v5.0

  maxv1=max(thk1(:,:));
  print(maxv1);
  printVarSummary(thk1)
  maxv50=max(thk50(:,:));
  print(maxv50);
  printVarSummary(thk50)

   thkdiff = new((/601,281/),float);
   thkdiff = thk50(:,:)-thk1(:,:);

;  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels 
;  res@cnLevels = (/ -1e-3,-7.5e-4,-5e-4,-2.5e-4,-1e-4,1e-4,2.5e-4,5e-4,7.5e-4,1e-3,2e-3,3e-3,4e-3,5e-3,7.5e-3,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.125,0.15,0.2,0.25,0.3,0.35,0.4,0.465/) ; set levels 

  res@tiMainString              = "" ; main title

  plot(0) = gsn_contour(wks,thk1,res)       ; contour the variable
  plot(1) = gsn_contour(wks,thk50,res)       ; contour the variable
  plot(2) = gsn_contour(wks,thkdiff,res)       ; contour the variable

  resP = True                          ; panel resources
;  resP@gsnPanelLaberBar  = True        ; common contour label bar
  resP@gsnMaximize = True              ; make as large as possible
  resP@txString    = "Thickness"  ; (optional) common title 
  gsn_panel(wks,plot,(/1,3/),resP)

  end
