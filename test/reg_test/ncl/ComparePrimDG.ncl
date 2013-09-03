; Simple script to generate a contour plot for the vorticity
; on the JW baroclinic test. It is assume that a netCDF file
; set to fnames contains a variable named zeta as such
;         double zeta(time, lev, lat, lon) ;
;                zeta:long_name = "Relative vorticity" ;
;                zeta:units = "1/s" ;
; Also it is required the hybrid coordinate variables
;         double hyam(lev) 
;                hyam:long_name = "hybrid A coefficiet at layer midpoints" ;
;                hyam:units = "dimensionless" ;
;         double hybm(lev) ;
;                hybm:long_name = "hybrid B coefficiet at layer midpoints" ;
;                hybm:units = "dimensionless" ;
; ******************************************************************

; MNL Note: I copied test/dg_jwbcl/movies/zeta.ncl and made two small changes
;           1) I plot days 0, 5, 10, 15, 20, 25, and 30 instead of 1 - 30
;           2) I include the day in the title of the plot

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;  fnames =  "jw_bcl1_diff.nc"
  fnames = new(2,"string")
  fnames(0) =  "../primdgx-trunk/baro-dg/movies/jw_bcl1.nc"
  fnames(1) =  "../primdgx-branch/baro-dg/movies/jw_bcl1.nc"
; fnames =  systemfunc("ls jw_baroclinic1.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	
  print ("files = "+fnames)
  fcnt = dimsizes(fnames)
  f=fall[0];      trunk (and file with dynamics)
  f2 = fall[1];   branch

  
; output_format = "ncgm"
   output_format = "pdf"
;  output_format@wkOrientation = "landscape"

  ntimes = dimsizes(f->time)
  nlev = dimsizes(f->lev)
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)
 ;end  
 print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)


; for HOMME dycore.  
  if (isatt(f,"ne")) then
    NE=f@ne
    NPTS = f@np
    deg=360.0 / (4*NE*(NPTS-1))  ; average grid spacing at equator
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    deg=360.0/nlon
  end if

  psname = "ps"
  if (isfilevar(f,"PS")) then
     psname = "PS"
  end if


; change the levels and or times to plot here
  firsttime=0
  lasttime=30
  timeinc=5


  time = f->time
  wks = gsn_open_wks(output_format,"zeta-diff")
  print ("Writing file zeta-diff"+"."+output_format)

  do n=firsttime,lasttime,timeinc
    print ("time = "+n)
    zeta1 = f->zeta(n,:,:,:)
    zeta2 = f2->zeta(n,:,:,:)
    ps1 = f->$psname$(n,:,:)
    ps2 = f2->$psname$(n,:,:)
	ps = ps1 - ps2

    psunits=stringtochar(ps1@units)
    if ( "Pa" .ne. psunits(0:1) ) then
      print("converting ps to Pascals")
      ps=ps*100
      ps@units="Pa"
    end if
    
    hyam1 = f->hyam
    hybm1 = f->hybm
    hyam2 = f2->hyam
    hybm2 = f2->hybm
    pnew = (/ 850 /)
    
    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = True
    P0mb = 1000
	lone = wgt_areaave(abs(ps),f->gw,1.0,0)/wgt_areaave(abs(ps1),f->gw,1.0,0)
	ltwo = sqrt(wgt_areaave(ps*ps,f->gw,1.0,0)/wgt_areaave(ps1*ps1,f->gw,1.0,0))
	linf = max(abs(ps))/max(abs(ps1))
    print ("ps1 min/max = "+min(ps1(:,:))+", "+max(ps1(:,:)) )
    print ("ps2 min/max = "+min(ps2(:,:))+", "+max(ps2(:,:)) )
    print ("diff min/max = "+min(ps(:,:))+", "+max(ps(:,:)) )
	print ("l1 / l2 / linf err = " + lone + ", " + ltwo + ", " + linf)
    var1 = vinth2p(zeta1,hyam1,hybm1,pnew,ps1,interp,P0mb,1,extrap)
    var2 = vinth2p(zeta2,hyam2,hybm2,pnew,ps2,interp,P0mb,1,extrap)
	var = var1 ; grab units
	var = var1 - var2
	lone = wgt_areaave(abs(var(0,:,:)),f->gw,1.0,0)/wgt_areaave(abs(var1(0,:,:)),f->gw,1.0,0)
	ltwo = sqrt(wgt_areaave(var(0,:,:)*var(0,:,:),f->gw,1.0,0)/wgt_areaave(var1(0,:,:)*var1(0,:,:),f->gw,1.0,0))
	linf = max(abs(var(0,:,:)))/max(abs(var1(0,:,:)))
    print ("zeta1 850 min/max = "+min(var1(0,:,:))+", "+max(var1(0,:,:)) )
    print ("zeta2 850 min/max = "+min(var2(0,:,:))+", "+max(var2(0,:,:)) )
    print ("diff 850 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
	print ("l1 / l2 / linf err = " + lone + ", " + ltwo + ", " + linf)
    var=var*1e5

    
    res                     = True
    res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnSpreadColors     = True         ; Use full colormap, but start
    res@gsnSpreadColorStart =  35           ; at color index 24.
    res@gsnSpreadColorEnd   = -35           ; at color index 24.
    
    res@cnFillOn            = True         ; Turn on contour fill
    res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLinesOn           = True        ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@cnLevelSelectionMode = "ManualLevels"
    
    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = 0
    res@mpShapeMode          = "FreeAspect"
    
    res@tmYLLabelFont = 21    ;  21 = helvetica
    res@tmXBLabelFont = 21    ;  22 = helvetica-bold 
    res@lbLabelFont = 21
    res@lbLabelBarOn = True
    
    thick=4
    res@tmXBMajorThicknessF = thick
    res@tmYLMajorThicknessF = thick
    res@tmXBMinorThicknessF = thick 
    res@tmYLMinorThicknessF = thick 
    res@tmBorderThicknessF = thick
;  res@tmXBLabelFontHeightF = .015
;
    res@mpCenterLonF         = 180

;    if ( time(n) .le. 7  ) then
;       res@cnMinLevelValF = -2.
;       res@cnMaxLevelValF = 5.
;       res@cnLevelSpacingF = 1.
;    else
       res@cnMinLevelValF = -30e-5
       res@cnMaxLevelValF = 15e-5
       res@cnLevelSpacingF = 3e-5
;    end if 

    gsn_define_colormap(wks,"BlWhRe")
	day_str = ":F21: (Day " + time(n) + ")"
    res@tiMainString        = ":F21:Diff between trunk and branch "+hdeg+":F34:0"+day_str
    print("drawing contour: "+res@tiMainString)
    contour = gsn_csm_contour_map(wks,var(0,:,:),res) 

    

end do 
end

