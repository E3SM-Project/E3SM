load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;***********************************************************************;
; function  || T || on 2D horizontal array
;***********************************************************************;
function norml2(varz[*][*]:double,gw[*]:double)
local i,output
begin
  s2 =  dimsizes(varz) 
  gs =  dimsizes(gw)

  if ( s2(0) .ne. gs(0) ) then
     print ("norml2: error: first dimension does not match Gauss weight dimension: " + s2(0) + " " + gs(0) )
  end if
  if ( 2*s2(0) .ne. s2(1) ) then
     print ("norml2: error: second dimsion is not 2x first dimension") 
  end if
  
  varl = ( gw # (varz^2) )/sum(gw)
  output = sum(varl)/s2(1)
  output = sqrt(output)
  return(output)
end


;***********************************************************************;
; detect type of lat/lon grid
;
; gauss_grid =-1    routine failed 
; gauss_grid = 1    Gauss grid (use shaeg routines)
; gauss_grid = 2    Equal angle grid with poles  (use shaec routines)
; gauss_grid = 3    Equal angle grid without poles  (must interpolate first)
;
;***********************************************************************;
function get_gridtype(lat)
local d1,d2
begin
  gauss_grid = -1
  if (lat(0) .lt. -89.99999 ) then
    print ("detected equal angle grid with poles")
    gauss_grid = 2    ; use shaec routines
  else
    d1=lat(1)-lat(0)
    d2=lat(2)-lat(1)
    ;print(d1-d2)
    if (abs(d1-d2) .lt. 1e-8 ) then
      print ("detected equal angle grid without poles")
      gauss_grid = 3    ; need to interpolate using fo2fsh
    else
      print("detected Gauss grid")
      gauss_grid = 1  
    end if
  end if
  return(gauss_grid)
end



begin
  fnames =  systemfunc("ls movies/baro2?-asp_baroclinic?.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)
  f=fall[0];      file with dynamics
  fQ=fall[1];     f has tracers

;  original HOMME names
  PS = f->ps
  Q1 = fQ->Q

; asp converted:  
;  PS = f->PS
;  Q1 = fQ->Q1


  Q2 = fQ->Q2
  Q3 = fQ->Q3

  
; output_format = "ncgm"
   output_format = "pdf"
  output_format@wkOrientation = "landscape"

  ntimes = dimsizes(f->time)
  nlev = dimsizes(f->lev)
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

  print("lat min/max: "+min(f->lat)+" "+max(f->lat))
  print("lon min/max: "+min(f->lon)+" "+max(f->lon))
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


; change the levels and or times to plot here
  firsttime=ntimes-1
  lasttime = ntimes-1
  timeinc=1


  time = f->time
  wks = gsn_open_wks(output_format,"q")
  print ("Writing file zeta"+"."+output_format)

  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  do n=firsttime,lasttime,timeinc

    res                     = True
    res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnSpreadColors     = True         ; Use full colormap, but start
;    res@gsnSpreadColorStart =  35           ; at color index 24.
;    res@gsnSpreadColorEnd   = -35           ; at color index 24.
    
    res@cnFillOn            = True         ; Turn on contour fill
    res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLinesOn           = True        ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@cnLevelSelectionMode = "ManualLevels"
    
    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = 180
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
    res@mpMinLonF = 90
    res@mpMaxLonF = 270
    res@mpMinLatF = 0
    res@mpMaxLatF = 90

    ps = PS(n,:,:)
    copyatt( ps, PS(n,:,:))
    if ( "Pa" .ne. ps@units) then
      print("converting ps to Pascals")
      ps=ps*100
      ps@units="Pa"
    end if
    
    hyam = f->hyam
    hybm = f->hybm
    pnew = (/ 850 /)
    
    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = False
    P0mb = 1000
    
    res@lbLabelBarOn = True

    print ("Q1 min/max = "+min(Q1(n,:,:,:))+", "+max(Q1(n,:,:,:)) )
    var = vinth2p(Q1(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("Q1 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = .0
    res@cnMaxLevelValF = .05
    res@cnLevelSpacingF = .002
    res@tiMainString        = "Q1"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    
    var = vinth2p(Q2(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("Q2 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = .25
    res@cnLevelSpacingF = .025
    res@tiMainString        = "Q2"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

    var = vinth2p(Q3(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("Q3 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = .05
    res@tiMainString        = "Q3"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    
    var(0,:,:)=PS(n,:,:)/100
    print ("PS min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = 940.
    res@cnMaxLevelValF = 1020.
    res@cnLevelSpacingF = 10.0
    res@tiMainString        = "PS"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

    var = 1e5*vinth2p(f->zeta(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    print ("zeta*1e5 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
    res@cnMinLevelValF = -5.
    res@cnMaxLevelValF = 15.
    res@cnLevelSpacingF = 1.
    res@tiMainString        = "zeta*1e5"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;


end do 

;
;  now compute l2 error against T340 reference solution 
;
  fref    = addfile("./T340ref.nc","r")       ; Open netCDF files.	

  if (get_gridtype(f->lat) .ne. 1 ) then
    print("Error: reference solution must be on a Gauss grid") 
    return
  end if

  ; interpolate output to 512x1024
  nlatref=dimsizes(fref->lat)
  nlonref=dimsizes(fref->lon)
  print("Interpolating solution to "+nlatref+" x "+nlonref)
  do n=0,lasttime,lasttime
     PSn = g2gsh_Wrap(PS(n,:,:),(/nlatref,nlonref/),nlat-1)
     error = PSn - fref->PS(n,:,:)
     error = error/100 ; ERROR computed in mb
     print("time="+fref->time(n)+" PS l2 error vs. T340 reference: "+norml2(error,fref->gw)+"mb")
   end do
end

