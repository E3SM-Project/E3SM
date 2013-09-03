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


;***********************************************************************;
; NCL script to plot PS
;***********************************************************************;
begin

fnames =  systemfunc("ls movies/asp_baroclinic1.nc; "+\
"ls rsplit*/movies/asp_baroclinic1.nc; " +\
"ls jwasp-ne30*/movies/asp_baroclinic1.nc; ")
 f      = addfiles(fnames,"r")       ; Open netCDF files.	
print("fnames = "+fnames)
plotsol=1
plotref=1

;
;  now compute l2 error against T340 reference solution 
;
  fref    = addfile("~/codes/homme/test/reg_test/ref_sol/T340ref.nc","r")       ; Open netCDF files.	
  if (ismissing(fref)) then
    plotref=0
  end if


  output_format = "pdf"
;  output_format = "X11"
;  output_format = "ncgm"


if (isatt(f[0],"ne")) then
   NE=f[0]@ne
   NPTS = f[0]@np
   deg=360.0 / (4*NE*(NPTS-1))
   hdeg=sprintf("%4.2f", deg)
   print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
end if 

  psname = "ps"
  if (isfilevar(f[0],"PS")) then
     psname = "PS"
  end if
  ntimes=dimsizes(f[0]->time)


; change the levels and or times to plot here
  firsttime=9
  lasttime=9
  timeinc=1



  if (lasttime .gt. ntimes-1) then
     lasttime=ntimes-1
  end if

sizevar = dimsizes(f[0]->$psname$)
ntime = sizevar(0)
nlat=sizevar(1)
nlon=sizevar(2)



lat = f[0]->lat
lon = f[0]->lon




  fcnt = dimsizes(fnames)
  print ("number of input files = " + fcnt )

  
  ntimes = dimsizes(f[0]->time)
  nlev = dimsizes(f[0]->lev)
  nlat = dimsizes(f[0]->lat)
  nlon = dimsizes(f[0]->lon)



  var = new((/nlat,nlon/),"double")
  diffl2 = new((/fcnt,lasttime+1/),"double")
  diffl2=0

   print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

  time = f[0]->time
  wks = gsn_open_wks(output_format,"ps")
  print ("Writing file ps."+output_format)


do n=firsttime,lasttime,timeinc
do kf = 0,fcnt-1
;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)
; Any time you don't want a plot drawn, just set the resources:
;  res@gsnDraw  = False
;  res@gsnFrame = False
;
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 20           ; at color index 24.


  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  res@mpGridAndLimbOn  = False
  res@mpGridSpacingF  = 10.0
  res@cnLevelSelectionMode = "ManualLevels"

  res@mpOutlineOn          = False
  res@mpFillOn             = False
  res@mpCenterLonF         = 0
  res@mpShapeMode          = "FreeAspect"



;  res@gsnAddCyclic = False                  ; This is not global data.
;    res@mpProjection = "CylindricalEquidistant"
;    res@mpLimitMode = "LatLon"
  res@mpCenterLonF         = 180
  res@mpMinLatF = -90
  res@mpMaxLatF = 90

  ;res@vpWidthF            = 1           ; change aspect ratio of plot
  ;res@vpHeightF           = .3

  varmin = 1.e9
  varmax = -1.
  var(:,:)= f[kf]->$psname$(n,:,:)
  copyatt( var, f[kf]->$psname$(n,:,:) )

  if ( "Pa" .eq. var@units) then
     var=var/100
     var@units="mb"
  end if

  if ( time(n) .ge. 7  ) then
     res@cnMinLevelValF = 940.
     res@cnMaxLevelValF = 1020.
     res@cnLevelSpacingF = 10.0
  else
     res@cnMinLevelValF = 993.
     res@cnMaxLevelValF = 1005.
     res@cnLevelSpacingF = 1.0
  end if 


  ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
  ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

  res@tiMainString        = "HOMME "+hdeg+":F34:0"+":F22:   Day "+time(n)

  ;gsn_define_colormap(wks,"rainbow")                ; choose colormap
  gsn_define_colormap(wks,"cosam")
  if (plotsol .eq. 1) then
     print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) + "  time="+time(n) )
     contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
     print("Contour complete: "+res@tiMainString)
   end if

  if (plotref .eq. 999) then
    ; now plot the CAM contour plot:
    res@tiMainString        = "CAM T340"
    var2= fref->PS(n,:,:)
    copyatt( var2, fref[0]->PS(n,:,:) )
    if ( "Pa" .eq. var2@units) then
       var2=var2/100
       var2@units="mb"
    end if
    print ("min/max = "+min(var2(:,:))+", "+max(var2(:,:)) + "  time="+time(n) )
    contour = gsn_csm_contour_map(wks,var2(:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
  end if

  if (plotref.eq.1) then
    ; compute l2 errors vs. ref solution
    if (get_gridtype(f[kf]->lat) .ne. 1 ) then
      print("Error: output must be on a Gauss grid") 
    else
      nref = ind( fref->time .eq. f[kf]->time(n) )
      if (ismissing(nref) ) then
        print("time="+f[kf]->time(n)+" not found in reference solution") 
      else
                                ; interpolate output to 512x1024
        nlatref=dimsizes(fref->lat)
        nlonref=dimsizes(fref->lon)
        print("Interpolating solution to "+nlatref+" x "+nlonref)
        PSn = g2gsh_Wrap(f[kf]->ps(n,:,:),(/nlatref,nlonref/),nlat-1)
        error = PSn - fref->PS(nref,:,:)
        error = error/100 ; ERROR computed in mb
        print("fname="+fnames(kf))
        print("time="+fref->time(nref)+" PS l2 error vs. T340 reference: "+norml2(error,fref->gw)+" mb")
      end if
    end if
  end if
  

end do  ; loop over file names
end do  ; loop over time level




end

