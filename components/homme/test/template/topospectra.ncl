load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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


;
;  NCL script to plot topography: contour plots and spherical harmonic spectra
;
;  spectral from multiple data sets are plotted together to show the 
;  relative smoothing
;
;  Edit fnames below to correct list of filenames
;  The first file in the list should be the unsmoothed reference topography
;
;  Original version: Mark Taylor, 2012
;
begin
setvalues NhlGetWorkspaceObjectId() 
   "wsMaximumSize" : 50000000
end setvalues 

;  note: for T85, use PHIS from history file, since smothing is
;  applied after model starts running.  (T85 input topo file is unsmoothed)
;

  fnames =  systemfunc(\
;+"ls ~/scratch1/topodata/fv0.9x1.25-default-PHIS-new-SGH-SGH30.nc; "\
; "ls ~/scratch1/topodata/USGS-gtopo30_ne120np4_latlon.nc; "\
 "ls ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4_latlon.nc; "\
+"ls ~/scratch1/topodata/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc;"\
+"ls ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4_7xdel2_latlon.nc; "\
+"ls ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4_8xdel2_latlon.nc; "\
+"ls ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4_16xdel2_latlon.nc; "\
;+"ls ~/scratch1/topodata/ne120/USGS-gtopo30_ne120np4_224xdel2_latlon.nc; "\
;+"ls ~/scratch1/topodata/ne120/USGS-gtopo30_ne120np4_384xdel2_latlon.nc; "\
);


;  output_format = "ncgm"
output_format = "pdf"
plot_spectra=False

f      = addfiles(fnames,"r")       ; Open netCDF files.	
print("f="+fnames)


pnew = 0.0  ;  set to > 0 to interpolate to pressure level
res                     = True

field = "PHIS"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -500
res@cnMaxLevelValF = 55000
res@cnLevelSpacingF = 1000
plot_spectra=True
kmax= 512  ; max wave number for spectra plot

; Open a file to send output to. You can use "ps", "pdf", "x11",
; or "ncgm".
;
wks = gsn_open_wks(output_format,field)
print ("Writing file "+field+"."+output_format)


fcnt = dimsizes(fnames)
nlat = dimsizes(f[0]->lat)
nlon = dimsizes(f[0]->lon)

nlat_max = kmax
do fn=0,fcnt-1
   if (dimsizes(f[fn]->lat) .gt. nlat_max ) then
     nlat_max = dimsizes(f[fn]->lat)
   end if
end do
print("nlat_max = "+nlat_max)
spc_all = new ( (/fcnt,nlat_max/), "double", 1d-99 )
spc_all = spc_all@_FillValue
spc_all!0 = "wavenumber"
x = new ( (/nlat_max/), "double", 1d-99 )
x = x@_FillValue

do fn=0,fcnt-1
  f1=f[fn]
  
  nlat = dimsizes(f1->lat)
  nlon = dimsizes(f1->lon)
  
  tmp1 = stringtochar(fnames(fn))
  l1 = 0
  l2 = dimsizes(tmp1)
  l2 = l2-1
  do while isStrSubset(chartostring(tmp1(l1:l2)),"/")
                                ;print("range: "+l1+" "+l2) 
    ii = indStrSubset(chartostring(tmp1(l1:l2)),"/")
                                ;print("ii = "+ii)
    l1 = l1+ii+1
                                ; print("substr "+tmp1(l1:l2))
  end do
                                ; remove .cam2* 
  if isStrSubset(chartostring(tmp1(l1:l2)),".cam2")  then
    iv = indStrSubset(chartostring(tmp1(l1:l2)),".cam2")
    ii = iv(0)
    l2 = l1+ii-1
  end if
  f1name = chartostring(tmp1(l1:l2))
  delete(tmp1)

;
; note: spherepack:  gauss grid or equally-spaced grid.
;  both grids do not include periodic point in longitidue (lon(0) <> lon(nlon-1))
;  equally spaced grid includes points at both poles
;
; CHECK FV PHIS data - what grid is it on?
;
  gauss_grid=get_gridtype(f1->lat)
  
    
  var = new((/nlat,nlon/),"double")
  ndim = dimsizes( getfilevardimsizes(f1,field))
  print ("nlat x nlon "+nlat+"x"+nlon+"  ndim="+ndim)
  
  if (ndim .eq. 1) then
    print("error computing ndim.  variable not in file?")
    print("file = "+fnames(fn))
    exit
  end if
  
;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;     gsn_define_colormap(wks,"BlWhRe")
  
;  gsn_define_colormap(wks,"tbr_var_0-500")
;  res@gsnSpreadColorStart = 60           ; at color index 24.
  
;  gsn_reverse_colormap(wks)
  
  
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "portrait"
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  
  
  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False     ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  
  
  
  res@mpOutlineOn          = True
  res@mpFillOn             = False
  
  varmin = 1.e9
  varmax = -1.
  var(:,:)= (f1->$field$(:,:))
  copyatt( var, f1->$field$(:,:) )

  print("name = "+f1name+" field="+field)
  print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) )
      
  res@tiMainString        = f1name
  contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
  print("Contour complete: "+res@tiMainString)
      
  if (plot_spectra) then
    print("computing spectra...")
    if (gauss_grid .eq. 1) then
      var_hat = shagC (var(:,:))
    else 
      if (gauss_grid .eq. 2) then
        var_hat = shaeC (var(:,:))
      else 
        if (gauss_grid .eq. 3) then
                                ; interpolate from fixed-offset to fixed grid
                                ; areg = fo2fshv(aoff)
                                ; interpolate from fixed-offset to fixed grid
                                ; for vectors: fo2fshv_Wrap(uoff,voff,ureg,vreg)
          var_hat = shaeC (fo2fsh(var(:,:)) )
        else
          print("ERROR: bad gauss_grid="+gauss_grid)
          var_hat = 0
        end if
      end if
    end if
        
    spc = new ( (/nlat/), "double", 1d-99 )
    spc = spc@_FillValue
    spc!0 = "wavenumber"

    tmpdim = dimsizes(var_hat)
    nmax = tmpdim(2)
    pwr = var_hat(0,:,:)*var_hat(0,:,:) 
    do n1=0,nmax-1
      x(n1) = n1
      spc(n1) = (/pwr(n1,0) /)
      do m=1,n1
        spc(n1) = spc(n1) + 2.*(pwr(n1,m))
      end do
    end do
    athick=4
    lthick=2
    
    resxy=True
                                ;res@gsLineThicknessF = 3.0
    resxy@xyXStyle ="Log"
    resxy@xyYStyle ="Log"
    resxy@tmXBMinorPerMajor = 8   
    resxy@tmXTMinorPerMajor = 8   
    resxy@tmYRMinorPerMajor = 8   
    resxy@tmYLMinorPerMajor = 8   
    resxy@tmXMajorGridThicknessF = athick
    resxy@tmYLLabelFont = 21    ;  21 = helvetica
    resxy@tmXBLabelFont = 21    ;  22 = helvetica-bold 
    
    resxy@tmXBMajorThicknessF = athick
    resxy@tmYLMajorThicknessF = athick
    resxy@tmXBMinorThicknessF = athick 
    resxy@tmYLMinorThicknessF = athick 
    resxy@tmBorderThicknessF = athick
    resxy@tiXAxisString          = "~F21~spherical wavenumber"
    resxy@tiYAxisString          = "~F21~Energy"
    resxy@gsnFrame = False
    resxy@xyLineThicknessF = lthick
    resxy@xyDashPattern = 0
    
    resxy@trXMinF = 1       
    resxy@trXMaxF = 1e3
    resxy@trYMinF = 1e-1      
    resxy@trYMaxF = 1e8
    spc_all(fn,0:nlat-1)=spc(0:nlat-1)
    delete(spc)
    delete(var_hat)
    delete(pwr)
  end if

  delete(var)
; loop over files
end do 

if ( plot_spectra) then
  resxy@xyLineColors = (/"black","purple","blue","green","orange","red3","cyan","brown"/)

  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = 1e-1      
  resxy@trYMaxF = 1e8
  plot  = gsn_xy (wks,x(0:kmax-1),spc_all(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance

  if (fcnt.gt.1)  then
  spc_r = spc_all(:,:)
  do n1=0,fcnt-1
    thresh=.50
    abovethresh=True
    do j=0,nlat_max-1
      spc_r(n1,j)=0
      if (.not. ismissing(spc_all(0,j)) .and. spc_all(0,j) .gt. 0) then
        spc_r(n1,j) = spc_all(n1,j)/spc_all(0,j)
        if ( .not. ismissing(spc_r(n1,j))) then
        if (abovethresh .and. ( spc_r(n1,j) .lt. thresh ) ) then
           abovethresh=False
           print("kf="+n1+" smoothing(thresh="+thresh+") for wavenumbers k>="+j)
        end if
        end if
      end if
    end do
    print("kf="+n1+" E("+x(0)+")="+spc_all(n1,0))
  end do
  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = .01
  resxy@trYMaxF = 2e0
  plot  = gsn_xy (wks,x(0:kmax-1),spc_r(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance
  end if 

  do n1=0,nlat_max-1
    spc_all(:,n1) = spc_all(:,n1)*(n1*(n1+1))
  end do
  resxy@trXMinF = 1       
  resxy@trXMaxF = 1e3
  resxy@trYMinF = 1e6      
  resxy@trYMaxF = 1e10
  plot  = gsn_xy (wks,x(0:kmax-1),spc_all(:,0:kmax-1),resxy) ; create plot
  frame(wks)                                     ; now advance

end if


delete(wks)
end


