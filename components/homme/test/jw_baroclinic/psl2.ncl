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


  athick=3
  lthick=5

  field = "ps"

if (True) then
labelnames = (/"preqx","theta-h","nh"/)
fnames =  systemfunc(\
"ls ~/scratch1/preqx/jw-ne16-nlev26-qsize4/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwthetah-ne16-nlev26-qsize4/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwtheta-ne16-nlev26-qsize?/movies/jw_baroclinic1.nc; " +\
"")
end if
if (False) then
labelnames = (/"preqx","theta-h","theta-h-eul","nh","nh-eul"/)
fnames =  systemfunc(\
"ls ~/scratch1/preqx/jw-ne16-nlev26-qsize4/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwthetah-ne16-nlev26-qsize4/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwthetah-eul-ne16-nlev26-qsize4/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwtheta-ne16-nlev26-qsize?/movies/jw_baroclinic1.nc; " +\
"ls ~/scratch1/preqx/jwtheta-eul-ne16-nlev26-qsize?/movies/jw_baroclinic1.nc; " +\
"")
end if



; USED FOR PDES On SPHERE:
; resolution comparison
if (False) then
labelnames = (/"3.75","1.88","1.00","0.50","0.25"/)
fnames =  systemfunc(\
"ls jw/jw-ne8-4t360l26-hnu2e16-1/movies/jw_baroclinic1.nc; " +\
"ls jw/jw-ne15-4t180l26-hnu1e16-1/movies/jw_baroclinic1.nc; " +\
;"ls jwnew/jw-ne30-4t90l26-hnu3e14-1/movies/jw_baroclinic1.nc; " +\
"ls jwnew/jw-ne30-4t90l26-hnu7e14-1/movies/jw_baroclinic1.nc; " +\
;"ls jwnew/jw-ne60-4t40l26-hnu3e13-1/movies/jw_baroclinic1.nc; " +\
"ls jwnew/jw-ne60-4t40l26-hnu8e13-1/movies/jw_baroclinic1.nc; " \
;"ls jwnew/jw-ne120-4t20l26-hnu5e12-1/movies/jw_baroclinic1.nc; " \
)
end if

; ne24 vs ne21 w/viscosity comparison
if (False) then
fnames =  systemfunc("ls "+\
"jw-ne24-4t100l26-hnu1.5e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne24-4t100l26-hnu2e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne24-4t100l26-hnu3e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne21-4t120l26-hnu1.5e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne21-4t120l26-hnu2e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne21-4t120l26-hnu3e15-1/movies/jw_baroclinic1.nc " \
)
end if

; NPTS = 2,3,4,8 comparison:
if (False) then
labelnames3 = (/"NE78-P2","NE39-P3","NE26-P4","NE11-P8"/)
fnames3 =  systemfunc("ls  "+\
"jw-ne78-2t120l26-hnu2e15-1/movies/jw_baroclinic1.nc " +\
"; ls "+\
"jw-ne39-3t100l26-hnu2e15-1/movies/jw_baroclinic1.nc " +\
"; ls "+\
"jw-ne26-4t100l26-hnu2e15-1/movies/jw_baroclinic1.nc " +\
"; ls "+\
"jw-ne11t50l26-hnu2e15-4/movies/jw_baroclinic1.nc " \
)
end if

 f      = addfiles(fnames,"r")       ; Open netCDF files.	


fcnt = dimsizes(fnames)
print ("number of input files = " + fcnt )
print("fnames="+fnames)

sizevar = dimsizes(f[0]->$field$)
ntime = sizevar(0)
nlat=sizevar(1)
nlon=sizevar(2)




if ( False ) then
  fnames =  systemfunc("ls ~/scratch1/preqx/jw-ne120-4t20l26-hnu1e13-1/movies/jw_baroclinic1.nc ")
  fref    = addfile(fnames,"r")       ; Open netCDF files.	
  PSref2048=fref->ps
  ; interpolate to 512x1024, store in PSref
  if (get_gridtype(fref->lat) .ne. 1 ) then
    print("Error: reference solution must be on a Gauss grid") 
    return
  end if
  ; spherepack interpolation to 
  print("Interpolating reference solution to "+nlat+" x "+nlon)
  PSref = g2gsh_Wrap(PSref2048(15:15,:,:),(/nlat,nlon/),nlat-1)
  print("Interpolation done")

  refname = "0.25"+":F34:0:F21:"
;  title = "HOMME vs HOMME 0.25:F34:0:F21:"
  title = "HOMME"
;  fcnt=fcnt-1
else
  fnames =  systemfunc("ls ~/scratch1/preqx/eaT340.ps.nc")
  ;fnames =  systemfunc("ls ~/scratch1/preqx/T340ref.nc")
  fref    = addfiles(fnames,"r")       ; Open netCDF files.	
  PSref=fref[0]->PS
  refname = "T340"
  title = "HOMME vs EUL SPECTRAL"
end if

plotsol=0
plotref=0
plotdiff=0
plotl2 = 1

firsttime=0
;lasttime=30
lasttime=9
timeinc=1




lat = f[0]->lat
lon = f[0]->lon
if ( 1.eq.1 ) then
   gw = f[0]->gw
   glat = gaus(nlat/2)
   print ("HOMME Gauss vs NCL Gauss latitudes:  max diff="+ max(abs(glat(:,0)- lat)) )
   print ("HOMME Gauss vs NCL Gauss weights:    max diff="+ max(abs(glat(:,1)- gw)) )
   d_eta = f[0]->ilev(1:26) - f[0]->ilev(0:25)
   print ("sum of gauss weights, eta weights: "+sum(gw)+" "+sum(d_eta))
end if


  output_format = "pdf"
;  output_format = "X11"

  vcnt=new((/fcnt,5/),"integer")
  do i=0,fcnt-1
    print ("i=" + i + " f[i]=" + getfiledimsizes(f[i]) )
    vcnt(i,:) = getfiledimsizes(f[i])
  end do	   	 
  
  ntimes = vcnt(0,4)
  nlev = vcnt(0,2)
  nlat = vcnt(0,1)
  nlon = vcnt(0,0)



  var = new((/nlat,nlon/),"double")
  diffl2 = new((/fcnt,lasttime+1/),"double")
  diffl2=diffl2@missing_value
  timesl2 = new((/lasttime+1/),"double")
  timesl2 = timesl2@missing_value


   print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

  wks = gsn_open_wks(output_format,"psl2")
  print ("Writing file psl2."+output_format)


do kf = 0,fcnt-1
  print ("kf="+kf+"   max time="+max(f[kf]->time))
end do

do n=firsttime,lasttime,timeinc
do kf = 0,fcnt-1

  nref = closest_val(n,fref[0]->time)
  nkf  = closest_val(n,f[kf]->time)
  if ( fref[0]->time(nref) .ne. f[kf]->time(nkf) ) 
     print("missing data - skipping time="+n)
     continue
  end if
  time = f[kf]->time

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
  res@mpGridAndLimbOn  = True
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
  res@mpMinLatF = 0
  res@mpMaxLatF = 90
  ;res@vpWidthF            = 1           ; change aspect ratio of plot
  ;res@vpHeightF           = .3

  varmin = 1.e9
  varmax = -1.
  var(:,:)= f[kf]->$field$(nkf,:,:)
  copyatt( var, f[kf]->$field$(nkf,:,:) )

  if ( "Pa" .eq. var@units) then
     var=var/100
     var@units="mb"
  end if
  if ( time(nkf) .ge. 7  ) then
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

  res@tiMainString        = "HOMME "

  ;gsn_define_colormap(wks,"rainbow")                ; choose colormap
  ;gsn_define_colormap(wks,"wh-bl-gr-ye-re")
  ;gsn_define_colormap(wks,"WhViBlGrYeOrRe")
  gsn_define_colormap(wks,"cosam")
  if (plotsol .eq. 1) then
     print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) + "  time="+time(nkf) )
     contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
     print("Contour complete: "+res@tiMainString)
   end if

  ; now plot the CAM contour plot:
  res@tiMainString        = "CAM T340"
  var2= PSref(nref,:,:)
  copyatt( var2, PSref(nref,:,:) )
  if ( "Pa" .eq. var2@units) then
     var2=var2/100
     var2@units="mb"
  end if
  if (plotref .eq. 1) then
    print ("min/max = "+min(var2(:,:))+", "+max(var2(:,:)) + "  time="+time(nkf) )
    contour = gsn_csm_contour_map(wks,var2(:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
  end if

  ; interpolate ref solution 
  nlatref=dimsizes(f[kf]->lat)
  nlonref=dimsizes(f[kf]->lon)
  print("Interpolating ref solution to "+nlatref+" x "+nlonref)
  var3 = g2gsh_Wrap(var2,(/nlatref,nlonref/),nlat-1)





  var=var-var3
  res@tiMainString        = "DIFF"
  res@cnMinLevelValF = -12.0
  res@cnMaxLevelValF = 12.0
  res@cnLevelSpacingF = 2.0
  res@mpGridAndLimbOn  = False
  if (plotdiff .eq. 1) then
    print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) + "  time="+time(nkf) )
    contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
  end if

  diffl2(kf,n) = norml2(var,gw)
  timesl2(n) = time(nkf)
  print ( kf + " time="+time(nkf) + "  || HOMME-REF || = " + diffl2(kf,n)  )
  delete(time)
end do  ; loop over file names
end do  ; loop over time level



if (plotl2 .eq. 1 ) then
  ; read in the uncertainty data
  unc = addfile("./reference_uncertainty.nc","r")
  l2_unc = unc->l2_max
  time_unc = unc->time

  resxy=True
  resxy@tmXBMajorThicknessF = athick
  resxy@tmYLMajorThicknessF = athick
  resxy@tmXBMinorThicknessF = athick 
  resxy@tmYLMinorThicknessF = athick 
  resxy@tmBorderThicknessF = athick

  resxy@tmXBMajorOutwardLengthF = 0
  resxy@tmXBMinorOutwardLengthF = 0
  resxy@tmYLMajorOutwardLengthF = 0
  resxy@tmYLMinorOutwardLengthF = 0

  resxy@gsnFrame = False
  resxy@trXMaxF   = lasttime
  resxy@trXMaxF   = 10
  resxy@trXMinF   = 0
  resxy@trYMaxF   = 1000
  resxy@trYMinF   = .001
  resxy@trYLog   = True
    resxy@tmXBMinorPerMajor = 4   
    resxy@tmXTMinorPerMajor = 4   
    resxy@tmYRMinorPerMajor = 8   
    resxy@tmYLMinorPerMajor = 8   
    resxy@tmYLLabelFont = 21    ;  21 = helvetica
    resxy@tmXBLabelFont = 21    ;  22 = helvetica-bold 

   resxy@tmXBLabelFontHeightF = .025  ; default = .02
   resxy@tmYLLabelFontHeightF = .025  ; default = .02

  colors = (/"white","black","Blue","green",\
"OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
"red","red1","red2","red3","red4",\
"DarkTurquoise","Turquoise","Cyan","gold","greenyellow",\
"khaki1","GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
/)
  gsn_define_colormap(wks, colors)



  resxy@tiMainString           = ":F21:"+title
  resxy@tiXAxisString          = ":F21:days"
  resxy@tiYAxisString          = ":F21:l:B1:2  SURFACE PRESSURE (hPa)"
  ; plot uncertainty as shaded yellow region between 0 and l2_unc
  resxy@gsnYRefLine           = 0.0         

;  resxy@gsnAboveYRefLineColor = "lightgoldenrod" 
  resxy@gsnAboveYRefLineColor = "khaki1"

  plot=gsn_csm_xy(wks,time_unc,l2_unc,resxy)



  ; add my curves:
  polyres                  = True
  polyres@gsLineThicknessF = lthick


  ystart=500.
  do kf = 0,fcnt-1
    if (kf .eq. 4 ) then
       polyres@gsLineColor = "purple"
       ypos = ystart * 2.0^-kf  
    end if
    if (kf .eq. 3 ) then
       polyres@gsLineColor = "black"
       ypos = ystart * 2.0^-kf  
    end if
    if (kf .eq. 2 ) then
       polyres@gsLineColor = "blue"
       ypos = ystart * 2.0^-kf  
    end if
    if (kf .eq. 1 ) then
       polyres@gsLineColor = "green"
       ypos = ystart * 2.0^-kf  
    end if
    if (kf .eq. 0 ) then
       polyres@gsLineColor = "red" 
       ypos = ystart * 2.0^-kf  
    end if
    dlabel = labelnames(kf)
    tmp=stringtocharacter(dlabel)
    if (tmp(1).eq.".") then
       ; for numbers like 0.25, add a degree symbol
       dlabel = dlabel+":F34:0:F21:"
    end if
    delete(tmp)
    gsn_polyline(wks,plot,timesl2(0:lasttime),diffl2(kf,:),polyres) 

    xx = (/2,4/)
    yy = (/ypos,ypos/)
    gsn_polyline(wks,plot,xx,yy,polyres)              ; add polyline
    polytext = True
    polytext@txFontHeightF      = 0.020                 ; change text size
    polytext@txJust             = "CenterLeft"          ; text justification
    gsn_text(wks,plot,":F21:"+dlabel+" - "+refname,5,ypos,polytext); add text

  end do
  frame(wks)                                     ; now advance

end if



end

