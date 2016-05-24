
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;  mergenc = systemfunc("cd movies; perl $HOMME_HOME/utils/mergenc.pl *.nc")
  fnames =  systemfunc("ls movies/baroclinic1.nc")
;  fnames =  systemfunc("ls combined/*1_.nc")
  f      = addfiles(fnames,"r")       ; Open netCDF files.	

  ; change the name of the field to plot here
  ; current field names are Q, ps, ke, u, v, T, zeta - names are case sensitive
  ; prec, accum

  
;  field = "T"
  field = "zeta"

; Open a file to send output to. You can use "ps", "pdf", "x11",
; or "ncgm".
;
  output_format = "ncgm"
;  output_format = "pdf"


;  lat1d=doubletofloat(f[:]->lat) *180.0/(4*atan(1.0))

;  lon1d=doubletofloat(f[:]->lon) *180.0/(4*atan(1.0)) + 90

  fcnt = dimsizes(fnames)
  vcnt=new((/fcnt,5/),"integer")
  do i=0,fcnt-1
    print ("i=" + i + " f[i]=" + getfiledimsizes(f[i]) )
    vcnt(i,:) = getfiledimsizes(f[i])
  end do	   	 
  
  ntimes = vcnt(0,4)
  nlev = vcnt(0,2)
  nlat = vcnt(0,1)
  nlon = vcnt(0,0)

; change the levels and or times to plot here
  firsttime=ntimes-1
  lasttime=ntimes-1
  timeinc=1


  nplev=1
  plev=(/nlev-1,10,0/)
  

  var = new((/nlat,nlon/),"float")
;  var!0="time"
;  var!0="lev"
;  var!0="lat"
;  var!1="lon"
;  var&lat = lat1d  
;  var&lon = lon1d
;  var&lat@units="degrees_north"
;  var&lon@units="degrees_east"



 print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

  fieldis2d = False
  if(field .eq. "ps" .or. field .eq. "prec" .or. field .eq. "accum") then
    fieldis2d = True
  end if
  if(field .eq. "udrag" .or. field .eq. "vdrag" .or. field .eq. "tsflx" .or. field .eq. "qsflx") then
    fieldis2d = True
  end if
  time = f[0]->time

  wks = gsn_open_wks(output_format,field)
  print ("Writing file "+field+"."+output_format)


  do n=firsttime,lasttime,timeinc



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
;  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.

    res@cnLevelSelectionMode = "ManualLevels"

    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = -90
    res@mpShapeMode          = "FreeAspect"


;   ZOOM into region in PolvaniMWR2004:
;  res@gsnAddCyclic = False                  ; This is not global data.
;    res@mpProjection = "CylindricalEquidistant"
;    res@mpLimitMode = "LatLon"
    res@mpMinLatF = 20
    res@mpMaxLatF = 80
    ;res@vpWidthF            = 1           ; change aspect ratio of plot
    res@vpHeightF           = .3

  do kk= 1,nplev
    k=plev(kk-1)

    if(fieldis2d) then
      varmin = 1.e9
      varmax = -1.
       var(:,:)= doubletofloat(f[:]->$field$(n,:,:))
       firstlev=0
       lastlev=0
    else
       var(:,:)=doubletofloat(f[0]->$field$(n,k,:,:))
       copyatt( var, f[0]->$field$(n,k,:,:) )
    end if

    print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) )

    res@cnMinLevelValF = -7.5e-5
    res@cnMaxLevelValF = 7.5e-5
    res@cnLevelSpacingF = 1e-5
    ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
    ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

    res@tiMainString        = "level = "+k+" time="+time(n)+" days"  
    contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;


    print("Contour complete: "+res@tiMainString)
 
  end do

; loop over time level
end do 


end

