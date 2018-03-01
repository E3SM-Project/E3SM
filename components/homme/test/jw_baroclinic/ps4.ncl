
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
fnames =  systemfunc("ls -r -S "+\
"jw-ne7-4t360l26-hpnu2.5e16-1/movies/jw_baroclinic1.nc " +\
"jw-ne13-4t200l26-hpnu1.2e16-1/movies/jw_baroclinic1.nc " +\
"jw-ne26-4t100l26-hpnu1e15-1/movies/jw_baroclinic1.nc " +\
"jw-ne52-4t50l26-hpnu1.5e14-1/movies/jw_baroclinic1.nc " +\
"jw-ne104-4t25l26-hpnu1.5e13-1/movies/jw_baroclinic1.nc " +\
"jw-ne120-4t20l26-hpnu1.5e13-1/movies/jw_baroclinic1.nc " )


 f      = addfiles(fnames,"r")       ; Open netCDF files.	

print(fnames)

  fcnt = dimsizes(fnames)
  vcnt=new((/fcnt,5/),"integer")
  do i=0,fcnt-1
    vcnt(i,:) = getfiledimsizes(f[i])
  end do	   	 

lat=0
lon=0
var=0

output_format = "pdf"
wks = gsn_open_wks(output_format,"ps4")


do n=9,9,1
do rf=0,fcnt-1

  delete(lat)  ; delete variables that we will resuse with different sizes
  delete(lon)  ; in this loop
  delete(var)


  ntimes = vcnt(rf,4)
  nlev = vcnt(rf,2)
  nlat = vcnt(rf,1)
  nlon = vcnt(rf,0)



sizevar = dimsizes(f[rf]->ps)
ntime = sizevar(0)
nlat=sizevar(1)
nlon=sizevar(2)

hdeg=0.0
if (nlat .eq. 1024 ) then
;  NE=104
    hdeg=0.3
end if
if (nlat .eq. 512 ) then
;  NE=52
    hdeg=0.6
end if
if (nlat .eq. 384 ) then
;  NE=42
   hdeg=0.7
end if
if (nlat .eq. 256 ) then
;  NE=26
    hdeg=1.2
end if
if (nlat .eq. 192 ) then
;  NE=21
    hdeg=1.4
end if
if (nlat .eq. 128 ) then
;  NE=13
    hdeg=2.3
end if
if (nlat .eq. 64 ) then
;  NE=13
    hdeg=4.3
end if


lat = f[rf]->lat
lon = f[rf]->lon

ntimes = vcnt(rf,4)
nlev = vcnt(rf,2)
nlat = vcnt(rf,1)
nlon = vcnt(rf,0)
print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon = "+nlat+"x"+nlon)

nplev=1
plev=(/nlev-1,10,0/)
var = new((/nlat,nlon/),"float")

time = f[rf]->time




;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; Any time you don't want a plot drawn, just set the resources:
;  res@gsnDraw  = False
;  res@gsnFrame = False
;
  res                     = True

  res@gsnFrame = False

;  res@lbAutoManage = False
  res@pmLabelBarDisplayMode = "NoCreate"
  if (rf .eq. fcnt-1 ) then
     res@gsnFrame = True
;     res@pmLabelBarDisplayMode = "Always" ; Turn on a label bar.
  end if

;  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 20           ; at color index 24.


  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.


; defaults:
;  res@vpXF      = 0.2
;  res@vpYF      = 0.8
;  res@vpWidthF  = 0.6
;  res@vpHeightF = 0.6

  res@vpXF      = 0.2
  res@vpYF      = .95 - .18*rf
  res@vpWidthF  = 0.5
  res@vpHeightF = 0.11


    res@cnLevelSelectionMode = "ManualLevels"

    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpShapeMode          = "FreeAspect"
    res@mpShapeMode          = "FixedAspectNoFitBB"

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
       var = doubletofloat(f[rf]->ps(n,:,:))
       copyatt( var, f[rf]->ps(n,:,:) )

    if ( "Pa" .eq. var@units) then
       var=var/100
       var@units="mb"
    end if

    print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) )

    if ( time(n) .ge. 7  ) then
       res@cnMinLevelValF = 940
       res@cnMaxLevelValF = 1020
       res@cnLevelSpacingF = 10
    else
       res@cnMinLevelValF = 993
       res@cnMaxLevelValF = 1005
       res@cnLevelSpacingF = 1
    end if 
    ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
    ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

    res@tiMainString        = "HOMME "+hdeg+":F34:0"

;    gsn_define_colormap(wks,"rainbow")                ; choose colormap
;    gsn_define_colormap(wks,"wh-bl-gr-ye-re")
;    gsn_define_colormap(wks,"WhViBlGrYeOrRe")
    gsn_define_colormap(wks,"cosam")

    contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;



    print("Contour complete: "+res@tiMainString)
 


; loop over file
end do 
; loop over time level
end do 


end

