
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
  fname =  "movies/baro1a-baroclinic1.nc"
  f      = addfile(fname,"r")       ; Open netCDF file.

  field = (/"T","zeta"/)

  output_format = "pdf"


  vcnt = getfiledimsizes(f)

  ntimes = dimsizes(f->time)
  nlev = dimsizes(f->lev)
  nlat = dimsizes(f->lat)
  nlon = dimsizes(f->lon)

; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1

  nplev=1
  plev=(/nlev-1,10,0/)

  var = new((/nlat,nlon/),"float")
  
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

  time = f->time

  wks = gsn_open_wks(output_format,"baro1a")
  print ("Writing file baro1a."+output_format)


  do n=0,1
; To see the various colormaps, go to:
; http://ngwww.ucar.edu/ncl/coltable.html
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start

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
    res@mpMinLatF = 20
    res@mpMaxLatF = 80
    res@vpHeightF           = .3

  do kk= 1,nplev
    k=plev(kk-1)
    var(:,:)=doubletofloat(f->$field(n)$(lasttime,k,:,:))
    copyatt( var, f->$field(n)$(lasttime,k,:,:) )

    print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) )

	if (n.eq.0) then
		; Levels for T
	    res@cnMinLevelValF = 250.
		res@cnMaxLevelValF = 295.
	    res@cnLevelSpacingF = 2.5
	else
		; Levels for zeta
		res@cnMinLevelValF = -7.5e-5
	    res@cnMaxLevelValF = 7.5e-5
		res@cnLevelSpacingF = 1e-5
	end if

    res@tiMainString        = field(n)+" at level = "+k+", time="+time(lasttime)+" days"  
    contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;


    print("Contour complete: "+res@tiMainString)
 
  end do

; loop over time level
end do 


end

