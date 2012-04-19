
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  fname =  "movies/baroclinic1.nc"
  f     = addfile(fname,"r")       ; Open netCDF file.	

  field = (/"T","zeta"/)

  output_format = "pdf"

  lat1d=doubletofloat(f->lat);*90./acos(0.)
  lon1d=doubletofloat(f->lon);*90./acos(0.)

  vcnt = getfiledimsizes(f)
  
  ntimes = vcnt(4)
  nlev = vcnt(1)
  nxy = vcnt(0)

; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1
  timeinc=12


  firstlev=nlev-1
  lastlev=nlev-1
  levinc=2
  var = new((/nlev,nxy/),"float")
  nxe=0
  
 print ("ntimes="+ntimes+" nlev="+nlev+" nxy="+nxy)


  time = f->time
 
  wks = gsn_open_wks(output_format,"baro1b")
  print ("Writing file baro1b."+output_format)

; To see the various colormaps, go to:
;
; http://ngwww.ucar.edu/ncl/coltable.html
;
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

    res@sfXArray            = lon1d        ; Required to tell NCL where to
    res@sfYArray            = lat1d        ; overlay data on globe.

    res@cnLevelSelectionMode = "ManualLevels"

    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = -90
    res@mpShapeMode          = "FreeAspect"

    res@mpMinLatF = 20
    res@mpMaxLatF = 80
    res@vpHeightF           = .3

  do k= firstlev,lastlev,levinc

    do n = 0,1
		var(k,:)=doubletofloat(f->$field(n)$(lasttime,k,:))
 
		if (n.eq.0) then
			res@cnMinLevelValF = 250.
			res@cnMaxLevelValF = 295.
			res@cnLevelSpacingF = 2.5
			var@long_name = "Temperature"
			var@units = "degrees kelvin"
		else
			res@cnMinLevelValF = -7.5e-5
		    res@cnMaxLevelValF = 7.5e-5
		    res@cnLevelSpacingF = 1e-5
			var@long_name = "Relative vorticity"
			var@units = "1/s"
		end if
		print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

        res@tiMainString        = field(n)+" at level = "+k+", time="+time(lasttime)+" days"    ;
        contour = gsn_csm_contour_map(wks,var(k,:),res)     ; Draw a contour plot.
        print("Contour complete: "+res@tiMainString)
    end do
  end do





end

