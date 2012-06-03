
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  fname =  "movies/asp_baroclinic2.nc"
  f     = addfile(fname,"r")       ; Open netCDF file.	
  fname1 =  "movies/asp_baroclinic1.nc"
  f1     = addfile(fname1,"r")       ; Open netCDF file.	

  ; field = (/"Q", "C4", "Q2", "C2", "Q3", "C3", "Q4", "C"/)
  ; field = (/"C4"/)
  field = (/"C4", "C2", "C3", "C"/)

  output_format = "ps"


  vcnt = getfiledimsizes(f)
  
  ntimes = vcnt(4)
  nlev = vcnt(1)
  nxy = vcnt(0)

; change the levels and or times to plot here
	firsttime=0
  lasttime=ntimes-1


  firstlev=nlev-1
  lastlev=nlev-1
  levinc=2
  nxe=0
  
 print ("ntimes="+ntimes+" nlev="+nlev+" nxy="+nxy)
 print ("plotting index = "+lasttime)

  time = f->time
 
  wks = gsn_open_wks(output_format,"baro")
  print ("Writing file baro."+output_format)

; To see the various colormaps, go to:
;
; http://ngwww.ucar.edu/ncl/coltable.html
;
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.



    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = -90
    res@mpShapeMode          = "FreeAspect"

    res@mpMinLatF = 20
    res@mpMaxLatF = 80
    res@vpHeightF           = .3

	do t=firsttime, lasttime
    do n = 0,dimsizes(field)-1
      print("field = "+field(n))
      var3d=(f->$field(n)$(t,:,:))
      print(dimsizes(var3d))

      if (f->$field(n)$!2 .eq. "nphys" ) then
        res@sfXArray = f->phys_lon
        res@sfYArray = f->phys_lat
      else
        res@sfXArray = f->lon
        res@sfYArray = f->lat
      end if
      
      
      ps = f1->ps(t,:)
      hyam = f->hyam
      hybm = f->hybm
      pnew = (/ 850 /)
                                ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
      interp = 2 
                                ; is extrapolation desired if data is outside the range of PS
      extrap = False
      P0mb = 1000
;    var = vinth2p(var3d,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      k=17
      var = var3d(k,:)
      
      
      if (n.eq.0) then
;    res@cnLevelSelectionMode = "ManualLevels"
;			res@cnMinLevelValF = 250.
;			res@cnMaxLevelValF = 295.
;			res@cnLevelSpacingF = 2.5
;			var@long_name = "Temperature"
;			var@units = "degrees kelvin"
      end if
;      print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 
      
      res@tiMainString        = field(n)+" at level = "+k+", time="+time(t)+" days"    ;
      contour = gsn_csm_contour_map(wks,var,res)     ; Draw a contour plot.
      print("Contour complete: "+res@tiMainString)
      delete(var3d)
      delete(var)
      delete(res@sfXArray)
      delete(res@sfYArray)
    end do
	end do


end

