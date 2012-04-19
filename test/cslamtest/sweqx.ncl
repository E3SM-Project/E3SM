; Christoph Erath, 15. February 2012
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  ; define problem, outputformat and output variables
  problem = "swtc5"
  output_format = "ps"	
  field = (/"c","c2"/)

  fname =  "movies/"+problem+"1.nc"
  f     = addfile(fname,"r")       ; Open netCDF file.
  vcnt = getfiledimsizes(f)
  ntimes = vcnt(3)
  nlev = vcnt(1)
  
  firsttime=0
  lasttime=ntimes-1
  timeinc=1

  firstlev=0
  lastlev=nlev-1
  levinc=1
  
 print ("ntimes="+ntimes+" nlev="+nlev)
 print ("plotting index = "+lasttime)

  time = f->time
 
  wks = gsn_open_wks(output_format,"sweqx")

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
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  res@cnMinLevelValF = 0.09
  res@cnMaxLevelValF = 2
  res@cnLevelSpacingF = 0.05
  
  do k = firstlev,lastlev,levinc
    do t = firsttime,lasttime,timeinc
      do n = 0,dimsizes(field)-1
        print("field = "+field(n))
        var=(f->$field(n)$(t,k,:))
        print(dimsizes(var))
        print("plotting var. min/max="+min(var)+" "+max(var))

        res@sfXArray = f->phys_lon
        res@sfYArray = f->phys_lat

        res@tiMainString        = field(n)+" at level = "+k+", time="+time(t)+" days"    ;
        contour = gsn_csm_contour_map(wks,var,res)     ; Draw a contour plot.
        print("Contour complete: "+res@tiMainString)
        delete(var)
        delete(res@sfXArray)
        delete(res@sfYArray)
      end do
    end do
  end do
end

