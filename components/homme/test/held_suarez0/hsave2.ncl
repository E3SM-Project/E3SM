load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  read and plot zonal-means computec by hsave.ncl
;

begin

tmin=360
tmax=38000


; NETCDF files to process
fnames =  systemfunc( \
  "ls hszonal0.nc;" \
)
fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
print("files = "+fnames)
fcnt = dimsizes(fnames)
print ("number of uv input files = " + fcnt )


if (isfilevar(fuv_all[0],"Omega")) then
  computeOmega=True
else
  computeOmega = False
end if

; check command line argument for pnew:
if (isvar("plvl")) then         ; pressure level on command line?
  pnew = plvl                   ; command line 
else
  pnew=250
end if


fuv=fuv_all[0];    file with u, v
nlev = dimsizes(fuv->lev)
nlat = dimsizes(fuv->lat)
print (" nlev="+nlev+" nlat="+nlat)

p0=1000*100;
;pressure = fuv->hyam*p0 + fuv->hybm*fuv->ps
eta = fuv->hyam+fuv->hybm



wks = gsn_open_wks("pdf","hszonal")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)

  res                     = True
;  res@gsnMaximize         = True         ; Maximize plot size [ps,eps,pdf]
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart =  2          ; subset of predefined  
;  res@gsnSpreadColorEnd   = 150          ; color table
    
  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True         ; contour lines [default=True]
  res@cnLineLabelsOn      = False        ; contour line labels [True]
  res@tmYLLabelFont       = 21    ;  21 = helvetica
  res@tmXBLabelFont       = 21    ;  22 = helvetica-bold 
  res@lbLabelFont         = 21
  res@lbLabelBarOn        = True
    
  thick             = 4
  res@tmXBMajorThicknessF = thick
  res@tmYLMajorThicknessF = thick
  res@tmXBMinorThicknessF = thick 
  res@tmYLMinorThicknessF = thick 
  res@tmBorderThicknessF  = thick
 ;res@tmXBLabelFontHeightF= .015
;
; res@lbLabelBarOn = True
;  res@trYReverse   = True
  res@cnLevelSelectionMode = "ManualLevels"
;  res@lbLabelAutoStride    = True
  res@trYReverse   = True
  res@sfYArray = eta
  res@gsnYAxisIrregular2Linear = True    ; plot with linear p axis

  res@vpWidthF            = .8           ; change aspect ratio of plot
  res@vpHeightF           = .42
;  res@tiMainString           = ":F21:U"
  res@tiYAxisString          = ":F21:eta"
;  res@tiXAxisString          = ":F21:Kinetic Energy"

   field=fuv->u
   field@long_name = "Zonal-mean Zonal Wind"
   res@cnMinLevelValF       = -16.
   res@cnMaxLevelValF       = 32.
   res@cnLevelSpacingF      = 4.
   print("U min/max "+min(field)+" "+max(field))
   c=gsn_csm_contour(wks,field,res)     ; Draw a contour plot.  ;

   field=fuv->T
   field@long_name = "Zonal-mean Temperature"
   res@cnMinLevelValF       = 180
   res@cnMaxLevelValF       = 340
   res@cnLevelSpacingF      = 10
   c=gsn_csm_contour(wks,field,res)     ; Draw a contour plot.  ;

   ; < (T - <T> )^2 > = < T^2 > - 2 < T <T> > + < <T>^2 > = 
   ;                   < T^2 > - < T >^2  
   field = fuv->T2 - fuv->T^2
   field@long_name = "Temperature Eddy Variance"
   field@units = "K^2"
   res@cnMinLevelValF       = 0.
   res@cnMaxLevelValF       = 50.
   res@cnLevelSpacingF      = 5.
   c=gsn_csm_contour(wks,field,res)     ; Draw a contour plot.  ;


   field=fuv->RH
   print("RH min/max = "+min(field)+" "+max(field))
   field@long_name = "Relative Humidity"
   field@units = ""
   res@cnMinLevelValF       = -.1
   res@cnMaxLevelValF       = 1.1  
   res@cnLevelSpacingF      = .10
   c=gsn_csm_contour(wks,field,res)     ; Draw a contour plot.  ;



end

