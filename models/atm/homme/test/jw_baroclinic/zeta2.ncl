
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$HOME/contrib/ncl/homme_get_resolution.ncl"




begin
  fnames =  systemfunc("ls movies/jw_baroclinic*.nc")
  fall   = addfiles(fnames,"r")       ; Open netCDF files.	

  fcnt = dimsizes(fnames)

  fps=fall[0];   file with ps
  f=fall[1];     f has zeta
;  fT=fall[2];    file with T
;  fQ=fall[3];    file with Q

if (isatt(f,"ne")) then
   NE=f@ne
   NPTS = f@nv
   deg=360.0 / (4*NE*(NPTS-1))
   hdeg=sprintf("%4.2f", deg)
   print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
else
   deg = get_resolution()
   hdeg=sprintf("%7.2f", deg)
end if 

;  field = "ps"
  field = "zeta"

  output_format = "ncgm"


  vcnt=new((/fcnt,5/),"integer")
  vcnt(0,:) = getfiledimsizes(f)
  
  ntimes = vcnt(0,4)
  nlev = vcnt(0,2)
  nlat = vcnt(0,1)
  nlon = vcnt(0,0)







; change the levels and or times to plot here
  firsttime=9
  lasttime=ntimes-1
  timeinc=1


;  var = new((/nlat,nlon/),"float")

 print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

  fieldis2d = False
  if(field .eq. "ps" .or. field .eq. "prec" .or. field .eq. "accum") then
    fieldis2d = True
  end if
  if(field .eq. "udrag" .or. field .eq. "vdrag" .or. field .eq. "tsflx" .or. field .eq. "qsflx") then
    fieldis2d = True
  end if
  time = f->time

  wks = gsn_open_wks(output_format,"zeta2")
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
    res@mpShapeMode          = "FreeAspect"


;   ZOOM into region in PolvaniMWR2004:
;  res@gsnAddCyclic = False                  ; This is not global data.

  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False

; default
;  res@mpCenterLonF      = 160.            ; choose center lon
;  res@mpCenterLatF      =  50.           ; choose center lat

; for 45 degree problem
  res@mpCenterLonF      = 160.  
  res@mpCenterLatF      =  15.  


res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpGridLatSpacingF = 30.               ; spacing for lat lines
  res@mpGridLonSpacingF = 30.               ; spacing for lon lines


    ;res@vpWidthF            = 1           ; change aspect ratio of plot
    ;res@vpHeightF           = .3


;    copyatt( var, f[0]->$field$(n,1,:,:) )

    ps = fps->ps(n,:,:)
    copyatt( ps, fps->ps(n,:,:))
    if ( "Pa" .ne. ps@units) then
       print("converting ps to Pascals")
       ps=ps*100
       ps@units="Pa"
    end if
  
    hyam = f->hyam
    hybm = f->hybm
    pnew = (/ 850 /)

    ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2 
    ; is extrapolation desired if data is outside the range of PS
    extrap = False
    P0mb = 1000

    var = vinth2p(f->zeta(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
    var=var*1e5
    print ("zeta 1e5 min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )


    if ( time(n) .le. 7  ) then
       res@cnMinLevelValF = -2
       res@cnMaxLevelValF = 5
       res@cnLevelSpacingF = 1
    else
       res@cnMinLevelValF = -5
       res@cnMaxLevelValF = 35
       res@cnLevelSpacingF = 5
    end if 
    ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
    ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

    gsn_define_colormap(wks,"BlWhRe")

;    res@tiMainString        = "level = "+pnew(0)+" time="+time(n)+" days"  
    res@tiMainString        = "HOMME "+hdeg+":F34:0"
    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)

;    var = vinth2p(fQ->Q(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
;    ;var=var*1e5
;    print ("Q min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
;    res@cnLevelSelectionMode = "EqualSpacedLevels"
;    res@tiMainString        = "Q"
;    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

;    var = vinth2p(fT->T(n,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
;    print ("T min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
;    res@cnLevelSelectionMode = "EqualSpacedLevels"
;    res@tiMainString        = "T"
;    contour = gsn_csm_contour_map(wks,var(0,:,:),res)     ; Draw a contour plot.  ;

 
 

; loop over time level
end do 


end

