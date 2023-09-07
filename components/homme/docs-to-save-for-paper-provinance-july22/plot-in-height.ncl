load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;run like this
; ncl plot-in-height.ncl fili=\"nh-dry-new.nc\" tt=45 ztop=1500 shortc=\"dry-new\" fff=\"Th\"

begin

;;; SET fili, tt, shortc, ztop, fff
;;; as file to read, time frame, shortc name for the plot suffix
;;; ztop is that same as in the run, fff as field

print("file is "+fili)
print("short cut  is "+shortc)
print("frame is "+tt)

gg=9.80616
P0=100000.0
;ztop=20000.0
numz=200
;fff="Th"


;-- define file name
  diri  = "./"
;  fili  = "planar_nonhydro_gravity_wave1.nc"
;  fili  = "planar_rising_bubble1.nc"
;  fili  = "hydro.nc"

;-- open file and read variable 
  f     =  addfile(diri+fili, "r")  

;  var   = f->Th(0,:,{0},{-150000:147487})     ;-- first time step, latitude=40N, longitude=0-60E.
  lon_t = f->lon            ;-- longitude=0-60E
  lat_t = f->lat
  lev_t = f->lev                    ;-- currently 17 levels

  hyam =   f->hyam         
  hybm =   f->hybm
  P0mb =   P0*0.01
  ps   =   f->ps  

  lev_p           = (/100, 200, 300, 400, 500, 600,700,800,900,1000/)
  lev_p@units     = "hPa"                        ; required for vinth2p

nlev = dimsizes(lev_t)
nlat = dimsizes(lat_t)
nlon = dimsizes(lon_t)

;print(lev_t)

;tt=19
field=f->$fff$(tt,:,:)
phi=f->geo(tt,:,:)

;  vp = vinth2p (field(:,:), hyam,hybm, lev_p ,ps(tt,:), 2, P0mb, 1, False )


selectedcols=ind(lat_t.eq.0)

;print(selectedcols)

;print(lat_t(selectedcols))
;print(lon_t(selectedcols))

nscols=dimsizes(selectedcols)
nplev=dimsizes(lev_p)

;print("num select cols="+nscols)

vals_z=new((/numz,nscols/),double)
zvals=new((/numz/),double)

;vals=new((/nplev,nscols/),double)
;psdummy=new((/nlat,1/),double)
;psdummy(:,0)=ps(tt,:)
;valsdummy=new((/nlev,nlat,1/),double)
;valsdummy(:,:,0)=field(:,:)
;vp = vinth2p (valsdummy, hyam,hybm, lev_p ,psdummy(:,:), 2, P0mb, 1, False )

;make z array
dz=(ztop - 0.0)/(numz-1)

do ii=0,numz-1
zvals(ii) = (ii)*dz
end do

;print("zvals="+zvals)

do ii=0,nscols-1
;print("ii="+ii+", aa="+aa(ii))
;print("vals(0:nlev-1,ii)"+vals(0:nlev-1,ii))
;print("f->Th(tt,0:nlev-1,aa(ii))"+f->Th(tt,0:nlev-1,aa(ii)))
;vals(0:nlev-1,ii)=f->Th(tt,0:nlev-1,aa(ii))
;sample=field(:,aa(i))
;vals(:,ii)=interpolated(:)

zvalues=phi(::-1,selectedcols(ii))/gg
values=field(::-1,selectedcols(ii))

;remapped_z = ftcurv(zvalues,values,zvals)
remapped_z = linint1(zvalues,values,False,zvals,0)

vals_z(:,ii) = remapped_z

;if ( ii .ge. 27 .and. ii .le. 31 ) then
;print ("ii="+ii)
;print ("original values "+values)
;print ("original geo "+zvalues/9.81)
;print ("remapped"+remapped_z)
;end if 


end do



;-- define workstation
  wks = gsn_open_wks("pdf",shortc+"-height-field-"+fff+"-tt-"+tostring(tt))

;-- set resources
  res                          =  True
  res@tiMainString             = "frame="+tostring(tt)+" "+fili   ;-- title string

;  cmap = read_colormap_file("MPL_GnBu")      ;-- read the color map
;  cmap = read_colormap_file("cmocean_haline")      ;-- read the color map
  cmap = read_colormap_file("MPL_Oranges")      ;-- read the color map
;  cmap = read_colormap_file("WhiteGreen")      ;-- read the color map

  res@cnFillOn                 =  True       ;-- turn on color fill
;  res@cnFillPalette            = cmap(::-1,:)
  res@cnFillPalette            = cmap(::1,:)
  res@cnLineLabelsOn           =  False      ;-- turns off contour line labels
  res@cnInfoLabelOn            =  False      ;-- turns off contour info label
  res@lbOrientation            = "vertical"  ;-- vertical label bar
;  res@tiYAxisString            =  var@long_name+"  [hPa]"
;  res@cnLinesOn                = False 


  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels    = 0.5+(/ 300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320/)   ; set levels


;;;;;;;;;; this will control appearance of color bar
;;;;;;;;;; disable this block first, then wrun witht his block to confirm that LabelStrings are consistent
  res@lbAutoManage          = False             ; we control label bar
  res@pmLabelBarDisplayMode = "Always"          ; turns on label bar
;;;  res@lbOrientation         = "Horizontal"      ; ncl default is vertical
;;;  res@pmLabelBarSide        = "Bottom"          ; default is right
;;; stride does not seem to work
  res@lbLabelStride         = 1                 ; skip every other label
;;;  res@lbLabelPosition  = "Center"             ; label position
  res@lbLabelAlignment = "BoxCenters"         ; label orientation
  res@lbLabelStrings =  (/"300","aa","302","aa","304","aa","306","aa",\
        "308","aaa","310","aa","312","aa","314","aa","316","aa","318","aa","320"/) 



;  res@cnMinLevelValF  = 300.			; set the minimum contour level
;  res@cnMaxLevelValF  = 320.			; set the maximum contour level
;  res@cnLevelSpacingF = 0.5			; set the interval between contours

;  res@cnMaxLevelCount      = 40
;  res@cnLevelSelectionMode = "EqualSpacedLevels" ; equally spaced contour levels

  res@sfXArray                 =  lon_t(selectedcols)      ;-- uses lon_t as plot x-axis
;  res@sfYArray                 =  lev_t/100  ;-- uses lev_t in hPa as plot y-axis
  res@sfYArray                 =  zvals  ;-- uses lev_t in hPa as plot y-axis

;  res@trYReverse               =  True       ;-- reverses y-axis

;-- generate the plot
  plot = gsn_csm_contour(wks,vals_z,res)

end
