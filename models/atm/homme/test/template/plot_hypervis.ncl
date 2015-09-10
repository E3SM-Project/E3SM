load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

fnames = systemfunc("ls movies/jw_baroclinic1.nc")
f      = addfiles(fnames,"r")       ; Open netCDF files.	

res                     = True
print("fnames = "+fnames)

;output_format = "ncgm"
output_format = "pdf"
;output_format = "X11"

    field = "hypervis"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 10.0
    res@cnLevelSpacingF = 10/50.


lat1d=doubletofloat(f[0]->lat);*180.0/3.14
lon1d=doubletofloat(f[0]->lon);*180.0/3.14
print("lat1d min/max = "+min(lat1d)+" "+max(lat1d))
print("lon1d min/max = "+min(lon1d)+" "+max(lon1d))


fcnt = dimsizes(fnames)
if (isfilevar(f[0],"time")) then
   ntimes = dimsizes(f[0]->time)
  time = f[0]->time
else
   ntimes = 1
   time = 0
end if
if (isfilevar(f[0],"lev")) then
   nlev = dimsizes(f[0]->lev)
else
   nlev =1
end if
nxy = dimsizes(f[0]->lat)


; loop over these timelevels
firsttime=0
lasttime=ntimes-1
timeinc=12

; contour plots of these levels:
firstlev=18
lastlev=18
;firstlev=0
;lastlev=nlev-1
levinc=1


var = new((/nlev,nxy/),"double")
  p = new((/nlev,nxy/),"double")
print ("ntimes="+ntimes+" nlev="+nlev+" nxy="+nxy)

 ndim = dimsizes(getfilevardimsizes(f[0],field))
; ndim = 1   assume ncols
; ndim = 2   assume time,ncols
; ndim = 3   assume time,nlev,ncols



wks = gsn_open_wks(output_format,field)
print ("Writing file "+field+"."+output_format)
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)



    res@sfXArray            = lon1d        ; Required to tell NCL where to
    res@sfYArray            = lat1d        ; overlay data on globe.
    res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "portrait"
    res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 24           ; at color index 24.
    res@cnFillOn            = True         ; Turn on contour fill
;res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
    res@cnFillMode          = "RasterFill"   ; Style of fill. You can also
    res@cnLinesOn           = False        ; Turn off contour lines
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
    res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
    res@mpOutlineOn          = True
    res@mpFillOn             = False
    res@mpShapeMode          = "FreeAspect"


    res@mpMinLatF = 0.
    res@mpMaxLatF = 90.
    res@mpMinLonF = -150.
    res@mpMaxLonF = -30.


do n=firsttime,lasttime,timeinc
  if(ndim .eq. 1) then
    varmin = 1.e9
    varmax = -1.
    var(0,:)= (f[0]->$field$(:))
    firstlev=0
    lastlev=0
    print("var min/max " + min(var(0,:)) + " " + max(var(0,:)) ) 
  end if
  if(ndim .eq. 2) then
    varmin = 1.e9
    varmax = -1.
    var(0,:)= (f[0]->$field$(n,:))
    firstlev=0
    lastlev=0
    print("var min/max " + min(var(0,:)) + " " + max(var(0,:)) ) 
  end if
  if(ndim .eq. 3) then
    if isStrSubset(field,"ICLDIWP") then
      var = 0
       do k=0,nlev-1
;          var(k,:) = f[0]->$field$(n,k,:)
;          sum over levels:  
         firstlev=0
         lastlev=0
          var(0,:) = var(0,:) +  f[0]->$field$(n,k,:)
       end do
    else
       do k=firstlev,lastlev,levinc
         var(k,:)=(f[0]->$field$(n,k,:))
       end do
    end if
  end if



  do k= firstlev,lastlev,levinc


    print("plotting min/max " + min(var(k,:)) + " " + max(var(k,:)) ) 
    res@tiMainString        = field+" at level = "+k+" time="+time(n)+" days"    ;
    contour = gsn_csm_contour_map(wks,var(k,:),res)     ; Draw a contour plot.
    print("Contour complete: "+res@tiMainString)
  end do

  
end do
end

