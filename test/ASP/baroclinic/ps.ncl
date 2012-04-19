load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"




;***********************************************************************;
; function  compute NE, NPTS from pwd
;***********************************************************************;
function get_resolution()
begin
  wdir =  systemfunc("pwd")
  print( "wdir = "+wdir)
  wname=stringtochar(wdir)
  ne=0
  npts=8
  ; find jw-ne
  do i=0,dimsizes(wname)-7
     if ( wname(i:i+4 ) .eq. "jw-ne" ) then  
         ;print ("i="+i+" char="+wname(i:i+4))
         istart=i+5
         ; find value of NE
         do j=istart+1,dimsizes(wname)-2
            if ( wname(j:j).eq."-" .or. wname(j:j).eq."t") then 
               ;print ("found a -"+wname(j:j) + "  string: "+wname(istart:j-1))
               ne=stringtointeger(chartostring(wname(istart:j-1)))
            
               ; find value of npts
               if ( wname(j:j).eq."-" ) then
                  do k=j+1,dimsizes(wname)-2
                     if ( wname(k:k).eq."t") then 
                        ;print ("found a t: "+wname(j:k) )
                        npts=stringtointeger(chartostring(wname(j+1:k-1)))
                        break
                     end if            
                  end do
               end if
            end if
            if (ne .ne. 0) then
               break
            end if
         end do
     end if 
     if (ne .ne. 0) then
        break
     end if
  end do
  if (ne .eq. 0) then
     deg=0
  else 
     deg=360.0 / (4*ne*(npts-1))
  end if 
  sdeg=sprintf("%7.2f", deg)
  print("Guessing resolution from filename: NE="+ne+" NPTS="+npts+" degree="+sdeg)
  return(sdeg)
end

;***********************************************************************;
; NCL script to plot PS
;***********************************************************************;
begin

fnames =  systemfunc("ls movies/asp_baroclinic1.nc")
 f      = addfiles(fnames,"r")       ; Open netCDF files.	

print("fnames = "+fnames)
plotsol=1
plotref=0

  output_format = "pdf"
;  output_format = "X11"
;  output_format = "ncgm"


if (isatt(f[0],"ne")) then
   NE=f[0]@ne
   NPTS = f[0]@nv
   deg=360.0 / (4*NE*(NPTS-1))
   hdeg=sprintf("%4.2f", deg)
   print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
else
   hdeg = get_resolution()
end if 

  psname = "ps"
  if (isfilevar(f[0],"PS")) then
     psname = "PS"
  end if
  ntimes=dimsizes(f[0]->time)


; change the levels and or times to plot here
  firsttime=9
  lasttime=9
  timeinc=1



  if (lasttime .gt. ntimes-1) then
     lasttime=ntimes-1
  end if

sizevar = dimsizes(f[0]->$psname$)
ntime = sizevar(0)
nlat=sizevar(1)
nlon=sizevar(2)



lat = f[0]->lat
lon = f[0]->lon




  fcnt = dimsizes(fnames)
  print ("number of input files = " + fcnt )
  vcnt=new((/fcnt,5/),"integer")
  do i=0,fcnt-1
    print ("i=" + i + " f[i]=" + getfiledimsizes(f[i]) )
    vcnt(i,:) = getfiledimsizes(f[i])
  end do	   	 
  
  ntimes = vcnt(0,4)
  nlev = vcnt(0,2)
  nlat = vcnt(0,1)
  nlon = vcnt(0,0)



  var = new((/nlat,nlon/),"double")
  diffl2 = new((/fcnt,lasttime+1/),"double")
  diffl2=0

   print ("ntimes="+ntimes+" nlat x nlon "+nlat+"x"+nlon)

  time = f[0]->time
  wks = gsn_open_wks(output_format,"ps")
  print ("Writing file ps."+output_format)


do n=firsttime,lasttime,timeinc
do kf = 0,fcnt-1
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
;  res@gsnSpreadColorStart = 20           ; at color index 24.


  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  res@mpGridAndLimbOn  = False
  res@mpGridSpacingF  = 10.0
  res@cnLevelSelectionMode = "ManualLevels"

  res@mpOutlineOn          = False
  res@mpFillOn             = False
  res@mpCenterLonF         = 0
  res@mpShapeMode          = "FreeAspect"



;  res@gsnAddCyclic = False                  ; This is not global data.
;    res@mpProjection = "CylindricalEquidistant"
;    res@mpLimitMode = "LatLon"
  res@mpCenterLonF         = 180
  res@mpMinLatF = -90
  res@mpMaxLatF = 90

  ;res@vpWidthF            = 1           ; change aspect ratio of plot
  ;res@vpHeightF           = .3

  varmin = 1.e9
  varmax = -1.
  var(:,:)= f[kf]->$psname$(n,:,:)
  copyatt( var, f[kf]->$psname$(n,:,:) )

  if ( "Pa" .eq. var@units) then
     var=var/100
     var@units="mb"
  end if

  if ( time(n) .ge. 7  ) then
     res@cnMinLevelValF = 940.
     res@cnMaxLevelValF = 1020.
     res@cnLevelSpacingF = 10.0
  else
     res@cnMinLevelValF = 993.
     res@cnMaxLevelValF = 1005.
     res@cnLevelSpacingF = 1.0
  end if 


  ;res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/16.
  ;print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 

  res@tiMainString        = "HOMME "+hdeg+":F34:0"+":F22:   Day "+time(n)

  ;gsn_define_colormap(wks,"rainbow")                ; choose colormap
  gsn_define_colormap(wks,"cosam")
  if (plotsol .eq. 1) then
     print ("min/max = "+min(var(:,:))+", "+max(var(:,:)) + "  time="+time(n) )
     contour = gsn_csm_contour_map(wks,var(:,:),res)     ; Draw a contour plot.  ;
     print("Contour complete: "+res@tiMainString)
   end if

  if (plotref .eq. 1) then
    ; now plot the CAM contour plot:
    res@tiMainString        = "CAM T340"
    var2= fref[0]->PS(n,:,:)
    copyatt( var2, fref[0]->PS(n,:,:) )
    if ( "Pa" .eq. var2@units) then
       var2=var2/100
       var2@units="mb"
    end if
    print ("min/max = "+min(var2(:,:))+", "+max(var2(:,:)) + "  time="+time(n) )
    contour = gsn_csm_contour_map(wks,var2(:,:),res)     ; Draw a contour plot.  ;
    print("Contour complete: "+res@tiMainString)
  end if


end do  ; loop over file names
end do  ; loop over time level




end

