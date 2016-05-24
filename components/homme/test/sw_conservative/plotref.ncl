load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"

;***********************************************************************;
; function  || T || on 2D horizontal array
;***********************************************************************;
function norml2(varz[*][*]:double,gw[*]:double)
local i,output
begin
  s2 =  dimsizes(varz) 
  gs =  dimsizes(gw)

  if ( s2(0) .ne. gs(0) ) then
     print ("norml2: error: first dimension does not match Gauss weight dimension: " + s2(0) + " " + gs(0) )
  end if
  if ( 2*s2(0) .ne. s2(1) ) then
     print ("norml2: error: second dimsion is not 2x first dimension") 
  end if
  
  varl = ( gw # (varz^2) )/sum(gw)
  output = sum(varl)/s2(1)
  output = sqrt(output)
  return(output)
end
function norml1(varz[*][*]:double,gw[*]:double)
local i,output
begin
  s2 =  dimsizes(varz) 
  gs =  dimsizes(gw)

  if ( s2(0) .ne. gs(0) ) then
     print ("norml1: error: first dimension does not match Gauss weight dimension: " + s2(0) + " " + gs(0) )
  end if
  if ( 2*s2(0) .ne. s2(1) ) then
     print ("norml1: error: second dimsion is not 2x first dimension") 
  end if
  
  varl = ( gw # (varz) )/sum(gw)
  output = sum(varl)/s2(1)
  return(output)
end

;***********************************************************************;
; NCL script to plot Jakob et al reference solutions
; Mark Taylor 2010-10-1
;
;***********************************************************************;
begin
; swtc6
;fref = addfile("../reg_test/ref_sol/REF0092.cdf","r")
;fnames =systemfunc("ls /tmp/sweqx/swtc6-NE?-*-LFT.nc ; ls /tmp/sweqx/swtc6-NE??-*-LFT.nc ")
;itime = 14
;case = 6

fref = addfile("../reg_test/ref_sol/REF0114.cdf","r")
;fnames =systemfunc("ls /tmp/sweqx/swtc5-NE?-*-LFT.nc ; ls /tmp/sweqx/swtc5-NE??-*-LFT.nc ")
fnames =systemfunc("ls ~/scratch1/sweqx/swtc5ref/movies/swtc51.nc")
itime = 15
case = 5

fhomme = addfiles(fnames,"r")
nf = dimsizes(fnames)
print("nf="+nf+" fnames = "+fnames)

; T106    160x320
tmp = dimsizes(fref->geopotential)
ntimes=tmp(0)
nwave=tmp(1)
mwave=tmp(2)




if (isfilevar(fhomme[0],"gw")) then
  print ("using gaussian grid")
  gauss_grid = 1    ; use shagc routines
else
  if (fhomme[0]->lat(0) .lt. -89.99999 ) then
    print ("using equally spaced grid")
    gauss_grid = 2    ; use shaec routines
  else
    print ("assuming shifted equally spaced grid")
    gauss_grid = 3     ; default, fixed-offset grid, need to interpolate using fo2fsh
  end if
end if



nlat = dimsizes(fhomme[0]->lat)
nlon = dimsizes(fhomme[0]->lon)
nlat@double = True
lat  = latGau(nlat, "lat", "latitude", "degrees_north")
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")
lon  = lonGlobeF(nlon, "lon", "longitude", "degrees_east")
print("HOMME & T213 reference solutions compared on grid: "+nlon+" x "+nlat)

geohatr = new ( (/nlat,nlat/), double)
geohati = new ( (/nlat,nlat/), double)
geo = new ( (/nlat,nlon/), double)
geo!0="lat"
geo!1="lon"
geo&lat = lat    
geo&lon = lon

geo_ref = geo
geo_topo = geo

g_homme = 9.80616D0


if (case .eq. 5 ) then
  ; need to subtract mountain from plot, but mountain is not in
  ; reference output
  h0 = 5960
  casename = "case5"
  pi = 4*atan(1d0)
  hm = 2000
  bigr=pi/9
  do i=0,nlon-1
    do j=0,nlat-1
      r=sqrt((lat(j)*(pi/180)-pi/6)^2 + (lon(i)*(pi/180) - 3*pi/2)^2)
      if ( r .gt. bigr) then 
         r = bigr
      end if
      geo_topo(j,i) = ( hm * (1 - r/bigr) ) * g_homme
    end do
  end do
  ; contour intervals for plots, as specified by test case
  cinterval= 100
  einterval = 20  

end if

if (case .eq. 6 ) then
  casename = "case6"
  geo_topo = 0
  h0 = 8000 
  ; contour intervals for plots, as specified by test case
  cinterval= 250
  einterval = 10  
end if




geohatr = 0
geohati = 0
geohatr(0:nwave-1,0:mwave-1) = 2*fref->geopotential(itime,:,:,0)
geohati(0:nwave-1,0:mwave-1) = 2*fref->geopotential(itime,:,:,1)
;geohatr(0,0)=geohatr(0,0)*2
;geohati(0,0)=geohati(0,0)*2
if (gauss_grid .eq. 1) then
   ; gauss grid
   shsgc(geohatr,geohati,geo_ref)
end if
if (gauss_grid .eq. 2) then
   ; equally spaced (fixed) grid
   shsec(geohatr,geohati,geo_ref)
   ; check this CODE!
   return
end if
if (gauss_grid .eq. 3) then
   ; we need output on offset grid
   ; interpolate to offset grid
   ; output nlat dimension will be 1 less than geo_ref_fixed
   geo_ref_fixed = new ( (/nlat+1,nlon/), double)
   shsec(geohatr,geohati,geo_ref_fixed)   ; output on equally spaced
   geo_ref = f2fosh(geo_ref_fixed)
   ;check this CODE!
   return
end if


geo_ref = (geo_ref   + geo_topo)  /g_homme

do n=0,nf-1
  time = fhomme[n]->time
  geo = g_homme  * (/fhomme[n]->geop(itime,0,:,:)/)
  geo = geo/g_homme

  refmax = max(geo_ref)
  diffmax = max(abs(geo-geo_ref))
  refl2 = norml2(geo_ref,gwt)
  diffl2 = norml2(abs(geo-geo_ref),gwt)
  refl1 = norml1(geo_ref,gwt)
  diffl1 = norml1(abs(geo-geo_ref),gwt)
  
;  print("ref mass = "+norml1(geo_ref,gwt) )
;  print("    mass = "+norml1(geo,gwt) )
  print("l1,l2,mx = "+diffl1/refl1+" "+diffl2/refl2+" "+diffmax/refmax)
  
end do
print("contouring at time="+itime)


;wks = gsn_open_wks("ncgm",casename)
wks = gsn_open_wks("pdf",casename)

res=True
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMinLevelValF = 1000. 
;res@cnMaxLevelValF = 20000.
res@cnLevelSpacingF = cinterval

res@tiMainString        = "NCAR T213"
contour = gsn_csm_contour_map(wks,geo_ref,res)     ; Draw a contour plot.  ;
res@tiMainString        = "HOMME"
contour = gsn_csm_contour_map(wks,geo,res)     ; Draw a contour plot.  ;

geodiff = geo
geodiff = geo-geo_ref
res@cnMinLevelValF = -200.
res@cnMaxLevelValF = 200.
res@cnLevelSpacingF = einterval
;contour = gsn_csm_contour_map(wks,geodiff,res)     ; Draw a contour plot.  ;


; make a color plot
  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True       ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.
  res@mpGridAndLimbOn  = False
;  res@mpGridSpacingF  = 30.0
  res@cnLevelSelectionMode = "ManualLevels"

  res@mpOutlineOn          = False
  res@mpFillOn             = False
  res@mpCenterLonF         = 0
  res@mpShapeMode          = "FreeAspect"

  res@cnMinLevelValF = 5000. 
  res@cnMaxLevelValF = 6000.
  res@cnLevelSpacingF = cinterval/2

    res@mpProjection = "CylindricalEquidistant"
    res@mpCenterLonF         = 0
contour = gsn_csm_contour_map(wks,geo,res)     ; Draw a contour plot.  ;

end

