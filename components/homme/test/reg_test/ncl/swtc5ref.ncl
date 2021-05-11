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
; edited by Mike Levy, 2011-01-19
;***********************************************************************;
begin
; swtc5
; -----
ref_dir = "~/codes/homme/test/reg_test"
fref = addfile(ref_dir+"/ref_sol/REF0114.cdf","r")
fhomme = addfile("./movies/swtc51.nc","r")
case = 5

if (isfilevar(fhomme,"gw")) then
  print ("using gaussian grid")
  gauss_grid = 1    ; use shagc routines
else
  if (fhomme->lat(0) .lt. -89.99999 ) then
    print ("using equally spaced grid")
    gauss_grid = 2    ; use shaec routines
  else
    print ("assuming shifted equally spaced grid")
    gauss_grid = 3     ; default, fixed-offset grid, need to interpolate using fo2fsh
  end if
end if

nlat = dimsizes(fhomme->lat)
nlon = dimsizes(fhomme->lon)
nlat@double = True
lat  = latGau(nlat, "lat", "latitude", "degrees_north")
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")
lon  = lonGlobeF(nlon, "lon", "longitude", "degrees_east")
print("HOMME & T213 reference solutions compared on grid: "+nlon+" x "+nlat)

; Convert NCAR geop reference solution to lat-lon
tmp = dimsizes(fref->geopotential)
ntimes=tmp(0)
nwave=tmp(1)
mwave=tmp(2)

geohatr = new ( (/nlat,nlat/), double)
geohati = new ( (/nlat,nlat/), double)
geo_topo = new ( (/nlat,nlon/), double)
geo_topo!0="lat"
geo_topo!1="lon"
geo_topo&lat = lat    
geo_topo&lon = lon
geo_ref = geo_topo
g_homme = 9.80616D0

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

geohatr = 0
geohati = 0
geohatr(0:nwave-1,0:mwave-1) = 2*fref->geopotential(15,:,:,0)
geohati(0:nwave-1,0:mwave-1) = 2*fref->geopotential(15,:,:,1)
shsgc(geohatr,geohati,geo_ref)
geo_ref = (geo_ref   + geo_topo)  /g_homme

wks = gsn_open_wks("pdf","./swtc5ref")
res=True
res@cnLevelSpacingF = 100.d 

res@tiMainString        = "NCAR T213 Reference (geop)"
contour = gsn_csm_contour_map(wks,geo_ref,res)     ; Draw a contour plot.  ;

res@tiMainOffsetYF = 0.025
print("geop error")
print("-----------")
geo = (/fhomme->geop(1,0,:,:)/)
geo!0 = "lat"
geo&lat = fhomme->lat
geo!1 = "lon"
geo&lon = fhomme->lon

refmax = max(geo_ref)
diffmax = max(abs(geo-geo_ref))
refl2 = norml2(geo_ref,gwt)
diffl2 = norml2(abs(geo-geo_ref),gwt)
refl1 = norml1(geo_ref,gwt)
diffl1 = norml1(abs(geo-geo_ref),gwt)
  
print("l1,l2,mx = "+diffl1/refl1+" "+diffl2/refl2+" "+diffmax/refmax)
  
res@tiMainString = "HOMME output (geop)"
contour = gsn_csm_contour_map(wks,geo,res)     ; Draw a contour plot

exit ; MNL Note: vorticity code doesn't work (yet?)

; Convert NCAR zeta reference solution to lat-lon
zetahatr = new ( (/nlat,nlat/), double)
zetahati = new ( (/nlat,nlat/), double)
zeta_ref = new ( (/nlat,nlon/), double)
zeta_ref!0="lat"
zeta_ref!1="lon"
zeta_ref&lat = lat    
zeta_ref&lon = lon

zetahatr = 0
zetahati = 0
zetahatr(0:nwave-1,0:mwave-1) = 2*fref->vorticity(15,:,:,0)
zetahati(0:nwave-1,0:mwave-1) = 2*fref->vorticity(15,:,:,1)
shsgc(zetahatr,zetahati,zeta_ref)
print((/min(zeta_ref),max(zeta_ref)/))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5e-5
res@cnMaxLevelValF = 5e-5
res@cnLevelSpacingF = 25d-6/3.
res@cnLineLabelFormat = "!*+^e"

res@tiMainString        = "NCAR T213 Reference (zeta)"
contour = gsn_csm_contour_map(wks,zeta_ref,res)     ; Draw a contour plot.  ;

res@tiMainOffsetYF = 0.025
print("zeta errors")
print("-----------")
zeta = (/fhomme->zeta(1,0,:,:)/)
zeta!0 = "lat"
zeta&lat = fhomme->lat
zeta!1 = "lon"
zeta&lon = fhomme->lon
print((/min(zeta),max(zeta)/))

refmax = max(zeta_ref)
diffmax = max(abs(zeta-zeta_ref))
refl2 = norml2(zeta_ref,gwt)
diffl2 = norml2(abs(zeta-zeta_ref),gwt)
refl1 = norml1(abs(zeta_ref),gwt)
diffl1 = norml1(abs(zeta-zeta_ref),gwt)
  
print("l1,l2,mx = "+diffl1/refl1+" "+diffl2/refl2+" "+diffmax/refmax)
  
res@tiMainString = "HOMME output (zeta)"
contour = gsn_csm_contour_map(wks,zeta,res)     ; Draw a contour plot.  ;

end

