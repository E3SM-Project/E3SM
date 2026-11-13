;*************************************************
; zonal_2.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;************************************************

;***********************************************************************;
; function  remove zonal mean
;***********************************************************************;
function remove_zonal_mean(var[*][*][*]:double)
begin
  varm=var         
  varz =dim_avg(var)
  sizevar = dimsizes(var)
  do i=0,sizevar(2)-1  ; lon
      varm(:,:,i)=varm(:,:,i)-varz
  end do     
  return(varm)
end

;***********************************************************************;
; function  || T || on zonally averaged (2D) array
;***********************************************************************;
function norml2_2d(varz[*][*]:double,gw[*]:double,etaw[*]:double)
local i,output
begin
  varl = ( (varz^2) # gw)/sum(gw)
  output = (varl # etaw)/sum(etaw)
  output = sqrt(output)
  return(output)
end
;***********************************************************************;
; function  || T || on 3D array
;***********************************************************************;
function norml2(var[*][*][*]:double,gw[*]:double,etaw[*]:double)
local i,output
begin
  varz =dim_avg(var^2)
  varl = (varz # gw)/sum(gw)
  ;print(varl)
  output = (varl # etaw)/sum(etaw)
  output = sqrt(output)
  ;dims = dimsizes(var)
  ;print(dims)
  return(output)
end



begin
;************************************************
; variable and file handling
;************************************************
wks  = gsn_open_wks("ncgm","zonal")                ; Opens a ps file 


in  = addfile("movies/jw_baroclinic1.nc","r")
time = in->time
sizevar = dimsizes(in->T)

ntime = sizevar(0)
nlev = sizevar(1)
nlat=sizevar(2)
nlon=sizevar(3)
lat = in->lat
lon = in->lon


glat = gaus(nlat/2)
print ("HOMME Gauss vs NCL Gauss latitudes:  max diff="+ max(abs(glat(:,0)- lat)) )
gw = glat(:,1)
d_eta = in->ilev(1:26) - in->ilev(0:25)
print ("sum of gauss weights, eta weights: "+sum(gw)+" "+sum(d_eta))


  print ("time = "+time(0) )
  ts  = in->T(0,:,:,:)
  nt = norml2(ts(:,:,:),gw,d_eta)
  print("  |T|| = "+nt)

;************************************************
; plotting
;************************************************

  res            = True                            ; plot mods desired
  res@trYReverse = True                            ; reverse Y axis
;  res@trXReverse = True                            ; reverse X axis
  res@cnLevelSpacingF = 5.0
  res@gsnYAxisIrregular2Linear = True  


; use wrapper function so that meta data is preserved
  zave = dim_avg_Wrap(ts)                            ; calculate zonal ave
  plot=gsn_csm_contour(wks,zave(:,:),res)             ; plot zonal ave


  ts  = in->zeta(0,:,:,:)
  zave = dim_avg_Wrap(ts)                            ; calculate zonal ave
  zave = zave*1e5
  res@cnLevelSpacingF = 0.25
  res@gsnYAxisIrregular2Linear = True  
  plot=gsn_csm_contour(wks,zave(:,:),res)             ; plot zonal ave


  ts  = in->u(0,:,:,:)
  zave = dim_avg_Wrap(ts)                            ; calculate zonal ave
  res@cnLevelSpacingF = 3.0
  res@gsnYAxisIrregular2Linear = True  
  plot=gsn_csm_contour(wks,zave(:,:),res)             ; plot zonal ave


; surface pressure
  phiz = dim_avg_Wrap(in->geos(:,:))                     
  plot=gsn_xy(wks,lat,phiz,False)


e1=time*0
e2=time*0

do t=0,ntime-1
  ts  = in->u(t,:,:,:)
  zave = dim_avg_Wrap(ts)
  if ( t .eq. 0 ) 
     zave0=zave
  end if

  ; compute l2 norm (3d) of u - uzonal
  um = remove_zonal_mean(ts(:,:,:))
  e1(t) = norml2(um,gw,d_eta)

  ; compute l2 norm (2d) of uzonal(t) - uzonal(0) 
  uz = zave(:,:) - zave0(:,:) 
  e2(t) = norml2_2d(uz,gw,d_eta)

  print ("time = "+time(t)+" "+e1(t)+" "+e2(t) )
end do


e1@long_name = ":F18:R:F25: u(t)-u:B:zonal:N:(t) :F18:R"
plot=gsn_xy(wks,time,e1,False)

e2@long_name = ":F18:R:F25: u(t)-u(0) :F18:R"
plot=gsn_xy(wks,time,e2,False)

end
