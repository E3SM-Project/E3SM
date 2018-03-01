load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$HOME/contrib/ncl/homme_get_resolution.ncl"
;
;  read u,v files, compute spectra, plot at t=30
;  save data to spectra.nc
;  plot from spectra.nc with ke2.ncl  (much faster)
;
;


begin

fnames =  systemfunc( \
;  "ls movies/jw_baroclinic-000086400-1.nc; " \
  "ls movies/jw_baroclinic1.nc; " \
)
fps_all   = addfiles(fnames,"r")       ; Open netCDF files.	
fcnt1 = dimsizes(fnames)

fnames =  systemfunc( \
;  "ls movies/jw_baroclinic-000086400-4.nc; " \
  "ls movies/jw_baroclinic4.nc; " \
)
fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
print("files = "+fnames)
fcnt = dimsizes(fnames)
print ("number of uv input files = " + fcnt + " ps input files = " + fcnt1 )

output_format = "ncgm"
wks = gsn_open_wks(output_format,"ke")
print ("Writing file ke."+output_format)


fps=fps_all[0];   file with ps
fuv=fuv_all[0];    file with u, v

ntimes = dimsizes(fps->time)
nlev = dimsizes(fps->lev)
nlat = dimsizes(fps->lat)
nlon = dimsizes(fps->lon)
print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)



spc = new ( (/nlat/), "double", 1d-99 )
;spc@_FillValue = 1d-99 ; has to be positive 
spc = spc@_FillValue
spc_c=spc
x=spc
n_avg = 0        



; change the levels and or times to plot here
firsttime=30
lasttime=30
firsttime=5
lasttime=5
timeinc=5
fk=0

do time=firsttime,lasttime,timeinc

  fps=fps_all[fk];
  fuv=fuv_all[fk];
  
  time_ps = fps->time
  time_uv = fuv->time
  nuv = ind(time_uv.eq.time)
  nps = ind(time_ps.eq.time)
  
  print("file = "+fnames(fk)+" record="+nps+" "+nuv  )
  if (isatt(fps,"ne") .and. fps@ne>0) then
    NE=fps@ne
    NPTS = fps@nv
    deg=360.0 / (4*NE*(NPTS-1))
    hdeg=sprintf("%4.2f", deg)
    print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
  else
    deg = get_resolution(systemfunc("pwd"))
    hdeg=sprintf("%7.2f", deg)
  end if 
  
  tmax = doubletoint(360/deg/2)
  tmax = nlat
  if (tmax .gt. nlat) then
    tmax=nlat
  end if
  print("tmax = "+tmax) 
  
  
  ps = fps->ps(nps,:,:)
  copyatt( ps, fps->ps(nps,:,:))
  if ( "Pa" .ne. ps@units) then
    print("converting ps to Pascals")
    ps=ps*100
    ps@units="Pa"
  end if
  
  hyam = fps->hyam
  hybm = fps->hybm
  pnew = (/ 850 /)
  
  interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
  extrap = True  ; is extrapolation desired if data is outside the range of PS
  P0mb = 1000    ; yes, this argument must be in mb
  var = vinth2p(fuv->u(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
  var2 = vinth2p(fuv->v(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
  print ("u min/max = "+min(var(0,:,:))+", "+max(var(0,:,:)) )
  print ("v min/max = "+min(var2(0,:,:))+", "+max(var2(0,:,:)) )
  
  print ("computing vector spherical harmonic transform...t="+time)
                                ;ab = vhaeC(var(0,:,:),var2(0,:,:));  ; for equally spaced grid, including both poles
  ab = vhagC(var(0,:,:),var2(0,:,:));   ; GAuss grid
  print ("done. ")
  
                                ; compute energy in u,v components
  cr = ab(0,:,:)                ; real coef  (nlat,nlat)
  ci = ab(1,:,:)                ; imaginary  (nlat,nlat)
  pwr = (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
  cr = ab(2,:,:)                ; real coef  (nlat,nlat)
  ci = ab(3,:,:)                ; imaginary  (nlat,nlat)
  pwr = pwr + (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
  ab=0
  
                                ; for clarity use do loops
  do n1=0,tmax-1
    spc(n1) = pwr(n1,0)
    do m=1,n1
      spc(n1) = spc(n1) + 2.*pwr(n1,m)
    end do
    spc(n1) = 0.25*spc(n1)
    x(n1) = n1
  end do
  
  ; delete variables that will be reused for next file
  delete(time_ps)
  delete(time_uv)


  lthick = 1.0

  res=True
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8   
  res@tmXTMinorPerMajor = 8   
  res@tmYRMinorPerMajor = 8   
  res@tmYLMinorPerMajor = 8   
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold 
  res@tmXBLabelFontHeightF = .10

  res@tiMainString           = ":F21:DAY "+time
  res@tiXAxisString          = ":F21:spherical wavenumber"
  res@tiYAxisString          = ":F21:Kinetic Energy"
  res@gsnFrame = False
  res@xyLineThicknessF = lthick

  colors = (/"white","black","OrangeRed4","Red3",\
  "red1","Red3","Red4","Blue",\
  "OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
  "red1","red2","red3","red4",\
  "DarkTurquoise","Turquoise","Cyan","blue",\
  "gold","greenyellow",\
  "GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
  /)
  gsn_define_colormap(wks, colors)
  
  res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
  res@xyDashPattern = 0
  
  
  res@trXMinF = 1       
  res@trXMaxF = 500       
  res@trYMinF = 1e-6      
  res@trYMaxF = 1e4   
  
  plot  = gsn_xy (wks,x,spc,res) ; create plot
  
  xref = (/15.,300./)
  spc_ref = 5*(xref/15)^(-3) ; make it 6 @ k=15
  polyres                  = True
  polyres@gsLineThicknessF = lthick
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
  frame(wks)                                     ; now advance
  
  
  do n1=0,nlat-1
    spc_c(n1) = spc(n1)*(n1^3)
  end do
  res@tiYAxisString          = ":F21:Compensated Kinetic Energy"
  res@trYMinF = 1
  res@trYMaxF = 1e6   
  plot  = gsn_csm_xy (wks,x,spc_c,res) ; create plot
  xref = (/15.,300./)
  spc_ref = (xref^3)*5*(xref/15)^(-3) ; make it 6 @ k=15
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
  frame(wks)                                     ; now advance 


end do  ; loop over time level


;save spectrum data
system("/bin/rm spectra.nc")           ; remove any pre-existing file
ncdf     = addfile("spectra.nc" ,"c")  ; open output netCDF file
; make time and UNLIMITED dimension      ; recommended  for most applications
;filedimdef(ncdf,"time",-1,True) 
ncdf->x = x
ncdf->spc = spc





end

