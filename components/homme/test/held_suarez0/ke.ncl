load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  read u,v files, compute spectra, plot at snapshots
;  save average to spectra.nc
;  plot from spectra.nc with ke2.ncl
;
;  set filename below, then run with:
;   ncl ke.ncl 'plvl=250'
;  (plvl = pressure level in mb.  default is 250)
;


;***********************************************************************;
; function  comput powerspectra of  U=(u1,u2)  V=(v1,v2)
;
;  vor = curl(U)
;  div = div(U)
;        VORTICITY COMPONENT                DIV COMPONENT
;  U = khat cross grad laplace^-1 (vor) +  grad laplace^-1 (div)
;
;
; 
;  Note: laplace^-1 = -1/(k(k+1))
;  Now compute via:  conj(V) dot DIV COMPONENT of U  
;                  = conj(V) dot grad laplace^-1 (div(U))
;                  = conj(V) dot [i k1 , i k2 ] div(U)(k) (-1) / (k(k+1)) 
;                  =  [i k1 conj(V1) + i k2 conj(V2) ]  div(U(k)) (-1)/(k(k+1))
;                  = conj[ -i k1 V1 +  -i k2 V2 ]  div(U(k)) (-1)/(k(k+1))
;                  = conj( div(V(k)) )   div(U(k))/(k(k+1))
;
;  Also:
;  V dot (khat cross grad X) = - vor(V) dot X ???
;                 conj(V) dot (khat cross grad laplace^-1(vor(U)) 
;              =  - conj(vor(V)) vor(U) (-1) / (k(k+1))
;              =  conj(vor(V)) vor(U) / (k(k+1)) 
; 
;  
;  so in the case U = V, we have:
;         VORTICITY COMPONNET                 DIV COMPONET
;  E(k) = |vor(U(k))|^2 / (k(k+1))             |div(U(k))|^2 / (k(k+1)) 
;
; in general:
;  E(k) = conj(vor(U(k))) * vor(V(k)) / (k(k+1)        VOR COMPONENT
;  E(k) = conj(div(U(k))) * div(V(k)) / (k(k+1)        DIV COMPONENT
;
;***********************************************************************;
function kespec_from_vordiv(p1[*][*][*]:double,p2[*][*][*]:double,nlat:integer)
local pwr,n1,x,spc,factor,m
;
; two scalars, p1 and p2.   
; p1 and p2 are the spherical harmonic coefficients of 
; the divergence of vectors U1 and U2
; or the vorticity of vectors U1 and U2
;
; power spectrum:  conj(p1)*p2
;                  (p1(0) - i p1(1) )*(p2(0) + i p2(1) ) 
;                = p1(0)*p2(0) + p1(1)*p2(1)  + Imag. parts.  
; but fields are real, so we know there is a -k wave number with:
;                  (p1(0) + i p1(1) )*(p2(0) - i p2(1) ) 
;                = p1(0)*p2(0) + p1(1)*p2(1)  - Imag. parts.  
;
; the imaginary parts will cancel, and the real parts are equal
; (hence the 2 below, since we only sum over positive wave numbers)
;
begin
spc = new ( (/nlat/), "double", 1d-99 )
spc = spc@_FillValue
spc!0 = "wavenumber"
x=spc

pwr = p1(0,:,:)*p2(0,:,:) + p1(1,:,:)*p2(1,:,:)
do n1=0,nlat-1
  x(n1) = n1
  spc(n1) = (/pwr(n1,0) /)
  factor = 4. * n1*(n1+1)
  do m=1,n1
    spc(n1) = spc(n1) + 2.*(pwr(n1,m))
  end do
  if (n1 .ne. 0) then
    spc(n1) = (6371229.^2) / factor * spc(n1)  ; scaling from NCAR/TN-388+STR
  end if
end do
; normalize so it agrees with Vector Spherical Harmonics code
spc = spc/2

return(spc)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


begin
plot_each_time = True  ; set to True to plot spectra from every napshot
                        ; in addition to time mean 

;dirname = "~/scratch1/preqx/hs-ne240-4t10l26-hnu1.6e12-1/movies/"
;dirname = "~/scratch1/preqx/hs-ne240-4t10l26-hnu5e11-1/movies/"
dirname = "~/scratch1/preqx/hs-ne240-4t10l26-hnu5e12-2/movies/"

if (isvar("fname")) then  
   ; command line specifed a single filename:
   cmd1  = "ls "+dirname+fname+"-1.nc;" 
   cmd2  = "ls "+dirname+fname+"-2.nc;" 
else
   fnamelist = "held_suarez0-*"
   cmd1  = "ls "+dirname+fnamelist+"-1.nc;" 
   cmd2  = "ls "+dirname+fnamelist+"-2.nc;" 
   fname = "held_suarez0-nhu1.6e12"
end if


; NETCDF file containing U,V
fnames =  systemfunc( cmd1 )
fuv_all   = addfiles(fnames,"r")       ; Open netCDF files.	
fcnt1 = dimsizes(fnames)
delete(fnames)

; NETCDF file containing FU,FV 
fnames =  systemfunc( cmd2 )
ff_all   = addfiles(fnames,"r")       ; Open netCDF files.	
print("files = "+fnames)
fcnt = dimsizes(fnames)
print ("number of FU input files = " + fcnt + " U input files = " + fcnt1 )


; check command line argument for pnew:
if (isvar("plvl")) then         ; pressure level on command line?
  pnew = plvl                   ; command line 
else
  ;pnew=250   ; CAM original
  ;pnew=500
  pnew=750
end if


output_format = "ncgm"
nclname = fname + "-" + pnew
wks = gsn_open_wks(output_format,nclname)
print ("computing at p="+pnew+" writing file ke."+output_format)


fuv=fuv_all[0];    file with u, v

ntimes = dimsizes(fuv->time)
nlev = dimsizes(fuv->lev)
nlat = dimsizes(fuv->lat)
nlon = dimsizes(fuv->lon)
print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)


;
; for output to spectra.nc
;
fnames_list = new ( (/fcnt,800/), "character")
fnames_list!0 = "run"
spc = new ( (/nlat/), "double", 1d-99 )
spc = spc@_FillValue
spc!0 = "wavenumber"
spc_vr = spc
spc_dv = spc

x=spc
do n1=0,nlat-1
  x(n1) = n1
end do

spc_avg = new ( (/nlat/), "double", 1d-99 )
spc_avg!0 = "wavenumber"
spc_avg=0
spc_avg_vr=spc_avg
spc_avg_dv=spc_avg

spc_avg_f = spc_avg
spc_avg_diff = spc_avg
spc_avg_conv = spc_avg

x_avg=spc_avg
n_avg = new ( (/fcnt/), "integer")
n_avg!0 = "run"
n_avg = 0        



; change the levels and or times to plot here
firsttime=240
lasttime=24600
timeinc=1

; loop over all times
do time=firsttime,lasttime,timeinc
  tfound=False
  ; for each time, loop over files looking for the one with this time
  do fk=0,fcnt-1        
    fuv=fuv_all[fk];
    nuv = ind(fuv->time.eq.time)

    ff =ff_all[fk];
    nff = ind(ff->time.eq.time)
    
    if ( ismissing(nuv) .or. ismissing(nff) ) then
    else
      tfound=True
      print("found time="+time+" file = "+fnames(fk)+" record="+nff+" "+nuv  )
      
      if (isatt(fuv,"ne")) then
        NE=fuv@ne
        NPTS = fuv@nv
        deg=360.0 / (4*NE*(NPTS-1))
        hdeg=sprintf("%4.2f", deg)
        ;print("NE="+NE+" NPTS="+NPTS+" degree="+hdeg)
      else
        ; lat/lon data
        deg = 360.0 / nlon
        NE = 1
      end if
      
      
        ps = fuv->ps(nuv,:,:)
        copyatt( ps, fuv->ps(nuv,:,:))
        if ( "Pa" .ne. ps@units) then
          print("converting ps to Pascals")
          ps=ps*100
          ps@units="Pa"
        end if
        
        hyam = fuv->hyam
        hybm = fuv->hybm
        
        interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
        extrap = True  ; is extrapolation desired if data is outside the range of PS
        P0mb = 1000    ; ps in Pa, but this argument and pnew must be in mb
        
        temp=vinth2p(fuv->u(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
	var=temp(0,:,:)
        temp=vinth2p(fuv->v(nuv,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
	var2=temp(0,:,:)

      ;print ("u min/max = "+min(var(:,:))+", "+max(var(:,:)) )
      ;print ("v min/max = "+min(var2(:,:))+", "+max(var2(:,:)) )
      print("computing UU spectra")
      vr = uv2vrG_Wrap (var,var2)
      u_vr = shagC (vr(:,:))
      vr = uv2dvG_Wrap (var,var2)
      u_dv = shagC (vr(:,:))

      spc_vr = kespec_from_vordiv(u_vr,u_vr,nlat)
      spc_dv = kespec_from_vordiv(u_dv,u_dv,nlat)
      spc = spc_vr + spc_dv
      print ("spc min/max = "+min(spc)+" "+max(spc))


      ; now spectra of F dot U
      print("computing U F spectra")
      temp=vinth2p(ff->FU(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var=temp(0,:,:)
      temp=vinth2p(ff->FV(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var2=temp(0,:,:)
      ;print ("fu min/max = "+min(var(:,:))+", "+max(var(:,:)) )
      ;print ("fv min/max = "+min(var2(:,:))+", "+max(var2(:,:)) )

      vr = uv2vrG_Wrap (var,var2)
      ff_vr = shagC (vr(:,:))
      vr = uv2dvG_Wrap (var,var2)
      ff_dv = shagC (vr(:,:))

      spc_f_vr = kespec_from_vordiv(u_vr,ff_vr,nlat)
      spc_f_dv = kespec_from_vordiv(u_dv,ff_dv,nlat)
      spc_f = spc_f_vr + spc_f_dv
      print ("spc_f min/max = "+min(spc_f)+" "+max(spc_f))

      ; now spectra of DIFF dot U
      print("computing U DIFF spectra")
      temp=vinth2p(ff->DIFFU(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var=temp(0,:,:)
      temp=vinth2p(ff->DIFFV(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var2=temp(0,:,:)

      vr = uv2vrG_Wrap (var,var2)
      ff_vr = shagC (vr(:,:))
      vr = uv2dvG_Wrap (var,var2)
      ff_dv = shagC (vr(:,:))

      spc_f_vr = kespec_from_vordiv(u_vr,ff_vr,nlat)
      spc_f_dv = kespec_from_vordiv(u_dv,ff_dv,nlat)
      spc_diff = spc_f_vr + spc_f_dv
      ;print ("spc_diff min/max = "+min(spc_diff)+" "+max(spc_diff))


      ; now spectra of CONV dot U
      print("computing CONV F spectra")
      temp=vinth2p(ff->CONVU(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var=temp(0,:,:)
      temp=vinth2p(ff->CONVV(nff,:,:,:),hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
      var2=temp(0,:,:)

      vr = uv2vrG_Wrap (var,var2)
      ff_vr = shagC (vr(:,:))
      vr = uv2dvG_Wrap (var,var2)
      ff_dv = shagC (vr(:,:))

      spc_f_vr = kespec_from_vordiv(u_vr,ff_vr,nlat)
      spc_f_dv = kespec_from_vordiv(u_dv,ff_dv,nlat)
      spc_conv = spc_f_vr + spc_f_dv
      ;print ("spc_conv min/max = "+min(spc_conv)+" "+max(spc_conv))



      spc_avg_f(:)= (n_avg(fk)*spc_avg_f(:) + spc_f(:) )/(n_avg(fk)+1)
      spc_avg_diff(:)= (n_avg(fk)*spc_avg_diff(:) + spc_diff(:) )/(n_avg(fk)+1)
      spc_avg_conv(:)= (n_avg(fk)*spc_avg_conv(:) + spc_conv(:) )/(n_avg(fk)+1)

      spc_avg(:)= (n_avg(fk)*spc_avg(:) + spc(:) )/(n_avg(fk)+1)
      spc_avg_vr(:)= (n_avg(fk)*spc_avg_vr(:) + spc_vr(:) )/(n_avg(fk)+1)
      spc_avg_dv(:)= (n_avg(fk)*spc_avg_dv(:) + spc_dv(:) )/(n_avg(fk)+1)
      x_avg(:) = x(:)        
      n_avg(fk) = n_avg(fk) + 1


      break    
    end if  ; compute specta, plot

  end do ; loop over files, looking for time

  res=True
                                ;res@gsLineThicknessF = 3.0
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8   
  res@tmXTMinorPerMajor = 8   
  res@tmYRMinorPerMajor = 8   
  res@tmYLMinorPerMajor = 8   
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold 
  
  res@tiMainString           = ":F21:DAY "+time
  res@tiXAxisString          = ":F21:spherical wavenumber"
  res@tiYAxisString          = ":F21:Kinetic Energy"
  res@gsnFrame = False
  res@xyLineThicknessF = 1.0
  
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

  polyres                  = True
  polyres@gsLineThicknessF = 1.0
  polyres@gsLineColor = "black"
  
  if ( plot_each_time .eq. True  .and. tfound ) then
    print("Plotting time="+time)
    
    plot  = gsn_xy (wks,x,spc,res) ; create plot
    
    xref = (/10.,300./)
    spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
    gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
    frame(wks)                                     ; now advance
    
    do n1=0,nlat-1
      spc(n1) = spc(n1)*(n1^(5./3.))
    end do
    res@tiYAxisString          = ":F21:Compensated Kinetic Energy"
    res@trYMinF = 1
    res@trYMaxF = 1e6   
    plot  = gsn_csm_xy (wks,x,spc,res) ; create plot
    xref = (/10.,300./)
    spc_ref = (xref^3)*100*(xref/10)^(-5./3.) ; make it 6 @ k=15
    gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
    frame(wks)                                     ; now advance 

    resL=True
    resL@xyXStyle ="Log"
    resL@trXMinF = 1       
    resL@trXMaxF = 300

    resL@trYMinF = -2.
    resL@trYMaxF =  2.
    ;print("Plotting spc_f time="+time)
    ;plot  = gsn_csm_xy (wks,x,1e5*spc_f,resL) ; create plot
    ;print("Plotting spc_diff time="+time)
    ;plot  = gsn_csm_xy (wks,x,1e7*spc_diff,resL) ; create plot
    ;print("Plotting spc_conv time="+time)
    ;plot  = gsn_csm_xy (wks,x,1e4*spc_conv,resL) ; create plot

  end if
  
end do  ; loop over time level
if (n_avg(0) .eq. 0 ) then
  print("no spectra computed")
  exit 
end if

res@trYMinF = 1e-6      
res@trYMaxF = 1e4   
res@tiMainString     = "Climatology"
res@xyLineColors = (/"blue","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg_dv,res) ; create plot
res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
plot  = gsn_csm_xy (wks,x_avg,spc_avg_vr,res) ; create plot

xref = (/10.,300./)
spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
frame(wks)



;save spectrum data
   fk=0
   name = fname+"-"+pnew+".nc"
   system("/bin/rm "+name)           ; remove any pre-existing file
   ncdf     = addfile(name ,"c")  ; open output netCDF file
   ncdf->x_avg = x_avg(:)
   ncdf->n_avg = n_avg(:)

   ncdf->spc_avg_conv = spc_avg_conv(:)
   ncdf->spc_avg_diff = spc_avg_diff(:)
   ncdf->spc_avg_f = spc_avg_f(:)

   ncdf->spc_avg = spc_avg(:)
   ncdf->spc_avg_vr = spc_avg_vr(:)
   ncdf->spc_avg_dv = spc_avg_dv(:)
   ncdf->ne = NE
   ncdf->deg = deg

do fk=0,fcnt-1
   c=sizeof(stringtochar(fnames(fk)));
   fnames_list(fk,0:c-1)=stringtochar(fnames(fk));
end do
   ncdf->fnames = fnames_list



end

