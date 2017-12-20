load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  read u,v,T,omega snapshots
;  output time-mean zonal-mean u,v,u^2,v^2,t,t^2,w,w^2 
;  plot zonal means with hsave2.ncl
;

begin

tmin=200
tmax=38000



; NETCDF files to process
fnames =  systemfunc( \
;  "ls ~/scratch1/preqx/hs-ne60-*/movies/held_suarez0*1.nc;"  \
  "ls movies/held_suarez0*1.nc;"  \
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
if (isfilevar(fuv_all[0],"DIFFT")) then
  computedifft=True
else
  computedifft = False
end if
if (isfilevar(fuv_all[0],"Q")) then
  computeQ=True
else
  computeQ = False
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
nlon = dimsizes(fuv->lon)
print (" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

ps = dim_avg_Wrap(fuv->ps(0,:,:))
ps=0*ps

u = dim_avg_Wrap(fuv->u(0,:,:,:))
u=0*u
v = dim_avg_Wrap(fuv->v(0,:,:,:))
v=0*v
T = dim_avg_Wrap(fuv->T(0,:,:,:))
T=0*T

u2 = dim_avg_Wrap(fuv->u(0,:,:,:)^2)
u2=0*u2

v2 = dim_avg_Wrap(fuv->v(0,:,:,:)^2)
v2=0*v2

T2 = dim_avg_Wrap(fuv->T(0,:,:,:)^2)
T2=0*T2

copy_VarMeta(u,u2)
u2@long_name=u@long_name + " squared"
u2@units=u@units + " squared"

copy_VarMeta(v,v2)
v2@long_name=v@long_name + " squared"
v2@units=v@units + " squared"

copy_VarMeta(T,T2)
T2@long_name=T@long_name + " squared"
T2@units=T@units + " squared"

if (computeOmega)  then
  omega = dim_avg_Wrap(fuv->omega(0,:,:,:))
  omega = 0*omega
  omega2 = dim_avg_Wrap(fuv->omega(0,:,:,:)^2)
  omega2=0*omega2

  copy_VarMeta(omega,omega2)
  omega2@long_name=omega@long_name + " squared"
  omega2@units=omega@units + " squared"

end if
if (computedifft)  then
  difft = (fuv->DIFFT(0,:,:,:))
  difft = 0*difft
end if
if (computeQ)  then
  Q = dim_avg_Wrap(fuv->Q(0,:,:,:))
  Q = 0*Q
  Q2 = dim_avg_Wrap(fuv->Q(0,:,:,:)^2)
  Q2=0*Q2
  RH = Q

  copy_VarMeta(Q,Q2)
  Q2@long_name=Q@long_name + " squared"
  Q2@units=Q@units + " squared"
  copy_VarMeta(Q,RH)
  RH@long_name="relative humidity"
  RH@units="percent"

  pressure=new( (/nlev,nlat,nlon/), "double")
  p0=1000.0*100
end if

um = new ( (/16,nlat/),"double")
um = 0
um!0 = "wave number"
um!1 = "lat"
um@name = "eddy kinetic energy"

pi = 4*atan(1d0)


; loop over all files
tmin_used=9999d0
tmax_used=0d0
count=0
do fk=0,fcnt-1        
  fuv=fuv_all[fk];    file with u, v
  ntimes = dimsizes(fuv->time)
  do t=0,ntimes-1
    time = fuv->time(t)
    if ( time .gt. tmin .and. time .le. tmax) then
      count=count+1
      if (time .gt. tmax_used) then
        tmax_used=time
      end if
      if (time .lt. tmin_used) then
        tmin_used=time
      end if
      if (mod(count,10).eq.1) then
        print(" processing time = "+time)
      end if

      ps = ps + dim_avg(fuv->ps(t,:,:))

      u = u + dim_avg(fuv->u(t,:,:,:))
      u2 = u2 + dim_avg(fuv->u(t,:,:,:)^2)
      v = v + dim_avg(fuv->v(t,:,:,:))
      v2 = v2 + dim_avg(fuv->v(t,:,:,:)^2)
      T = T + dim_avg(fuv->T(t,:,:,:))
      T2 = T2 + dim_avg(fuv->T(t,:,:,:)^2)
      if (computeOmega)  then
        omega = omega + dim_avg(fuv->omega(t,:,:,:))
        omega2 = omega2 + dim_avg(fuv->omega(t,:,:,:)^2)
      end if
      if (computedifft)  then
        difft = difft + (fuv->DIFFT(t,:,:,:))
      end if
      if (computeQ)  then
        ;Q = Q + dim_avg(fuv->Q(t,:,:,:))
        ;Q2 = Q2 + dim_avg(fuv->Q(t,:,:,:)^2)
       
        do l=0,nlev-1
           pressure(l,:,:) = fuv->hyam(l)*p0 + fuv->hybm(l)*fuv->ps(t,:,:)
        end do
        ; vapor pressure  Q p / (Rd/Rv)
        VP_f = fuv->Q(t,:,:,:)*pressure/( 287.04/461.5 )
        ; saturation vapor pressure (mb)
        SVP_f = exp(53.67957 - 6743.769/fuv->T(t,:,:,:) - 4.8451 * log( fuv->T(t,:,:,:) ) )
        SVP_f = 100*SVP_f

        RH = RH + dim_avg(VP_f/SVP_f)
        ;print("min/max RH "+min(VP_f/SVP_f)+" "+max(VP_f/SVP_f) ) 
      end if


      ; compute first 15 Fourier modes of u:  VERY SLOW. rewrite as matrix multiply
;      do m=1,15
;        print("computing m="+m)
;        do l=0,nlev-1
;          deta = (fuv->hyai(l+1) + fuv->hybi(l+1)) - (fuv->hyai(l) + fuv->hybi(l))
;          do j=0,nlat-1
;            ; integrate against sine and cosine
;            a=0d0
;            b=0d0
;            do i=0,nlon-1
;              theta = m*2*pi*i/nlon
;              a=a+ fuv->u(t,l,j,i)*cos(theta)/nlon
;              b=b+ fuv->u(t,l,j,i)*sin(theta)/nlon
;            end do
;            a=2*a
;            b=2*b
;            um(m,j)=um(m,j) + .5*(a*a + b*b)*deta
;          end do
;        end do
;      end do
        


    end if
  end do
end do
ps=ps/count
u=u/count
u2=u2/count
v=v/count
v2=v2/count
T=T/count
T2=T2/count

delete(ps@time) 
ps@time = tmin_used+":"+tmax_used

delete(u@time) 
u@time = tmin_used+":"+tmax_used
delete(u2@time) 
u2@time = tmin_used+":"+tmax_used

delete(v@time) 
v@time = tmin_used+":"+tmax_used
delete(v2@time) 
v2@time = tmin_used+":"+tmax_used

delete(T@time) 
T@time = tmin_used+":"+tmax_used
delete(T2@time) 
T2@time = tmin_used+":"+tmax_used

if (computeOmega)  then
   omega=omega/count
   omega2=omega2/count
   delete(omega@time) 
   omega@time = tmin_used+":"+tmax_used
   delete(omega2@time) 
   omega2@time = tmin_used+":"+tmax_used
end if
if (computedifft)  then
   difft=difft/count
   delete(difft@time) 
   difft@time = tmin_used+":"+tmax_used
end if
if (computeQ)  then
   Q=Q/count
   Q2=Q2/count
   RH=RH/count

   delete(Q@time) 
   Q@time = tmin_used+":"+tmax_used
   delete(Q2@time) 
   Q2@time = tmin_used+":"+tmax_used
   delete(RH@time) 
   RH@time = tmin_used+":"+tmax_used
end if


name = "hszonal.nc"
system("/bin/rm "+name)           ; remove any pre-existing file
ncdf     = addfile(name ,"c")  ; open output netCDF file

ncdf->ps = ps
ncdf->u = u
ncdf->u2 = u2
ncdf->v = v
ncdf->v2 = v2
ncdf->T = T
ncdf->T2 = T2

if (computeOmega)  then
  ncdf->omega = omega
  ncdf->omega2 = omega2
end if
if (computedifft)  then
  ncdf->difft = difft
end if
if (computeQ)  then
  ;ncdf->Q = Q
  ;ncdf->Q2 = Q2
  ncdf->RH = RH
end if

ncdf->gw = fuv->gw
ncdf->ilev = fuv->ilev
ncdf->hyam = fuv->hyam
ncdf->hyai = fuv->hyai
ncdf->hybm = fuv->hybm
ncdf->hybi = fuv->hybi

ncdf@ne = fuv@ne
ncdf@nv = fuv@nv

count!0 = "ntimes"
ncdf->num_snapshots = count
fnames_list = new ( (/fcnt,800/), "character")
fnames_list!0 = "nfiles"
fnames_list!1 = "MAX_LENGTH"
do fk=0,fcnt-1
   c=sizeof(stringtochar(fnames(fk)));
   fnames_list(fk,0:c-1)=stringtochar(fnames(fk));
end do
   ncdf->fnames = fnames_list

end

