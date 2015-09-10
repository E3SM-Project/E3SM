module crcp_mod

  implicit none
  real, private :: taunud
  integer, private :: nz_wr
  integer, private :: lord, nplt
  real :: zab, t_abs, epp
  integer :: itp
#ifdef TESTMODE
  real :: dtout, dttap
#endif
contains
  subroutine crcp_input_parameters
    use grid_init_mod, only : ex1
    use gcrk_mod, only :  itmn, iprc
    use rad_mod, only : hmx, hmn
    use noise_mod, only : nz_noise, nfr_noise,thn_a, qvn_a
    use prof_init_mod, only : kaver, eps_wind


#ifdef TESTMODE
#define DTTERMS ,dtout,dttap
#else
#define DTTERMS
#endif

    namelist /crcp_nl/ ex1, taunud, nz_wr, nz_noise, nfr_noise, &
         thn_a, qvn_a, zab, t_abs, epp, itp, lord, nplt, itmn, iprc, &
         hmx, hmn, kaver, eps_wind DTTERMS


#ifdef TESTMODE
#undef DTTERMS    
    dtout=360. ! in min
    dttap=360. ! in min
#endif

    ex1 = 4.1/2.9
    taunud = 3600.
    nz_wr = 1

    zab = 17000.     ! height at which absorber begins
    t_abs=600.      ! time scale

!cc calculate pressure gradient force:
    epp=1.e-6
    itp=100

    lord = 3
    nplt = 100

    itmn=5
!cc iprc 0-no preconditioning, 1-with precon
    iprc=1


!cc INITIAL noise amplitude, top level; will be modified in BLM
!cc frequency is prescribed for entire simulation...
    thn_a     = 0.01
    qvn_a     = 0.02e-3
    nz_noise  = 3  

!   jpe not a valid value, reset later in subroutine crcp
    nfr_noise = -1

! from rad_mod
    hmx = 15.e3
    hmn = 12.e3

! from prof_init_mod
    kaver = 5
    eps_wind = 1.e-5


    open(unit=7,file="crcp.nl",status="OLD")
    read(unit=7,nml=crcp_nl)
    close(7)

#ifdef TESTMODE
    write(*,crcp_nl)
#endif

  end subroutine crcp_input_parameters
  subroutine crcp_init(dt,dx,dz,nx, nz)
    use surfflux_mod, only : qcmx,rhmx,th,qv,uu,vv
    use gcrk_mod, only : r,qr,ar
    use grid_init_mod,  only : grid_init
    use mpdat_mod, only : v1,v2,f1,f2,cp,cn,mx,mn
    use moist_init_mod, only : moist_init
    
    
    integer, intent(in):: nx, nz
    real, intent(in) :: dt,dx,dz
    ! from surfflux
    allocate(qcmx(nz),rhmx(nz),th(nz),qv(nz),uu(nz),vv(nz))
    ! from gcrk_1
    allocate(r(nx*nz), qr(nx*nz), ar(nx*nz))
    ! from mpdat_2d
    allocate(v1(nx+1, nz),v2(nx, nz+1),f1(nx+1, nz), f2(nx, nz+1),   & 
       cp(nx,nz),cn(nx, nz),mx(nx, nz),mn(nx, nz))


    call grid_init(dt,dx,dz,nx,nz)

!cc initialize moisture parameters:

    call moist_init


  end subroutine crcp_init

  subroutine crcp(press,temp,zin,vap,uu,vv,st,dt,dx,dz,nzls, dux, duy,   &
       dth, dqv, dqc, dqr, sst, ssqv, t_strt, t_end, istep, nx, nz)
    use moist_init_mod, only : hlatv, hlats, tt0, ee0, gg, rv, rg, cp
    use grid_init_mod,  only : xx, zz, height, time, gac
    use mpdat_mod,      only : mpdat_2d
    use noise_mod,      only : noise,nz_noise, nfr_noise,thn_a, qvn_a
    use surfflux_mod,   only : surfflux, qcmx, rhmx
    use absor_mod,      only : absor
    use velprd_mod,     only : velprd_1
    use rain_mod,       only : rain_fall
    use thermo_mod,     only : thermo
    use rad_mod,        only : radcool
    use zero_mean_mod,  only : zero_mean
    use gcrk_mod,       only : gcrk_1, prforc_1
    use integxz_mod,    only : integxz, integxz_noise
#ifdef TESTMODE
     use diag_mod,       only : diagno_1, diagno_2
#endif
     use prof_init_mod, only  : prof_init, tm_e , rho0, th0, th_e, &
          ux_e, uy_e, qv_e

    implicit none
    integer, intent(in) :: nx, nz, nzls

    integer, intent(inout) :: istep
    real, intent(in) :: t_strt, t_end, sst, ssqv, dt, dx, dz, st
    real, intent(in) :: press(nzls),temp(nzls),zin(nzls),vap(nzls),&
         uu(nzls), vv(nzls)
    real, intent(out) :: dux(nzls),duy(nzls),dth(nzls),dqv(nzls),&
         dqc(nzls),dqr(nzls)

    real :: ft(nx, nz),fx(nx, nz),fy(nx, nz),fz(nx, nz),     &
         fqv(nx, nz), fqc(nx, nz), fqr(nx, nz), p(nx, nz) 
    
    REAL :: theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz) 
    ! current time level 
    REAL :: ux(nx,nz),uy(nx,nz),uz(nx,nz) 
    ! previous time level
    REAL :: uxp(nx,nz),uyp(nx,nz),uzp (nx,nz)
    !c  advective velocities:                                               
    REAL ::  uxa(nx+1,nz),uza(nx,nz+1) 
    integer :: i,k
#ifdef TESTMODE
    !c latent heating:                                                      
    REAL ::  heat_av (nz)
    !c surface fluxes:                                                      
    REAL ::  slat_av (nx)
    REAL ::  ssen_av (nx)
    !c surface precip                                                       
    REAL ::  sprec_av (nx)
    !c radiative cooling                                                    
    REAL ::  radc_av (nz)
#endif
    real :: scr1(nx, nz),scr2(nx, nz) 
    !c latent heating:                                                      
    REAL ::  heat(nz) 
    !c surface fluxes:                                                      
    REAL ::  slat(nx) 
    REAL ::  ssen(nx) 
    REAL ::  radc(nz) 
    REAL ::  thn(nx,nz), qvn(nx,nz) 
    real :: tau(nz)
    real :: sprecip(nx)
    !c required variables:                                                  
    real :: anor_inv, a, b, c, d, e, rh, aver, aveg, sums, suml, sump
    real :: qvsw,  thetme
    real :: pre, delt, tt, sum, sum1, towi, epsb, esw
    integer :: itime, liner, ntime, nxz
!    real, parameter :: taunud=3600.
    real :: den(nx,nz)

!BGL  additional variables for optimization
      real recth0, recdt, xscale, expd(nx)


      call crcp_init(dt,dx,dz,nx,nz)

      call prof_init(press,temp,zin,vap,uu,vv,st,nz,nzls)


      nxz = nx*nz


!cc frequency of noise application (every 10 min below):
      if(nfr_noise<=0) nfr_noise = 10*nint(60./dt)

!cc absorber (tau=1/time scale)

      towi=1./t_abs
      do k=1,nz
      tau(k)=towi*max(0.,height(k)-zab)/(height(nz)-zab)
      enddo


    do k=1,nz
       do i=1,nx
          theta(i,k)=th_e(k) 
          qv(i,k)=qv_e(k)   
          qc(i,k)=0.
          qr(i,k)=0.

          ux(i,k)=ux_e(k)
          uy(i,k)=uy_e(k)
          uz(i,k)=0.
          uxp(i,k)=ux_e(k)
          uzp(i,k)=0.

          ft(i,k)=0.
          fx(i,k)=0.
          fy(i,k)=0.
          fz(i,k)=0.
          
          fqv(i,k)=0.
          fqc(i,k)=0.
          fqr(i,k)=0.

          p(i,k)=0.
       enddo
    enddo

      do k=1,nz
         do i=1,nx
            den(i,k)=rho0(k)*gac(k)
         enddo
      enddo


#ifdef TESTMODE
!cc plot initial fields:
       call diagno_1(ux,uz,theta,scr1,scr2,den, nx, nz, time)
       call diagno_2(qv,qc,qr,scr1,scr2,den, nx, nz)
!      call plot_1(ux,uy,uz,theta)
!      call plot_2(theta,qv,qc,qr)


       do k=1,nz
       heat_av(k)=0.
       radc_av(k)=0.
       enddo
       do i=1,nx
       slat_av(i)=0.
       ssen_av(i)=0.
       sprec_av(i)=0.
       enddo

!cc save initial data:
       write(17) time,nx,nz
       write(17) ux,uz,theta,qv,qc,qr
       write(17) heat_av,radc_av,slat_av,ssen_av,sprec_av
#endif

!cCCC MARCH FORWARD IN TIME:
              ntime=nint ( (t_end-t_strt) * 60. / dt) 
              do itime=1,ntime   ! TIME LOOP
!               print*,'*** itime, time: ',itime,time

!cc extrapolate in time to get advective momentums:
                 call velprd_1(ux,uxp,uxa,uz,uzp,uza,den, nx, nz)
          
!cc save previous velocities:
       do i=1,nxz
        uxp(i,1)=ux(i,1)
        uzp(i,1)=uz(i,1)
       enddo


!c surface flux:
       call surfflux(theta,qv,ux,uy,ft,fqv,fx,fy,ssen,slat,sst,ssqv,nx,nz)
#ifdef TESTMODE
       do i=1,nx
       slat_av(i)=slat_av(i)+slat(i)
       ssen_av(i)=ssen_av(i)+ssen(i)
       enddo
#endif
!cc
      sum=0.
      do i=1,nx-1
      sum=sum+slat(i)/real(nx-1)
      enddo

!cc radiative cooling:
       call radcool(ft,radc, nx, nz)
#ifdef TESTMODE
       do k=1,nz
       radc_av(k)=radc_av(k)+radc(k)
       enddo
#endif

!cc add half of the force:
       do i=1,nxz
        theta(i,1)=theta(i,1)+.5*dt*ft(i,1)
        ux(i,1)   =   ux(i,1)+.5*dt*fx(i,1)
        uy(i,1)   =   uy(i,1)+.5*dt*fy(i,1)
        uz(i,1)   =   uz(i,1)+.5*dt*fz(i,1)
        qv(i,1)   =   qv(i,1)+.5*dt*fqv(i,1)
        qc(i,1)   =   qc(i,1)+.5*dt*fqc(i,1)
        qr(i,1)   =   qr(i,1)+.5*dt*fqr(i,1)
       enddo

!ccc maintain wind:
!       nz_wr=1
!c       nz_wr=12  ! above 2 km only
!       taunud=3600.   ! 60 min time scale
       xscale = 1.0/real(nx-1)

       do k=nz_wr,nz
        aver=0.
        aveg=0.
        do i=1,nx-1
          aver=aver+ux(i,k)*xscale
          aveg=aveg+uy(i,k)*xscale
        enddo
        do i=1,nx
          ux(i,k)=ux(i,k) -(aver-ux_e(k))*dt/taunud
          uy(i,k)=uy(i,k) -(aveg-uy_e(k))*dt/taunud
!          ux(i,k)=ux(i,k) -(aver-uxls(k))*dt/taunud
!          uy(i,k)=uy(i,k) -(aveg-uyls(k))*dt/taunud
        enddo
       enddo

!cC ADVECTION:
!c liner: 1-iord=1, 0-iord prescribed inside mpdata
        liner=0
        if(itime/6*6.eq.itime) liner=1

!cc advect velocities:
        call mpdat_2d(uxa,uza,   ux,den,1,liner,sprecip, nx, nz)
        call mpdat_2d(uxa,uza,   uy,den,2,liner,sprecip, nx, nz)
        call mpdat_2d(uxa,uza,   uz,den,2,liner,sprecip, nx, nz)
!cc advect thermodynamic variables:
        call mpdat_2d(uxa,uza,theta,den,3,liner,sprecip, nx, nz)
        call mpdat_2d(uxa,uza,   qv,den,4,liner,sprecip, nx, nz)
        call mpdat_2d(uxa,uza,   qc,den,5,liner,sprecip, nx, nz)
!cc
!cc
          call rain_fall(qr,tm_e,rho0,uza, nx, nz)

          call mpdat_2d(uxa,uza,   qr,den,6,liner, sprecip, nx, nz)
!cc
          do i=1,nx-1
          sprecip(i)=sprecip(i)*dz/dt
          enddo
          sprecip(nx)=sprecip(1)

#ifdef TESTMODE
           do i=1,nx
           sprec_av(i)=sprec_av(i)+sprecip(i) 
           enddo
#endif
!cc

!cc save velocities after advection into advective velocities:
!cc (used as scratch)
       do k=1,nz
       do i=1,nx
       uxa(i,k)=ux(i,k)
       uza(i,k)=uz(i,k)
       enddo
       enddo

!cc finish thermodynamics 
       call thermo(theta,qv,qc,qr,ft,fqv,fqc,fqr,heat, nx, nz)
#ifdef TESTMODE
       do k=1,nz
       heat_av(k)=heat_av(k)+heat(k)
       enddo
#endif       
!cc add absorber:
       if(zab.lt.zz(nz))  &
             call absor(ux,uz,theta,qv,qc,qr,ft,fqv,fqc,fqr,tau, nx, nz)

!cc add buoyancy
       epsb=rv/rg-1.
       do k=1,nz
         recth0 = 1.0/th0(k)
         do i=1,nx
            scr1(i,k)=gg*( (theta(i,k)-th_e(k))*recth0 &
                 + epsb*(qv(i,k)-qv_e(k))-qc(i,k)-qr(i,k) )
!            scr1(i,k)=gg*( (theta(i,k)-thls(k))*recth0 &
!                 + epsb*(qv(i,k)-qvls(k))-qc(i,k)-qr(i,k) )
         enddo
       enddo

!cc filter in vertical
!cc       call integz(scr1,scr2,nx,nz)
       call integxz(scr1,scr2,nx,nz)

!cc apply
       do k=1,nz
       do i=1,nx
       uz(i,k) = uz(i,k)+.5*dt*scr1(i,k)
       enddo
       enddo


      call gcrk_1(p,scr1,scr2,ux,uz,itp,epp, nx, nz, lord, nplt)
      call prforc_1(p,scr1,scr2,ux,uz, nx, nz)
      do k=1,nz
      do i=1,nx
      ux(i,k)=scr1(i,k)
      uz(i,k)=scr2(i,k)
      enddo
      enddo

!cc calculate velocity forces (using saved velocities after advection):
       recdt = 1.0/dt

       do k=1,nz
       do i=1,nx
       fx(i,k)=(ux(i,k)-uxa(i,k))  *2.*recdt
       fz(i,k)=(uz(i,k)-uza(i,k))  *2.*recdt
       fy(i,k)=0.
       enddo
       enddo

!cc add noise if needed:
        if(itime/nfr_noise*nfr_noise.eq.itime) then
        print*,'+++ rndm noise: ta,qa,z: ',thn_a,qvn_a,nz_noise
        call noise(thn,nx,nz,nz_noise)
        call integxz_noise(thn,scr2,nx,nz,nz_noise)
        call zero_mean(thn,nx,nz,nz_noise,height)
        do k=1,nz
        do i=1,nx
        theta(i,k)=theta(i,k) + thn_a*thn(i,k)
        qv(i,k)   =   qv(i,k) + qvn_a*thn(i,k)
        enddo
        enddo
        endif
!cc max cloud water for noise enhancement when clouds present...
        do k=1,nz
        qcmx(k)=0.
        do i=1,nx
        qcmx(k)=max(qcmx(k),qc(i,k))
        enddo
        enddo
!cc max RH for noise enhancement...
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do k=1,nz
        rhmx(k)=0.
        thetme=th_e(k)/tm_e(k)
        pre=1.e5*thetme**e
        do i=1,nx
          tt=theta(i,k)/thetme
          delt=(tt-tt0)/(tt*tt0)
          expd(i) = d * delt
        end do

        call vexp(expd, expd, nx)

        do i=1,nx
          esw=ee0*expd(i)
          qvsw=a * esw /(pre-esw)
          rh=qv(i,k)/qvsw * 100.
          rhmx(k)=max(rhmx(k),rh)
        enddo
      enddo

!cc update clock (in minutes...)
       time = real(itime)*dt/60. 
#ifdef TESTMODE
!cc output and plot:
!cc
!cc output mean precip every timestep:
!c      sum=0.
!c      do i=1,nx-1
!c      sum=sum+sprecip(i)/real(nx-1)
!c      enddo
!c      write(65,*) time,sum

       if(amod(time+.1*dt/60.,dtout).lt.0.5*dt/60.) then
!ccc plot selected fields:
!      call plot_1(ux,uy,uz,theta)
!      call plot_2(theta,qv,qc,qr)
!cc analysis of output:
       print*,'   '
       call diagno_1(ux,uz,theta,scr1,scr2,den,nx,nz,time)
       call diagno_2(qv,qc,qr,scr1,scr2,den,nx,nz)
!cc 
       endif
       if(amod(time+.1*dt/60.,dttap).lt.0.5*dt/60.) then
       anor_inv=dt/(dttap*60.)
       do k=1,nz
       heat_av(k)=heat_av(k)*anor_inv
       radc_av(k)=radc_av(k)*anor_inv
       enddo
       do i=1,nx
       slat_av(i)=slat_av(i)*anor_inv
       ssen_av(i)=ssen_av(i)*anor_inv
       sprec_av(i)=sprec_av(i)*anor_inv
       enddo

       write(17) time,nx,nz
       write(17) ux,uz,theta,qv,qc,qr
       write(17) heat_av,radc_av,slat_av,ssen_av,sprec_av

!cc print mean surface fluxes and precip:
      sums=0.
      suml=0.
      sump=0.
      do i=1,nx-1
      sums=sums+ ssen_av(i)/real(nx-1)
      suml=suml+ slat_av(i)/real(nx-1)
      sump=sump+sprec_av(i)/real(nx-1)
      enddo
      sump=sump*3600.*24.
      print*,'mean sensible flux (W/m2): ',sums*1005
      print*,'mean latent flux (W/m2):   ',suml*2.5e6
      print*,'mean prec (mm/day):        ',sump
      print*,'            (W/m2):        ',sump/3.46*100.
!cc print mean surface fluxes and precip:

      do k=1,nz
       heat_av(k)=0.
       radc_av(k)=0.
       enddo
       do i=1,nx
       slat_av(i)=0.
       ssen_av(i)=0.
       sprec_av(i)=0.
       enddo
       endif
#endif
      enddo                     ! TIME LOOP

      call crcp_finalize

  end subroutine crcp

  subroutine crcp_finalize
    use surfflux_mod, only : qcmx,rhmx,th,qv,uu,vv
    use gcrk_mod, only : r,qr,ar
    use mpdat_mod, only : v1,v2,f1,f2,cp,cn,mx,mn


    ! from surfflux
    deallocate(qcmx,rhmx,th,qv,uu,vv)
    ! from gcrk_1
    deallocate(r,qr,ar)
    ! from mpdat_2d
    deallocate(v1,v2,f1, f2,cp,cn,mx,mn)


  end subroutine crcp_finalize
end module crcp_mod
