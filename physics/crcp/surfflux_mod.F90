module surfflux_mod
  use grid_init_mod, only : dtl=>dt, dti, time, height
  use moist_init_mod, only : gg, cp, pi, hlatv
  use prof_init_mod, only : rho0, ux_e, uy_e
  use noise_mod, only : tpert=>thn_a, qpert=>qvn_a, nz_limit=> nz_noise
  implicit none

  real :: thsrf, qvsrf 
  real, private :: blh
  real, allocatable :: qcmx(:), rhmx(:) 
! used in bl_1d
  real, allocatable :: th(:), qv(:), uu(:), vv(:) 

contains          
      subroutine surfflux(theta,qv,ux,uy,fth,fqv,fx,fy,hfsen,hflat,&
           sst,ssqv,nx,nz)
!cc
!cc this a routine that calculates surface flux and distributes
!cc the flux in the vertical
!cc    spatial fluctuations are allowed to exist, see the logic below...
!cc
    integer, intent(in) :: nx, nz
    real, intent(in) :: theta (nx, nz), qv (nx, nz) 
    real, intent(in) :: ux (nx, nz)
    real, intent(in) :: uy (nx, nz)
    real, intent(inout) ::  fth (nx, nz), fqv (nx, nz) , fx (nx, nz), fy (nx, nz) 
    real, intent(out) ::  hfsen (nx), hflat (nx) 
    real, intent(in) :: sst, ssqv
    real :: ss1 (nx, nz), ss2 (nx, nz), scr (nx, nz) 
    real :: ss3 (nx, nz), ss4 (nx, nz) 
    real :: thsfl (nx), qvsfl (nx), usfl (nx), vsfl (nx) 
    real :: pthsfl (nx), pqvsfl (nx), pusfl (nx), pvsfl (nx) 

    real :: ventfc (2) , qvsfla, usfla, vsfla, thsfla
!    COMMON / surface / thsrf, qvsrf 
    real :: sigmau, ustar, aml, h_dis, thvsm, speed1, thvm, sstv

    integer :: i, k, k1
    real :: zr, qva, tha, dz1, expfun, zz, zrough, rho, uua, vva, anorm

!BGL  additional variables for optimization
      real recdz1

      double precision rand
!cc
      if(nz.gt.300) stop 'surfflux'

            rho=rho0(1)   !   surface density
            thsrf=sst     !   surface temperature
            qvsrf=ssqv    !   surface water vapor mixing ratio
            zrough=0.0002 !   surface roughness 

!cc NOTE: surface layer similarity applied only to mean profiles:
!cc get the mean surface air values:
        anorm=real(nx-1)
        vva=0.
        uua=0.
        tha=0.
        qva=0.
        do i=1,nx-1
        uua=uua+ux(i,1)/anorm
        vva=vva+uy(i,1)/anorm
        tha=tha+theta(i,1)/anorm
        qva=qva+qv(i,1)/anorm
        enddo

        zr=2.
        sstv=thsrf*(1.+.61*qvsrf)
        thvm=tha*(1.+.61*qva)
        thvsm=sstv-thvm
        speed1=sqrt(uua**2 + vva**2)
        sigmau=1.0 ! amplitude of wind gustiness

!c      print*,'in sflux2 sigmau,thvm,thvsm,speed1,zr: ',
!c     1    sigmau,thvm,thvsm,speed1,zr 

      call sflux2(sigmau,thvm,thvsm,speed1,zr,zrough,ustar,ventfc,aml)

!c      print*,'out sflux2 ustar,ventfc,aml: ',ustar,ventfc,aml

      do i=1,nx
!create surface momentum fluxes
      usfl(i) = rho0(1)*ventfc(1)*(0.-uua-.5*sigmau)
      vsfl(i) = rho0(1)*ventfc(1)*(0.-vva-.1*sigmau)
!create surface heat fluxes
      thsfl(i)=rho0(1)*ventfc(2)*(thsrf-tha)
      qvsfl(i)=rho0(1)*ventfc(2)*(qvsrf-qva)

      hfsen(i)=thsfl(i)
      hflat(i)=qvsfl(i)
      enddo

!cc get mean and perturbation of surface fluxes:
       thsfla=0.
       qvsfla=0.
        usfla=0.
        vsfla=0.
       do i=1,nx-1
       thsfla=thsfla+thsfl(i)/real(nx-1)
       qvsfla=qvsfla+qvsfl(i)/real(nx-1)
       usfla = usfla+ usfl(i)/real(nx-1)
       vsfla = vsfla+ vsfl(i)/real(nx-1)
       enddo

!c       print*,' sflux: s,l,u,v: ',thsfla*1.e3,qvsfla*2.5e6,
!c     1          usfla,vsfla
!ccc 
       do i=1,nx
       pthsfl(i)  = thsfl(i)   -thsfla
       pqvsfl(i)  = qvsfl(i)   -qvsfla
       pusfl(i)   =  usfl(i)   - usfla
       pvsfl(i)   =  vsfl(i)   - vsfla
       enddo

!ccc call 1D BL routine to deal with mean surface flux distribution:
       call bl_1d(thsfla,qvsfla,usfla,vsfla,theta,qv,fth,fqv, &
            ux,uy,fx,fy,nx,nz)

!ccc distribute perturbations from the mean:
      h_dis=.5*blh       !<-- e-folding scale
!c     h_dis=100.
      do k=1,nz
!cc note shift of flux positions:
        k1=min(k+1,nz)
        zz=.5*(height(k)+height(k1))
        expfun=exp(-zz/h_dis)
        do i=1,nx
          ss1(i,k)=pthsfl(i)*expfun
          ss2(i,k)=pqvsfl(i)*expfun
          ss3(i,k)= pusfl(i)*expfun
          ss4(i,k)= pvsfl(i)*expfun
        enddo
      enddo

      recdz1=2.0/(height(2)*rho0(1))
      do i=1,nx
!cc first level above the ground:
        fth(i,1)=fth(i,1)-(ss1(i,1)-pthsfl(i))*recdz1 * 2.
        fqv(i,1)=fqv(i,1)-(ss2(i,1)-pqvsfl(i))*recdz1 * 2.
         fx(i,1)= fx(i,1)-(ss3(i,1)- pusfl(i))*recdz1 * 2.
         fy(i,1)= fy(i,1)-(ss4(i,1)- pvsfl(i))*recdz1 * 2.
      end do
!cc higher levels:
      do k=2,nz
        k1=min(k+1,nz)
        recdz1=2.0/((height(k1)-height(k-1))*rho0(k))
        do i=1,nx
          fth(i,k)=fth(i,k)-(ss1(i,k)-ss1(i,k-1))*recdz1 *   2.
          fqv(i,k)=fqv(i,k)-(ss2(i,k)-ss2(i,k-1))*recdz1 *   2.
           fx(i,k)= fx(i,k)-(ss3(i,k)-ss3(i,k-1))*recdz1 *   2.
           fy(i,k)= fy(i,k)-(ss4(i,k)-ss4(i,k-1))*recdz1 *   2.
        enddo
      enddo

      return
    end subroutine surfflux

    subroutine bl_1d(tsfl,qsfl,usfl,vsfl,thl,qvl,fthl,fqvl, &
         ul,vl,ful,fvl, nx, nz)
!cc CRM arrays:
    integer, intent(in) :: nx, nz
    real, intent(in) :: tsfl, qsfl, usfl, vsfl
    real, intent(in) :: thl (nx, nz), qvl (nx, nz)  
    real, intent(inout) :: fthl (nx, nz), fqvl (nx, nz) 
    real, intent(in) ::ul (nx, nz), vl (nx, nz) 
    real, intent(inout) ::ful (nx, nz), fvl (nx, nz) 

    integer :: nzf
    real ::  fth (nz+1), fqv (nz+1), fuu (nz+1), fvv (nz+1) 
    real ::  hgt (nz), hgtf (nz+1) 

    !c arrays for communication                                             
    real ::  thf (nz), qvf (nz), uuf (nz), vvf (nz) 

    real :: qcmm, gammau, wstar, ri, gammat, rhmm, vel2, hbl,deltathv
    real :: thv, zz, gammaq, rho, fi_h, wt, dz1, qstar, capl, akm, akz
    real :: sum1, sum2, sum3, sum4, epsil, dt, tstar, ustar, ustar2, beta
    real :: zsrf, akt, fi_m, coe, wm, prand, t0, vonkar, thvsrf, tvsfl, ric
    integer :: kk, i, k, iter, ntime

!BGL  additional variables for optimization
      real xscale
      
      nzf=nz+1

      if(time.le.1.e-5) then
!cccccccccccccccccccc initial fields....
!cc initial BL height:
         blh=50.

         xscale = 1.0/real(nx-1)

         do k=1,nz
         sum1=0.
         sum2=0.
          do i=1,nx-1
          sum1=sum1+thl(i,k)*xscale
          sum2=sum2+qvl(i,k)*xscale
          enddo
         th(k)=sum1
         qv(k)=sum2
         uu(k)=ux_e(k)
         vv(k)=uy_e(k)
         enddo
!cccccccccccccccccccc end initial fields....
      endif

      do k=1,nz
      hgt(k)=height(k)
      enddo
      hgtf(1)=hgt(1)
      hgtf(nz+1)=hgt(nz)
      do k=2,nz
      hgtf(k)=.5*(hgt(k)+hgt(k-1))
      enddo

!ccc calculations in each time step of the CRM:
      vonkar=0.4
!      grav=9.8
      t0=300.
      epsil=1.e-8

      ric=0.25
      beta=gg/t0

      dt = dtl

!cc number of time steps:
      ntime=nint(dtl/dt)
      if(abs(ntime*dt-dtl).gt.0.0001) stop 'time steps do not match'

!cc derive CRM forcings on 1d BLM profiles:
!cc   CRM levels:

      xscale = 1.0/real(nx-1)

      do k=1,nz
        sum1=0.
        sum2=0.
        sum3=0.
        sum4=0.
        do i=1,nx-1
          sum1=sum1+thl(i,k)*xscale
          sum2=sum2+qvl(i,k)*xscale
          sum3=sum3+ ul(i,k)*xscale
          sum4=sum4+ vl(i,k)*xscale
        enddo
        thf(k)=(sum1-th(k))/dtl
        qvf(k)=(sum2-qv(k))/dtl
        uuf(k)=(sum3-uu(k))/dtl
        vvf(k)=(sum4-vv(k))/dtl
      enddo

!cCCCC  march forward in time in BLM:
           do iter=1,ntime

!c      print*,'t q sfl(W/m**2): ',tsfl*1.17*1005.,qsfl*1.17*2.5e6
!c      print*,'usfl,vsfl :',usfl,vsfl

!cc u*, T* ,q*, L:
      ustar2=sqrt(usfl**2 + vsfl**2)
      ustar=sqrt(ustar2)

      tvsfl=tsfl + .61*t0*qsfl
      tstar=-tvsfl/ustar
      qstar= -qsfl/ustar

      capL=ustar2/(vonkar*beta*(tstar+epsil))

!c       print*,'ustar,tstar,qstar: ',ustar,tstar,qstar
!c       print*,'MO l, blh: ',capL,blh

!cc  BL height:
!c below is when surface temp and qv is known:
      thvsrf=thsrf*(1.+.61*qvsrf)
      deltathv=thvsrf-th(1)*(1.+.61*qv(1))
!c end: below is when surface temp and qv is known:

!cc below is when only fluxes are known:
!c      if(capL.lt.0.) then
!c      wstar=(beta*blh*max(0.,tvsfl))**(1./3.)
!c      wm=(ustar**3 + .6*wstar**3)**(1./3.)
!c      deltathv = 8.5*tvsfl/wm
!c      deltathv=min(8.,deltathv)
!c      else
!c      deltathv=0.
!c      endif
!c      thvsrf=th(1)*(1.+.61*qv(1)) + deltathv
!cc end: below is when only fluxes are known:

!ccc
      deltathv = .5*deltathv
!ccc

       hbl=50.
       do k=2,nz
       kk=k
       zz=hgt(k)
       thv=th(k)*(1.+.61*qv(k))
       vel2=uu(k)**2 + vv(k)**2
       ri=zz*beta*(thv - thvsrf) / (vel2 + 0.1)
       if(ri.gt.ric) go to 100
       enddo
!c         stop 'blh .gt. model top'
       go to 101
100    hbl=hgt(kk-1)
101    hbl=max(50.,hbl)
!c this is diagnostic boundary layer height:
       blh=hbl

!c       print*,'++ time,blh,deltathv: ',time,blh,deltathv

!cc  w*:
      wstar=(beta*blh*max(0.,tvsfl))**(1./3.)

!cc temp and qv perturbations for CRM:
        tpert = deltathv*max(0.,tsfl) / max(1.e-3,tvsfl)
        qpert = deltathv*max(0.,qsfl) / max(1.e-3,tvsfl)

!c          print*,'tpert,qpert (tbad): ',tpert,qpert

        do k=2,nz
        nz_limit=k
        if(height(k).gt.blh) go to 3011
        enddo
        nz_limit=4
3011    nz_limit = max(4,nz_limit)

!cc enhance pert when cloud water present:
        qcmm=0.
        rhmm=0.
        do k=1,nz_limit
        qcmm=max(qcmm,qcmx(k))
        rhmm=max(rhmm,rhmx(k))
        enddo

        if(qcmm.lt.1.e-5 .and. rhmm.gt.70.) qcmm=3.e-4

        tpert=tpert+hlatv/cp*qcmm
        qpert=qpert+         qcmm

        tpert = max(0.02,min(1.  ,tpert))
        qpert = max(0.04e-3,min(2.e-3,qpert))

!cc fluxes:
      fuu(1)=usfl
      fvv(1)=vsfl
      fth(1)=tsfl
      fqv(1)=qsfl
      fuu(nzf)=0.
      fvv(nzf)=0.
      fth(nzf)=0.
      fqv(nzf)=0.

      if(capL.ge.0.) then
!cc stable case:
      gammau=0.
      gammat=0.
      gammaq=0.
      do k=2,nzf-1
      dz1=hgt(k)-hgt(k-1)
      rho=.5*(rho0(k)+rho0(k-1))
!cc velocity scale wt and wm for stable case:
      zz=min(blh,hgtf(k))
      fi_h=1.+5.*zz/capL
      if(zz.ge.capL) fi_h=5.+zz/capL
      wt=ustar/fi_h
      wm=wt
      akm=vonkar*wm*zz*(1. - zz/blh)**2
      akt=vonkar*wt*zz*(1. - zz/blh)**2
      fuu(k)=-rho*akm*((uu(k)-uu(k-1))/dz1 - gammau)
      fvv(k)=-rho*akm*((vv(k)-vv(k-1))/dz1 - gammau)
      fth(k)=-rho*akt*((th(k)-th(k-1))/dz1 - gammat)
      fqv(k)=-rho*akt*((qv(k)-qv(k-1))/dz1 - gammaq)
      enddo

         else

!cc unstable case
      zsrf=.1*blh
      do k=2,nzf-1
      dz1=hgt(k)-hgt(k-1)
      rho=.5*(rho0(k)+rho0(k-1))
      zz=min(blh,hgtf(k))
      if(zz.le.zsrf) then
      fi_h=(1. - 15.*zz/capL)**(-0.5)
      wt=ustar/fi_h
      fi_m=(1. - 15.*zz/capL)**(-0.33333)
      wm=ustar/fi_m
      gammau=0.
      gammat=0.
      gammaq=0.
      else
      wm=(ustar**3 + .6*wstar**3)**(1./3.)
      prand=(1-15.*.1*blh/capL)**(-.16667) + 7.2*vonkar*.1*wstar/wm
      wt=wm/prand
      gammau=0.
      gammat=7.2*wstar*tsfl/(wm**2*blh)
      gammaq=7.2*wstar*qsfl/(wm**2*blh)
      endif

      akm=vonkar*wm*zz*(1. - zz/blh)**2
      akt=vonkar*wt*zz*(1. - zz/blh)**2
      fuu(k)=-rho*akm*((uu(k)-uu(k-1))/dz1 - gammau)
      fvv(k)=-rho*akm*((vv(k)-vv(k-1))/dz1 - gammau)
      fth(k)=-rho*akt*((th(k)-th(k-1))/dz1 - gammat)
      fqv(k)=-rho*akt*((qv(k)-qv(k-1))/dz1 - gammaq)
      enddo

         endif

!cccc flux divergence:
      do k=1,nz
      dz1=hgtf(k+1)-hgtf(k)
      coe=dt/dz1/rho0(k)

      uu(k)=uu(k) - coe * (fuu(k+1)-fuu(k))    + uuf(k)*dt
      vv(k)=vv(k) - coe * (fvv(k+1)-fvv(k))    + vvf(k)*dt
      th(k)=th(k) - coe * (fth(k+1)-fth(k))    + thf(k)*dt
      qv(k)=qv(k) - coe * (fqv(k+1)-fqv(k))    + qvf(k)*dt

      enddo

           enddo
!cCCCC  END of march forward in time in BLM:

!ccc look at LBM levels only 
      do k=1,nz
         sum1=0.
         sum2=0.
         sum3=0.
         sum4=0.
          do i=1,nx-1
          sum1=sum1+thl(i,k)*xscale
          sum2=sum2+qvl(i,k)*xscale
          sum3=sum3+ ul(i,k)*xscale
          sum4=sum4+ vl(i,k)*xscale
          enddo
      do i=1,nx
      fthl(i,k)=fthl(i,k) + 2.*(th(k)-sum1)/dtl
      fqvl(i,k)=fqvl(i,k) + 2.*(qv(k)-sum2)/dtl
      ful(i,k)= ful(i,k)  + 2.*(uu(k)-sum3)/dtl
      fvl(i,k)= fvl(i,k)  + 2.*(vv(k)-sum4)/dtl
      enddo
      enddo

      return
    end subroutine bl_1d

      SUBROUTINE SFLUX2 (SIGMAU,THVM,THVSM,SPEED1,ZR,ZROUGH,USTAR, &
           VENTFC, MOLEN )
!c
!c     *** REVISED 13 DEC 1985: IMPROVED ITERATION SCHEME ***
!c     *** REVISED 23 SEP 1988: CHANGE STOPIT TO AVOID NONCONVERGENCE ***
!c     *** REVISED 10 NOV 1988: SPEED1 NO LONGER CHANGED BY ROUTINE ***
!c     *** IC(40) => IC(60) -- 8 SEP 89 ***
!c
!c     INPUT VARIABLES :
!c
!c     SIGMAU -- HORIZONTAL VELOCITY FLUCTUATION ( RMS )
!c     THVM -- VIRTUAL POTENTIAL TEMPERATURE AT ANEMOMETER LEVEL
!c     THVSM -- SURFACE-ANEMOMETER LEVEL DEFICIT OF VIRTUAL POTENTIAL
!c              TEMPERATURE
!c     SPEED1 -- HORIZONTAL WIND SPEED
!c     ZR -- HEIGHT OF ANEMOMETER LEVEL ABOVE SURFACE
!c     ZROUGH -- ROUGHNESS LENGTH
!c
!c     OUTPUT VARIABLES :
!c
!c     USTAR -- FRICTION VELOCITY
!c     VENTFC -- VENTILATION FACTORS
!c     MOLEN -- MONIN-OBUKHOV LENGTH
!c
!c     HOW TO CALCULATE SURFACE FLUXES OF U, V, T, Q:
!c
!c     U(2), V(2): values of U, V at height ZR
!c     T=dry static energy (or potential temperature)
!c     TSFC=surface value of T, T(2)=T at height ZR
!c     Q=water vapor mixing ratio
!c     QSFC=saturation mixing ratio for TSFC,PSFC; Q(2)=Q at height ZR
!c     GWET = ground wetness (range: 0 to 1)
!c
!c     UW(1) = - VENTFC(1) * U(2)
!c     VW(1) = - VENTFC(1) * V(2)
!c     WT(1) = VENTFC(2) * ( TSFC - T(2) )
!c     WQ(1) = VENTFC(2) * ( GWET * QSFC - Q(2) )
!c
      LOGICAL STABLE,STOPIT

    real, intent(in) :: zrough, sigmau, zr, thvm, thvsm, speed1
    real, intent(out) :: ustar,  MOLEN , VENTFC (2) 
    integer, parameter :: maxit=5
    real, parameter :: bus=0.74, crit=0.003, vk=0.35, delta=0.609
    real :: tem1, tem2, tem3, x, y, custar, ctstar, cui, cti, cuni, ctni, speedm
    real :: ct, zeta, cu
    integer :: it

!      DATA BUS,CRIT,MAXIT/0.74,0.003,5/
!      DATA PI,VK,GRAV,DELTA /3.1415927,0.35,9.806,0.609/
!c
!c
      STOPIT = .FALSE.
!c
!c*    SPEEDM = AMAX1 ( SPEEDM, 1.E-03 )
      SPEEDM = AMAX1 ( SPEED1, 1.E-03 )
!c
!c     NEUTRAL VALUES OF CU AND CT : CUN AND CTN
!c
      TEM1 = ALOG ( ZR / ZROUGH )
      CUNI = TEM1 / VK
      CTNI = CUNI * BUS
!c
!c     SURFACE - AIR DEFICIT OF VIRTUAL POTENTIAL TEMPERATURE : THVSM
!c
      STABLE = THVSM .LT. 0.
!c
!c     START ITERATION WITH NEUTRAL VALUES FOR CU AND CT
!c
      IT = 0
      CU = 1. / CUNI
      CT = 1. / CTNI
      IF ( .NOT. STABLE ) SPEEDM = AMAX1 ( SPEEDM, SIGMAU )
!c
    1 CONTINUE
!c
      IT = IT + 1
!c
      ZETA = - ZR * CT * VK * GG * THVSM / ( THVM * CU **2 &
           * SPEEDM **2 )
      IF ( STABLE ) GO TO 2
!c
!c     UNSTABLE OR NEUTRAL CASE
!c
      X = ( 1. - 15. * ZETA ) ** ( 1. / 4. )
      Y = ( 1. -  9. * ZETA ) ** ( 1. / 4. )
!c
      TEM2 = TEM1 - ( ALOG ( ( 1. + X **2 ) / 2. ) &
           + 2. * ALOG ( ( 1. + X ) / 2. ) - 2. * ATAN ( X ) + PI / 2. )
!c
      TEM3 = TEM1 - 2. * ALOG ( ( 1. + Y **2 ) / 2. )
!c
      CUI = TEM2 / VK
      CUI = AMAX1 ( CUI, 0.5 * CUNI )
      CTI = BUS * TEM3 / VK
      CTI =  AMAX1 ( CTI, 0.3 * CTNI )
!c
      GO TO 3
!c
    2 CONTINUE
!c
!c     STABLE CASE
!c
!c     ENFORCE ZETA LESS THAN 2.45 ( EQUIVALENT TO RICHARDSON NUMBER LESS
!c     THAN 0.9 * CRITICAL RICHARDSON NUMBER ).
!c
      IF ( ZETA .LT. 2.45 ) GO TO 4
      STOPIT = .TRUE.
      ZETA = 2.45
!c
    4 CONTINUE
!c
      TEM2 = TEM1 + 4.7 * ZETA
!c
      TEM3 = TEM1 + 4.7 / BUS * ZETA
!c
      CUI = TEM2 / VK
      CTI = BUS * TEM3 / VK
!c
    3 CONTINUE
      STOPIT = STOPIT .OR. IT .EQ. MAXIT
      IF ( STOPIT ) GO TO 5
!c
!c      CHECK FOR CONVERGENCE
!c
      CUSTAR = CU
      CTSTAR = CT
      CU = 1. / CUI
      CT = 1. / CTI
      STOPIT = ABS ( CU / CUSTAR - 1. ) .LE. CRIT &
           .AND.   ABS ( CT / CTSTAR - 1. ) .LE. CRIT
      GO TO 6
!c
    5 CU = 1. / CUI
      CT = 1. / CTI
    6 CONTINUE
!c
      IF ( .NOT. STOPIT ) GO TO 1
!c
!c     ITERATION COMPLETED. CALCULATE USTAR AND VENTFC
!c
      IF ( STABLE ) GO TO 7
!c
!c     UNSTABLE OR NEUTRAL CASE ( ALGORITHM REVISED 9/3/85 )
!c
      USTAR = CU * SPEEDM
      VENTFC(1) = CU * USTAR
      VENTFC(2) = CT * USTAR
!c
!c     CHECK THAT VENTFC EXCEEDS TOWNSEND'S (1964) FREE CONVECTION VALUE.
!c
      IF ( CTI .LT. 0.3 * CTNI ) &
           VENTFC(2) = AMAX1 ( VENTFC(2), 0.0019 * THVSM ** ( 1. / 3. ) )
      GO TO 8
!c
!c     STABLE CASE
!c
    7 USTAR = CU * SPEEDM
      VENTFC(1) = CU * USTAR
      VENTFC(2) = CT * USTAR
    8 CONTINUE
!c
!c     MONIN-OBUKHOV LENGTH
!c 
      ZETA = - ZR * CT * VK * GG * THVSM / ( THVM * CU **2 &
           * SPEEDM **2 )
      ZETA = AMAX1 ( ABS ( ZETA ), 1.E-06 ) * SIGN ( 1., ZETA )
      MOLEN = ZR / AMIN1 ( ZETA, 2.45 )

      RETURN
    END SUBROUTINE SFLUX2

  end module surfflux_mod
