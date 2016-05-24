      module physicsUtest
        use kinds, only : real_kind
        use physical_constants, only : gg=>g, cp, rg=>rgas, rv=>rwater_vapor

        use control_mod, only : dt, np, ntime

        real (kind=real_kind) :: time

        real (kind=real_kind), allocatable :: pre0(:),rho0(:),th0(:)
        real (kind=real_kind), allocatable :: th_e(:),ux_e(:),uy_e(:)
        real (kind=real_kind), allocatable :: qv_e(:),tm_e(:)
        

        real (kind=real_kind), parameter :: hlatv=2.53e6, hlats=2.84e6
        real (kind=real_kind), parameter :: tt0=273.16, ee0=611.
        
        real (kind=real_kind), parameter:: thsrf=302.4706,qvsrf=2.8133571e-2

      contains


      subroutine prof_init
        implicit none
        
!c initial profiles; note that lowest model level (975 mb) is taken
!c as within boundary layer using surface sounding data...

!c NOTE: since model will develop its own equilibrium sounding,
!c the initial sounding is not important, thus only a simple
!c interpolation is used here...
      
      integer, parameter :: npin=23
      real(kind=real_kind) ::  press(npin),temp(npin),thin(npin),&
               zin(npin), vap(npin),uu(npin),vv(npin),rhoin(npin)


      real(kind=real_kind) :: thetme, coe2, a, b, c, d, e
      integer :: k, kk, iisn
       data zin /npin*0./
       data press  / &
       1008.00, 991.25, 945.50, 893.79, 836.06, 772.82, 705.22, &
        635.05, 564.48, 495.73, 430.71, 370.78, 316.72, 268.82, &
        226.98, 190.82, 159.87, 133.55, 111.29,  92.56,  52.31, &
         22.08,   9.32/
       data temp  / &
         25.26,  24.13,  21.04,  18.66,  16.50,  13.41,   9.06, &
          3.73,  -1.51,  -6.97, -14.09, -22.44, -30.57, -39.60, &
        -48.69, -57.40, -65.21, -72.58, -76.71, -74.98, -74.98, &
        -74.98, -74.98/
       data vap  / &
        0.178E+02, 0.172E+02, 0.156E+02, 0.134E+02, 0.111E+02, &
        0.888E+01, 0.631E+01, 0.487E+01, 0.396E+01, 0.200E+01, &
        0.984E+00, 0.806E+00, 0.370E+00, 0.135E+00, 0.599E-01, &
        0.258E-01, 0.123E-01, 0.582E-02, 0.367E-02, 0.589E-02, &
        0.104E-02, 0.247E-02, 0.585E-02/
       data uu / npin*5./
!       data uu / 3.,3.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,
!     1   14.,15.,16.,17.,18.,19.,20.,21.,22.,23./
       data vv / npin*5./

!c surface data:

      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

       do k=1,npin
       thetme=(1.e3/press(k))**(rg/cp)
       temp(k)=temp(k)+273.16
       vap(k)=vap(k)*1.e-3
       thin(k)=temp(k)*thetme
       rhoin(k)=(press(k)*1.e2)/(rg*temp(k)*(1.+.62*vap(k)))
       enddo

      print*,'INPUT SOUNDING'
      do k=1,npin
      print 921,press(k),temp(k),thin(k),vap(k),rhoin(k)
921   format(1x,'p,t,th,qv,rho: ',3f8.1,2e13.4)
      enddo

!c lowest level: assume it is within mixed layer:
       th0(np)=thin(1)
       th_e(np)=thin(1)       
       qv_e(np)=vap(1)
       thetme=(1.e3/pre0(np))**(rg/cp)
       tm_e(np)=th_e(np)/thetme
       rho0(np)=(pre0(np)*1.e2)/(rg*tm_e(np)*(1.+.62*qv_e(np)))
       ux_e(np)=uu(1)
       uy_e(np)=vv(1)

      do k=np-1,1,-1
        do kk=2,npin
          iisn=kk-1
          if(press(kk).le.pre0(k)) go to 665
        enddo
!       print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
        stop 'SOUNDING'
 665    continue
!       print*,'iisn=',iisn
        coe2=(pre0(k)-press(iisn))/(press(iisn+1)-press(iisn))
        th_e(k)=coe2*thin(iisn+1) + (1.-coe2)*thin(iisn)
        th0(k)=th_e(k)
        qv_e(k)=coe2*vap(iisn+1) + (1.-coe2)*vap(iisn)
        ux_e(k)=coe2*uu(iisn+1) + (1.-coe2)*uu(iisn)
        uy_e(k)=coe2*vv(iisn+1) + (1.-coe2)*vv(iisn)
        thetme=(1.e3/pre0(k))**(rg/cp)
        tm_e(k)=th_e(k)/thetme
        rho0(k)=(pre0(k)*1.e2)/(rg*tm_e(k)*(1.+.62*qv_e(k)))
        enddo

      print*,'PROFILES'
      do k=np,1,-1
       print 200,pre0(k),th0(k),rho0(k),th_e(k), &
                  tm_e(k),qv_e(k)*1.e3,ux_e(k)
 200    format(1x,'p,th0,rho0,the,tme,qve,ue:', &
               2f7.1,f5.2,2f7.1,e10.3,f6.1)
      enddo

      return
      end subroutine prof_init

      subroutine radcool(ft)
        implicit none
      real(kind=real_kind) :: ft(np)
!
      real(kind=real_kind) :: day, coe
      integer :: k

       day=24.*3600.
       do k=1,np
       coe=th_e(k)/tm_e(k)
       if(pre0(k).lt.200.) coe=0.
       ft(k)=ft(k) - 1.5/day*coe
       enddo

      return
      end subroutine radcool

      subroutine surfflux(ux,uy,fx,fy,theta,qv,fth,fqv,np1)
        implicit none
        integer :: np1
        real (kind=real_kind), dimension(np1) :: fth(np1),fqv(np1)      
        real (kind=real_kind), dimension(np1) :: theta(np1),qv(np1)      
        real (kind=real_kind), dimension(np1) :: ux(np1),uy(np1)      
        real (kind=real_kind), dimension(np1) :: fx(np1),fy(np1)      
        real (kind=real_kind) :: wind, w_star, ssen, smom, slat, &
             drag, p_dis, coeth, coemo,coeqv
!c
      p_dis=50.       ! depth of the mixed layer (hPa)
      p_dis=p_dis*1.e2

!c simple surface flux formula (can be made more sophisticated)
      drag=1.3e-3
      coeth=drag
      coeqv=drag
      coemo=drag
      w_star=1.
      wind=sqrt(ux(np)**2 + uy(np)**2 + w_star**2)
      ssen=coeth*rho0(np)*wind*(thsrf-theta(np))  ! sensible flux
      slat=coeqv*rho0(np)*wind*(qvsrf-   qv(np))  ! latent flux
      smom=coemo*rho0(np)*wind*(   0.-   wind  )  ! momentum flux

                      print*,gg/p_dis*ssen,gg,p_dis,ssen

           print*,' ++ ssen, slat: ',ssen*1005., slat*2.5e6
           print*,' ++ mom fl: ',smom

!c we assume that the temperature and moisture at the first model level 
!c represent values within the boundary layer; the temperature
!c change corresponds to the divergence of the surface flux:
      fth(np)=fth(np) + gg/p_dis * ssen 
      fqv(np)=fqv(np) + gg/p_dis * slat
!c surface momentum fluxes set to zero for SCM tests...
!      fx(np) = fx(np) + gg/p_dis * ux(np)/(wind+1.e-6) *smom
!      fy(np) = fy(np) + gg/p_dis * uy(np)/(wind+1.e-6) *smom
      return
      end subroutine surfflux


      subroutine m_adjustment(l,ux,uy,fx,fy,th,qv,fth,fqv,con_precip,cbmf)
      use convect43c
      implicit none
!c profiles  
      integer, intent(in) :: l
      real (kind=real_kind) :: q(l),t(l),p(l), &
           tp(l),qs(l),peh1(l), &
            thetme(l),ph(l+1),FT(l),FQ(l),ph1(l+1) &
           ,u(l),v(l),fu(l),fv(l),tra(l,1),ftra(l,1) 
      real (kind=real_kind) fth(l),fqv(l),th(l),qv(l)
      real (kind=real_kind) ux(l),uy(l),fx(l),fy(l)
      real (kind=real_kind) con_precip, wd, delt, tprime, &
                 qprime, precip, cbmf1, cbmf
      real (kind=real_kind) :: thi,del,ees,pg,t00,hlat,&
           a,b,c,d,e,cap,capi
      integer :: nl,nd,k,iconv_2,iconv_3,iconv_4,iconv_5,ntra,iflag
      integer :: k1,lk,iconv_add

!ccccccccccccccccccccc
!c comments: 1. we do not distinguish between mixing ratios and
!c           specific humidity (see comments to Kerry's code).

      iconv_add=0
      iconv_4=0
      iconv_3=0
      iconv_2=0

      t00=tt0
      hlat=hlatv
      cap=rg/cp
      capi=1./cap
      a=rg/rv
      b=hlat/(rv*t00)
      c=hlat/cp
      d=hlat/rv
      e=-cp/rg
      lk=18 + 1 ! highest level reached by convection + 1

        cbmf1=CBMF

       do k1=1,np
       k=np-k1+1
        thetme(k1)=th_e(k1)/tm_e(k1)
        p(k)=pre0(k1)
        t(k)=th(k1)/thetme(k1)
        q(k)=qv(k1)
        thi=1./th(k1)
       del=b*thetme(k1)*t00*thi
       ees=ee0*exp(b-del)
       qs(k)=a*ees/(p(k)*1.e2-ees)
       u(k)=ux(k1)
       v(k)=uy(k1)
       FT(k)=0.
       FQ(k)=0.
       enddo

       do k=1,l-1
        ph(k+1)=.5*(p(k)+p(k+1))
       enddo
        ph(1)=2.*p(1)-ph(2)
        ph(l+1)=max(1._real_kind, 2._real_kind*p(l)-ph(l))
        pg=p(1)


       ND=l
       NL=lk
       DELT=dt
       NTRA=1
       TRA=0.
#ifdef DEBUG
!ccccccccccccccccc
          print*,'      ------------- B E F O R E -----------------'
          print*,' -- T: ',T(1),th(np),thetme(np),th_e(np),tm_e(np)
          print*,' -- Q: ',Q
          print*,' -- QS: ',QS
          print*,' -- u: ',u
          print*,' -- v: ',v
          print*,' -- TRA: ',TRA
          print*,' -- P: ',P
          print*,' -- PH: ',PH
          print*,' -- ND: ',ND
          print*,' -- NL: ',NL
          print*,' -- NTRA: ',NTRA
          print*,' -- DELT: ',DELT
          print*,' -- IFLAG: ',IFLAG
          print*,' -- FT: ',FT
          print*,' -- FQ: ',FQ
          print*,' -- Fu: ',Fu
          print*,' -- Fv: ',Fv
          print*,' -- FTRA: ',FTRA
          print*,' -- PRECIP: ',PRECIP
          print*,' -- TPRIME: ',TPRIME
          print*,' -- QPRIME: ',QPRIME
          print*,' -- CBMF1: ',CBMF1
          print*,'      ------------- B E F O R E -----------------'
          print*,'      ------------------------- -----------------'
#endif
!
        CALL CONVECT  &
         (T,   Q,    QS,     U,    V,      TRA,    P,    PH, &
          ND,  NL,   NTRA,   DELT, IFLAG,  FT,     FQ,   FU, &
          FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, CBMF1   )
#ifdef DEBUG
          print*,'      ------------- A F T E R  -----------------'
          print*,' -- T: ',T
          print*,' -- Q: ',Q
          print*,' -- QS: ',QS
          print*,' -- u: ',u
          print*,' -- v: ',v
          print*,' -- TRA: ',TRA
          print*,' -- P: ',P
          print*,' -- PH: ',PH
          print*,' -- ND: ',ND
          print*,' -- NL: ',NL
          print*,' -- NTRA: ',NTRA
          print*,' -- DELT: ',DELT
          print*,' -- IFLAG: ',IFLAG
          print*,' -- FT: ',FT
          print*,' -- FQ: ',FQ
          print*,' -- Fu: ',Fu
          print*,' -- Fv: ',Fv
          print*,' -- FTRA: ',FTRA
          print*,' -- PRECIP: ',PRECIP
          print*,' -- TPRIME: ',TPRIME
          print*,' -- QPRIME: ',QPRIME
          print*,' -- CBMF1: ',CBMF1
          print*,'      ------------- A F T E R  -----------------'
          print*,'      ------------------------------------------'
#endif
       con_precip=PRECIP
       CBMF=cbmf1
       
       print *,'FTH: ',fth(np),ft(1),thetme(np)

       do k=1,np
       k1=np-k+1
       thetme(k1)=th_e(k1)/tm_e(k1)
       fth(k1)=fth(k1)+FT(k)*thetme(k1)
       fqv(k1)=fqv(k1)+FQ(k)
       fx(k1) = fx(k1)+Fu(k)
       fy(k1) = fy(k1)+Fv(k)
       th(k1)=t(k)*thetme(k1)
       qv(k1)=q(k)
       enddo



       IF(IFLAG.eq.1) iconv_add=iconv_add+1
       IF(IFLAG.eq.2) iconv_2=iconv_2+1
       IF(IFLAG.eq.3) iconv_3=iconv_3+1
       IF(IFLAG.eq.4) iconv_4=iconv_4+1
       IF(IFLAG.eq.5) iconv_5=iconv_5+1

       print *,'** ADJUSTMENT DONE IN ' &
     ,iconv_add,' OF 1  COLUMN**'
       print *,'CFL condition violated in               ' &
                                 ,iconv_4,'COLUMN'
       print *,'NO convection cloud base to high        ' &
                                 ,iconv_3,'COLUMN'
       print *,'NO convection lfted condensation to high' &
                                 ,iconv_2,'COLUMN'
       print *,'************************************************'

      RETURN
      END subroutine m_adjustment
#ifdef NCARG
      subroutine plot(ux,uy,th,qv)
      real(kind=real_kind) ::  ux(np),uy(np),th(np),qv(np) 
      real ::   pnd(np) ! crazy nondim pressure: 0 surface, 1 top
      character*50 lhead

      real :: qvlocal(np)

      do k=1,np
        pnd(k)=(1000.-pre0(k))/1000. 
      enddo

!c ux
      call set(.2,.5,.15,.9, -20.,20.,0.,1.,1)
      call labmod('(f5.0)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,1,10)
      call curved(real(ux),pnd,np)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,951) time/(24.*60.)
  951 format('  ux and uy      time (days)  ',f6.1)
      CALL plchhq(.55,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.35,0.05, 'ux (m/s)', 0.016,0.,0)
      CALL plchhq(.07,0.525, 'nondim pre', 0.016,90.,0)
!c uy
      call set(.6,.9,.15,.9, -20.,20.,0.,1.,1)
      call labmod('(f5.0)','(f5.0)',5,5,2,2,20,20,0)
      call periml(2,5,1,10)
      call curved(real(uy),pnd,np)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.75,0.05, 'uy (m/s)', 0.016,0.,0)
      call frame

!c theta
      call set(.2,.5,.15,.9, 250.,400.,0.,1.,1)
      call labmod('(f5.0)','(f5.1)',5,5,2,2,20,20,0)
      call periml(3,5,1,10)
      call curved(real(th(2:)),pnd(2:),np-1)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time/(24.*60.)
  952 format('  theta and qv      time (days)  ',f6.1)
      CALL plchhq(.55,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.35,0.05, 'theta (K)', 0.016,0.,0)
      CALL plchhq(.07,0.525, 'nondim pre', 0.016,90.,0)
!c qv
      do k=1,np
      qvlocal(k)=qv(k)*1.e3
      enddo
      call set(.6,.9,.15,.9, 0.,25.,0.,1.,1)
      call labmod('(f5.0)','(f5.0)',5,5,2,2,20,20,0)
      call periml(5,5,1,10)
      call curved(qvlocal,pnd,np)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.75,0.05, 'qv (g/kg)', 0.016,0.,0)
      call frame

!      do k=1,np
!      qv(k)=qv(k)*1.e-3
!      enddo

      return
      end subroutine plot
#endif

      end module physicsUtest


      program moist_convec
      use physicsUtest
      use control_mod
      implicit none
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Single Column Model using Emanuel's convection scheme...
!  -- pressure as vertival coordinate --
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c

      
      real (kind=real_kind), allocatable:: tt(:),prec(:),sw(:)

!C  MODEL VARIABLES AND FORCES
      real (kind=real_kind), allocatable:: ux(:),uy(:),fx(:),fy(:)
      real (kind=real_kind), allocatable:: theta(:),qv(:),ft(:),fqv(:)

      integer :: itime, it, k
      real (kind=real_kind) ::  cbmf, con_precip
      real :: eps, dp


      call readnl()
      print *, 'np= ',np,' ntime=',ntime
      allocate(tt(ntime),prec(ntime),sw(ntime))
      allocate(ux(np),uy(np),fx(np),fy(np))
      allocate(theta(np),qv(np),ft(np),fqv(np))

      allocate(pre0(np),rho0(np),th0(np))
      allocate(th_e(np),ux_e(np),uy_e(np))
      allocate(qv_e(np),tm_e(np))

#ifdef NCARG
      call opngks
      call gsclip(0)
#endif

!c grid:
      time=0.

      pre0(1)=press_toa
      dp = (srfpress-press_toa)/float(np-1)

      do k=2,np
         pre0(k)=pre0(k-1)+dp
      enddo

       do k=np,1,-1
       print*,' -- k,pre: ',k,pre0(k)
       enddo

!c initialize model profiles:
      call prof_init

       do k=1,np
        theta(k)=th_e(k)
        qv(k)=qv_e(k)
        ux(k)=ux_e(k)
        uy(k)=uy_e(k)

        fx(k)=0.
        fy(k)=0.
        ft(k)=0.
        fqv(k)=0.
       enddo
      
!cc plot initial fields:
#ifdef NCARG
       call plot(ux,uy,theta,qv)
#else
       write(17) np
       write(17) 0.0,ux,uy,theta,qv
#endif


       CBMF=0.


!CCC MARCH FORWARD IN TIME:

       do itime=1,ntime   ! TIME LOOP

          print*,'*** itime, time: ',itime,time

          do k=1,np
             fx(k)=0.
             fy(k)=0.
             ft(k)=0.
             fqv(k)=0.
          enddo

!c surface flux, radiative cooling:
          call surfflux(ux,uy,fx,fy,theta,qv,ft,fqv,np)
          print *, 'FT1 ',ft(np)
          call radcool(ft)
          print *, 'FT2 ',ft(np)
!c convective parameterization
       call m_adjustment(np,ux,uy,fx,fy,theta,qv,ft,fqv,con_precip,CBMF)

       print*,'--- time, con_prec: ',time,con_precip

       do k=1,np
          theta(k)=theta(k)+dt*ft(k)
          qv(k)   =   qv(k)+dt*fqv(k)
          ux(k)   =   ux(k)+dt* fx(k)
          uy(k)   =   uy(k)+dt* fy(k)
       enddo

       prec(itime)=con_precip
       sw(itime)=sqrt(ux(np)**2 + uy(np)**2)

!c update clock (in minutes...)
       time=itime*dt/60. 

!c output and plot:
#ifdef NCARG
       if(itime/nplot*nplot.eq.itime) then
         call plot(ux,uy,theta,qv)
       endif
#else
       if(itime/nplot*nplot.eq.itime) then
       write(17) time,ux,uy,theta,qv
       endif
#endif

       enddo      ! TIME LOOP

!c plot precip vs time:
 
      eps=1.e-2

      do it=1,ntime
        tt(it)=float(it)*dt/(3600.*24)
        prec(it)=max(real(prec(it)),(eps))
        print*,'-- t,pr(mm/d),sw(m/s): ',tt(it),prec(it),sw(it)
      enddo

      print *,'precip     (min,max,sum): ', &
             minval(prec),maxval(prec),sum(prec)
      print *,'wind speed (min,max,sum): ', &
             minval(sw),maxval(sw),sum(sw)
#ifdef NCARG
!      call set(.2,.8,.2,.6,0.,100.,eps,10.,2)
      call set(.2,.8,.2,.6,0.,100.,0.,20.,1)
      call labmod('(f5.0)','(f5.2)',5,5,2,2,20,20,0)
!      call periml(5,4,1,9)
      call periml(5,4,4,5)
      call curved(real(tt),real(prec),ntime)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.50,0.67,' precip versus time',0.016,0.,0)
      CALL plchhq(.50,0.05, 'time (days)', 0.016,0.,0)
      CALL plchhq(.07,0.40, 'precip (mm/day)', 0.016,90.,0)
      call frame

!c    finished...

        call clsgks
#endif
        stop
        end

