module thermo_mod
  use grid_init_mod, only : dti
  use moist_init_mod
  use prof_init_mod, only : th_e, tm_e, rho0
#if !defined(_AIX) && !defined(_BGL)
  use mpdat_mod, only : cvmgm
#endif
contains   
  subroutine thermo(th,qv,qc,qr,fth,fqv,fqc,fqr,heat0, nx, nz)
    implicit none
    integer, intent(in) :: nx, nz
    real ::  th(nx, nz), qv(nx, nz), qc(nx, nz), qr(nx, nz), fth(nx, nz), &
         fqv(nx, nz), fqc(nx, nz), fqr(nx, nz)                                
    real :: heat0(nz), heat(nx, nz)
    integer :: l,n,nl

    integer i, k
    real :: cons, conr, aut, auti, a, b, c, d, e
    real :: thn, delta, fff, fffp, esi, esw, delt, qvs, qvsi, qvsw, qps, qpr
    real :: del2, autc, qcc, qci, tt, ventr, vents, ssw, ssi, tc, coe_l, thetme
    real :: cf1, dep_r, acc, acc_s, acc_r, dep_s, devp, dcol, dep, thfun, times, pre
    real :: g_acc_r, g_acc_s, g_dep_r, g_dep_s

    real :: lambdr,lambds,massr,masss

!BGL  additional variables for optimization
      real rectt0, time1, time2, rtc
      real rfactor, sfactor, recthetme, xscale
      real recth(nx), expd(nx), expb(nx), recw(nx), reci(nx)
      real tmpr(nx), tmps(nx), recthn(nx), recr(nx), recs(nx), fi, tm, td, tu, ad, au
      real alim01, comb,diamr(nx), diams(nx), rootr(nx), roots(nx), rer(nx), res(nx)
      integer, save :: call_count
      data call_count /0/

!c statement functions:
      alim01(fi)=max(0.,min(1.,fi))
      comb(tm,td,tu,ad,au)= &
           alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

!ondensation/evaporation
      l = nz
      n = nx
      nl = nx*nz
      call_count = call_count + 1

      time1 = rtc()
#ifdef TESTMODE
      call tbeg('thermo')
#endif

      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      rectt0 = 1.0/tt0

      do 100 k=1,l
        thetme = th_e(k)/tm_e(k)
        coe_l = comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution
        pre = 1.e5*thetme**e

        call vrec(recth, th(1,k), n)

        do i=1,n
          delt = rectt0 - thetme*recth(i)
          expd(i) = d * delt
          expb(i) = b * delt
        end do

        call vexp(expd, expd, n)
        call vexp(expb, expb, n)

        do i = 1, n
          recw(i) = pre - ee0*expd(i)
          reci(i) = pre - ee0*expb(i)
        end do

        call vrec(recw, recw, n)
        call vrec(reci, reci, n)

        do i = 1, n
          esw = ee0*expd(i)
          esi = ee0*expb(i)
          qvsw = a * esw * recw(i)
          qvsi = a * esi * reci(i)
          qvs = coe_l*qvsw + (1.-coe_l)*qvsi
!cc linearized condensation rate is next:
          cf1 = thetme*recth(i)
          cf1 = cf1*cf1
          cf1 = c*cf1*pre*recw(i)*d
          tmpr(i) = (qv(i,k)-qvs)/(1.+qvs*cf1)
        end do
!--->
!cc one Newton-Raphson iteration is next:
        do i = 1, n
          recthn(i) = th(i,k) + c*thetme*tmpr(i)
        end do

        call vrec(recthn, recthn, n)

        do i = 1, n
          delt = rectt0 - thetme*recthn(i)
          expd(i)  = d * delt
          expb(i)  = b * delt
        end do

        call vexp(expd, expd, n)
        call vexp(expb, expb, n)

        do i = 1, n
          recw(i) = pre - ee0*expd(i)
          reci(i) = pre - ee0*expb(i)
        end do

        call vrec(recw, recw, n)
        call vrec(reci, reci, n)

        do i = 1, n
          esw = ee0*expd(i)
          esi = ee0*expb(i)
          qvsw = a * esw * recw(i)
          qvsi = a * esi * reci(i)
          qvs = coe_l*qvsw + (1.-coe_l)*qvsi
          delta = tmpr(i)
          fff = qv(i,k) - delta - qvs
          cf1 = thetme*recthn(i)
          cf1 = cf1*cf1
          cf1 = c*cf1*pre*recw(i)*d
          fffp = -1. -qvs*cf1
          delta = delta - fff/fffp
!cc end of the iteration; if required, it can be repeated
!--->
          delta = min( qv(i,k), max(-qc(i,k),delta) )
          qv(i,k) = qv(i,k) - delta
          qc(i,k) = qc(i,k) + delta
          th(i,k) = th(i,k) + c*thetme*delta

          heat(i,k) = c*thetme*delta*dti

          delta = min( qv(i,k), max(-qc(i,k),delta) )
          fqv(i,k) = -delta*2.*dti
          fth(i,k) = -c*thetme*fqv(i,k)
          fqc(i,k) = -fqv(i,k)

          heat(i,k) = heat(i,k) + c*thetme*delta*dti

        end do
  100 end do

!cc remove trace of water variables:
!      nl=n*l
      do i=1,nl
        qc(i,1)=cvmgm(0.,qc(i,1),qc(i,1)-1.e-9)
        qr(i,1)=cvmgm(0.,qr(i,1),qr(i,1)-1.e-10)
      enddo


!ompute moist forces update
      do 300 k=1,l
        thetme=th_e(k)/tm_e(k)
        recthetme = 1.0/thetme
        pre=1.e5*thetme**e
        coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution

        call vrec(recth, th(1,k), n)

        do i=1,n
          delt = rectt0 - thetme*recth(i)
          expd(i) = d * delt
          expb(i) = b * delt
        end do

        call vexp(expd, expd, n)
        call vexp(expb, expb, n)

        do i = 1, n
          recw(i) = pre - ee0*expd(i)
          reci(i) = pre - ee0*expb(i)
        end do

        call vrec(recw, recw, n)
        call vrec(reci, reci, n)

        rfactor = rho0(k)/(ar*anor*gamb1r)
        sfactor = rho0(k)/(as*anos*gamb1s)

        do i=1,n
          qpr = qr(i,k)*coe_l   ! divide between rain and snow
          qps = qr(i,k) - qpr
          tmpr(i) = rfactor*(qpr + 1.e-6)
          tmps(i) = sfactor*(qps + 1.e-6)
        end do

#ifdef _BGL
        call vsqrt(tmpr, tmpr, n)
        call vsqrt(tmpr, tmpr, n)
#else
        call vlog(tmpr, tmpr, n)
        tmpr(1:n) = 0.25*tmpr(1:n)
        call vexp(tmpr, tmpr, n)
#endif

#ifdef _BGL
        call vcbrt(tmps, tmps, n)
#else
        call vlog(tmps, tmps, n)
        tmps(1:n) = (1./3.)*tmps(1:n)
        call vexp(tmps, tmps, n)
#endif

        do i = 1, n
          recr(i) = anor*tmpr(i)
          recs(i) = anos*tmps(i)
        end do
   
        call vrec(recr, recr, n)
        call vrec(recs, recs, n)

        do i = 1, n
          qpr = qr(i,k)*coe_l                ! divide between rain and snow
          qps = qr(i,k) - qpr                ! divide between rain and snow

          massr = rho0(k)*(qpr+1.e-7) * recr(i)  ! mass
          masss = rho0(k)*(qps+1.e-7) * recs(i)  ! mass

          diamr(i) = massr/ar ! diameter
          diams(i) = masss/as ! diameter
        end do

#ifdef _BGL
        call vcbrt(diamr, diamr, n)
#else
        call vlog(diamr, diamr, n)
        diamr(1:n) = (1./3.)*diamr(1:n)
        call vexp(diamr, diamr, n)
#endif

        call vsqrt(rootr, diamr, n)

        call vsqrt(diams, diams, n)

        call vsqrt(roots, diams, n)
        call vsqrt(roots, roots, n)

        do i = 1, n
          rer(i) = cr*diamr(i)*rootr(i)/2.e-5   ! Reynolds number
          res(i) = cs*diams(i)*roots(i)/2.e-5   ! Reynolds number
        end do

        call vsqrt(rer, rer, n)
        call vsqrt(res, res, n)

        do i = 1, n
          esw = ee0*expd(i)
          esi = ee0*expb(i)
          qvsw = a * esw * recw(i)
          qvsi = a * esi * reci(i)

          ssw = qv(i,k) / qvsw      ! saturation ratio
          ssi = qv(i,k) / qvsi      ! saturation ratio

          qcc = qc(i,k)*coe_l                ! divide between ice and water
          qci = qc(i,k) - qcc                ! divide between ice and water

!C AUTOCONVERSION:
!c rain - Berry:
          del2 = 1.e3*rho0(k)*qcc
          autc = 1./rho0(k)*1.67e-5*del2*del2 / &
               (5. + .0366*dconc/(ddisp*(del2+1.E-6)))
!c snow:
          tt = th(i,k)*recthetme
          tc = tt - tt0
          times = min(1.e3,(3.56*tc+106.7)*tc+1.e3) ! time scale for
          auti = qci/times
          AUT = autc + auti

!C GROWTH:
          conr = anor*tmpr(i)! concentration
          cons = anos*tmps(i)! concentration

          qpr = qr(i,k)*coe_l                ! divide between rain and snow
          qps = qr(i,k) - qpr                ! divide between rain and snow

          ventr = max(1.,.78+.27*rer(i))  ! ventilation factor
          vents = max(1.,.65+.39*res(i))  ! ventilation factor

          thfun = 1.e-7/(2.2*tm_e(k)/esw+2.2e-2/tm_e(k))  ! thermodynamic fun.

          rfactor = (pi/4.)*cr*er*alphr*rho0(k)*qc(i,k)
          sfactor = (pi/4.)*cs*es*alphs*rho0(k)*qc(i,k)

          g_acc_r = rfactor*diamr(i)**2*rootr(i) ! growth
          g_acc_s = sfactor*diams(i)**2*roots(i) ! growth

          g_dep_r = 4.*pi*diamr(i)/betr*(ssw-1.)*ventr*thfun   ! growth/evap
          g_dep_s = 4.*pi*diams(i)/bets*(ssi-1.)*vents*thfun   ! growth/evap

          acc_r = conr * g_acc_r * qpr / (qpr + 1.e-9)
          acc_s = cons * g_acc_s * qps / (qps + 1.e-9)

          ACC = acc_r + acc_s  ! growth by accretion

          dep_r = conr * g_dep_r * qpr / (qpr + 1.e-9)
          dep_s = cons * g_dep_s * qps / (qps + 1.e-9)

          DEP = dep_r + dep_s  ! growth by deposition

          dcol = 2.*(AUT + ACC)
          dcol = min(dcol,  2.*dti*qc(i,k)+fqc(i,k))
          devp = 2.*DEP
          devp = max(devp, -2.*dti*qr(i,k)-dcol)
!c
          fqr(i,k) = devp + dcol
          fqc(i,k) = fqc(i,k) - dcol
          fqv(i,k) = fqv(i,k) - devp
          fth(i,k) = fth(i,k) + c*devp*thetme

          heat(i,k) = heat(i,k) + 0.5*c*devp*thetme
        end do

  300 end do

      xscale = 1.0/real(nx-1)

      do k=1,nz
        heat0(k)=0.
        do i=1,nx-1
          heat0(k)=heat0(k)+heat(i,k)*xscale
        enddo
      enddo

      time2 = rtc()

!     if (call_count .eq. 100) then
!     open(10, file='thermo.opt', form='unformatted')
!     write(10) qv
!     write(10) qc
!     write(10) qr
!     write(10) th
!     write(10) heat
!     write(10) fqv
!     write(10) fth
!     write(10) fqc
!     write(10) fqr
!     write(10) heat0
!     close(10)
!     print *, 'msec per call =', 1.0e3*(time2 - time1)
!     call tend('thermo')
!     call tend('main')
!     call tprt()
!     stop
!     endif
#ifdef TESTMODE 
      call tend('thermo')
#endif
      return
    end subroutine thermo
  end module thermo_mod
