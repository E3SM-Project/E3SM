module rain_mod                                                                        
  use grid_init_mod
  use moist_init_mod
  implicit none
contains
      subroutine rain_fall(qr,tm_e,rho,uza, nx, nz)
      integer, intent(in) :: nx, nz
!cc modify vertical advective velocity for rain fallout
      real qr (nx, nz), uza (nx, nz + 1), tm_e (nz), rho (nz) 


      real lambdr,lambds,massr,masss
      real :: alim01, comb, fi, tm, td, tu, ad, au, vts, vtr, qpr, qps, qrv
      real :: vtf, coe_l, gc3, dens 
      real rfactor, sfactor
      real rec8, rec12
      real tmpr(nx), tmps(nx) 
      integer :: i, k

!cc statement functions:
      alim01(fi)=max(0.,min(1.,fi))
      comb(tm,td,tu,ad,au)= &
       alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad
#ifdef TESTMODE
      call tbeg('rain')
#endif
      gc3=dt*dzi

      rec8 = 1.0/8.0
      rec12 = 1.0/12.0

      rfactor = 1.0/(ar*anor*gamb1r)
      sfactor = 1.0/(as*anos*gamb1s)

      do k=2,nz
        dens = 0.5*(rho(k) + rho(k-1))
        coe_l = comb(tm_e(k),tdn,tup,0.,1.)   ! liquid part

        do i=1,nx
          qrv = 0.5*(qr(i,k)+ qr(i,k-1))
          qpr = qrv*coe_l         ! divide between rain and snow
          qps = qrv - qpr         ! divide between rain and snow

          tmpr(i) = rfactor * dens * (qpr+1.e-6)
          tmps(i) = sfactor * dens * (qps+1.e-6)
        end do

        call vlog(tmpr, tmpr, nx)
        call vlog(tmps, tmps, nx)
        tmpr(1:nx) = rec8  * tmpr(1:nx)
        tmps(1:nx) = rec12 * tmps(1:nx)
        call vexp(tmpr, tmpr, nx)
        call vexp(tmps, tmps, nx)

        do i=1,nx
          vtr = (cr*gambd1r/gamb1r) * tmpr(i)  ! terminal velocity
          vts = (cs*gambd1s/gamb1s) * tmps(i)  ! terminal velocity

          vtf = coe_l*vtr + (1. - coe_l)*vts   ! TERMINAL VELOCITY
          vtf = min(6.,vtf)
          uza(i,k) = uza(i,k) - vtf*dens*gc3
        end do
      end do
!ccc
!cC LOWER AND UPPER BOUNDARIES:
!cc lower:
        coe_l=comb(tm_e(1),tdn,tup,0.,1.)   ! liquid part
        dens=1.5*rho(1) - 0.5*rho(2)
        do i=1,nx
          qrv = max(0.,1.5*qr(i,1) - 0.5*qr(i,2))
          qpr = qrv*coe_l         ! divide between rain and snow
          qps = qrv - qpr           ! divide between rain and snow
          tmpr(i) = rfactor * dens * (qpr+1.e-6)
          tmps(i) = sfactor * dens * (qps+1.e-6)
        end do
        
        call vlog(tmpr, tmpr, nx)
        call vlog(tmps, tmps, nx)
        tmpr(1:nx) = rec8  * tmpr(1:nx)
        tmps(1:nx) = rec12 * tmps(1:nx)
        call vexp(tmpr, tmpr, nx)
        call vexp(tmps, tmps, nx)

        do i=1,nx
          vtr = (cr*gambd1r/gamb1r) * tmpr(i)  ! terminal velocity
          vts = (cs*gambd1s/gamb1s) * tmps(i)  ! terminal velocity

          vtf = coe_l*vtr + (1. - coe_l)*vts   ! TERMINAL VELOCITY
          vtf = min(6.,vtf)
          uza(i,1) = uza(i,1) - vtf*dens*gc3
        end do
!cc upper:
        coe_l = comb(tm_e(nz),tdn,tup,0.,1.)   ! liquid part
        dens = 1.5*rho(nz) - 0.5*rho(nz-1)
        do i=1,nx
          qrv = max(0.,1.5*qr(nz,1) - 0.5*qr(nz-1,2))
          qpr = qrv*coe_l         ! divide between rain and snow
          qps = qrv-qpr           ! divide between rain and snow
          tmpr(i) = rfactor * dens * (qpr+1.e-6)
          tmps(i) = sfactor * dens * (qps+1.e-6)
        end do

        call vlog(tmpr, tmpr, nx)
        call vlog(tmps, tmps, nx)
        tmpr(1:nx) = rec8  * tmpr(1:nx)
        tmps(1:nx) = rec12 * tmps(1:nx)
        call vexp(tmpr, tmpr, nx)
        call vexp(tmps, tmps, nx)

        do i=1,nx
          vtr = (cr*gambd1r/gamb1r) * tmpr(i)  ! terminal velocity
          vts = (cs*gambd1s/gamb1s) * tmps(i)  ! terminal velocity

          vtf = coe_l*vtr + (1. - coe_l)*vts   ! TERMINAL VELOCITY

          vtf = min(6.,vtf)

          uza(i,nz+1) = uza(i,nz+1) - vtf*dens*gc3
        end do

!       open(10, file='uza.opt', form='unformatted')
!       write(10) uza
!       close(10)
!       stop
#ifdef TESTMODE
      call tend('rain')
#endif
      return
    end subroutine rain_fall
  end module rain_mod
