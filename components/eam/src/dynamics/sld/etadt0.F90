subroutine etadt0(lat     ,nlon    ,                            &
                  rcoslat ,d       ,u       ,v       ,dpsl    , &
                  dpsm    ,pdel    ,pstar   ,etadot  )
!
!-----------------------------------------------------------------------
!
! Purpose:
! to compute (1/ps)etadot(dp/deta) for time step 0.
! NOTE:  Computing only the bottom "plev" levels of
!        a half-level field.
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use hycoef, only : hybd, hybi, nprlev
  implicit none

!------------------------------Arguments--------------------------------
!
  integer, intent(in) :: lat                     ! latitude
  integer, intent(in) :: nlon                    ! number of longitudes

  real(r8), intent(in)   :: rcoslat(nlon)        ! 1 / cos(lat)
  real(r8), intent(in)   :: d     (plon,plev)    ! divergence
  real(r8), intent(in)   :: u     (plon,plev)    ! zonal wind * cos(lat)
  real(r8), intent(in)   :: v     (plon,plev)    ! meridional wind * cos(lat)
  real(r8), intent(in)   :: dpsl  (plon)         ! longitudinal component of grad ln(ps)
  real(r8), intent(in)   :: dpsm  (plon)         ! latitudinal  component of grad ln(ps)
  real(r8), intent(in)   :: pdel  (plon,plev)    ! layer thicknesses (pressure)
  real(r8), intent(in)   :: pstar (plon)         ! surface pressure
  real(r8), intent(out)  :: etadot(plon,plev)    ! vertical velocity in eta coordinates
!
!---------------------------Local workspace-----------------------------
!
  integer i,k                   ! indices
  real(r8) ddpk (plon)          ! parial sum of   d  *delta p
  real(r8) ddpn (plon)          ! complete sum of d  *delta p
  real(r8) vkdp (plon)          ! V dot grad(ps)
  real(r8) vpdsk(plon)          ! partial sum of V dot grad(ps) delta b
  real(r8) vpdsn(plon)          ! complete sum of V dot grad(ps) delta
  real(r8) pspsl(plon)          ! Ps*d(lnPs)/d(long.)
  real(r8) pspsm(plon)          ! Ps*d(lnPs)/d(lat. )
!
!-----------------------------------------------------------------------
!
! zero auxiliary fields in blank common
!
  do i = 1,nlon
     ddpk    (i) = 0.0_r8
     ddpn    (i) = 0.0_r8
     vpdsk   (i) = 0.0_r8
     vpdsn   (i) = 0.0_r8
     pspsl   (i) = pstar(i)*dpsl(i)
     pspsm   (i) = pstar(i)*dpsm(i)
     etadot(i,plev) = 0.0_r8
  end do
!
! calculate some auxiliary quantities first
!
  do k = 1,plev
!
! sum(l = 1,k) [ d(l)*dp(l) ]
!
     do i = 1,nlon
        ddpn (i) = ddpn (i) + d  (i,k)*pdel(i,k)
     end do
     if(k.ge.nprlev) then
        do i = 1,nlon
!
! sum(l = 1,k) [ v(k)*grad(ps)*dB(l) ]
!
           vkdp(i) = rcoslat(i)*(u(i,k)*pspsl(i) + v(i,k)*pspsm(i))
           vpdsn(i) = vpdsn(i) + vkdp(i)*hybd(k)
        end do
     end if
  end do
!
! Compute (1/ps)*etadot(dp/deta)
!
  do k = 1,plev-1
!
! sum(l = 1,k) [ d(l)*dp(l) ]
!
     do i = 1,nlon
        ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
        etadot(i,k) = -ddpk(i)
     end do
!
! sum(l = 1,k) [ v(l)*grad(pstar)*dB(l) ] +
! B(k+1/2) * sum(j = 1,K) [ d(j)*dp(j) + v(j)*grad(pstar)*dB(j) ]
!
     if(k.ge.nprlev) then
        do i = 1,nlon
           vkdp(i)  = rcoslat(i)*(u(i,k)*pspsl(i) + v(i,k)*pspsm(i))
           vpdsk(i) = vpdsk(i) + vkdp(i)*hybd(k)
           etadot(i,k) = etadot(i,k) - vpdsk(i) + hybi(k+1)*(ddpn(i)+vpdsn(i))
        end do
     end if
  end do
!
  do k = 1,plev-1
     do i = 1,nlon
        etadot(i,k) = etadot(i,k)/pstar(i)
     end do
  end do
!
  return
end subroutine etadt0

