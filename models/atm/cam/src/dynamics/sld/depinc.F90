subroutine depinc(jcen    ,dt      ,ra      ,locgeo  ,lvert   , &
                  lam     ,phi     ,ump     ,vmp     ,wmp     , &
                  upr     ,vpr     ,lamdp   ,phidp   ,lampr   , &
                  phipr   ,sigpr   ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Compute departure point increments based upon interpolated wind components.
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
  use pmgrid,       only: plon, plev
  use scanslt,      only: i1, plond, platd
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: jcen              ! index of lat slice(extend)
  real(r8), intent(in)   :: dt                ! time step (seconds)
  real(r8), intent(in)   :: ra                ! 1./(radius of earth)
  logical , intent(in)   :: locgeo            ! local geodesic flag
  logical , intent(in)   :: lvert             ! flag to compute vertical trajectory
  real(r8), intent(in)   :: lam  (plond)      ! long. coord of model grid
  real(r8), intent(in)   :: phi  (platd)      ! lat.  coord of model grid
  real(r8), intent(inout):: ump  (plon,plev)  ! interpolated u field
  real(r8), intent(inout):: vmp  (plon,plev)  ! interpolated v field
  real(r8), intent(in)   :: wmp  (plon,plev)  ! interpolated w field
  real(r8), intent(out)  :: upr  (plon,plev)  ! interpolated u field (local geodesic)
  real(r8), intent(out)  :: vpr  (plon,plev)  ! interpolated v field (local geodesic)
  real(r8), intent(in)   :: lamdp(plon,plev)  ! zonal      dep pt. coord.
  real(r8), intent(in)   :: phidp(plon,plev)  ! meridional dep pt. coord.
  real(r8), intent(out)  :: lampr(plon,plev)  ! trajectory increment (x-direction)
  real(r8), intent(inout):: phipr(plon,plev)  ! trajectory increment (y-direction)
  real(r8), intent(out)  :: sigpr(plon,plev)  ! trajectory increment (vertical)
  integer , intent(in)   :: nlon              ! number of longitudes for this lat
!
!---------------------------Local workspace-----------------------------
!
  real(r8) phi0    ! Current latitude (radians)           
  real(r8) cphi0   ! cosine of latitude
  real(r8) sphi0   ! sine   of latitude

  real(r8) cdlam   ! |
  real(r8) clamgc  ! |
  real(r8) cphid   ! |
  real(r8) cphigc  ! |
  real(r8) dlamgc  ! | -- temporary variables
  real(r8) sdlam   ! |
  real(r8) slamgc  ! |
  real(r8) sphid   ! |
  real(r8) sphigc  ! |

  integer i        ! |
  integer ii       ! | -- indices
  integer k        ! |
!
!-----------------------------------------------------------------------
!
  phi0    = phi(jcen)
  cphi0   = cos(phi0)
  sphi0   = sin(phi0)
!
! Compute trajectory increment
!
  do k = 1,plev
!
! Place u/v on unit sphere
!
     do i = 1,nlon
        ump(i,k) = ump(i,k)*ra
        vmp(i,k) = vmp(i,k)*ra
        if (lvert) sigpr(i,k) = -dt*wmp(i,k)
     end do
!
! If near the pole, convert to g.c.
!
     if (locgeo) then
!
! Compute winds at departure points in g.c.
!
        do i = 1,nlon
           ii = i + i1 - 1
           dlamgc = lam(ii) - lamdp(i,k)
           sdlam  = sin( dlamgc )
           cdlam  = cos( dlamgc )
           sphid  = sin( phidp(i,k) )
           cphid  = cos( phidp(i,k) )
           sphigc = sphid*cphi0 - cphid*sphi0*cdlam
           cphigc = cos( asin( sphigc ) )
           slamgc = -sdlam*cphid/cphigc
           clamgc = cos( asin( slamgc ) )
           vpr(i,k) = (vmp(i,k)*(cphid*cphi0 + sphid*sphi0*cdlam) &
                - ump(i,k)*sphi0*sdlam)/cphigc
           upr(i,k) = (ump(i,k)*cdlam + vmp(i,k)*sphid*sdlam + &
                vpr(i,k)*slamgc*sphigc)/clamgc
!
! Compute trajectory increment between arrival and departure points.
!
           lampr(i,k) = -dt*upr(i,k)/cos(phipr(i,k))
           phipr(i,k) = -dt*vpr(i,k)
        end do
     else
        do i = 1,nlon
           lampr(i,k) = -dt*ump(i,k)/cos(phidp(i,k))
           phipr(i,k) = -dt*vmp(i,k)
        end do
     end if
  end do
!
  return
end subroutine depinc

