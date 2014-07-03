
subroutine grdxy(dlam    ,lam     ,phi     ,w       ,sinlam  , &
                 coslam  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define the "extended" grid used in the semi-Lagrangian transport
! scheme.  The longitudes are equally spaced and the latitudes are
! Gaussian.  The global grid is extended to include "wraparound" points
! on all sides.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat
  use scanslt,      only: nxpt, jintmx, plond, platd, nlonex
  use gauaw_mod, only: gauaw
  implicit none

!------------------------------Parameters-------------------------------
  integer, parameter :: istart = nxpt+1         ! index for first model long.
  integer, parameter :: jstart = nxpt+jintmx+1  ! index for first model lat.
  integer, parameter :: jstop  = jstart-1+plat  ! index for last  model lat.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  real(r8), intent(out) :: dlam(platd)          ! longitudinal increment
  real(r8), intent(out) :: lam   (plond,platd)  ! long. coords. in extended grid
  real(r8), intent(out) :: phi   (platd)        ! lat.  coords. in extended grid
  real(r8), intent(out) :: w     (plat)         ! Gaussian weights
  real(r8), intent(out) :: sinlam(plond,platd)  ! sin(lam)
  real(r8), intent(out) :: coslam(plond,platd)  ! cos(lam)
!
! dlam    Length of increment in longitude grid.
! lam     Longitude values in the extended grid.
! phi     Latitude values in the extended grid.
! w       Gauss weights for latitudes in the global grid.  (These sum
!         to 2.0 like the ones in CCM1.)
! sinlam  Sine of longitudes in global grid (no extension points).
! coslam  Cosine of longitudes in global grid (no extension points).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,j,ig            ! indices
  integer nlond             ! extended long dim
  real(r8) lam0             ! lamda = 0
  real(r8) pi               ! 3.14...
  real(r8) wrk(platd)       ! work space
!-----------------------------------------------------------------------
!
  lam0 = 0.0_r8
  pi = 4._r8*atan(1._r8)
!
! Interval length in equally spaced longitude grid.
!
  do j=1,platd
     dlam(j) = 2._r8*pi/real(nlonex(j),r8)
!
! Longitude values on extended grid.
!
     nlond = nlonex(j) + 1 + 2*nxpt
     do i = 1,nlond
        lam(i,j) = real(i-istart,r8)*dlam(j) + lam0
     end do
  end do
!
! Compute Gauss latitudes and weights.  On return; phi contains the
! sine of the latitudes starting closest to the north pole and going
! toward the south; w contains the corresponding Gauss weights.
!
  call gauaw(phi     ,w       ,plat    )
!
! Reorder and compute latitude values.
!
  do j = jstart,jstop
     wrk(j) = asin( phi(jstop-j+1) )
  end do
  phi(jstart:jstop) = wrk(jstart:jstop)
!
! North and south poles.
!
  phi(jstart-1) = -pi/2.0_r8
  phi(jstop +1) =  pi/2.0_r8
!
! Extend Gauss latitudes below south pole so that the spacing above
! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
!
  if( jstart > 2 )then
     do j = 1,jstart-2
        phi(j) = -pi - phi(2*jstart-2-j)
     end do
  end if
!
! Analogously for Northern Hemisphere
!
  if( platd > jstop+1 )then
     do j = jstop+2,platd
        phi(j) = pi - phi(2*jstop+2-j)
     end do
  end if
!
! Sine and cosine of longitude.
!
  do j=1,platd
     ig = 0
     do i = istart,nlonex(j)+nxpt
        ig = ig + 1
        sinlam(ig,j) = sin( lam(i,j) )
        coslam(ig,j) = cos( lam(i,j) )
     end do
  end do

  return
end subroutine grdxy
