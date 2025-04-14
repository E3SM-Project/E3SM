module gw_oro

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use gw_utils, only: r8
use gw_common, only: pver

implicit none
private
save

! Public interface
public :: gw_oro_init
public :: gw_oro_src

! 1/2 * horizontal wavenumber
real(r8) :: oroko2 = huge(1._r8)

contains

!==========================================================================

subroutine gw_oro_init(errstring)

  use gw_common, only: kwv
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  errstring = ""

  oroko2 = 0.5_r8 * kwv

end subroutine gw_oro_init

!==========================================================================

subroutine gw_oro_src(ncol, &
     u, v, t, sgh, pmid, pint, dpm, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_common, only: pgwv, fcrit2, rair
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(ncol)
  ! Midpoint and interface pressures.
  real(r8), intent(in) :: pmid(ncol,pver), pint(ncol,0:pver)
  ! Layer thickness (pressure delta).
  real(r8), intent(in) :: dpm(ncol,pver)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)
  ! Midpoint Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-pgwv:pgwv,0:pver)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,0:pver)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-pgwv:pgwv)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real(r8) :: hdsp(ncol)
  ! Max orographic standard deviation to use.
  real(r8) :: sghmax
  ! c=0 stress from orography.
  real(r8) :: tauoro(ncol)
  ! Averages over source region.
  real(r8) :: nsrc(ncol) ! B-V frequency.
  real(r8) :: rsrc(ncol) ! Density.
  real(r8) :: usrc(ncol) ! Zonal wind.
  real(r8) :: vsrc(ncol) ! Meridional wind.

  ! Difference in interface pressure across source region.
  real(r8) :: dpsrc(ncol)

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves
  real(r8), parameter :: orohmin = 10._r8
  ! min wind speed for orographic waves
  real(r8), parameter :: orovmin = 2._r8

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp = 2.0_r8 * sgh

  k = pver
  src_level = k-1
  rsrc = pmid(:,k)/(rair*t(:,k)) * dpm(:,k)
  usrc = u(:,k) * dpm(:,k)
  vsrc = v(:,k) * dpm(:,k)
  nsrc = nm(:,k)* dpm(:,k)

  do k = pver-1, pver/2, -1
     do i = 1, ncol
        if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
           src_level(i) = k-1
           rsrc(i) = rsrc(i) + pmid(i,k) / (rair*t(i,k))* dpm(i,k)
           usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
        end if
     end do
  end do

  do i = 1, ncol
     dpsrc(i) = pint(i,pver) - pint(i,src_level(i))
  end do

  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,pver))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,0) = ubm(:,1)

  ubi(:,1:pver-1) = midpoint_interp(ubm)

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  do i = 1, ncol
     if ((ubi(i,pver) .gt. orovmin) .and. (hdsp(i) .gt. orohmin)) then
        sghmax = fcrit2 * (ubi(i,pver) / nsrc(i))**2
        tauoro(i) = oroko2 * min(hdsp(i)**2, sghmax) * rsrc(i) * nsrc(i) &
             * ubi(i,pver)
     else
        tauoro(i) = 0._r8
        src_level(i) = pver
     end if
  end do

  ! Set the phase speeds and wave numbers in the direction of the source
  ! wind. Set the source stress magnitude (positive only, note that the
  ! sign of the stress is the same as (c-u).
  do k = pver, minval(src_level), -1
     where (src_level <= k) tau(:,0,k) = tauoro
  end do

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0._r8

end subroutine gw_oro_src

end module gw_oro
