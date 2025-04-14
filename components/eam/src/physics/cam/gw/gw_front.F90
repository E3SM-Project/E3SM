module gw_front

!
! This module handles gravity waves from frontal sources, and was extracted
! from gw_drag in May 2013.
!

use gw_utils, only: r8
use gw_common, only: pver, pgwv, cref

implicit none
private
save

public :: gw_front_init
public :: gw_cm_src

! Tuneable settings.

! Frontogenesis function critical threshold.
real(r8) :: frontgfc = huge(1._r8)

! Level at which to check the frontogenesis function to determine when
! waves will be launched.
integer :: kfront = 0

! Average value of gaussian over gravity wave spectrum bins, multiplied by
! background source strength (taubgnd).
real(r8), allocatable :: fav(:)

contains

! Initialize module parameters.
subroutine gw_front_init(taubgnd, frontgfc_in, kfront_in, errstring)
  use gw_common, only: dc

  ! The following are used to set the module data.
  real(r8), intent(in) :: taubgnd
  real(r8), intent(in) :: frontgfc_in
  integer, intent(in) :: kfront_in
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  ! Parameters to calculate fav (average value of gaussian over bin).

  ! Integration interval to get bin average.
  real(r8), parameter :: dca  = 0.1_r8
  ! Width of gaussian in phase speed.
  real(r8), parameter :: c0   = 30._r8

  ! Integration bounds and indices.
  integer :: l, n
  real(r8) :: cmn, cmx

  integer :: ierr

  errstring = ""

  frontgfc = frontgfc_in
  kfront = kfront_in

  ! Allocate and calculate fav.
  allocate(fav(-pgwv:pgwv), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return

  if (pgwv > 0) then
     do l = -pgwv,pgwv
        ! Lower bound of bin.
        cmn = cref(l) - 0.5_r8*dc
        cmx = cref(l) + 0.5_r8*dc
        ! Loop over integration intervals in bin.
        fav(l) = 0.5_r8 * dca * (exp(-(cmn/c0)**2) + exp(-(cmx/c0)**2))
        do n = 1,nint(dc/dca)-1
           fav(l) = fav(l) + dca * exp(-((cmn+n*dca)/c0)**2)
        end do
     end do
     fav = fav/dc
  else
     fav(0) = 1._r8
  end if

  ! Multiply by source strength.
  fav = taubgnd * fav
  ! Prohibit wavenumber 0.
  fav(0) = 0._r8

end subroutine gw_front_init

!==========================================================================
subroutine gw_cm_src(ncol, ngwv, kbot, u, v, frontgf, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Driver for multiple gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
  ! Column and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, ngwv

  ! Index of source interface.
  integer, intent(in) :: kbot

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Frontogenesis function.
  real(r8), intent(in) :: frontgf(:,:)

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
  ! Column and wavenumber indices.
  integer :: k, l

  ! Whether or not to launch waves in this column.
  logical :: launch_wave(ncol)

  ! Zonal/meridional wind averaged over source region.
  real(r8) :: usrc(ncol), vsrc(ncol)

  !------------------------------------------------------------------------
  ! Determine the source layer wind and unit vectors, then project winds.
  !------------------------------------------------------------------------

  ! Just use the source level interface values for the source wind speed
  ! and direction (unit vector).
  src_level = kbot
  tend_level = kbot
  usrc = 0.5_r8*(u(:,kbot+1)+u(:,kbot))
  vsrc = 0.5_r8*(v(:,kbot+1)+v(:,kbot))

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,kbot))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, kbot
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,0) = ubm(:,1)

  ubi(:,1:kbot-1) = midpoint_interp(ubm(:,1:kbot))

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------

  tau = 0._r8

  ! GW generation depends on frontogenesis at specified level (may be below
  ! actual source level).
  launch_wave = (frontgf(:ncol,kfront) > frontgfc)

  do l = 0, ngwv
     where (launch_wave)
        tau(:,l,kbot) = fav(l)
        tau(:,-l,kbot) = fav(l)
     end where
  end do

  ! Set phase speeds as reference speeds plus the wind speed at the source
  ! level.
  c = spread(cref, 1, ncol) + spread(abs(ubi(:,kbot)),2,2*ngwv+1)

end subroutine gw_cm_src

end module gw_front
