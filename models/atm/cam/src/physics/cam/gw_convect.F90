module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use gw_utils, only: r8

use gw_common, only: pver, pgwv

implicit none
private
save

public :: gw_convect_init
public :: gw_beres_src

! Dimension for heating depth.
integer :: maxh
! Dimension for mean wind in heating.
integer :: maxuh

! Index for level at 700 mb.
integer :: k700

! Table of source spectra.
real(r8), allocatable :: mfcc(:,:,:)

contains

!==========================================================================

subroutine gw_convect_init(k700_in, mfcc_in, errstring)
  ! Index at 700 mb.
  integer, intent(in) :: k700_in
  ! Source spectra to keep as table.
  real(r8), intent(in) :: mfcc_in(:,:,:)
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  integer :: ierr

  errstring = ""

  k700 = k700_in

  ! First dimension is maxh.
  maxh = size(mfcc_in,1)
  ! Second dimension is -maxuh:maxuh (size 2*maxuh+1).
  maxuh = (size(mfcc_in,2)-1)/2

  allocate(mfcc(maxh,-maxuh:maxuh,-pgwv:pgwv), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return
  mfcc = mfcc_in

end subroutine gw_convect_init

!==========================================================================

subroutine gw_beres_src(ncol, ngwv, lat, u, v, netdt, &
     zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, &
     hdepth, maxq0)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  use gw_common, only: dc, cref
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
!------------------------------Arguments--------------------------------
  ! Column and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, ngwv

  ! Column latitudes [rad].
  real(r8), intent(in) :: lat(ncol)

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Heating rate due to convection.
  real(r8), intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)

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

  ! Heating depth and maximum heating in each column.
  real(r8), intent(out) :: hdepth(ncol), maxq0(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at 700mb.
  real(r8) :: u700(ncol), v700(ncol)
  ! 3.14...
  real(r8), parameter :: pi = 4._r8*atan(1._r8)

  ! Maximum heating rate.
  real(r8) :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: mini(ncol), maxi(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer  :: Umini,Umaxi
  ! Mean wind in heating region.
  real(r8) :: uh(ncol)
  ! Min/max projected wind value in each column.
  real(r8) :: Umin(ncol), Umax(ncol)
  ! Source level tau for a column.
  real(r8) :: tau0(-PGWV:PGWV)
  ! Speed of convective cells relative to storm.
  integer :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Heating rate conversion factor.
  real(r8), parameter :: CF = 20._r8
  ! Averaging length.
  real(r8), parameter :: AL = 1.0e5_r8

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

  tau = 0.0_r8
  hdepth = 0.0_r8
  q0 = 0.0_r8
  tau0 = 0.0_r8

  !------------------------------------------------------------------------
  ! Determine 700 mb layer wind and unit vectors, then project winds.
  !------------------------------------------------------------------------

  ! Just use the 700 mb interface values for the source wind speed and
  ! direction (unit vector).

  u700 = u(:,k700)
  v700 = v(:,k700)

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(u700, v700, xv, yv, ubi(:,k700))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,0) = ubm(:,1)

  ubi(:,1:pver-1) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  mini = 0
  maxi = 0
  do k = pver, 1, -1
     do i = 1, ncol
        if (mini(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              mini(i) = k
              maxi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0_r8) mini(i) = k
           end if
        else if (maxi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              maxi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (.not. (netdt(i,k) > 0.0_r8)) maxi(i) = k
           end if
        end if
     end do
     ! When all done, exit
     if (all(maxi /= 0)) exit
  end do

  ! Heating depth in km.
  hdepth = [ ( (zm(i,maxi(i))-zm(i,mini(i)))/1000._r8, i = 1, ncol ) ]
  ! Confine depth to table range.
  hdepth = min(hdepth, real(maxh, r8))

  ! Maximum heating rate.
  do k = minval(maxi), maxval(mini)
     where (k >= maxi .and. k <= mini)
        q0 = max(q0, netdt(:ncol,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0 = q0*24._r8*3600._r8

  ! Multipy by conversion factor
  q0 = q0 * CF

  ! Taking ubm at 700 mb to be the storm speed, find the cell speed where
  ! the storm speed is > 10 m/s.
  CS = int(sign(max(abs(ubm(:,k700))-10._r8, 0._r8), ubm(:,k700)))

  uh = 0._r8
  do k = minval(maxi), maxval(mini)
     where (k >= maxi .and. k <= mini)
        uh = uh + ubm(:,k)/(mini-maxi+1)
     end where
  end do

  uh = uh - real(CS, r8)

  ! Limit uh to table range.
  uh = min(uh, real(maxuh, r8))
  uh = max(uh, -real(maxuh, r8))

  ! Speeds for critical level filtering.
  Umin =  pgwv*dc
  Umax = -pgwv*dc
  do k = minval(maxi), maxval(mini)
     where (k >= maxi .and. k <= mini)
        Umin = min(Umin, ubm(:,k))
        Umax = max(Umax, ubm(:,k))
     end where
  end do

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if depth >= 2.5 km, else set tau0 = 0.
     !---------------------------------------------------------------------

     if ((hdepth(i) >= 2.5_r8) .and. (abs(lat(i)) < (pi/2._r8))) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = mfcc(NINT(hdepth(i)),NINT(uh(i)),:)

        ! Shift spectrum so that it is relative to the ground.
        shift = -nint(real(CS(i), r8)/dc)
        tau0 = cshift(tau0,shift)

        ! Adjust magnitude.
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        Umini = max(nint(Umin(i)/dc),-PGWV)
        Umaxi = min(nint(Umax(i)/dc),PGWV)

        if (Umaxi > Umini) then
           tau0(Umini:Umaxi) = 0.0_r8
        end if

        tau(i,-ngwv:ngwv,maxi(i)) = tau0(-ngwv:ngwv)

     end if ! depth >= 2.5

  enddo

  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = maxi
  tend_level = maxi

  ! Set phase speeds; just use reference speeds.
  c = spread(cref, 1, ncol)

end subroutine gw_beres_src

end module gw_convect
