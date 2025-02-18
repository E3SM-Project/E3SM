module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!
use cam_logfile, only: iulog
use spmd_utils,  only: masterproc
use gw_utils,    only: r8
use gw_common,   only: pver, pgwv

implicit none
private
save

public :: gw_convect_init
public :: gw_beres_src

! Dimension for heating depth.
integer :: maxh
! Dimension for mean wind in heating.
integer :: maxuh

! Index for level for storm/steering flow (usually 700 mb)
integer :: k_src_wind

! Table of source spectra.
real(r8), allocatable :: mfcc(:,:,:)

contains

!==========================================================================

subroutine gw_convect_init( plev_src_wind, mfcc_in, errstring)
  use ref_pres, only: pref_edge
  real(r8), intent(in) :: plev_src_wind        ! reference pressure value [Pa] to set k_src_wind (previously hardcoded to 70000._r8)
  real(r8), intent(in) :: mfcc_in(:,:,:)       ! Source spectra to keep as table
  character(len=*), intent(out) :: errstring   ! Report any errors from this routine
  integer :: ierr
  integer :: k

  errstring = ""

  do k = 0, pver
    if ( pref_edge(k+1) < plev_src_wind ) k_src_wind = k+1
  end do

  if (masterproc) write (iulog,*) 'gw_convect: steering flow level = ',k_src_wind

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
     hdepth, maxq0_out, maxq0_conversion_factor, hdepth_scaling_factor, &
     hdepth_min, storm_speed_min, &
     use_gw_convect_old)
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

  ! Heating conversion factor
  real(r8), intent(in) :: maxq0_conversion_factor

  ! Scaling factor for the heating depth
  real(r8), intent(in) :: hdepth_scaling_factor

  ! minimum hdepth for for spectrum lookup table
  real(r8), intent(in) :: hdepth_min

  ! minimum convective storm speed
  real(r8), intent(in) :: storm_speed_min

  ! switch for restoring legacy method
  logical, intent(in) :: use_gw_convect_old

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
  real(r8), intent(out) :: hdepth(ncol), maxq0_out(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional source wind
  real(r8) :: u_src(ncol), v_src(ncol)
  ! 3.14...
  real(r8), parameter :: pi = 4._r8*atan(1._r8)

  ! Maximum heating rate.
  real(r8) :: maxq0(ncol)

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
  integer :: storm_speed(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! fixed parameters (we may want to expose these in the namelist for tuning)
  real(r8), parameter :: tau_avg_length       = 100e3 ! spectrum averaging length [m]
  real(r8), parameter :: heating_altitude_max = 20e3  ! max altitude [m] to check for max heating

  ! note: the heating_altitude_max is probably not needed because there is
  ! rarely any convective heating above this level and the performance impact
  ! of skipping the iteration over higher levels is likely negilible.

  integer :: ndepth_pos
  integer :: ndepth_tot

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

  tau    = 0.0_r8
  hdepth = 0.0_r8
  maxq0  = 0.0_r8
  tau0   = 0.0_r8

  !------------------------------------------------------------------------
  ! Determine source layer wind and unit vectors, then project winds.
  !------------------------------------------------------------------------

  ! source wind speed and direction
  u_src = u(:,k_src_wind)
  v_src = v(:,k_src_wind)

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(u_src, v_src, xv, yv, ubi(:,k_src_wind))

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

  ! Find indices for the top and bottom of the heating range.
  mini = 0
  maxi = 0

  if (use_gw_convect_old) then
    !---------------------------------------------------------------------
    ! original version used in CAM4/5/6 and EAMv1/2/3
    do k = pver, 1, -1
      do i = 1, ncol
        if (mini(i) == 0) then
          ! Detect if we are outside the maximum range (where z = 20 km).
          if (zm(i,k) >= heating_altitude_max) then
            mini(i) = k
            maxi(i) = k
          else
            ! First spot where heating rate is positive.
            if (netdt(i,k) > 0.0_r8) mini(i) = k
          end if
        else if (maxi(i) == 0) then
          ! Detect if we are outside the maximum range (z = 20 km).
          if (zm(i,k) >= heating_altitude_max) then
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
    !---------------------------------------------------------------------
  else
    !---------------------------------------------------------------------
    ! cleaner version that addresses bug in original where heating max and
    ! depth were too low whenever heating <=0 occurred in the middle of
    ! the heating profile (ex. at the melting level)
    do i = 1, ncol
      do k = pver, 1, -1
        if ( zm(i,k) < heating_altitude_max ) then
          if ( netdt(i,k) > 0.0_r8 ) then
            ! Set mini as first spot where heating rate is positive
            if ( mini(i)==0 ) mini(i) = k
            ! Set maxi to current level
            maxi(i) = k
          end if
        else
          ! above the max check if indices were found
          if ( mini(i)==0 ) mini(i) = k
          if ( maxi(i)==0 ) maxi(i) = k
        end if
      end do
    end do
    !---------------------------------------------------------------------
   end if


  ! Heating depth in km.
  hdepth = [ ( (zm(i,maxi(i))-zm(i,mini(i)))/1000._r8, i = 1, ncol ) ]
  ! Confine depth to table range.
  hdepth = min(hdepth, real(maxh, r8))

  ! apply tunable scaling factor for the heating depth
  hdepth = hdepth * hdepth_scaling_factor

  ! Maximum heating rate.
  do k = minval(maxi), maxval(mini)
     where (k >= maxi .and. k <= mini)
        maxq0 = max(maxq0, netdt(:ncol,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0_out = maxq0*24._r8*3600._r8

  ! Multipy by conversion factor
  maxq0 = maxq0 * maxq0_conversion_factor

  ! Taking ubm at assumed source level to be the storm speed, 
  ! find the cell speed where the storm speed is > storm_speed_min
  storm_speed = int(sign(max(abs(ubm(:,k_src_wind))-storm_speed_min, 0._r8), ubm(:,k_src_wind)))

  uh = 0._r8
  do k = minval(maxi), maxval(mini)
     where (k >= maxi .and. k <= mini)
        uh = uh + ubm(:,k)/(mini-maxi+1)
     end where
  end do

  uh = uh - real(storm_speed, r8)

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
      if ((hdepth(i) >= hdepth_min)) then
         ndepth_pos = 0
         ndepth_tot = 0
         do k = 1,pver
            if ( k>=maxi(i).and.k<=mini(i) )  then
               ndepth_tot = ndepth_tot + 1
               if ( netdt(i,k)>0 )  ndepth_pos = ndepth_pos + 1
            end if
         end do
         ! write (iulog,*) 'WHDEBUG - i: ',i,' Hd: ',hdepth(i),' Hn: ',ndepth_pos,' N: ',ndepth_tot
      end if

     if ((hdepth(i) >= hdepth_min) .and. (abs(lat(i)) < (pi/2._r8))) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = mfcc(NINT(hdepth(i)),NINT(uh(i)),:)

        ! Shift spectrum so that it is relative to the ground.
        shift = -nint(storm_speed(i)/dc)
        tau0 = cshift(tau0,shift)

        ! Adjust magnitude.
        tau0 = tau0*maxq0(i)*maxq0(i)/tau_avg_length

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
