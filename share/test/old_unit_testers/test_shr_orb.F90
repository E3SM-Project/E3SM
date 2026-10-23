 program test_shr_orb
!
! Simple unit-test program for the shr_orb_mod module.
!
! Erik Kluzek
!
! $Id: test_shr_orb.F90 7482 2007-11-07 20:54:58Z erik $
!
 use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN
 use shr_orb_mod, only: shr_orb_cosz, shr_orb_params, shr_orb_decl, shr_orb_print
 implicit none
 integer, parameter :: nyears = 5
 integer, parameter :: ndays = 5
 real   (SHR_KIND_R8), parameter :: jday(ndays) = &
           (/ 0.0_SHR_KIND_R8, 0.25_SHR_KIND_R8, 0.5_SHR_KIND_R8, 180.0_SHR_KIND_R8, 365.0_SHR_KIND_R8 /)   ! Julian cal day (1.xx to 365.xx)
 real   (SHR_KIND_R8) :: lat = 42.0_SHR_KIND_R8 ! Centered latitude (radians)
 real   (SHR_KIND_R8) :: lon =  0.0_SHR_KIND_R8 ! Centered longitude (radians)
 real   (SHR_KIND_R8) :: declin     ! Solar declination (radians)
 real   (SHR_KIND_R8) :: eccen      ! orbital eccentricity
 real   (SHR_KIND_R8) :: obliq      ! obliquity in degrees
 real   (SHR_KIND_R8) :: mvelp      ! moving vernal equinox long
 integer(SHR_KIND_IN), parameter :: iyear_AD(nyears) = &
                    (/-900000, -1650, 1950, 3600, 1000000/)
 logical              :: log_print = .true. ! Flags print of status/error
 real   (SHR_KIND_R8) :: obliqr     ! Earths obliquity in rad
 real   (SHR_KIND_R8) :: lambm0     ! Mean long of perihelion at
                                    ! vernal equinox (radians)
 real   (SHR_KIND_R8) :: mvelpp     ! moving vernal equinox long
                                    ! of perihelion plus pi (rad)
 real   (SHR_KIND_R8) :: cosz       ! cosine of solar zenith angle
 real   (SHR_KIND_R8) :: eccf       ! Earth-sun distance factor
 integer i, j   ! Indices

 print *, 'Test orbit calculation for ', nyears, ' years and ', ndays, ' days '
 do i = 1, nyears
    call shr_orb_params( iyear_AD(i) , eccen  , obliq , mvelp     ,     &
              &               obliqr   , lambm0 , mvelpp, log_print )
    call shr_orb_print( iyear_AD(i), eccen, obliq, mvelp )
    do j = 1, ndays
      call shr_orb_decl(jday(j),eccen ,mvelpp ,lambm0 ,obliqr ,declin,eccf)
      cosz = shr_orb_cosz(jday(j),lat,lon,declin)
      print *, 'jday = ', jday(j), ' declin = ', declin, ' cosz = ', cosz
    end do
 end do
 call test_equation_of_time()
 print *, 'PASS'

 contains

 subroutine test_equation_of_time()
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use shr_kind_mod, only: SHR_KIND_I8
   use shr_const_mod, only: pi => SHR_CONST_PI
   use shr_orb_mod, only: shr_orb_azimuth, shr_orb_decl_gregorian
   real(SHR_KIND_R8) :: day, eot, delta, factor, reference_delta, reference_factor
   real(SHR_KIND_R8) :: actual, expected, phase, previous_eot
  real(SHR_KIND_R8) :: latitude
   integer :: sample
  integer :: latitude_index, time_index, season

   previous_eot = 0.0_SHR_KIND_R8
   do sample = 1, 365
   day = real(sample, SHR_KIND_R8)
   call shr_orb_decl(day, 0.0_SHR_KIND_R8, 0.0_SHR_KIND_R8, 0.0_SHR_KIND_R8, &
     0.0_SHR_KIND_R8, delta, factor, eot_seconds=eot)
   if (abs(eot) > 1.e-8_SHR_KIND_R8) stop 1

   call shr_orb_decl(day, 0.0167_SHR_KIND_R8, 1.8_SHR_KIND_R8, 0.0_SHR_KIND_R8, &
     0.409_SHR_KIND_R8, reference_delta, reference_factor)
   call shr_orb_decl(day, 0.0167_SHR_KIND_R8, 1.8_SHR_KIND_R8, 0.0_SHR_KIND_R8, &
     0.409_SHR_KIND_R8, delta, factor, eot_seconds=eot)
   if (transfer(delta, 0_SHR_KIND_I8) /= transfer(reference_delta, 0_SHR_KIND_I8)) stop 2
   if (transfer(factor, 0_SHR_KIND_I8) /= transfer(reference_factor, 0_SHR_KIND_I8)) stop 3
   if (abs(eot) > 1200.0_SHR_KIND_R8) stop 4
   if (sample > 1 .and. abs(eot - previous_eot) > 60.0_SHR_KIND_R8) stop 5
   previous_eot = eot

   actual = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, 0.2_SHR_KIND_R8, delta)
   expected = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, 0.2_SHR_KIND_R8, &
     delta, eot_seconds=0.0_SHR_KIND_R8)
   if (transfer(actual, 0_SHR_KIND_I8) /= transfer(expected, 0_SHR_KIND_I8)) stop 6
   actual = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, 0.2_SHR_KIND_R8, &
     delta, 1800.0_SHR_KIND_R8)
   expected = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, 0.2_SHR_KIND_R8, &
     delta, 1800.0_SHR_KIND_R8, eot_seconds=0.0_SHR_KIND_R8)
   if (transfer(actual, 0_SHR_KIND_I8) /= transfer(expected, 0_SHR_KIND_I8)) stop 7

   phase = 2.0_SHR_KIND_R8*pi*eot/86400.0_SHR_KIND_R8
   actual = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, 0.2_SHR_KIND_R8, &
     delta, 1800.0_SHR_KIND_R8, eot_seconds=eot)
   expected = shr_orb_cosz(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, &
     0.2_SHR_KIND_R8 + phase, delta, 1800.0_SHR_KIND_R8)
   if (abs(actual - expected) > 1.e-12_SHR_KIND_R8) stop 8
   actual = shr_orb_azimuth(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, &
     0.2_SHR_KIND_R8, delta, 1.0_SHR_KIND_R8, eot_seconds=eot)
   expected = shr_orb_azimuth(day + 0.25_SHR_KIND_R8, 0.4_SHR_KIND_R8, &
     0.2_SHR_KIND_R8 + phase, delta, 1.0_SHR_KIND_R8)
   if (abs(actual - expected) > 1.e-12_SHR_KIND_R8) stop 9
   end do

   do sample = -1, 1, 2
   call shr_orb_decl(80.5_SHR_KIND_R8, 0.0167_SHR_KIND_R8, &
     real(sample, SHR_KIND_R8)*pi/2, 0.0_SHR_KIND_R8, 0.0_SHR_KIND_R8, &
     delta, factor, eot_seconds=eot)
   if (real(sample, SHR_KIND_R8)*eot <= 0.0_SHR_KIND_R8) stop 10
   actual = shr_orb_cosz(80.5_SHR_KIND_R8 - eot/86400.0_SHR_KIND_R8, &
     0.0_SHR_KIND_R8, 0.0_SHR_KIND_R8, 0.0_SHR_KIND_R8, eot_seconds=eot)
   if (abs(actual - 1.0_SHR_KIND_R8) > 1.e-12_SHR_KIND_R8) stop 11
   end do

   call shr_orb_decl_gregorian(-0.5_SHR_KIND_R8, 0.0167_SHR_KIND_R8, 1.8_SHR_KIND_R8, &
     0.0_SHR_KIND_R8, 0.409_SHR_KIND_R8, reference_delta, reference_factor)
   call shr_orb_decl_gregorian(-0.5_SHR_KIND_R8, 0.0167_SHR_KIND_R8, 1.8_SHR_KIND_R8, &
     0.0_SHR_KIND_R8, 0.409_SHR_KIND_R8, delta, factor, eot_seconds=eot)
   if (transfer(delta, 0_SHR_KIND_I8) /= transfer(reference_delta, 0_SHR_KIND_I8)) stop 12
   if (transfer(factor, 0_SHR_KIND_I8) /= transfer(reference_factor, 0_SHR_KIND_I8)) stop 13
  do sample = -1, 1, 2
    eot = real(sample, SHR_KIND_R8)*1200.0_SHR_KIND_R8
    phase = 2.0_SHR_KIND_R8*pi*eot/86400.0_SHR_KIND_R8
    do season = -1, 1, 2
      delta = real(season, SHR_KIND_R8)*0.409_SHR_KIND_R8
      do latitude_index = -1, 1
        latitude = real(latitude_index, SHR_KIND_R8)*pi/2.0_SHR_KIND_R8
        do time_index = 0, 8
          day = 80.0_SHR_KIND_R8 + real(time_index, SHR_KIND_R8)/8.0_SHR_KIND_R8
          actual = shr_orb_cosz(day,latitude,0.2_SHR_KIND_R8,delta, &
              1800.0_SHR_KIND_R8,eot_seconds=eot)
          expected = shr_orb_cosz(day,latitude,0.2_SHR_KIND_R8+phase,delta,1800.0_SHR_KIND_R8)
          if (.not. ieee_is_finite(actual) .or. abs(actual-expected) > 1.e-12_SHR_KIND_R8) stop 14
          actual = shr_orb_cosz(day,latitude,0.2_SHR_KIND_R8,delta,eot_seconds=eot)
          expected = shr_orb_cosz(day,latitude,0.2_SHR_KIND_R8+phase,delta)
          if (.not. ieee_is_finite(actual) .or. abs(actual-expected) > 1.e-12_SHR_KIND_R8) stop 15
        end do
      end do
    end do
  end do
 end subroutine test_equation_of_time

 end program test_shr_orb
