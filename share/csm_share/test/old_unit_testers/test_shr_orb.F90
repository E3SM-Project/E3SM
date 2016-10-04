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
 print *, 'PASS'

 end program test_shr_orb
