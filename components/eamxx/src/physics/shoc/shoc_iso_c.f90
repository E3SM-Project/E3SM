module shoc_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to shoc fortran.
!

contains

  subroutine shoc_init_c(nlev, gravit, rair, rh2o, cpair, &
                         zvir, latvap, latice, karman, p0) bind(c)
    use shoc, only: shoc_init, npbl

    integer(kind=c_int), value, intent(in) :: nlev ! number of levels

    real(kind=c_real), value, intent(in)  :: gravit ! gravity
    real(kind=c_real), value, intent(in)  :: rair   ! dry air gas constant
    real(kind=c_real), value, intent(in)  :: rh2o   ! water vapor gas constant
    real(kind=c_real), value, intent(in)  :: cpair  ! specific heat of dry air
    real(kind=c_real), value, intent(in)  :: zvir   ! rh2o/rair - 1
    real(kind=c_real), value, intent(in)  :: latvap ! latent heat of vaporization
    real(kind=c_real), value, intent(in)  :: latice ! latent heat of fusion
    real(kind=c_real), value, intent(in)  :: karman ! Von Karman's constant
    real(kind=c_real), value, intent(in)  :: p0     ! Reference pressure

    real(kind=c_real) :: pref_mid(nlev) ! unused values

    pref_mid = 0
    call shoc_init(nlev, gravit, rair, rh2o, cpair, &
                   zvir, latvap, latice, karman, p0, &
                   pref_mid, nlev, 1)
    npbl = nlev ! set pbl layer explicitly so we don't need pref_mid.
  end subroutine shoc_init_c

  ! shoc_init for shoc_main_bfb testing
  subroutine shoc_init_for_main_bfb_c(nlev, gravit, rair, rh2o, cpair, &
                                      zvir, latvap, latice, karman, p0, &
                                      pref_mid, nbot_shoc, ntop_shoc) bind(c)
    use shoc, only: shoc_init

    integer(kind=c_int), value, intent(in) :: nlev ! number of levels
    integer(kind=c_int), value, intent(in) :: nbot_shoc ! Bottom level to which SHOC is applied
    integer(kind=c_int), value, intent(in) :: ntop_shoc ! Top level to which SHOC is applied

    real(kind=c_real), value, intent(in)  :: gravit ! gravity
    real(kind=c_real), value, intent(in)  :: rair   ! dry air gas constant
    real(kind=c_real), value, intent(in)  :: rh2o   ! water vapor gas constant
    real(kind=c_real), value, intent(in)  :: cpair  ! specific heat of dry air
    real(kind=c_real), value, intent(in)  :: zvir   ! rh2o/rair - 1
    real(kind=c_real), value, intent(in)  :: latvap ! latent heat of vaporization
    real(kind=c_real), value, intent(in)  :: latice ! latent heat of fusion
    real(kind=c_real), value, intent(in)  :: karman ! Von Karman's constant
    real(kind=c_real), value, intent(in)  :: p0     ! Reference pressure

    real(kind=c_real), intent(in), dimension(nlev) :: pref_mid ! reference pressures at midpoints
    call shoc_init(nlev, gravit, rair, rh2o, cpair, &
                   zvir, latvap, latice, karman, p0, &
                   pref_mid, nbot_shoc, ntop_shoc)
  end subroutine shoc_init_for_main_bfb_c


end module shoc_iso_c

