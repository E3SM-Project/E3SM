module cam_control_mod
!----------------------------------------------------------------------- 
! 
! Purpose: Model control variables
! replaces comctl.h
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8=>shr_kind_r8
  public

  integer :: nsrest            ! run type flag (0=initial, 1=restart-continuation, 3=restart-branch)

  logical :: ideal_phys        ! true => run "idealized" model configuration
  logical :: adiabatic         ! true => no physics
  logical :: moist_physics     ! true => moist physics enabled, i.e. ((.not. ideal_phys) .and. (.not. adiabatic))
  logical :: aqua_planet       ! Flag to run model in "aqua planet" mode

  logical :: print_step_cost   ! true => print per-timestep cost info

  logical :: indirect          ! True => include indirect radiative effects of sulfate aerosols

  logical :: use_64bit_nc = .true.   ! true => use new 64-bit netCDF format for cam history files

! Earth's orbital characteristics
!	
! Orbital information after processed by orbit_params
!
      real(r8) :: eccen       ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
      real(r8) :: obliqr      ! Earth's obliquity in radians
      real(r8) :: lambm0      ! Mean longitude of perihelion at the 
                              ! vernal equinox (radians)
      real(r8) :: mvelpp      ! Earth's moving vernal equinox longitude
!                             ! of perihelion plus pi (radians)
!
!------------------------------------------------------------

! from perturb.h
      real(r8) :: pertlim = 0.0_r8

! from comadj.h
      integer :: nlvdry = 3

! from comtsc.h


      real(r8) :: latice     ! Latent heat of fusion
      real(r8) :: tmelt      ! Melting temperature of snow and ice
      real(r8) :: latvap     ! Latent heat of vaporization
      real(r8) :: rair       ! Gas constant for dry air
      real(r8) :: stebol     ! Stefan-Boltzmann constant
      real(r8) :: snwedp     ! Snow equivalent depth factor

      integer :: magfield_fix_year = 1995

end module cam_control_mod
