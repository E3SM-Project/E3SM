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
      real(r8) :: pertlim     = 0.0_r8
      logical  :: new_random  = .false.
      integer  :: seed_custom = 0
      logical  :: seed_clock  = .false.

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

contains

subroutine cam_ctrl_set_physics_type(phys_package)
  use cam_abortutils, only : endrun
  use spmd_utils,     only: masterproc

  ! Dummy argument
  character(len=*), intent(in) :: phys_package
  ! Local variable
  character(len=*), parameter :: subname = 'cam_ctrl_set_physics_type'

  adiabatic = trim(phys_package) == 'adiabatic'
  ideal_phys = trim(phys_package) == 'held_suarez'
  moist_physics = .not. (adiabatic .or. ideal_phys)
  if (adiabatic .and. ideal_phys) then
    call endrun (subname//': FATAL: Only one of ADIABATIC or HELD_SUAREZ can be .true.')
  end if

  if ((.not. moist_physics) .and. aqua_planet) then
    call endrun (subname//': FATAL: AQUA_PLANET not compatible with dry physics package, ('//trim(phys_package)//')')
  end if

  if (masterproc) then
    if (adiabatic) then
      write(iulog,*) 'Run model ADIABATICALLY (i.e. no physics)'
    end if
    if (ideal_phys) then
      write(iulog,*) 'Run model with Held-Suarez physics forcing'
    end if
  end if

end subroutine cam_ctrl_set_physics_type

end module cam_control_mod
