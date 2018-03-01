module mpp_varcon

#ifdef USE_PETSC_LIB

  !
  ! !PUBLIC TYPES:
  implicit none
  save

#include <petsc/finclude/petsc.h>

  PetscReal :: grav        = 9.80616d0   ! gravity constant [m/s2]

  PetscReal :: cpliq       = 4.188d3     ! Specific heat of water [J/kg-K]
  PetscReal :: cpice       = 2.11727d3   ! Specific heat of ice [J/kg-K]

  PetscReal :: denh2o      = 1.000d3     ! density of liquid water [kg/m3]
  PetscReal :: denice      = 0.917d3     ! density of ice [kg/m3]

  PetscReal :: tkair       = 0.023d0     ! thermal conductivity of air   [W/m/K]
  PetscReal :: tkice       = 2.290d0     ! thermal conductivity of ice   [W/m/K]
  PetscReal :: tkwat       = 0.57d0      ! thermal conductivity of water [W/m/K]
  PetscReal :: thk_bedrock = 3.0d0       ! thermal conductivity of 'typical' saturated granitic rock 
                                         ! (Clauser and Huenges, 1995)(W/m/K)

  PetscReal :: tfrz        = 273.15d0    ! freezing temperature [K]

  PetscReal :: cnfac       = 0.5d0       ! Crank Nicholson factor between 0 and 1

  PetscReal :: capr        = 0.34d0      ! Tuning factor to turn first layer T into surface T

#endif


  ! Initialize landunit type constants

  integer :: istsoil    !soil         landunit type (natural vegetation)
  integer :: istcrop    !crop         landunit type
  integer :: istice     !land ice     landunit type (glacier)
  integer :: istice_mec !land ice     (multiple elevation classes) landunit type
  integer :: istdlak    !deep lake    landunit type (now used for all lakes)
  integer :: istwet     !wetland      landunit type (swamp, marsh, etc.)

  integer :: max_lunit  !!maximum value that lun%itype can have

  integer, public :: icol_roof
  integer, public :: icol_sunwall
  integer, public :: icol_shadewall
  integer, public :: icol_road_imperv
  integer, public :: icol_road_perv

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public mpp_varcon_init_landunit
  public mpp_varcon_init_column
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine mpp_varcon_init_landunit(istsoil_val, istcrop_val, istice_val, istice_mec_val, &
      istdlak_val, istwet_val, max_lunit_val)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent (in) :: istsoil_val
    integer, intent (in) :: istcrop_val
    integer, intent (in) :: istice_val
    integer, intent (in) :: istice_mec_val
    integer, intent (in) :: istdlak_val
    integer, intent (in) :: istwet_val
    integer, intent (in) :: max_lunit_val

    istsoil    = istsoil_val
    istcrop    = istcrop_val
    istice     = istice_val
    istice_mec = istice_mec_val
    istdlak    = istdlak_val
    istwet     = istwet_val
    max_lunit  = max_lunit_val

  end subroutine mpp_varcon_init_landunit

  !------------------------------------------------------------------------------
  subroutine mpp_varcon_init_column(icol_roof_val, icol_sunwall_val, icol_shadewall_val, &
      icol_road_imperv_val, icol_road_perv_val)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent (in) :: icol_roof_val
    integer, intent (in) :: icol_sunwall_val
    integer, intent (in) :: icol_shadewall_val
    integer, intent (in) :: icol_road_imperv_val
    integer, intent (in) :: icol_road_perv_val

    icol_roof        = icol_roof_val
    icol_sunwall     = icol_sunwall_val
    icol_shadewall   = icol_shadewall_val
    icol_road_imperv = icol_road_imperv_val
    icol_road_perv   = icol_road_perv_val

  end subroutine mpp_varcon_init_column

end module mpp_varcon
