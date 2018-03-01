 module mpp_varctl
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6

  !----------------------------------------------------------
  ! VSFM switches
  !----------------------------------------------------------
  logical          , public :: use_vsfm                    = .false.
  logical          , public :: vsfm_use_dynamic_linesearch = .false.
  logical          , public :: vsfm_include_seepage_bc     = .false.
  logical          , public :: lateral_connectivity        = .false.
  logical          , public :: restart_vsfm                = .false.
  character(len=32), public :: vsfm_satfunc_type           = 'smooth_brooks_corey_bz3'
  character(len=32), public :: vsfm_lateral_model_type     = 'none'

  !----------------------------------------------------------
  ! PETSc-based thermal model switches
  !----------------------------------------------------------
  logical, public :: use_petsc_thermal_model = .false.

  ! !PUBLIC MEMBER FUNCTIONS:
  public mpp_varctl_init_vsfm
  public mpp_varctl_init_petsc_thermal
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine mpp_varctl_init_vsfm(use_vsfm_val, vsfm_use_dynamic_linesearch_val, &
      vsfm_include_seepage_bc_val, lateral_connectivity_val, restart_vsfm_val, &
      vsfm_satfunc_type_val, vsfm_lateral_model_type_val)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    logical         , intent (in) :: use_vsfm_val
    logical         , intent (in) :: vsfm_use_dynamic_linesearch_val
    logical         , intent (in) :: vsfm_include_seepage_bc_val
    logical         , intent (in) :: lateral_connectivity_val
    logical         , intent (in) :: restart_vsfm_val
    character(len=*), intent (in) :: vsfm_satfunc_type_val
    character(len=*), intent (in) :: vsfm_lateral_model_type_val

    use_vsfm                    = use_vsfm_val
    vsfm_use_dynamic_linesearch = vsfm_use_dynamic_linesearch_val
    vsfm_include_seepage_bc     = vsfm_include_seepage_bc_val
    lateral_connectivity        = lateral_connectivity_val
    restart_vsfm                = restart_vsfm_val
    vsfm_satfunc_type           = trim(vsfm_satfunc_type_val)
    vsfm_lateral_model_type     = trim(vsfm_lateral_model_type_val)

  end subroutine mpp_varctl_init_vsfm

  !------------------------------------------------------------------------------
  subroutine mpp_varctl_init_petsc_thermal(use_petsc_thermal_model_val)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    logical, intent (in) :: use_petsc_thermal_model_val

    use_petsc_thermal_model = use_petsc_thermal_model_val

  end subroutine mpp_varctl_init_petsc_thermal


end module mpp_varctl
