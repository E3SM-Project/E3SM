module med_constants_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use ESMF

  implicit none

  public

  integer           , parameter :: med_constants_dbug_flag = 6
 !logical           , parameter :: med_constants_statewrite_flag = .true.
  logical           , parameter :: med_constants_statewrite_flag = .false.
  real(ESMF_KIND_R8), parameter :: med_constants_spval_init = 0.0_ESMF_KIND_R8  ! spval for initialization
  real(ESMF_KIND_R8), parameter :: med_constants_spval = 0.0_ESMF_KIND_R8  ! spval
  real(ESMF_KIND_R8), parameter :: med_constants_czero = 0.0_ESMF_KIND_R8  ! spval
  integer           , parameter :: med_constants_ispval_mask = -987987     ! spval for RH mask values
  character(len=*)  , parameter :: med_constants_spval_rhfile = 'SkiP_FilE'

  !-----------------------------------------------------------------------------

end module med_constants_mod
