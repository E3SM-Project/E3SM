module med_constants_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use shr_kind_mod, only : R8=>shr_kind_r8

  implicit none

  private

  integer,          public, parameter :: med_constants_dbug_flag = 0
  logical,          public, parameter :: med_constants_statewrite_flag = .false.
  real(R8),         public, parameter :: med_constants_spval_init = 0.0_R8  ! spval for initialization
  real(R8),         public, parameter :: med_constants_spval = 0.0_R8  ! spval
  real(R8),         public, parameter :: med_constants_czero = 0.0_R8  ! spval
  integer,          public, parameter :: med_constants_ispval_mask = -987987     ! spval for RH mask values

  !-----------------------------------------------------------------------------

end module med_constants_mod
