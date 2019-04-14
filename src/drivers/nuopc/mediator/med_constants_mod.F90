module med_constants_mod

  use shr_kind_mod   , only : R8=>SHR_KIND_R8
  use shr_kind_mod   , only : R4=>SHR_KIND_R4
  use shr_kind_mod   , only : IN=>SHR_KIND_IN
  use shr_kind_mod   , only : I8=>SHR_KIND_I8
  use shr_kind_mod   , only : CL=>SHR_KIND_CL
  use shr_kind_mod   , only : CS=>SHR_KIND_CS
  use shr_kind_mod   , only : CX=>SHR_KIND_CX
  use shr_kind_mod   , only : CXX=>SHR_KIND_CXX

  use shr_cal_mod    , only : med_constants_noleap => shr_cal_noleap
  use shr_cal_mod    , only : med_constants_gregorian => shr_cal_gregorian
  use shr_cal_mod    , only : shr_cal_ymd2date
  use shr_cal_mod    , only : shr_cal_noleap
  use shr_cal_mod    , only : shr_cal_gregorian

  implicit none

  logical,  parameter :: med_constants_statewrite_flag = .false.
  real(R8), parameter :: med_constants_spval_init = 0.0_R8  ! spval for initialization
  real(R8), parameter :: med_constants_spval = 0.0_R8  ! spval
  real(R8), parameter :: med_constants_czero = 0.0_R8  ! spval
  integer,  parameter :: med_constants_ispval_mask = -987987     ! spval for RH mask values
  integer,  parameter :: med_constants_SecPerDay = 86400    ! Seconds per day

  integer :: med_constants_dbug_flag = 5

end module med_constants_mod
