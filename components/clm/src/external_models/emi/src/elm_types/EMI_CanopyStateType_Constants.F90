module EMI_CanopyStateType_Constants
  !
  implicit none
  private
  !
  integer, parameter, public :: L2E_STATE_ALTMAX           = 0601
  integer, parameter, public :: L2E_STATE_ALTMAX_LASTYEAR  = 0602
  integer, parameter, public :: L2E_STATE_LBL_RSC_H2O      = 0603
  integer, parameter, public :: L2E_STATE_ELAI             = 0604

  integer, parameter, public :: E2L_STATE_FSUN             = 0605
  integer, parameter, public :: E2L_STATE_LAISUN           = 0606
  integer, parameter, public :: E2L_STATE_LAISHA           = 0607

end module EMI_CanopyStateType_Constants
