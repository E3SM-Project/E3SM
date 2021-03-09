module EMI_SoilHydrologyType_Constants
  !
  implicit none
  private
  !
  integer, parameter, public :: L2E_STATE_WTD              = 0201
  integer, parameter, public :: L2E_STATE_QCHARGE          = 0202
  integer, parameter, public :: L2E_STATE_FRACICE          = 0203

  integer, parameter, public :: E2L_STATE_WTD              = 0204

  integer, parameter, public :: E2L_FLUX_AQUIFER_RECHARGE  = 0205

end module EMI_SoilHydrologyType_Constants
