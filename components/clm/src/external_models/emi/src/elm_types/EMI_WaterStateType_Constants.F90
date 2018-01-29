module EMI_WaterStateType_Constants
  !
  implicit none
  private
  !
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVGRND            = 0101
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVGRND            = 0102
  integer, parameter, public :: L2E_STATE_VSFM_PROGNOSTIC_SOILP          = 0103
  integer, parameter, public :: L2E_STATE_FRAC_H2OSFC                    = 0104
  integer, parameter, public :: L2E_STATE_FRAC_INUNDATED                 = 0105
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI         = 0106
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI         = 0107
  integer, parameter, public :: L2E_STATE_H2OSOI_VOL_NLEVSOI             = 0108
  integer, parameter, public :: L2E_STATE_AIR_VOL_NLEVSOI                = 0109
  integer, parameter, public :: L2E_STATE_RHO_VAP_NLEVSOI                = 0110
  integer, parameter, public :: L2E_STATE_RHVAP_SOI_NLEVSOI              = 0111
  integer, parameter, public :: L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI  = 0112
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVSOI             = 0113
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVSOI             = 0114
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVSNOW            = 0115
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVSNOW            = 0116
  integer, parameter, public :: L2E_STATE_H2OSNOW                        = 0117
  integer, parameter, public :: L2E_STATE_H2OSFC                         = 0118
  integer, parameter, public :: L2E_STATE_FRAC_SNOW_EFFECTIVE            = 0119

  integer, parameter, public :: E2L_STATE_H2OSOI_LIQ                     = 0120
  integer, parameter, public :: E2L_STATE_H2OSOI_ICE                     = 0121
  integer, parameter, public :: E2L_STATE_VSFM_PROGNOSTIC_SOILP          = 0122

end module EMI_WaterStateType_Constants
