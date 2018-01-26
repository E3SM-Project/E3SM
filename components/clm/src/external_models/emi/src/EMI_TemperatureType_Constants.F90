module EMI_TemperatureType_Constants
  
  implicit none
  private
    
  integer, parameter, public :: L2E_STATE_TSOIL_NLEVGRND                         = 0001
  integer, parameter, public :: L2E_STATE_TSNOW                                  = 0002
  integer, parameter, public :: L2E_STATE_TH2OSFC                                = 0003
  integer, parameter, public :: L2E_STATE_T_SOI10CM                              = 0004
  integer, parameter, public :: L2E_STATE_T_SOIL_NLEVSOI                         = 0005
  integer, parameter, public :: L2E_STATE_T_VEG                                  = 0006

  integer, parameter, public :: E2L_STATE_TSOIL_NLEVGRND                         = 0701
  integer, parameter, public :: E2L_STATE_TSNOW_NLEVSNOW                         = 0702
  integer, parameter, public :: E2L_STATE_TH2OSFC                                = 0703

end module EMI_TemperatureType_Constants
