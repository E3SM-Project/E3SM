module EMI_WaterFluxType_Constants

  implicit none
  private

  integer, parameter, public :: L2E_FLUX_INFIL_MASS_FLUX                         = 0801
  integer, parameter, public :: L2E_FLUX_VERTICAL_ET_MASS_FLUX                   = 0802
  integer, parameter, public :: L2E_FLUX_DEW_MASS_FLUX                           = 0803
  integer, parameter, public :: L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX              = 0804
  integer, parameter, public :: L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX         = 0805
  integer, parameter, public :: L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX = 0806
  integer, parameter, public :: L2E_FLUX_DRAINAGE_MASS_FLUX                      = 0807

  integer, parameter, public :: E2L_FLUX_AQUIFER_RECHARGE                        = 1101
  integer, parameter, public :: E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX         = 1102

end module EMI_WaterFluxType_Constants
