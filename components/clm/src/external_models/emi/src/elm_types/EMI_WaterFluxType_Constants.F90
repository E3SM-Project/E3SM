module EMI_WaterFluxType_Constants
  !
  implicit none
  private
  !
  integer, parameter, public :: L2E_FLUX_INFIL_MASS_FLUX                          = 0801
  integer, parameter, public :: L2E_FLUX_VERTICAL_ET_MASS_FLUX                    = 0802
  integer, parameter, public :: L2E_FLUX_DEW_MASS_FLUX                            = 0803
  integer, parameter, public :: L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX               = 0804
  integer, parameter, public :: L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX          = 0805
  integer, parameter, public :: L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX  = 0806
  integer, parameter, public :: L2E_FLUX_DRAINAGE_MASS_FLUX                       = 0807
  integer, parameter, public :: L2E_FLUX_INFL                                     = 0808
  integer, parameter, public :: L2E_FLUX_TOTDRAIN                                 = 0809
  integer, parameter, public :: L2E_FLUX_GROSS_EVAP_SOIL                          = 0810
  integer, parameter, public :: L2E_FLUX_GROSS_INFL_SOIL                          = 0811
  integer, parameter, public :: L2E_FLUX_SURF                                     = 0812
  integer, parameter, public :: L2E_FLUX_DEW_GRND                                 = 0813
  integer, parameter, public :: L2E_FLUX_DEW_SNOW                                 = 0814
  integer, parameter, public :: L2E_FLUX_SUB_SNOW_VOL                             = 0815
  integer, parameter, public :: L2E_FLUX_SUB_SNOW                                 = 0816
  integer, parameter, public :: L2E_FLUX_H2OSFC2TOPSOI                            = 0817
  integer, parameter, public :: L2E_FLUX_SNOW2TOPSOI                              = 0818
  integer, parameter, public :: L2E_FLUX_ROOTSOI                                  = 0819
  integer, parameter, public :: L2E_FLUX_ADV                                      = 0820
  integer, parameter, public :: L2E_FLUX_DRAIN_VR                                 = 0821
  integer, parameter, public :: L2E_FLUX_TRAN_VEG                                 = 0822
  integer, parameter, public :: L2E_FLUX_ROOTSOI_FRAC                             = 0823

  integer, parameter, public :: E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX          = 0824

end module EMI_WaterFluxType_Constants
