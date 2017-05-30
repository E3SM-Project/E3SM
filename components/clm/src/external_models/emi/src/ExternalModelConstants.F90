module ExternalModelConstants

  implicit none
  private

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ID for various external models
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer, parameter, public :: EM_INITIALIZATION_STAGE                          = 000

  integer, public, parameter :: EM_ID_BETR                                       = 001
  integer, parameter, public :: EM_BETR_BEGIN_MASS_BALANCE_STAGE                 = 002
  integer, parameter, public :: EM_BETR_PRE_DIAG_WATER_FLUX_STAGE                = 003

  integer, public, parameter :: EM_ID_FATES                                      = 101
  integer, parameter, public :: EM_FATES_SUNFRAC_STAGE                           = 102

  integer, public, parameter :: EM_ID_PFLOTRAN                                   = 200

  integer, public, parameter :: EM_ID_VSFM                                       = 300
  integer, parameter, public :: EM_VSFM_SOIL_HYDRO_STAGE                         = 301

  integer, public, parameter :: EM_ID_PTM                                        = 400
  integer, parameter, public :: EM_PTM_TBASED_SOLVE_STAGE                        = 401

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for state variables sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! temeprature_type
  integer, parameter, public :: L2E_STATE_TSOIL_NLEVGRND                         = 0001
  integer, parameter, public :: L2E_STATE_TSNOW                                  = 0002
  integer, parameter, public :: L2E_STATE_TH2OSFC                                = 0003

  ! waterstate_type
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVGRND                    = 0101
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVGRND                    = 0102
  integer, parameter, public :: L2E_STATE_VSFM_PROGNOSTIC_SOILP                  = 0103
  integer, parameter, public :: L2E_STATE_FRAC_H2OSFC                            = 0104
  integer, parameter, public :: L2E_STATE_FRAC_INUNDATED                         = 0105
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI                 = 0106
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI                 = 0107
  integer, parameter, public :: L2E_STATE_H2OSOI_VOL_NLEVSOI                     = 0108
  integer, parameter, public :: L2E_STATE_AIR_VOL_NLEVSOI                        = 0109
  integer, parameter, public :: L2E_STATE_RHO_VAP_NLEVSOI                        = 0110
  integer, parameter, public :: L2E_STATE_RHVAP_SOI_NLEVSOI                      = 0111
  integer, parameter, public :: L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI          = 0112
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVSOI                     = 0113
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVSOI                     = 0114
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVSNOW                    = 0115
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVSNOW                    = 0116
  integer, parameter, public :: L2E_STATE_H2OSNOW                                = 0117
  integer, parameter, public :: L2E_STATE_H2OSFC                                 = 0118
  integer, parameter, public :: L2E_STATE_FRAC_SNOW_EFFECTIVE                    = 0119

  ! soilhydrology_type
  integer, parameter, public :: L2E_STATE_WTD                                    = 0201

  ! IDs for states sent from External Model to ALM

  ! waterstate_type
  integer, parameter, public :: E2L_STATE_H2OSOI_LIQ                             = 0301
  integer, parameter, public :: E2L_STATE_H2OSOI_ICE                             = 0302
  integer, parameter, public :: E2L_STATE_VSFM_PROGNOSTIC_SOILP                  = 0303

  ! soilhydrology_type
  integer, parameter, public :: E2L_STATE_WTD                                    = 0401

  ! soilstate_type
  integer, parameter, public :: E2L_STATE_SOIL_MATRIC_POTENTIAL                  = 0501

  ! canopystate_type
  integer, parameter, public :: E2L_STATE_FSUN                                   = 0601
  integer, parameter, public :: E2L_STATE_LAISUN                                 = 0602
  integer, parameter, public :: E2L_STATE_LAISHA                                 = 0603

  ! temperature_type
  integer, parameter, public :: E2L_STATE_TSOIL_NLEVGRND                         = 0701
  integer, parameter, public :: E2L_STATE_TSNOW_NLEVSNOW                         = 0702
  integer, parameter, public :: E2L_STATE_TH2OSFC                                = 0703

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for fluxes sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! waterflux_type
  integer, parameter, public :: L2E_FLUX_INFIL_MASS_FLUX                         = 0801
  integer, parameter, public :: L2E_FLUX_VERTICAL_ET_MASS_FLUX                   = 0802
  integer, parameter, public :: L2E_FLUX_DEW_MASS_FLUX                           = 0803
  integer, parameter, public :: L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX              = 0804
  integer, parameter, public :: L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX         = 0805
  integer, parameter, public :: L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX = 0806
  integer, parameter, public :: L2E_FLUX_DRAINAGE_MASS_FLUX                      = 0807

  ! atm2lnd_type
  integer, parameter, public :: L2E_FLUX_SOLAR_DIRECT_RADDIATION                 = 0901
  integer, parameter, public :: L2E_FLUX_SOLAR_DIFFUSE_RADDIATION                = 0902

  ! energyflux_type
  integer, parameter, public :: L2E_FLUX_ABSORBED_SOLAR_RADIATION                = 1001
  integer, parameter, public :: L2E_FLUX_SOIL_HEAT_FLUX                          = 1002
  integer, parameter, public :: L2E_FLUX_SNOW_HEAT_FLUX                          = 1003
  integer, parameter, public :: L2E_FLUX_H2OSFC_HEAT_FLUX                        = 1004
  integer, parameter, public :: L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX                 = 1005

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for fluxes sent from External Model to ALM
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! waterflux_type
  integer, parameter, public :: E2L_FLUX_AQUIFER_RECHARGE                        = 1101
  integer, parameter, public :: E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX         = 1102

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for filter variables sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: L2E_FILTER_HYDROLOGYC                            = 1201
  integer, parameter, public :: L2E_FILTER_NUM_HYDROLOGYC                        = 1202
  integer, parameter, public :: L2E_FILTER_NOLAKEC                               = 1203
  integer, parameter, public :: L2E_FILTER_NUM_NOLAKEC                           = 1204
  integer, parameter, public :: L2E_FILTER_NOLAKEC_AND_NOURBANC                  = 1205
  integer, parameter, public :: L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC              = 1206

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for column-level attributes sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: L2E_COLUMN_ACTIVE                                = 1301
  integer, parameter, public :: L2E_COLUMN_TYPE                                  = 1302
  integer, parameter, public :: L2E_COLUMN_LANDUNIT_INDEX                        = 1303
  integer, parameter, public :: L2E_COLUMN_ZI                                    = 1304
  integer, parameter, public :: L2E_COLUMN_DZ                                    = 1305
  integer, parameter, public :: L2E_COLUMN_Z                                     = 1306
  integer, parameter, public :: L2E_COLUMN_AREA                                  = 1307
  integer, parameter, public :: L2E_COLUMN_GRIDCELL_INDEX                        = 1308
  integer, parameter, public :: L2E_COLUMN_PATCH_INDEX                           = 1309
  integer, parameter, public :: L2E_COLUMN_NUM_SNOW_LAYERS                       = 1310
  integer, parameter, public :: L2E_COLUMN_ZI_SNOW_AND_SOIL                      = 1311
  integer, parameter, public :: L2E_COLUMN_DZ_SNOW_AND_SOIL                      = 1312
  integer, parameter, public :: L2E_COLUMN_Z_SNOW_AND_SOIL                       = 1313

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for landunit-level attributes sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: L2E_LANDUNIT_TYPE                                = 1401
  integer, parameter, public :: L2E_LANDUNIT_LAKEPOINT                           = 1402
  integer, parameter, public :: L2E_LANDUNIT_URBANPOINT                          = 1403

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for parameters sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! soilstate_type
  integer, parameter, public :: L2E_PARAMETER_WATSATC                            = 1501
  integer, parameter, public :: L2E_PARAMETER_HKSATC                             = 1502
  integer, parameter, public :: L2E_PARAMETER_BSWC                               = 1503
  integer, parameter, public :: L2E_PARAMETER_SUCSATC                            = 1504
  integer, parameter, public :: L2E_PARAMETER_EFFPOROSITYC                       = 1505
  integer, parameter, public :: L2E_PARAMETER_CSOL                               = 1506
  integer, parameter, public :: L2E_PARAMETER_TKMG                               = 1507
  integer, parameter, public :: L2E_PARAMETER_TKDRY                              = 1508

end module ExternalModelConstants
