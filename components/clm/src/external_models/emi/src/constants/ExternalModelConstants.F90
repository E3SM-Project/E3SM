module ExternalModelConstants
  !
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

  integer, public, parameter :: EM_ID_STUB                                       = 500
  integer, parameter, public :: EM_STUB_SOIL_HYDRO_STAGE                         = 501
  integer, parameter, public :: EM_STUB_SOIL_THERMAL_STAGE                       = 502

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for variable sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: L2E_VAR_MAX_PATCH_PER_COL                        = 1601

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! IDs for patch-level attributes sent from ALM to External Model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: L2E_PATCH_ACTIVE                                 = 1701
  integer, parameter, public :: L2E_PATCH_TYPE                                   = 1702
  integer, parameter, public :: L2E_PATCH_WT_COL                                 = 1703

end module ExternalModelConstants
