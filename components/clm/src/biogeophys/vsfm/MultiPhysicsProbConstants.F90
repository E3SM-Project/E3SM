#ifdef USE_PETSC_LIB


module MultiPhysicsProbConstants

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Constants for multi-physcis problems
  !-----------------------------------------------------------------------
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  ! mpp_itype
  PetscInt, parameter, public :: MPP_VSFM_SNES_CLM                 = 11

  ! soe_itype
  PetscInt, parameter, public :: SOE_RE_ODE                        = 101

  ! ge_itype
  PetscInt, parameter, public :: GE_RE                             = 201

  ! mesh_itype
  PetscInt, parameter, public :: MESH_CLM_SOIL_COL                 = 301
  PetscInt, parameter, public :: MESH_ALONG_GRAVITY                = 311
  PetscInt, parameter, public :: MESH_AGAINST_GRAVITY              = 312

  ! region_itype
  PetscInt, parameter, public :: SOIL_TOP_CELLS                    = 401
  PetscInt, parameter, public :: SOIL_BOTTOM_CELLS                 = 402
  PetscInt, parameter, public :: SOIL_CELLS                        = 403
  PetscInt, parameter, public :: SOIL_CELLS_OF_NEIGH_MESH          = 404
  PetscInt, parameter, public :: SNOW_TOP_CELLS                    = 405
  PetscInt, parameter, public :: SSW_TOP_CELLS                     = 406
  PetscInt, parameter, public :: ALL_CELLS                         = 407
  PetscInt, parameter, public :: DEFINED_BY_CELL_ID                = 408
  PetscInt, parameter, public :: FACE_TOP                          = 409
  PetscInt, parameter, public :: FACE_BOTTOM                       = 410

  ! condition%itype
  PetscInt, parameter, public :: COND_NULL                         = 500
  PetscInt, parameter, public :: COND_BC                           = 501
  PetscInt, parameter, public :: COND_SS                           = 502
  PetscInt, parameter, public :: COND_MASS_RATE                    = 503
  PetscInt, parameter, public :: COND_MASS_FLUX                    = 504
  PetscInt, parameter, public :: COND_DIRICHLET                    = 505
  PetscInt, parameter, public :: COND_DIRICHLET_FRM_OTR_GOVEQ      = 506
  PetscInt, parameter, public :: COND_HEAT_FLUX                    = 507
  PetscInt, parameter, public :: COND_DARCY_RATE                   = 508

  !
  PetscInt, parameter, public :: VAR_XI                            = 601
  PetscInt, parameter, public :: VAR_DXI_DP                        = 602
  PetscInt, parameter, public :: VAR_DENxSATxPOR                   = 611
  PetscInt, parameter, public :: VAR_DDENxSATxPOR_DP               = 612
  PetscInt, parameter, public :: VAR_DXI_DTIME                     = 603
  PetscInt, parameter, public :: VAR_PRESSURE                      = 604
  PetscInt, parameter, public :: VAR_TEMPERATURE                   = 605
  PetscInt, parameter, public :: VAR_PRESSURE_PREV                 = 606
  PetscInt, parameter, public :: VAR_BC_SS_CONDITION               = 607
  PetscInt, parameter, public :: VAR_LIQ_SAT                       = 608
  PetscInt, parameter, public :: VAR_DENSITY_TYPE                  = 609
  PetscInt, parameter, public :: VAR_MASS                          = 610
  PetscInt, parameter, public :: VAR_SOIL_MATRIX_POT               = 611
  PetscInt, parameter, public :: VAR_FRAC_LIQ_SAT                  = 612

  !
  PetscInt, parameter, public :: AUXVAR_INTERNAL                   = 701
  PetscInt, parameter, public :: AUXVAR_BC                         = 702
  PetscInt, parameter, public :: AUXVAR_SS                         = 703

  PetscInt, parameter, public :: PETSC_TS                          = 801
  PetscInt, parameter, public :: PETSC_SNES                        = 802
  PetscInt, parameter, public :: PETSC_KSP                         = 803

  !
  PetscReal, parameter, public :: PRESSURE_REF                     = 101325.d0     ! [Pa]
  PetscReal, parameter, public :: GRAVITY_CONSTANT                 = 9.8068d0      ! [m s^{-2}]
  PetscReal, parameter, public :: FMWH2O                           = 18.01534d0    ! [kg kmol^{-1}]


end module MultiPhysicsProbConstants

#endif
