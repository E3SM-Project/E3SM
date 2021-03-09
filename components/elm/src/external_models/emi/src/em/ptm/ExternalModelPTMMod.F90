module ExternalModelPTMMod

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod, only : emi_data_list, emi_data
  use decompMod                    , only : bounds_type
  use mpp_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use MultiPhysicsProbThermal      , only : mpp_thermal_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  !
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_ptm_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsat
     integer :: index_l2e_init_parameter_csol
     integer :: index_l2e_init_parameter_tkmg
     integer :: index_l2e_init_parameter_tkdry

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------
     integer :: index_l2e_col_num_snow_lyrs
     integer :: index_l2e_col_zi
     integer :: index_l2e_col_dz
     integer :: index_l2e_col_z
     integer :: index_l2e_col_active
     integer :: index_l2e_col_landunit_index

     integer :: index_l2e_landunit_lakepoint
     integer :: index_l2e_landunit_urbanpoint

     integer :: index_l2e_filter_nolakec_and_nourbanc
     integer :: index_l2e_filter_num_nolakec_and_nourbanc

     integer :: index_l2e_state_frac_snow_eff
     integer :: index_l2e_state_frac_h2osfc
     integer :: index_l2e_state_h2osno
     integer :: index_l2e_state_h2osfc
     integer :: index_l2e_state_h2osoi_liq_nlevgrnd
     integer :: index_l2e_state_h2osoi_ice_nlevgrnd
     integer :: index_l2e_state_h2osoi_liq_nlevsnow
     integer :: index_l2e_state_h2osoi_ice_nlevsnow

     integer :: index_l2e_state_temperature_soil_nlevgrnd
     integer :: index_l2e_state_temperature_snow
     integer :: index_l2e_state_temperature_h2osfc

     integer :: index_l2e_flux_hs_soil
     integer :: index_l2e_flux_hs_top_snow
     integer :: index_l2e_flux_hs_h2osfc
     integer :: index_l2e_flux_dhsdT
     integer :: index_l2e_flux_sabg_lyr

     integer :: index_e2l_state_temperature_soil_nlevgrnd
     integer :: index_e2l_state_temperature_snow
     integer :: index_e2l_state_temperature_h2osfc

     type(mpp_thermal_type) :: thermal_mpp

   contains
     procedure, public :: Populate_L2E_Init_List  => EM_PTM_Populate_L2E_Init_List
     procedure, public :: Populate_L2E_List       => EM_PTM_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_PTM_Populate_E2L_List
     procedure, public :: Init                    => EM_PTM_Init
     procedure, public :: Solve                   => EM_PTM_Solve
  end type em_ptm_type

contains

  !------------------------------------------------------------------------
  subroutine EM_PTM_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PTM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ptm_type)                  :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                      = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active          = index

    id                                      = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index  = index

    id                                      = L2E_COLUMN_ZI_SNOW_AND_SOIL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi              = index

    id                                      = L2E_COLUMN_DZ_SNOW_AND_SOIL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz              = index

    id                                      = L2E_COLUMN_Z_SNOW_AND_SOIL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z               = index

    id                                      = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type       = index

    id                                      = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint  = index

    id                                      = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint = index

    id                                      = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsat    = index

    id                                      = L2E_PARAMETER_CSOL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_csol      = index

    id                                      = L2E_PARAMETER_TKMG
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_tkmg      = index

    id                                      = L2E_PARAMETER_TKDRY
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_tkdry     = index

    deallocate(em_stages)

  end subroutine EM_PTM_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PTM_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PETSc-based Thermal Model (PTM) from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ptm_type)                  :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages                               = 1
    allocate(em_stages(number_em_stages))
    em_stages(1)                                   = EM_PTM_TBASED_SOLVE_STAGE

    id                                             = L2E_COLUMN_NUM_SNOW_LAYERS
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_num_snow_lyrs               = index

    id                                             = L2E_COLUMN_ZI_SNOW_AND_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_zi                          = index

    id                                             = L2E_COLUMN_DZ_SNOW_AND_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz                          = index

    id                                             = L2E_COLUMN_Z_SNOW_AND_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_z                           = index

    id                                             = L2E_COLUMN_ACTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_active                      = index

    id                                             = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_landunit_index              = index

    id                                             = L2E_LANDUNIT_LAKEPOINT
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_landunit_lakepoint              = index

    id                                             = L2E_LANDUNIT_URBANPOINT
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_landunit_urbanpoint             = index

    id                                             = L2E_FILTER_NOLAKEC_AND_NOURBANC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_nolakec_and_nourbanc     = index

    id                                             = L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_nolakec_and_nourbanc = index

    id                                             = L2E_STATE_FRAC_SNOW_EFFECTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_frac_snow_eff             = index

    id                                             = L2E_STATE_FRAC_H2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_frac_h2osfc               = index

    id                                             = L2E_STATE_H2OSNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osno                    = index

    id                                             = L2E_STATE_H2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osfc                    = index

    id                                             = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq_nlevgrnd       = index

    id                                             = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice_nlevgrnd       = index

    id                                             = L2E_STATE_H2OSOI_LIQ_NLEVSNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq_nlevsnow       = index

    id                                             = L2E_STATE_H2OSOI_ICE_NLEVSNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice_nlevsnow       = index

    id                                             = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_temperature_soil_nlevgrnd = index

    id                                             = L2E_STATE_TSNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_temperature_snow          = index

    id                                             = L2E_STATE_TH2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_temperature_h2osfc        = index

    id                                             = L2E_FLUX_ABSORBED_SOLAR_RADIATION
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_sabg_lyr                   = index

    id                                             = L2E_FLUX_SOIL_HEAT_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_hs_soil                    = index

    id                                             = L2E_FLUX_SNOW_HEAT_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_hs_top_snow                = index

    id                                             = L2E_FLUX_H2OSFC_HEAT_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_hs_h2osfc                  = index

    id                                             = L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dhsdT                      = index

    deallocate(em_stages)

  end subroutine EM_PTM_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_PTM_Populate_E2L_List(this, e2l_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by PTM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ptm_type)                   :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PTM_TBASED_SOLVE_STAGE

    id                                             = E2L_STATE_TSOIL_NLEVGRND
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_temperature_soil_nlevgrnd = index

    id                                             = E2L_STATE_TSNOW_NLEVSNOW
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_temperature_snow          = index

    id                                             = E2L_STATE_TH2OSFC
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_temperature_h2osfc        = index

    deallocate(em_stages)

  end subroutine EM_PTM_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_PTM_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use mpp_varctl                , only : use_petsc_thermal_model
    !
    implicit none
    ! !ARGUMENTS
    class(em_ptm_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent (in)   :: bounds_clump

    if (.not.use_petsc_thermal_model) return

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(this, iam)

    ! 2. Add all meshes needed for the MPP
    call add_meshes(this, l2e_init_list, bounds_clump)

    ! 3. Add all governing equations
    call add_goveqns(this)

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(this, l2e_init_list, bounds_clump)

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars(this)

    ! 6. Setup the MPP
    call this%thermal_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call add_material_properties(this, l2e_init_list, bounds_clump)

  end subroutine EM_PTM_Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(this, iam)
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    !
    implicit none
    ! !ARGUMENTS
    class(em_ptm_type)              :: this
    integer            , intent(in) :: iam

    !
    ! Set up the multi-physics problem
    !
    call this%thermal_mpp%Init       ()
    call this%thermal_mpp%SetName    ('Snow + Standing water + Soil thermal model using temperature')
    call this%thermal_mpp%SetID      (MPP_THERMAL_TBASED_KSP_CLM)
    call this%thermal_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Add meshes to the thermal MPP problem
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use filterMod                 , only : filter
    use mpp_varpar                , only : nlevgrnd, nlevsno
    use shr_infnan_mod            , only : shr_infnan_isnan
    use abortutils                , only : endrun
    use shr_infnan_mod            , only : isnan => shr_infnan_isnan
    use elm_varcon                , only : spval
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : GE_THERM_SNOW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SOIL_TBASED
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use petscsys
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ptm_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer           :: c,g,fc,j,l               ! do loop indices
    real(r8), pointer :: z(:,:)                   ! centroid at "z" level [m]
    real(r8), pointer :: zi(:,:)                  ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)                  ! layer thickness at "z" level (m)

    integer           :: imesh
    integer           :: ncells_local
    integer           :: ncells_ghost
    integer           :: nlev
    integer           :: first_active_soil_col_id
    integer           :: col_id
    integer           :: icell
    integer           :: iconn
    integer           :: snow_nconn
    integer           :: soil_nconn
    integer           :: ieqn
    integer           :: ieqn_1
    integer           :: ieqn_2
    integer           :: icond
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)

    real(r8), pointer :: snow_xc(:)               ! x-position of grid cell [m]
    real(r8), pointer :: snow_yc(:)               ! y-position of grid cell [m]
    real(r8), pointer :: snow_zc(:)               ! z-position of grid cell [m]
    real(r8), pointer :: snow_dx(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_dy(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_dz(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_area(:)             ! area of grid cell [m^2]
    integer , pointer :: snow_filter(:)           !
    integer , pointer :: snow_conn_id_up(:)       !
    integer , pointer :: snow_conn_id_dn(:)       !
    real(r8), pointer :: snow_conn_dist_up(:)     !
    real(r8), pointer :: snow_conn_dist_dn(:)     !
    real(r8), pointer :: snow_conn_area(:)        !
    integer , pointer :: snow_conn_type(:)        !

    real(r8), pointer :: ssw_xc(:)                ! x-position of grid cell [m]
    real(r8), pointer :: ssw_yc(:)                ! y-position of grid cell [m]
    real(r8), pointer :: ssw_zc(:)                ! z-position of grid cell [m]
    real(r8), pointer :: ssw_dx(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_dy(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_dz(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_area(:)              ! area of grid cell [m^2]
    integer , pointer :: ssw_filter(:)            !

    real(r8), pointer :: soil_xc(:)               ! x-position of grid cell [m]
    real(r8), pointer :: soil_yc(:)               ! y-position of grid cell [m]
    real(r8), pointer :: soil_zc(:)               ! z-position of grid cell [m]
    real(r8), pointer :: soil_dx(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_dy(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_dz(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_area(:)             ! area of grid cell [m^2]
    integer , pointer :: soil_conn_type(:)        !
    integer , pointer :: soil_filter(:)           !
    integer, pointer  :: soil_conn_id_up(:)       !
    integer, pointer  :: soil_conn_id_dn(:)       !
    real(r8), pointer :: soil_conn_dist_up(:)     !
    real(r8), pointer :: soil_conn_dist_dn(:)     !
    real(r8), pointer :: soil_conn_area(:)        !

    integer  , pointer :: col_active(:)
    integer  , pointer :: col_landunit(:)
    integer  , pointer :: lun_lakpoi(:)
    integer  , pointer :: lun_urbpoi(:)

    integer :: bounds_proc_begc, bounds_proc_endc

    !----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi             , zi           )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz             , dz           )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_z              , z            )

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active          , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index  , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_lakepoint  , lun_lakpoi   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_urbanpoint , lun_urbpoi   )

    ! Allocate memory and setup data structure for VSFM-MPP

    allocate (snow_xc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_yc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_zc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_dx           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_dy           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_dz           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_area         ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_filter       ((bounds_proc_endc-bounds_proc_begc+1 )*nlevsno      ))
    allocate (snow_conn_id_up   ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))
    allocate (snow_conn_id_dn   ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))
    allocate (snow_conn_dist_up ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))
    allocate (snow_conn_dist_dn ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))
    allocate (snow_conn_area    ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))
    allocate (snow_conn_type    ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevsno-1)  ))

    allocate (ssw_xc            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_yc            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_zc            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_dx            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_dy            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_dz            ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_area          ((bounds_proc_endc-bounds_proc_begc+1                )))
    allocate (ssw_filter        ((bounds_proc_endc-bounds_proc_begc+1                )))

    allocate (soil_xc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_yc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_zc           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_dx           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_dy           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_dz           ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_area         ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_filter       ((bounds_proc_endc-bounds_proc_begc+1 )*nlevgrnd     ))
    allocate (soil_conn_id_up   ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_id_dn   ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_dist_up ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_dist_dn ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_area    ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_type    ((bounds_proc_endc-bounds_proc_begc+1 )*(nlevgrnd-1) ))

    first_active_soil_col_id = -1
    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       if ((col_active(c) == 1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
          first_active_soil_col_id = c
          exit
       endif
    end do

    if (first_active_soil_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Save geometric attributes for snow mesh
    icell = 0
    do c = bounds_proc_begc, bounds_proc_endc
       do j = -nlevsno+1, 0
          icell = icell + 1

          if ((col_active(c) == 1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
             col_id = c
             snow_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             snow_filter(icell) = 0
          end if

          snow_xc(icell)   = 0.d0
          snow_yc(icell)   = 0.d0
          snow_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j))
          snow_dx(icell)   = 1.d0
          snow_dy(icell)   = 1.d0
          snow_dz(icell)   = dz(col_id,j)
          snow_area(icell) = 1.d0

       end do
    end do

    ! Save geometric attributes for standing surface water mesh
    icell = 0
    j     = 0
    do c = bounds_proc_begc, bounds_proc_endc

       icell = icell + 1

       if ((col_active(c) == 1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
          col_id = c
          ssw_filter(icell) = 1
       else
          col_id = first_active_soil_col_id
          ssw_filter(icell) = 0
       end if

       ssw_xc(icell)   = 0.d0
       ssw_yc(icell)   = 0.d0
       ssw_zc(icell)   = -0.5d0*(1.d-6 + zi(col_id,j))
       ssw_dx(icell)   = 1.d0
       ssw_dy(icell)   = 1.d0
       ssw_dz(icell)   = 1.d-6
       ssw_area(icell) = 1.d0

    end do

    ! Save geometric attributes for soil mesh
    icell = 0
    do c = bounds_proc_begc, bounds_proc_endc
       do j = 1, nlevgrnd
          icell = icell + 1

          if ((col_active(c) == 1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
             col_id = c
             soil_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             soil_filter(icell) = 0
          end if

          soil_xc(icell)   = 0.d0
          soil_yc(icell)   = 0.d0
          soil_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j))
          soil_dx(icell)   = 1.d0
          soil_dy(icell)   = 1.d0
          soil_dz(icell)   = dz(col_id,j)
          soil_area(icell) = 1.d0

       end do
    end do

    ! Save information about internal connections for snow mesh
    iconn = 0
    do c = bounds_proc_begc, bounds_proc_endc
       do j = -nlevsno+1, -1

          iconn = iconn + 1
          snow_conn_id_up(iconn)   = (c-bounds_proc_begc)*nlevsno + j  + nlevsno
          snow_conn_id_dn(iconn)   = snow_conn_id_up(iconn) + 1
          snow_conn_dist_up(iconn) = zi(c,j)   - z(c,j)
          snow_conn_dist_dn(iconn) = z( c,j+1) - zi(c,j)
          snow_conn_area(iconn)    = 1.d0
          snow_conn_type(iconn)    = CONN_VERTICAL

       end do
    end do
    snow_nconn = iconn

    ! Save information about internal connections for soil mesh
    iconn = 0
    do c = bounds_proc_begc, bounds_proc_endc

       do j = 1, nlevgrnd-1

          iconn = iconn + 1
          soil_conn_id_up(iconn)   = (c-bounds_proc_begc)*nlevgrnd + j
          soil_conn_id_dn(iconn)   = soil_conn_id_up(iconn) + 1
          soil_conn_dist_up(iconn) = zi(c,j)   - z(c,j)
          soil_conn_dist_dn(iconn) = z( c,j+1) - zi(c,j)
          soil_conn_area(iconn)    = 1.d0
          soil_conn_type(iconn)    = CONN_VERTICAL

       end do
    end do
    soil_nconn = iconn

    !
    ! Set up the meshes
    !

    call this%thermal_mpp%SetNumMeshes (3)

    !
    ! Set mesh for snow
    !
    imesh        = 1
    nlev         = nlevsno
    ncells_local = (bounds_proc_endc - bounds_proc_begc + 1)*nlev
    ncells_ghost = 0

    call this%thermal_mpp%MeshSetName        (imesh, 'CLM snow thermal mesh')
    call this%thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call this%thermal_mpp%MeshSetID          (imesh, MESH_CLM_SNOW_COL)
    call this%thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call this%thermal_mpp%MeshSetGridCellFilter      (imesh, snow_filter)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , snow_xc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , snow_yc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , snow_zc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , snow_dx)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , snow_dy)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , snow_dz)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , snow_area)
    call this%thermal_mpp%MeshComputeVolume          (imesh)

    call this%thermal_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, snow_nconn,  &
         snow_conn_id_up, snow_conn_id_dn, snow_conn_dist_up, snow_conn_dist_dn, &
         snow_conn_area, snow_conn_type)

    !
    ! Set mesh for standing water
    !
    imesh        = 2
    nlev         = 1
    ncells_local = (bounds_proc_endc - bounds_proc_begc + 1)*nlev
    ncells_ghost = 0
    call this%thermal_mpp%MeshSetName(imesh, 'CLM standing water thermal snow mesh')
    call this%thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call this%thermal_mpp%MeshSetID          (imesh, MESH_CLM_SSW_COL)
    call this%thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call this%thermal_mpp%MeshSetGridCellFilter      (imesh, ssw_filter)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , ssw_xc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , ssw_yc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , ssw_zc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , ssw_dx)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , ssw_dy)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , ssw_dz)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , ssw_area)
    call this%thermal_mpp%MeshComputeVolume          (imesh)

    !
    ! Set mesh for standing water
    !
    imesh        = 3
    nlev         = nlevgrnd
    ncells_local = (bounds_proc_endc - bounds_proc_begc + 1)*nlev
    ncells_ghost = 0

    call this%thermal_mpp%MeshSetName        (imesh, 'CLM soil thermal mesh')
    call this%thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call this%thermal_mpp%MeshSetID          (imesh, MESH_CLM_THERMAL_SOIL_COL)
    call this%thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call this%thermal_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call this%thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call this%thermal_mpp%MeshComputeVolume          (imesh)

    call this%thermal_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, soil_nconn,  &
         soil_conn_id_up, soil_conn_id_dn, soil_conn_dist_up, soil_conn_dist_dn, &
         soil_conn_area, soil_conn_type)

    ! Free up memory
    deallocate (snow_xc           )
    deallocate (snow_yc           )
    deallocate (snow_zc           )
    deallocate (snow_dz           )
    deallocate (snow_area         )
    deallocate (snow_filter       )
    deallocate (snow_conn_id_up   )
    deallocate (snow_conn_id_dn   )
    deallocate (snow_conn_dist_up )
    deallocate (snow_conn_dist_dn )
    deallocate (snow_conn_area    )
    deallocate (snow_conn_type    )

    deallocate (ssw_xc            )
    deallocate (ssw_yc            )
    deallocate (ssw_zc            )
    deallocate (ssw_dz            )
    deallocate (ssw_area          )
    deallocate (ssw_filter        )

    deallocate (soil_xc           )
    deallocate (soil_yc           )
    deallocate (soil_zc           )
    deallocate (soil_dz           )
    deallocate (soil_area         )
    deallocate (soil_filter       )
    deallocate (soil_conn_id_up   )
    deallocate (soil_conn_id_dn   )
    deallocate (soil_conn_dist_up )
    deallocate (soil_conn_dist_dn )
    deallocate (soil_conn_area    )
    deallocate (soil_conn_type    )

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : GE_THERM_SNOW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SOIL_TBASED
    use MultiPhysicsProbConstants , only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    !
    implicit none
    !
    class(em_ptm_type)                   :: this

    call this%thermal_mpp%AddGovEqn(GE_THERM_SNOW_TBASED,                                                  &
         'Thermal equation using temprature formulation in snow (KSP formulation)',                   &
         MESH_CLM_SNOW_COL)

    call this%thermal_mpp%AddGovEqn(GE_THERM_SSW_TBASED,                                                   &
         'Thermal equation using temprature formulation in standing surface water (KSP formulation)', &
         MESH_CLM_SSW_COL)

    call this%thermal_mpp%AddGovEqn(GE_THERM_SOIL_TBASED,                                                  &
         'Thermal equation using temprature formulation in soil (KSP formulation)',                   &
         MESH_CLM_THERMAL_SOIL_COL)

    call this%thermal_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ptm_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    integer           :: c                        ! do loop indices
    integer           :: ieqn, ieqn_1, ieqn_2
    integer           :: icond
    real(r8), pointer :: z(:,:)                   ! centroid at "z" level [m]
    real(r8), pointer :: zi(:,:)                  ! interface level below a "z" level (m)
    real(r8), pointer :: soil_top_conn_dist_dn(:) !

    integer  , pointer :: col_active(:)
    integer  , pointer :: col_landunit(:)
    integer  , pointer :: lun_lakpoi(:)
    integer  , pointer :: lun_urbpoi(:)
    integer            :: bounds_proc_begc, bounds_proc_endc

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi, zi)
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_z , z )

    allocate (soil_top_conn_dist_dn(bounds_proc_endc-bounds_proc_begc+1))

    do c = bounds_proc_begc, bounds_proc_endc
       soil_top_conn_dist_dn(c - bounds_proc_begc + 1) = z(c,1) - zi(c,0)
    end do

    !
    ! Add BC/SS
    !
    ieqn = 1
    call this%thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_snow', 'W/m^2', COND_HEAT_FLUX, &
         SNOW_TOP_CELLS)
    call this%thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,           &
         'Absorbed_solar_radiation', 'W/m^2', COND_HEAT_RATE, &
         ALL_CELLS)

    ieqn = 2
    call this%thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_standing_surface_water', 'W/m^2', COND_HEAT_FLUX, &
         SSW_TOP_CELLS)

    ieqn = 3
    call this%thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_soil', 'W/m^2', COND_HEAT_FLUX, &
         SOIL_TOP_CELLS)
    call this%thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,           &
         'Absorbed_solar_radiation', 'W/m^2', COND_HEAT_RATE, &
         ALL_CELLS)

    !
    ! Add coupling boundary condition
    !
    ieqn_1 = 1
    ieqn_2 = 3
    call this%thermal_mpp%GovEqnAddCouplingCondition(ieqn_1, ieqn_2, &
         SNOW_BOTTOM_CELLS, SOIL_TOP_CELLS)

    ieqn_1 = 2
    ieqn_2 = 3
    call this%thermal_mpp%GovEqnAddCouplingCondition(ieqn_1, ieqn_2, &
         SSW_TOP_CELLS, SOIL_TOP_CELLS)

    !
    ! Hack to update the interface distance for soil governing equation because
    ! the centroid of a soil layer is not located at a depth which is
    ! average of the bounding interface depths
    !
    !    Location of       Centroid
    !  CLM's centroid          that is equidistant
    !  ---------             ---------
    !      |                     |
    !      o                     |
    !      |                     o
    !      |                     |
    !      |                     |
    !  ---------             ---------
    !
    ieqn_2 = 3
    icond  = 2
    call this%thermal_mpp%GovEqnUpdateBCConnectionSet(ieqn_2, icond, VAR_DIST_DN, &
         bounds_proc_endc - bounds_proc_begc + 1, &
         soil_top_conn_dist_dn)
    icond  = 3
    call this%thermal_mpp%GovEqnUpdateBCConnectionSet(ieqn_2, icond, VAR_DIST_DN, &
         bounds_proc_endc - bounds_proc_begc + 1, &
         soil_top_conn_dist_dn)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : VAR_DZ
    !
    implicit none
    !
    class(em_ptm_type)                   :: this
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)

    !
    ! Allocate auxvars
    !
    call this%thermal_mpp%AllocateAuxVars()

    !
    ! Set variables to be exchanged among the coupled equations
    !

    ieqn = 1
    nvars_for_coupling = 2
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_TEMPERATURE
    var_ids_for_coupling    (2) = VAR_THERMAL_COND
    goveqn_ids_for_coupling (1) = 3
    goveqn_ids_for_coupling (2) = 3

    call this%thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    ieqn = 2
    call this%thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    ieqn = 3
    nvars_for_coupling = 10
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    ( 1) = VAR_TEMPERATURE
    var_ids_for_coupling    ( 2) = VAR_THERMAL_COND
    var_ids_for_coupling    ( 3) = VAR_FRAC
    var_ids_for_coupling    ( 4) = VAR_ACTIVE
    var_ids_for_coupling    ( 5) = VAR_DIST_UP
    var_ids_for_coupling    ( 6) = VAR_TEMPERATURE
    var_ids_for_coupling    ( 7) = VAR_THERMAL_COND
    var_ids_for_coupling    ( 8) = VAR_FRAC
    var_ids_for_coupling    ( 9) = VAR_ACTIVE
    var_ids_for_coupling    (10) = VAR_DZ

    goveqn_ids_for_coupling ( 1) = 1
    goveqn_ids_for_coupling ( 2) = 1
    goveqn_ids_for_coupling ( 3) = 1
    goveqn_ids_for_coupling ( 4) = 1
    goveqn_ids_for_coupling ( 5) = 1
    goveqn_ids_for_coupling ( 6) = 2
    goveqn_ids_for_coupling ( 7) = 2
    goveqn_ids_for_coupling ( 8) = 2
    goveqn_ids_for_coupling ( 9) = 2
    goveqn_ids_for_coupling (10) = 2

    call this%thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine add_material_properties(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use elm_varcon                , only : spval
    use mpp_varpar                , only : nlevgrnd, nlevsno
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ptm_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    integer           :: c,g,fc,j,l               ! do loop indices
    integer, pointer  :: thermal_filter(:)
    integer, pointer  :: thermal_lun_type(:)
    real(r8), pointer :: thermal_watsat(:,:)
    real(r8), pointer :: thermal_csol(:,:)
    real(r8), pointer :: thermal_tkmg(:,:)
    real(r8), pointer :: thermal_tkdry(:,:)
    real(r8), pointer :: clm_watsat(:,:)
    real(r8), pointer :: clm_csol(:,:)
    real(r8), pointer :: clm_tkmg(:,:)
    real(r8), pointer :: clm_tkdry(:,:)

    integer  , pointer :: col_active(:)
    integer  , pointer :: col_landunit(:)

    integer  , pointer :: lun_type(:)
    integer  , pointer :: lun_lakpoi(:)
    integer  , pointer :: lun_urbpoi(:)

    integer            :: bounds_proc_begc, bounds_proc_endc

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active          , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index  , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type       , lun_type     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_lakepoint  , lun_lakpoi   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_urbanpoint , lun_urbpoi   )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_watsat, clm_watsat)
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_csol  , clm_csol  )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_tkmg  , clm_tkmg  )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_tkdry , clm_tkdry )

    ! Allocate memory and setup data structure for VSFM-MPP

    allocate (thermal_filter   (bounds_proc_begc:bounds_proc_endc      ))
    allocate (thermal_lun_type (bounds_proc_begc:bounds_proc_endc      ))
    allocate (thermal_watsat   (bounds_proc_begc:bounds_proc_endc,nlevgrnd ))
    allocate (thermal_csol     (bounds_proc_begc:bounds_proc_endc,nlevgrnd ))
    allocate (thermal_tkmg     (bounds_proc_begc:bounds_proc_endc,nlevgrnd ))
    allocate (thermal_tkdry    (bounds_proc_begc:bounds_proc_endc,nlevgrnd ))

    thermal_filter(:)   = 0
    thermal_lun_type(:) = 0
    thermal_watsat(:,:) = 0._r8
    thermal_csol(:,:)   = 0._r8
    thermal_tkmg(:,:)   = 0._r8
    thermal_tkdry(:,:)  = 0._r8

    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       thermal_lun_type(c)  = lun_type(l)

       if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
          thermal_filter(c) = 1
          do j = 1,nlevgrnd
             if (clm_watsat(c,j) /= spval) thermal_watsat(c,j) = clm_watsat(c,j)
             if (clm_csol(  c,j) /= spval) thermal_csol  (c,j) = clm_csol  (c,j)
             if (clm_tkmg(  c,j) /= spval) thermal_tkmg  (c,j) = clm_tkmg  (c,j)
             if (clm_tkdry( c,j) /= spval) thermal_tkdry (c,j) = clm_tkdry (c,j)
          enddo
          j = 1;

       endif
    enddo

    call MPPThermalSetSoils(this%thermal_mpp, bounds_proc_begc, &
                           bounds_proc_endc,               &
                           thermal_filter,                 &
                           thermal_lun_type,               &
                           thermal_watsat,                 &
                           thermal_csol,                   &
                           thermal_tkmg,                   &
                           thermal_tkdry                   &
                           )

  end subroutine add_material_properties

    !------------------------------------------------------------------------
  subroutine EM_PTM_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, &
       e2l_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! The PTM dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ptm_type)                   :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump
    !

    select case (em_stage)
    case (EM_PTM_TBASED_SOLVE_STAGE)
       call EM_PTM_TBased_Solve(this, dt, nstep, l2e_list, e2l_list, bounds_clump)
    case default
       write(iulog,*)'EM_PTM_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_PTM_Solve

    !------------------------------------------------------------------------
  subroutine EM_PTM_TBased_Solve(this, dt, nstep, l2e_list, e2l_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! The PETSc-based Thermal Model dirver
    !
#include <petsc/finclude/petsc.h>
    !
    use mpp_varpar                , only : nlevgrnd, nlevsno
    use mpp_varcon                , only : capr
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_AREAL_DEN
    use MultiPhysicsProbConstants , only : VAR_ICE_AREAL_DEN
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_SNOW_WATER
    use MultiPhysicsProbConstants , only : VAR_NUM_SNOW_LYR
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : VAR_TUNING_FACTOR
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_DHS_DT
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use SystemOfEquationsBaseType   , only : sysofeqns_base_type
    use SystemOfEquationsThermalType, only : sysofeqns_thermal_type
    
    use petscsys
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ptm_type)                   :: this
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: j,c,l,idx
    integer                              :: offset
    integer                              :: fc
    integer                              :: soe_auxvar_id
    integer                              :: begc, endc

    ! internal auxvars
    real(r8) , pointer                   :: temperature_1d         (:)
    real(r8) , pointer                   :: liq_areal_den_1d       (:)
    real(r8) , pointer                   :: ice_areal_den_1d       (:)
    real(r8) , pointer                   :: snow_water_1d          (:)
    real(r8) , pointer                   :: dz_1d                  (:)
    real(r8) , pointer                   :: dist_up_1d             (:)
    real(r8) , pointer                   :: dist_dn_1d             (:)
    real(r8) , pointer                   :: frac_1d                (:)
    integer  , pointer                   :: num_snow_layer_1d      (:)
    logical  , pointer                   :: is_active_1d           (:)

    ! boundary auxvars
    real(r8) , pointer                   :: hs_snow_1d             (:)
    real(r8) , pointer                   :: hs_sh2o_1d             (:)
    real(r8) , pointer                   :: hs_soil_1d             (:)
    real(r8) , pointer                   :: dhsdT_snow_1d          (:)
    real(r8) , pointer                   :: dhsdT_sh2o_1d          (:)
    real(r8) , pointer                   :: dhsdT_soil_1d          (:)
    real(r8) , pointer                   :: frac_soil_1d           (:)

    ! source sink
    real(r8) , pointer                   :: sabg_snow_1d           (:)
    real(r8) , pointer                   :: sabg_soil_1d           (:)

    real(r8) , pointer                   :: tsurf_tuning_factor_1d (:)

    real(r8) , pointer                   :: l2e_frac_snow_eff(:)
    real(r8) , pointer                   :: l2e_frac_h2osfc(:)
    real(r8) , pointer                   :: l2e_h2osno(:)
    real(r8) , pointer                   :: l2e_h2osfc(:)
    real(r8) , pointer                   :: l2e_th2osfc(:)
    real(r8) , pointer                   :: l2e_hs_soil(:)
    real(r8) , pointer                   :: l2e_hs_top_snow(:)
    real(r8) , pointer                   :: l2e_hs_h2osfc(:)
    real(r8) , pointer                   :: l2e_dhsdT(:)

    real(r8) , pointer                   :: l2e_h2osoi_liq_nlevgrnd(:,:)
    real(r8) , pointer                   :: l2e_h2osoi_ice_nlevgrnd(:,:)
    real(r8) , pointer                   :: l2e_h2osoi_liq_nlevsnow(:,:)
    real(r8) , pointer                   :: l2e_h2osoi_ice_nlevsnow(:,:)
    real(r8) , pointer                   :: l2e_tsoil(:,:)
    real(r8) , pointer                   :: l2e_tsnow(:,:)
    real(r8) , pointer                   :: l2e_sabg_lyr(:,:)

    real(r8) , pointer                   :: e2l_th2osfc(:)
    real(r8) , pointer                   :: e2l_tsoil(:,:)
    real(r8) , pointer                   :: e2l_tsnow(:,:)

    real(r8) , pointer                   :: col_dz(:,:)
    real(r8) , pointer                   :: col_zi(:,:)
    real(r8) , pointer                   :: col_z(:,:)
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_type(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: col_snl(:)

    integer  , pointer                   :: lun_lakpoi(:)
    integer  , pointer                   :: lun_urbpoi(:)

    integer  , pointer                   :: l2e_filter(:)
    integer                              :: l2e_num_filter

    PetscErrorCode                       :: ierr
    PetscBool                            :: converged
    PetscInt                             :: converged_reason

    class(sysofeqns_base_type), pointer  :: soe

    begc = bounds_clump%begc
    endc = bounds_clump%endc

    allocate(temperature_1d         ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(liq_areal_den_1d       ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(ice_areal_den_1d       ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(snow_water_1d          ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(dz_1d                  ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(dist_up_1d             ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(dist_dn_1d             ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(frac_1d                ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))

    allocate(num_snow_layer_1d      ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))
    allocate(is_active_1d           ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))

    allocate(hs_snow_1d             ((endc-begc+1                      )))
    allocate(hs_sh2o_1d             ((endc-begc+1                      )))
    allocate(hs_soil_1d             ((endc-begc+1                      )))
    allocate(dhsdT_snow_1d          ((endc-begc+1                      )))
    allocate(dhsdT_sh2o_1d          ((endc-begc+1                      )))
    allocate(dhsdT_soil_1d          ((endc-begc+1                      )))
    allocate(frac_soil_1d           ((endc-begc+1                      )))

    allocate(sabg_snow_1d           ((endc-begc+1)*nlevsno ))
    allocate(sabg_soil_1d           ((endc-begc+1)*nlevgrnd))

    allocate(tsurf_tuning_factor_1d ((endc-begc+1)*(nlevgrnd+nlevsno+1 )))

    call l2e_list%GetPointerToReal1D(this%index_l2e_state_frac_snow_eff            , l2e_frac_snow_eff)
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_frac_h2osfc              , l2e_frac_h2osfc)
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_h2osno                   , l2e_h2osno)
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_h2osfc                   , l2e_h2osfc)
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_temperature_h2osfc       , l2e_th2osfc)
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_hs_soil                   , l2e_hs_soil)
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_hs_top_snow               , l2e_hs_top_snow)
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_hs_h2osfc                 , l2e_hs_h2osfc)
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dhsdT                     , l2e_dhsdT)

    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq_nlevgrnd      , l2e_h2osoi_liq_nlevgrnd)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice_nlevgrnd      , l2e_h2osoi_ice_nlevgrnd)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq_nlevsnow      , l2e_h2osoi_liq_nlevsnow)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice_nlevsnow      , l2e_h2osoi_ice_nlevsnow)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_temperature_soil_nlevgrnd, l2e_tsoil)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_temperature_snow         , l2e_tsnow)
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_sabg_lyr                  , l2e_sabg_lyr)

    call l2e_list%GetPointerToReal2D(this%index_l2e_col_zi             , col_zi           )
    call l2e_list%GetPointerToReal2D(this%index_l2e_col_dz             , col_dz           )
    call l2e_list%GetPointerToReal2D(this%index_l2e_col_z              , col_z            )

    call l2e_list%GetPointerToInt1D(this%index_l2e_col_num_snow_lyrs    , col_snl   )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_active          , col_active   )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_landunit_index  , col_landunit )
    call l2e_list%GetPointerToInt1D(this%index_l2e_landunit_lakepoint  , lun_lakpoi   )
    call l2e_list%GetPointerToInt1D(this%index_l2e_landunit_urbanpoint , lun_urbpoi   )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_nolakec_and_nourbanc , l2e_filter )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_nolakec_and_nourbanc   , l2e_num_filter    )

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_temperature_soil_nlevgrnd, e2l_tsoil)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_temperature_snow         , e2l_tsnow)
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_temperature_h2osfc       , e2l_th2osfc)

    ! Initialize
    temperature_1d(:)         = 273.15_r8
    liq_areal_den_1d(:)       = 0._r8
    ice_areal_den_1d(:)       = 0._r8
    frac_soil_1d(:)           = 1._r8
    frac_1d(:)                = 1._r8
    num_snow_layer_1d(:)      = 0
    is_active_1d(:)           = .false.
    hs_snow_1d(:)             = 0._r8
    hs_sh2o_1d(:)             = 0._r8
    hs_soil_1d(:)             = 0._r8
    dhsdT_snow_1d(:)          = 0._r8
    dhsdT_sh2o_1d(:)          = 0._r8
    dhsdT_soil_1d(:)          = 0._r8
    sabg_snow_1d(:)           = 0._r8
    sabg_soil_1d(:)           = 0._r8
    snow_water_1d(:)          = 0._r8
    tsurf_tuning_factor_1d(:) = 1._r8
    dz_1d(:)                  = 0._r8
    dist_up_1d(:)             = 0._r8
    dist_dn_1d(:)             = 0._r8

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Save data for snow
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    offset = 0

    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)
       do j = -nlevsno+1, 0

          l = col_landunit(c)

          ! Is this a soil column on which PETSc based thermal solver works?
          if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then

             if (j >= col_snl(c)+1) then

                ! Index for internal SoE auxvars
                idx = (c-begc)*nlevsno + j + nlevsno + offset

                ! Save data for internal SoE auxvars
                temperature_1d(idx)    = l2e_tsnow(c,j)
                dz_1d(idx)             = col_dz(c,j)
                liq_areal_den_1d(idx)  = l2e_h2osoi_liq_nlevsnow(c,j)
                ice_areal_den_1d(idx)  = l2e_h2osoi_ice_nlevsnow(c,j)
                num_snow_layer_1d(idx) = -col_snl(c)
                is_active_1d(idx)      = .true.
                dist_up_1d(idx)        = col_zi(c,j  ) - col_z(c,j)
                dist_dn_1d(idx)        = col_z(c,j) - col_zi(c,j-1)
                frac_1d(idx)           = l2e_frac_snow_eff(c)

                ! If not the top snow layer, save amount of absorbed solar
                ! radiation
                if (j /= col_snl(c) +  1) then
                   sabg_snow_1d(idx)      = l2e_sabg_lyr(c,j)
                endif

                ! Save follow data only for the top snow layer
                if (j == col_snl(c)+1) then

                   ! Save tuning_factor for internal SoE auxvars
                   tsurf_tuning_factor_1d(idx) = col_dz(c,j) / &
                        (0.5_r8*(col_z(c,j)-col_zi(c,j-1)+capr*(col_z(c,j+1)-col_zi(c,j-1))))

                   ! Index for boundary SoE auxvars
                   idx = (c-begc)+1

                   ! Save data for boundary SoE auxvars
                   hs_snow_1d(idx)    = l2e_hs_top_snow(c)
                   dhsdT_snow_1d(idx) = l2e_dhsdT(c)
                   frac_soil_1d(idx)  = frac_soil_1d(idx) - l2e_frac_snow_eff(c)
                endif

             endif
          endif
       enddo
    enddo

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Save data for h2osfc
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    offset = (endc - begc + 1)*nlevsno ! Number of data for snow

    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)
       l = col_landunit(c)

       if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then

          if (l2e_frac_h2osfc(c) > 0._r8) then

             ! Index for internal SoE auxvars
             idx = (c-begc) + 1 + offset

             ! Save data for internal SoE auxvars
             temperature_1d(idx) = l2e_th2osfc(c)
             dz_1d(idx)          = 1.0e-3*l2e_h2osfc(c)/l2e_frac_h2osfc(c)
             is_active_1d(idx)   = .true.
             frac_1d(idx)        = l2e_frac_h2osfc(c)
             dist_up_1d(idx)     = dz_1d(idx)/2.d0/l2e_frac_h2osfc(c)
             dist_dn_1d(idx)     = dz_1d(idx)/2.d0/l2e_frac_h2osfc(c)

             ! Index for boundary SoE auxvars
             idx                 = (c-begc) + 1

             ! Save data for boundary SoE auxvars
             frac_soil_1d(idx)   = frac_soil_1d(idx) - l2e_frac_h2osfc(c)
             dhsdT_sh2o_1d(idx)  = l2e_dhsdT(c)
             hs_sh2o_1d(idx)     = l2e_hs_h2osfc(c)
          endif
       endif
    enddo

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Save data for soil
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    offset = (endc - begc + 1)*(nlevsno + 1) ! Number of data for snow + sh2o

    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)
       l = col_landunit(c)

       do j = 1,nlevgrnd

          ! Index for internal SoE auxvars
          idx = (c-begc)*nlevgrnd + j + offset

          if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then

             ! Save data for internal SoE auxvars
             temperature_1d(idx)    = l2e_tsoil(c,j)
             dz_1d(idx)             = col_dz(c,j)
             is_active_1d(idx)      = .true.
             liq_areal_den_1d(idx)  = l2e_h2osoi_liq_nlevgrnd(c,j)
             ice_areal_den_1d(idx)  = l2e_h2osoi_ice_nlevgrnd(c,j)
             frac_1d(idx)           = 1.0_r8
             dist_up_1d(idx)        = col_zi(c,j) - col_z(c,j)
             dist_dn_1d(idx)        = col_zi(c,j) - col_z(c,j)

             ! Is this the top soil layer?
             if (j == 1) then

                dz_1d(idx)             = col_z(c,j)*2.d0
                num_snow_layer_1d(idx) = -col_snl(c)

                ! Save data for internal SoE auxvars
                snow_water_1d(idx)                         = l2e_h2osno(c)

                if (col_snl(c) /= 0) then
                   ! Save data for boundary SoE auxvars
                   sabg_soil_1d((c-begc)*nlevgrnd + j) = l2e_frac_snow_eff(c)*l2e_sabg_lyr(c,j)
                else
                   ! Save data for internal SoE auxvars
                   tsurf_tuning_factor_1d(idx) = col_dz(c,j) / &
                        (0.5_r8*(col_z(c,j)-col_zi(c,j-1)+capr*(col_z(c,j+1)-col_zi(c,j-1))))
                endif

                ! Save data for boundary SoE auxvars
                hs_soil_1d(   (c-begc) + 1) = l2e_hs_soil(c)
                dhsdT_soil_1d((c-begc) + 1) = l2e_dhsdT(c)

             endif
          endif
       enddo
    enddo

    soe => this%thermal_mpp%soe
    select type(soe)
    class is(sysofeqns_thermal_type)

       ! Set temperature
       call soe%SetSolnPrevCLM(temperature_1d)

       ! Set h2soi_liq
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_LIQ_AREAL_DEN, soe_auxvar_id, liq_areal_den_1d)

       ! Set h2osi_ice
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_ICE_AREAL_DEN, soe_auxvar_id, ice_areal_den_1d)

       ! Set snow water
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_SNOW_WATER, soe_auxvar_id, snow_water_1d)

       ! Set dz
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_DZ, soe_auxvar_id, dz_1d)

       ! Set dist_up
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_DIST_UP, soe_auxvar_id, dist_up_1d)

       ! Set dist_dn
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_DIST_DN, soe_auxvar_id, dist_dn_1d)

       ! Set number of snow layers
       soe_auxvar_id = 1;
       call soe%SetIDataFromCLM(AUXVAR_INTERNAL, &
            VAR_NUM_SNOW_LYR, soe_auxvar_id, num_snow_layer_1d)

       ! Set if cell is active
       call soe%SetBDataFromCLM(AUXVAR_INTERNAL, &
            VAR_ACTIVE, is_active_1d)

       ! Set tuning factor
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_TUNING_FACTOR, soe_auxvar_id, tsurf_tuning_factor_1d)

       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_FRAC, soe_auxvar_id, frac_1d)

       !
       ! Set heat flux for:
       !

       ! 1) top snow layer
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, hs_snow_1d)

       ! 2) top standing water layer
       soe_auxvar_id = 2;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, hs_sh2o_1d)

       ! 3) soil
       soe_auxvar_id = 3;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, hs_soil_1d)

       !
       ! Set derivative of heat flux w.r.t temperature for:
       !

       ! 1) top snow layer
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_DHS_DT, soe_auxvar_id, dhsdT_snow_1d)

       ! 2) top standing water layer
       soe_auxvar_id = 2;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_DHS_DT, soe_auxvar_id, dhsdT_sh2o_1d)

       ! 3) soil
       soe_auxvar_id = 3;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_DHS_DT, soe_auxvar_id, dhsdT_soil_1d)

       !
       ! Set fraction of soil not covered by snow and standing water
       !
       soe_auxvar_id = 3;
       call soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_FRAC, soe_auxvar_id, frac_soil_1d)


       ! Set absorbed solar radiation
       soe_auxvar_id = 1;
       call soe%SetRDataFromCLM(AUXVAR_SS, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, sabg_snow_1d)

       ! Set absorbed solar radiation
       soe_auxvar_id = 2;
       call soe%SetRDataFromCLM(AUXVAR_SS, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, sabg_soil_1d)
    end select


    ! Preform Pre-StepDT operations
    call this%thermal_mpp%soe%PreStepDT()

    ! Solve
    call this%thermal_mpp%soe%StepDT(dt, nstep, &
         converged, converged_reason, ierr); CHKERRQ(ierr)

    ! Did the model converge
    if (.not. converged) then
       call endrun(msg=' ERROR: PETSc thermal model failed to converge '//&
            errMsg(__FILE__, __LINE__))
    endif

    soe => this%thermal_mpp%soe
    select type(soe)
    class is(sysofeqns_thermal_type)
       ! Get the updated soil tempreature
       call soe%GetSoln(temperature_1d)
    end select

    ! Put temperature back in ALM structure for snow
    offset = 0

    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)

       do j = -nlevsno+1, 0
          idx = (c-begc)*nlevsno + j + nlevsno + offset
          l = col_landunit(c)
          if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
             if (j >= col_snl(c)+1) then
                e2l_tsnow(c,j) = temperature_1d(idx)
             endif
          endif
       enddo
    enddo

    ! Put temperature back in ALM structure for soil
    offset = (endc - begc + 1)*(nlevsno + 1)
    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)
       do j = 1,nlevgrnd
          idx = (c-begc)*nlevgrnd + j + offset
          l = col_landunit(c)
          if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
             e2l_tsoil(c,j) = temperature_1d(idx)
          endif
       enddo
    enddo

    ! Put temperature back in ALM structure for standing water
    ! NOTE: Soil temperature needs to be updated first
    offset = (endc - begc + 1)*nlevsno
    do fc = 1, l2e_num_filter
       c = l2e_filter(fc)
       idx = (c-begc) + 1 + offset
       l = col_landunit(c)
       if ((col_active(c)==1) .and. (lun_lakpoi(l) == 0) .and. (lun_urbpoi(l) == 0)) then
          if (l2e_frac_h2osfc(c) > 0._r8) then
             e2l_th2osfc(c) = temperature_1d(idx)
          else
             e2l_th2osfc(c) = e2l_tsoil(c,1)
          endif
       endif
    enddo

    ! Free up memory
    deallocate (temperature_1d         )
    deallocate (liq_areal_den_1d       )
    deallocate (ice_areal_den_1d       )
    deallocate (snow_water_1d          )
    deallocate (dz_1d                  )
    deallocate (dist_up_1d             )
    deallocate (dist_dn_1d             )
    deallocate (frac_1d                )
    deallocate (num_snow_layer_1d      )
    deallocate (is_active_1d           )
    deallocate (hs_snow_1d             )
    deallocate (hs_sh2o_1d             )
    deallocate (hs_soil_1d             )
    deallocate (dhsdT_snow_1d          )
    deallocate (dhsdT_sh2o_1d          )
    deallocate (dhsdT_soil_1d          )
    deallocate (frac_soil_1d           )
    deallocate (sabg_snow_1d           )
    deallocate (sabg_soil_1d           )
    deallocate (tsurf_tuning_factor_1d )

  end subroutine EM_PTM_TBased_Solve

#endif

end module ExternalModelPTMMod
