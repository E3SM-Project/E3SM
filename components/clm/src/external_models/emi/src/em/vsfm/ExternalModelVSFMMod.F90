module ExternalModelVSFMMod

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use ExternalModelInterfaceDataMod, only : emi_data_list, emi_data
  use mpp_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use MultiPhysicsProbVSFM         , only : mpp_vsfm_type
  use decompMod                    , only : bounds_type
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
  implicit none
  !

  type, public, extends(em_base_type) :: em_vsfm_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     integer :: index_e2l_init_state_smp
     integer :: index_e2l_init_state_wtd

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_smp
     integer :: index_e2l_state_wtd
     integer :: index_e2l_state_soilp

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi

     ! IDs to indentify the conditions for VSFM
     integer :: vsfm_cond_id_for_infil
     integer :: vsfm_cond_id_for_et
     integer :: vsfm_cond_id_for_dew
     integer :: vsfm_cond_id_for_drainage
     integer :: vsfm_cond_id_for_snow
     integer :: vsfm_cond_id_for_sublimation
     integer :: vsfm_cond_id_for_lateral_flux

     type(mpp_vsfm_type) :: vsfm_mpp

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_VSFM_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_VSFM_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_VSFM_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_VSFM_Populate_E2L_List
     procedure, public :: Init                    => EM_VSFM_Init
     procedure, public :: Solve                   => EM_VSFM_Solve
  end type em_vsfm_type

contains

  !------------------------------------------------------------------------
  subroutine EM_VSFM_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by VSFM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp            = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

  end subroutine EM_VSFM_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_VSFM_Populate_E2L_Init_List(this, e2l_init_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by VSFM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                  :: this
    class(emi_data_list) , intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data)      , pointer       :: data
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_smp             = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

    deallocate(em_stages)

  end subroutine EM_VSFM_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_VSFM_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by VSFM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_VSFM_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    deallocate(em_stages)

  end subroutine EM_VSFM_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_VSFM_Populate_E2L_List(this, e2l_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by VSFM from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                  :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_VSFM_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_smp        = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_FLUX_AQUIFER_RECHARGE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrecharge   = index

    deallocate(em_stages)

end subroutine EM_VSFM_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_VSFM_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
#include <petsc/finclude/petsc.h>
    !
    !
    ! !USES:
    use mpp_varctl                , only : vsfm_use_dynamic_linesearch
    use petscsnes
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump

    !
    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(this, iam)

    ! 2. Add all meshes needed for the MPP
    call add_meshes(this, l2e_init_list, bounds_clump)

    ! 3. Add all governing equations
    call add_goveqns(this)

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(this)

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars(this)

    ! 6. Setup the MPP
    call this%vsfm_mpp%SetupProblem(vsfm_use_dynamic_linesearch)

    ! 7. Add material properities associated with all governing equations
    call set_material_properties(this, l2e_init_list, bounds_clump)

    ! 8. Set initial conditions
    call set_initial_conditions(this, l2e_init_list, bounds_clump)

    ! 9. Determine IDs for various source-sink condition
    call determine_condition_ids(this)

    !10.
    call extract_data_for_alm(this, l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine EM_VSFM_Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(em_vsfm, iam)
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type) :: em_vsfm
    integer, intent(in) :: iam

    !
    ! Set up the multi-physics problem
    !
    call em_vsfm%vsfm_mpp%Init       ()
    call em_vsfm%vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model')
    call em_vsfm%vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call em_vsfm%vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Add meshes used in the VSFM MPP
    !
    ! !USES:
    use mpp_varcon                , only : istcrop, istsoil
    use mpp_varcon                , only : icol_road_perv
    use mpp_varpar                , only : nlevgrnd
    use mpp_varctl                , only : lateral_connectivity
    use mpp_varcon                , only : max_lunit
    use mpp_varctl                , only : vsfm_lateral_model_type
    !use mpp_bounds                , only : bounds_proc_begg_all, bounds_proc_endg_all
    !use mpp_bounds                , only : bounds_proc_begc_all, bounds_proc_endc_all
    !use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    !use mpp_bounds                , only : nclumps
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : DISCRETIZATION_VERTICAL_ONLY
    use MultiPhysicsProbConstants , only : DISCRETIZATION_VERTICAL_WITH_SS
    use MultiPhysicsProbConstants , only : DISCRETIZATION_THREE_DIM
    !
    implicit none
    !
    class(em_vsfm_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g,fc,j,l           ! do loop indices
    integer            :: imesh
    integer            :: nlev
    integer            :: first_active_soil_col_id
    integer            :: col_id
    integer            :: icell
    integer            :: iconn
    integer            :: horz_nconn
    integer            :: vert_nconn
    integer            :: comb_nconn
    integer            :: ieqn
    integer            :: ieqn_1
    integer            :: ieqn_2
    integer            :: icond
    integer            :: ncells_local
    integer            :: mpi_rank

    real(r8), pointer  :: z(:,:)               ! centroid at "z" level [m]
    real(r8), pointer  :: zi(:,:)              ! interface level below a "z" level (m)
    real(r8), pointer  :: dz(:,:)              ! layer thickness at "z" level (m)

    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt , pointer :: soil_filter(:)       ! 

    PetscInt, pointer  :: vert_conn_id_up(:)   !
    PetscInt, pointer  :: vert_conn_id_dn(:)   !
    PetscReal, pointer :: vert_conn_dist_up(:) !
    PetscReal, pointer :: vert_conn_dist_dn(:) !
    PetscReal, pointer :: vert_conn_area(:)    !
    PetscInt , pointer :: vert_conn_type(:)    !

    PetscInt, pointer  :: horz_conn_id_up(:)   !
    PetscInt, pointer  :: horz_conn_id_dn(:)   !
    PetscReal, pointer :: horz_conn_dist_up(:) !
    PetscReal, pointer :: horz_conn_dist_dn(:) !
    PetscReal, pointer :: horz_conn_area(:)    !
    PetscInt, pointer  :: horz_conn_type(:)    !

    PetscInt, pointer  :: comb_conn_id_up(:)   !
    PetscInt, pointer  :: comb_conn_id_dn(:)   !
    PetscReal, pointer :: comb_conn_dist_up(:) !
    PetscReal, pointer :: comb_conn_dist_dn(:) !
    PetscReal, pointer :: comb_conn_area(:)    !
    PetscInt, pointer  :: comb_conn_type(:)    !

    real(r8), pointer  :: xc_col(:)            ! x-position of grid cell [m]
    real(r8), pointer  :: yc_col(:)            ! y-position of grid cell [m]
    real(r8), pointer  :: zc_col(:)            ! z-position of grid cell [m]
    real(r8), pointer  :: area_col(:)          ! area of grid cell [m^2]
    integer, pointer   :: grid_owner(:)        ! MPI rank owner of grid cell

    integer            :: ncells_ghost         ! total number of ghost gridcells on the processor
    integer            :: ncols_ghost

    integer, pointer   :: vsfm_landunit_ind(:,:)
    integer, pointer   :: vsfm_coli(:)
    integer, pointer   :: vsfm_colf(:)

    PetscInt           :: discretization_type  !

    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_type(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)
    integer  , pointer                   :: lun_lakpoi(:)
    integer  , pointer                   :: lun_urbpoi(:)

    integer :: bounds_proc_begg_all, bounds_proc_endg_all
    integer :: bounds_proc_begc_all, bounds_proc_endc_all
    integer :: bounds_proc_begc, bounds_proc_endc

    !----------------------------------------------------------------------

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi             , zi           )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz             , dz           )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_z              , z            )

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active          , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_type            , col_type     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index  , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type       , lun_type     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_lakepoint  , lun_lakpoi   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_urbanpoint , lun_urbpoi   )

    !if (nclumps /= 1) then
    !   call endrun(msg='ERROR VSFM only supported for clumps = 1')
    !endif

    bounds_proc_begg_all = bounds_clump%begg_all
    bounds_proc_endg_all = bounds_clump%endg_all
    bounds_proc_begc_all = bounds_clump%begc_all
    bounds_proc_endc_all = bounds_clump%endc_all
    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    allocate(xc_col            (bounds_proc_begc_all:bounds_proc_endc_all                    ))
    allocate(yc_col            (bounds_proc_begc_all:bounds_proc_endc_all                    ))
    allocate(zc_col            (bounds_proc_begc_all:bounds_proc_endc_all                    ))
    allocate(area_col          (bounds_proc_begc_all:bounds_proc_endc_all                    ))
    allocate(grid_owner        (bounds_proc_begg_all:bounds_proc_endg_all                    ))

    allocate (soil_xc           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_yc           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_zc           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_dx           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_dy           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_dz           ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_area         ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))
    allocate (soil_filter       ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*nlevgrnd     ))

    allocate (vert_conn_id_up   ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_id_dn   ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_dist_up ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_dist_dn ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_area    ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_type    ((bounds_proc_endc_all-bounds_proc_begc_all+1 )*(nlevgrnd-1) ))

    xc_col(:)     = 0.d0
    yc_col(:)     = 0.d0
    zc_col(:)     = 0.d0
    area_col(:)   = 0.d0
    grid_owner(:) = 0

    first_active_soil_col_id = -1
    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       if ((col_active(c) == 1) .and. .not.(lun_lakpoi(l) == 1) .and. .not. (lun_urbpoi(l) == 1)) then
          first_active_soil_col_id = c
          exit
       endif
    end do

    if (first_active_soil_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    xc_col(:)     = 0._r8
    yc_col(:)     = 0._r8
    zc_col(:)     = 0._r8
    grid_owner(:) = 0
    area_col(:)   = 1._r8
    ncols_ghost   = 0
    horz_nconn    = 0

#if 0
    if (lateral_connectivity) then

       call this%vsfm_mpp%GetMPIRank(mpi_rank)

       call update_mesh_information (ncells_ghost, &
            ncols_ghost, xc_col, yc_col, zc_col,   &
            area_col, grid_owner, mpi_rank )

       call setup_lateral_connections (grid_owner, &
            bounds_proc_begg_all, bounds_proc_endg_all, &
            bounds_proc_begc_all, bounds_proc_endc_all, &
            zc_col,                                     &
            horz_nconn,                            &
            horz_conn_id_up, horz_conn_id_dn,      &
            horz_conn_dist_up, horz_conn_dist_dn,  &
            horz_conn_area, horz_conn_type)

    endif
#endif

    ! Save geometric attributes for soil mesh
    icell = 0
    do c = bounds_proc_begc, bounds_proc_endc
       do j = 1, nlevgrnd
          icell = icell + 1
          l = col_landunit(c)

          if ((col_active(c) == 1).and. &
               (lun_type(l) == istsoil .or. col_type(c) == icol_road_perv .or. &
                lun_type(l) == istcrop)) then
             col_id = c
             soil_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             soil_filter(icell) = 0
          end if

          soil_xc(icell)   = xc_col(c)
          soil_yc(icell)   = yc_col(c)
          soil_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j)) + zc_col(c)
          soil_dx(icell)   = 1.d0
          soil_dy(icell)   = 1.d0
          soil_dz(icell)   = dz(col_id,j)
          soil_area(icell) = area_col(c)

       end do
    end do

    ! Save information about internal connections for soil mesh
    iconn = 0
    do c = bounds_proc_begc, bounds_proc_endc

       do j = 1, nlevgrnd-1

          iconn = iconn + 1
          vert_conn_id_up(iconn)   = (c-bounds_proc_begc)*nlevgrnd + j
          vert_conn_id_dn(iconn)   = vert_conn_id_up(iconn) + 1
          vert_conn_dist_up(iconn) = 0.5d0*dz(c,j  ) !zi(c,j)   - z(c,j)
          vert_conn_dist_dn(iconn) = 0.5d0*dz(c,j+1) !z( c,j+1) - zi(c,j)
          vert_conn_area(iconn)    = area_col(c)
          vert_conn_type(iconn)    = CONN_VERTICAL

       end do
    end do
    vert_nconn = iconn

    !
    ! Set up the meshes
    !
    call this%vsfm_mpp%SetNumMeshes(1)

    !
    ! Set mesh value
    !
    imesh        = 1
    nlev         = nlevgrnd
    ncells_local = (bounds_proc_endc - bounds_proc_begc + 1)*nlev
    ncells_ghost = ncols_ghost*nlev

    call this%vsfm_mpp%MeshSetName        (imesh, 'Soil mesh')
    call this%vsfm_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call this%vsfm_mpp%MeshSetID          (imesh, MESH_CLM_SOIL_COL)
    call this%vsfm_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call this%vsfm_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call this%vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call this%vsfm_mpp%MeshComputeVolume          (imesh)

    !
    ! Set connection
    !
    if (vsfm_lateral_model_type == 'none') then
       discretization_type = DISCRETIZATION_VERTICAL_ONLY

       call this%vsfm_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, &
            vert_nconn,  vert_conn_id_up, vert_conn_id_dn, &
            vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area, vert_conn_type)

    else if (vsfm_lateral_model_type == 'source_sink') then
       discretization_type = DISCRETIZATION_VERTICAL_WITH_SS

       call this%vsfm_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, &
            vert_nconn,  vert_conn_id_up, vert_conn_id_dn, &
            vert_conn_dist_up, vert_conn_dist_dn, vert_conn_area, vert_conn_type)

       call this%vsfm_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_LATERAL, &
            horz_nconn,  horz_conn_id_up, horz_conn_id_dn, &
            horz_conn_dist_up, horz_conn_dist_dn, horz_conn_area, horz_conn_type)

    else if (vsfm_lateral_model_type == 'three_dimensional') then
       discretization_type = DISCRETIZATION_THREE_DIM

       comb_nconn = vert_nconn + horz_nconn

       allocate (comb_conn_id_up   (comb_nconn))
       allocate (comb_conn_id_dn   (comb_nconn))
       allocate (comb_conn_dist_up (comb_nconn))
       allocate (comb_conn_dist_dn (comb_nconn))
       allocate (comb_conn_area    (comb_nconn))
       allocate (comb_conn_type    (comb_nconn))

       comb_conn_id_up   (1:vert_nconn) = vert_conn_id_up   (1:vert_nconn)
       comb_conn_id_dn   (1:vert_nconn) = vert_conn_id_dn   (1:vert_nconn)
       comb_conn_dist_up (1:vert_nconn) = vert_conn_dist_up (1:vert_nconn)
       comb_conn_dist_dn (1:vert_nconn) = vert_conn_dist_dn (1:vert_nconn)
       comb_conn_area    (1:vert_nconn) = vert_conn_area    (1:vert_nconn)
       comb_conn_type    (1:vert_nconn) = vert_conn_type    (1:vert_nconn)

       comb_conn_id_up   (vert_nconn+1:comb_nconn) = horz_conn_id_up   (1:horz_nconn)
       comb_conn_id_dn   (vert_nconn+1:comb_nconn) = horz_conn_id_dn   (1:horz_nconn)
       comb_conn_dist_up (vert_nconn+1:comb_nconn) = horz_conn_dist_up (1:horz_nconn)
       comb_conn_dist_dn (vert_nconn+1:comb_nconn) = horz_conn_dist_dn (1:horz_nconn)
       comb_conn_area    (vert_nconn+1:comb_nconn) = horz_conn_area    (1:horz_nconn)
       comb_conn_type    (vert_nconn+1:comb_nconn) = horz_conn_type    (1:horz_nconn)

       call this%vsfm_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, &
            comb_nconn,  comb_conn_id_up, comb_conn_id_dn, &
            comb_conn_dist_up, comb_conn_dist_dn, comb_conn_area, comb_conn_type)

       deallocate (comb_conn_id_up  )
       deallocate (comb_conn_id_dn  )
       deallocate (comb_conn_dist_up)
       deallocate (comb_conn_dist_dn)
       deallocate (comb_conn_area   )
       deallocate (comb_conn_type   )

    else
       call endrun(msg='ERROR: ' // &
            'Unknown vsfm_lateral_model_type = ' // trim(vsfm_lateral_model_type) // &
            errMsg(__FILE__, __LINE__))
    endif

    ! Free up memory
    deallocate(xc_col            )
    deallocate(yc_col            )
    deallocate(zc_col            )
    deallocate(area_col          )
    deallocate(grid_owner        )

    deallocate (soil_xc          )
    deallocate (soil_yc          )
    deallocate (soil_zc          )
    deallocate (soil_dx          )
    deallocate (soil_dy          )
    deallocate (soil_dz          )
    deallocate (soil_area        )
    deallocate (soil_filter      )

    deallocate (vert_conn_id_up  )
    deallocate (vert_conn_id_dn  )
    deallocate (vert_conn_dist_up)
    deallocate (vert_conn_dist_dn)
    deallocate (vert_conn_area   )
    deallocate (vert_conn_type   )

    if (lateral_connectivity) then
       deallocate (horz_conn_id_up   )
       deallocate (horz_conn_id_dn   )
       deallocate (horz_conn_dist_up )
       deallocate (horz_conn_dist_dn )
       deallocate (horz_conn_area    )
       deallocate (horz_conn_type    )
    endif

  end subroutine add_meshes

  !------------------------------------------------------------------------

#if 0
  subroutine update_mesh_information (ncells_ghost, ncols_ghost,  &
       xc_col, yc_col, zc_col, area_col, &
       grid_owner, mpi_rank &
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use decompMod                 , only : bounds_type, get_proc_bounds
    use domainMod                 , only : ldomain
    use domainLateralMod          , only : ldomain_lateral
    use mpp_varpar                , only : nlevgrnd
    use mpp_varcon                , only : istsoil
    use initGridCellsMod          , only : initGhostGridCells
    use decompMod                 , only : get_proc_total_ghosts
    use UnstructuredGridType      , only : ScatterDataG2L
    use clm_instMod               , only : soilstate_vars
    use mpp_bounds                , only : bounds_proc_begc_all, bounds_proc_endc_all
    use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    use mpp_bounds                , only : bounds_proc_begg_all, bounds_proc_endg_all
    use mpp_bounds                , only : bounds_proc_begg, bounds_proc_endg
    !
    implicit none
    !
    integer            :: ncols_ghost    ! number of ghost columns
    real(r8), pointer  :: xc_col(:)      ! x-position of grid cell [m]
    real(r8), pointer  :: yc_col(:)      ! y-position of grid cell [m]
    real(r8), pointer  :: zc_col(:)      ! z-position of grid cell [m]
    real(r8), pointer  :: area_col(:)    ! area of grid cell [m^2]
    integer, pointer   :: grid_owner(:)  ! MPI rank owner of grid cell
    integer, intent(in):: mpi_rank       ! MPI rank
    !
    type(bounds_type)  :: bounds_proc
    integer            :: c,g,fc,j,l     ! do loop indices
    integer            :: nblocks
    integer            :: beg_idx
    integer            :: ndata_send     ! number of data sent by local mpi rank
    integer            :: ndata_recv     ! number of data received by local mpi rank
    real(r8), pointer  :: data_send(:)   ! data sent by local mpi rank
    real(r8), pointer  :: data_recv(:)   ! data received by local mpi rank

    integer            :: ncells_ghost   ! total number of ghost gridcells on the processor
    integer            :: nlunits_ghost  ! total number of ghost landunits on the processor
    integer            :: npfts_ghost    ! total number of ghost pfts on the processor
    integer            :: nCohorts_ghost ! total number of ghost cohorts on the processor
    !----------------------------------------------------------------------

    call initGhostGridCells()

    call get_proc_bounds(bounds_proc)
    call soilstate_vars%InitColdGhost(bounds_proc)

    nblocks    = 4
    ndata_send = nblocks*ldomain_lateral%ugrid%ngrid_local
    ndata_recv = nblocks*ldomain_lateral%ugrid%ngrid_ghosted

    allocate(data_send(ndata_send))
    allocate(data_recv(ndata_recv))

    ! Aggregate the data to send
    do g = bounds_proc_begg, bounds_proc_endg

       beg_idx = (g-bounds_proc_begg)*nblocks

       if (isnan(ldomain%xCell(g))) then
          call endrun(msg='ERROR initialize3: xCell = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%xCell(g)

       if (isnan(ldomain%yCell(g))) then
          call endrun(msg='ERROR initialize3: yCell = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%yCell(g)

       if (isnan(ldomain%topo(g))) then
          call endrun(msg='ERROR initialize3: topo = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%topo(g)

       beg_idx = beg_idx + 1;
       data_send(beg_idx) = real(mpi_rank)

    enddo

    ! Scatter: Global-to-Local
    call ScatterDataG2L(ldomain_lateral%ugrid, nblocks, &
         ndata_send, data_send, ndata_recv, data_recv)

    ! Save data for ghost subgrid category
    do c = bounds_proc_begc_all, bounds_proc_endc_all

       g       = col%gridcell(c)
       beg_idx = (g-bounds_proc_begg)*nblocks

       beg_idx = beg_idx + 1; xc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; yc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; zc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; grid_owner(g) = data_recv(beg_idx)

       area_col(c) = ldomain_lateral%ugrid%areaGrid_ghosted(g-bounds_proc_begg + 1)

    enddo

    deallocate(data_send)
    deallocate(data_recv)

    call get_proc_total_ghosts(ncells_ghost, nlunits_ghost, &
         ncols_ghost, npfts_ghost, nCohorts_ghost)

  end subroutine update_mesh_information

  !------------------------------------------------------------------------

  subroutine setup_lateral_connections (grid_owner, &
       begg, endg, begc, endc, zc_col,              &
       nconn_horz,                                  &
       horz_conn_id_up, horz_conn_id_dn,            &
       horz_conn_dist_up, horz_conn_dist_dn,        &
       horz_conn_area, horz_conn_type               &
       )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    use domainLateralMod          , only : ldomain_lateral
    use mpp_varpar                , only : nlevgrnd
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use mpp_varcon                , only : istsoil
    !
    implicit none
    !
    integer, pointer   :: grid_owner(:)        ! MPI rank owner of grid cell
    integer            :: begg,endg
    integer            :: begc,endc
    real(r8), pointer  :: zc_col(:)            ! z-position of grid cell [m]
    PetscInt           :: nconn_horz
    PetscInt, pointer  :: horz_conn_type(:)    !
    PetscInt, pointer  :: horz_conn_id_up(:)   !
    PetscInt, pointer  :: horz_conn_id_dn(:)   !
    PetscReal, pointer :: horz_conn_dist_up(:) !
    PetscReal, pointer :: horz_conn_dist_dn(:) !
    PetscReal, pointer :: horz_conn_area(:)    !
    !
    ! !LOCAL VARIABLES:
    PetscInt           :: c,j,l                !indices
    PetscInt           :: icell
    PetscInt           :: iconn
    PetscInt           :: nconn
    PetscInt           :: id_up, id_dn
    PetscInt           :: first_active_hydro_col_id
    PetscInt           :: col_id
    PetscInt           :: g_up, g_dn, iedge
    PetscInt           :: c_idx_up, c_idx_dn
    PetscInt           :: l_idx_up, l_idx_dn
    PetscInt           :: ltype, ctype
    PetscInt           :: tmp
    PetscReal          :: dist_x, dist_y, dist_z, dist
    PetscReal          :: dc, dv
    PetscErrorCode     :: ierr
    real(r8), pointer  :: dz(:,:)              ! layer thickness at "z" level (m)
    !----------------------------------------------------------------------

    dz         =>    col%dz

    !
    ! Sets up lateral connection between columns of type 'istsoil'
    !
    ! Assumptions:
    ! - There is only ONE 'istsoil' column per landunit per grid cell.
    ! - Grid cells that are laterally connected via cellsOnCell
    !   field defined in domain netcdf file have at least ONE
    !   column of 'istsoil' type

    nconn_horz = 0

    ltype = istsoil
    ctype = istsoil

    ! Determine number of lateral connections

    do icell = 1, ldomain_lateral%ugrid%ngrid_local
       do iedge = 1, ldomain_lateral%ugrid%maxEdges

          if (ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) > icell) then
             g_up = icell + begg - 1
             g_dn = ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

             l_idx_up = grc%landunit_indices(ltype, g_up)
             l_idx_dn = grc%landunit_indices(ltype, g_dn)

             c_idx_up = -1
             c_idx_dn = -1

             do c = lun_pp%coli(l_idx_up), lun_pp%colf(l_idx_up)
                if (col%itype(c) == ctype) then
                   if (c_idx_up /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_up = c
                endif
             enddo

             do c = lun_pp%coli(l_idx_dn), lun_pp%colf(l_idx_dn)
                if (col%itype(c) == ctype) then
                   if (c_idx_dn /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_dn = c
                endif
             enddo

             if (c_idx_up > -1 .and. c_idx_dn > -1) then
                nconn_horz = nconn_horz + 1
             else
                write(iulog,*)'CreateFromCLMCols: No column of ctype = ', ctype, &
                     ' found between following grid cells: ',g_up,g_dn
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif
          endif
       enddo
    enddo

    nconn_horz = nconn_horz * nlevgrnd

    allocate (horz_conn_id_up   (nconn_horz))
    allocate (horz_conn_id_dn   (nconn_horz))
    allocate (horz_conn_dist_up (nconn_horz))
    allocate (horz_conn_dist_dn (nconn_horz))
    allocate (horz_conn_area    (nconn_horz))
    allocate (horz_conn_type    (nconn_horz))

    iconn = 0
    do icell = 1, ldomain_lateral%ugrid%ngrid_local

       do iedge = 1, ldomain_lateral%ugrid%maxEdges
          if (ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) > icell) then
             g_up = icell + begg - 1
             g_dn = ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

             l_idx_up = grc%landunit_indices(ltype, g_up)
             l_idx_dn = grc%landunit_indices(ltype, g_dn)

             c_idx_up = -1
             c_idx_dn = -1

             do c = lun_pp%coli(l_idx_up), lun_pp%colf(l_idx_up)
                if (col%itype(c) == ctype) then
                   if (c_idx_up /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_up = c
                endif
             enddo

             do c = lun_pp%coli(l_idx_dn), lun_pp%colf(l_idx_dn)
                if (col%itype(c) == ctype) then
                   if (c_idx_dn /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_dn = c
                endif
             enddo

             if (c_idx_up > -1 .and. c_idx_dn > -1) then

                if (grid_owner(g_up) > grid_owner(g_dn)) then
                   tmp      = g_up;
                   g_up     = g_dn
                   g_dn     = tmp

                   tmp      = l_idx_up;
                   l_idx_up = l_idx_dn
                   l_idx_dn = tmp

                   tmp      = c_idx_up;
                   c_idx_up = c_idx_dn
                   c_idx_dn = tmp
                endif

                dc = ldomain_lateral%ugrid%dcOnGrid_local(iedge, icell)
                dv = ldomain_lateral%ugrid%dvOnGrid_local(iedge, icell)

                do j = 1, nlevgrnd
                   iconn = iconn + 1

                   id_up = (c_idx_up - begc)*nlevgrnd + j
                   id_dn = (c_idx_dn - begc)*nlevgrnd + j

                   horz_conn_type         = CONN_HORIZONTAL
                   horz_conn_id_up(iconn) = id_up
                   horz_conn_id_dn(iconn) = id_dn

                   !this%is_active(id_up) = PETSC_TRUE
                   !this%is_active(id_dn) = PETSC_TRUE

                   horz_conn_area(iconn) = dz(c_idx_up,j)*dv
                   dist = (dc**2.d0 + (zc_col(c_idx_up) - zc_col(c_idx_dn))**2.d0)**0.5d0

                   horz_conn_dist_up(iconn) = 0.5d0*dist
                   horz_conn_dist_dn(iconn) = 0.5d0*dist

                enddo
             endif ! if (c_idx_up > -1 .and. c_idx_dn > -1)
          endif ! if (ugrid%gridsOnGrid_local(iedge,icell) > icell)
       enddo
    enddo

  end subroutine setup_lateral_connections
#endif

  !------------------------------------------------------------------------

  subroutine add_goveqns(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none
    !
    class(em_vsfm_type)                  :: this

    call this%vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

    call this%vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use mpp_varctl                , only : vsfm_lateral_model_type
    use mpp_varctl                , only : vsfm_include_seepage_bc
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : COND_SEEPAGE_BC
    !
    ! !ARGUMENTS
    implicit none
    !
    class(em_vsfm_type)                  :: this
    PetscInt :: ieqn

    ieqn = 1

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Infiltration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Evapotranspiration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Dew_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Drainage_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Snow_Disappearance_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Sublimation_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    if (vsfm_lateral_model_type == 'source_sink' ) then

       call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
            'Lateral_flux', 'kg/s', COND_MASS_RATE, &
            SOIL_CELLS)

       if (vsfm_include_seepage_bc) then
          call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
               'Seepage_Flux', 'kg/s', COND_SEEPAGE_BC, &
               SOIL_TOP_CELLS)
       endif

    else if (vsfm_lateral_model_type == 'three_dimensional') then

       if (vsfm_include_seepage_bc) then
          call this%vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
               'Seepage_Flux', 'kg/s', COND_SEEPAGE_BC, &
               SOIL_TOP_CELLS)
       endif

    endif

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars(this)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    class(em_vsfm_type)              :: this

    !
    ! Allocate auxvars
    !
    call this%vsfm_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    use clm_instMod               , only : soilstate_vars
    use clm_instMod               , only : soilhydrology_vars
    use mpp_varcon                , only : istcrop
    use mpp_varcon                , only : istsoil
    use mpp_varcon                , only : icol_road_perv
    use mpp_varpar                , only : nlevgrnd
    use mpp_varctl                , only : vsfm_satfunc_type
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use EOSWaterMod               , only : DENSITY_TGDPB01
    !use mpp_bounds                , only : bounds_proc_begc_all, bounds_proc_endc_all
    !use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type)              :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer    :: clm_watsat(:,:)
    real(r8), pointer    :: clm_hksat(:,:)
    real(r8), pointer    :: clm_bsw(:,:)
    real(r8), pointer    :: clm_sucsat(:,:)
    real(r8), pointer    :: clm_eff_porosity(:,:)
    real(r8), pointer    :: clm_zwt(:)
    real(r8), pointer    :: vsfm_watsat(:,:)
    real(r8), pointer    :: vsfm_hksat(:,:)
    real(r8), pointer    :: vsfm_bsw(:,:)
    real(r8), pointer    :: vsfm_sucsat(:,:)
    real(r8), pointer    :: vsfm_eff_porosity(:,:)
    real(r8), pointer    :: vsfm_residual_sat(:,:)
    integer, pointer     :: vsfm_filter(:)
    !
    integer              :: c,g,fc,j,l            ! do loop indices
    integer              :: ncells_ghost          ! total number of ghost gridcells on the processor
    integer              :: nlunits_ghost         ! total number of ghost landunits on the processor
    integer              :: ncols_ghost           ! total number of ghost columns on the processor
    integer              :: npfts_ghost           ! total number of ghost pfts on the processor
    integer              :: nCohorts_ghost        ! total number of ghost cohorts on the processor
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_type(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)

    integer :: bounds_proc_begc_all, bounds_proc_endc_all
    integer :: bounds_proc_begc, bounds_proc_endc

    !-----------------------------------------------------------------------

    bounds_proc_begc_all = bounds_clump%begc_all
    bounds_proc_endc_all = bounds_clump%endc_all
    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_watsatc      , clm_watsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_hksatc       , clm_hksat        )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_bswc         , clm_bsw          )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_sucsatc      , clm_sucsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_effporosityc , clm_eff_porosity )
    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_state_wtd         , clm_zwt          )

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active              , col_active       )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_type                , col_type         )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index      , col_landunit     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type           , lun_type         )

    ! Allocate memory
    allocate(vsfm_filter       (bounds_proc_begc_all:bounds_proc_endc_all           ))
    allocate(vsfm_watsat       (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(vsfm_hksat        (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(vsfm_bsw          (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(vsfm_sucsat       (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(vsfm_eff_porosity (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(vsfm_residual_sat (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))

    ! Initialize
    vsfm_filter       (:)   = 0
    vsfm_watsat       (:,:) = 0._r8
    vsfm_hksat        (:,:) = 0._r8
    vsfm_bsw          (:,:) = 0._r8
    vsfm_sucsat       (:,:) = 0._r8
    vsfm_residual_sat (:,:) = 0._r8

    ! Save data to initialize VSFM
    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       if ((col_active(c) == 1) .and. &
           (lun_type(l) == istsoil .or. col_type(c) == icol_road_perv .or. &
           lun_type(l) == istcrop)) then

          vsfm_filter    (c) = 1

          do j = 1 ,nlevgrnd
             vsfm_watsat(c,j)       = clm_watsat(c,j)
             vsfm_hksat(c,j)        = clm_hksat(c,j)
             vsfm_bsw(c,j)          = clm_bsw(c,j)
             vsfm_sucsat(c,j)       = clm_sucsat(c,j)
             vsfm_eff_porosity(c,j) = clm_eff_porosity(c,j)
          enddo

       endif
    enddo

    ncols_ghost = 0

    call VSFMMPPSetSoils(this%vsfm_mpp, bounds_proc_begc, bounds_proc_endc, &
         ncols_ghost, vsfm_filter, &
         vsfm_watsat, vsfm_hksat, vsfm_bsw, vsfm_sucsat, &
         vsfm_residual_sat, vsfm_satfunc_type, DENSITY_TGDPB01)

    ! Free up memory
    deallocate(vsfm_filter       )
    deallocate(vsfm_watsat       )
    deallocate(vsfm_hksat        )
    deallocate(vsfm_bsw          )
    deallocate(vsfm_sucsat       )
    deallocate(vsfm_eff_porosity )
    deallocate(vsfm_residual_sat )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    use clm_instMod               , only : soilstate_vars
    use clm_instMod               , only : soilhydrology_vars
    use mpp_varcon                , only : istcrop
    use mpp_varcon                , only : istsoil
    use mpp_varcon                , only : icol_road_perv
    use mpp_varpar                , only : nlevgrnd
    use mpp_varctl                , only : vsfm_satfunc_type
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use MultiPhysicsProbConstants , only : PRESSURE_REF
    !use mpp_bounds                , only : bounds_proc_begc_all, bounds_proc_endc_all
    !use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type)              :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer    :: clm_zi(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer    :: clm_zwt(:)            !
    integer              :: c,g,fc,j,l            ! do loop indices
    integer              :: ncells_ghost          ! total number of ghost gridcells on the processor
    integer              :: nlunits_ghost         ! total number of ghost landunits on the processor
    integer              :: ncols_ghost           ! total number of ghost columns on the processor
    integer              :: npfts_ghost           ! total number of ghost pfts on the processor
    integer              :: nCohorts_ghost        ! total number of ghost cohorts on the processor
    real(r8), pointer    :: press_ic_1d(:)        ! pressure initial condition (m)
    integer :: icell
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_type(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)
    integer :: bounds_proc_begc_all, bounds_proc_endc_all
    integer :: bounds_proc_begc, bounds_proc_endc
    !-----------------------------------------------------------------------

    bounds_proc_begc_all = bounds_clump%begc_all
    bounds_proc_endc_all = bounds_clump%endc_all
    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_state_wtd    , clm_zwt      )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi            , clm_zi       )

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active         , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_type           , col_type     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type      , lun_type     )

    ! Allocate memory
    allocate(press_ic_1d ((bounds_proc_endc_all - bounds_proc_begc_all + 1)*nlevgrnd))

    ! Initialize
    press_ic_1d(:) = 101325.d0

    ! Save data to initialize VSFM
    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       do j = 1, nlevgrnd
          icell = (c - bounds_proc_begc)*nlevgrnd + j
          if ((col_active(c) == 1) .and. &
               (lun_type(l) == istsoil .or. col_type(c) == icol_road_perv .or. &
               lun_type(l) == istcrop)) then

             press_ic_1d(icell) = PRESSURE_REF + &
                  997.16d0*GRAVITY_CONSTANT * &
                  (-clm_zwt(c) - (-0.5d0*(clm_zi(c,j-1) + clm_zi(c,j))))
          endif
       enddo
    enddo

    !call VSFMMPPSetICs(vsfm_mpp, bounds_proc_begc, bounds_proc_endc, &
    !zc_col, vsfm_filter, vsfm_zwt)
    call this%vsfm_mpp%Restart(press_ic_1d)

    ! Free up memory
    deallocate(press_ic_1d)

  end subroutine set_initial_conditions

  !-----------------------------------------------------------------------
  subroutine determine_condition_ids(this)
    !
    !DESCRIPTION
    !  Determines the IDs of various source-sink conditions in VSFM
    !
    use mpp_varctl                       , only : vsfm_lateral_model_type
    use MultiPhysicsProbConstants        , only : COND_SS
    use MultiPhysicsProbConstants        , only : COND_NULL
    use mpp_varctl                       , only : iulog
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
    use SystemOfEquationsBaseType        , only : sysofeqns_base_type
    use SystemOfEquationsVSFMType        , only : sysofeqns_vsfm_type
    ! !ARGUMENTS:
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type)          :: this
    integer                      :: ier ! error status
    !
    ! !LOCAL VARIABLES:
    character (len=256), pointer :: cond_names(:)
    integer                      :: num_conds
    integer                      :: num_conds_expected
    integer                      :: nn
    integer                      :: kk
    character (len=256)          :: cond_name
    class(sysofeqns_base_type), pointer  :: soe
    !------------------------------------------------------------------------------

    this%vsfm_cond_id_for_infil        = -1
    this%vsfm_cond_id_for_et           = -1
    this%vsfm_cond_id_for_dew          = -1
    this%vsfm_cond_id_for_drainage     = -1
    this%vsfm_cond_id_for_snow         = -1
    this%vsfm_cond_id_for_sublimation  = -1
    this%vsfm_cond_id_for_lateral_flux = -1

    num_conds_expected = 6

    if (vsfm_lateral_model_type == 'source_sink' ) then
       num_conds_expected = num_conds_expected + 1
    end if

    soe => this%vsfm_mpp%soe
    select type(soe)
    class is(sysofeqns_vsfm_type)
       ! Get the number of conditions
       call soe%GetConditionNames(COND_SS, COND_NULL, num_conds, cond_names)
    end select

    if (num_conds /= num_conds_expected) then
      write(iulog,*)'In init_vsfm_condition_ids: Source-sink conditions /= ', num_conds_expected
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do nn = 1, num_conds
       select case(trim(cond_names(nn)))

       case ("Infiltration_Flux")
          this%vsfm_cond_id_for_infil        = nn

       case ("Evapotranspiration_Flux")
          this%vsfm_cond_id_for_et           = nn

       case ("Dew_Flux")
          this%vsfm_cond_id_for_dew          = nn

       case ("Drainage_Flux")
          this%vsfm_cond_id_for_drainage     = nn

       case ("Snow_Disappearance_Flux")
          this%vsfm_cond_id_for_snow         = nn

       case ("Sublimation_Flux")
          this%vsfm_cond_id_for_sublimation  = nn

       case ("Lateral_flux")
          this%vsfm_cond_id_for_lateral_flux = nn

       case default
          write(iulog,*) trim(cond_names(nn))
          write(iulog,*)'In init_vsfm_condition_ids: Unknown flux.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    enddo

    if (this%vsfm_cond_id_for_infil == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_infil not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_et == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_et not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_dew == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_dew not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_drainage == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_drainage not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_snow == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_snow not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_sublimation == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_sublimation not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_lateral_model_type == 'source_sink') then
       if (this%vsfm_cond_id_for_lateral_flux == -1) then
          write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_lateral_flux not defined.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

    endif

    deallocate(cond_names)

  end subroutine determine_condition_ids

  !-----------------------------------------------------------------------
  subroutine extract_data_for_alm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use mpp_varctl                , only : restart_vsfm
    !use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    use mpp_varpar                , only : nlevgrnd
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g  ! do loop indices

    real(r8)  , pointer                  :: l2e_col_zi(:,:)

    real(r8)  , pointer                  :: l2e_soilp(:,:)
    real(r8)  , pointer                  :: l2e_mflx_snowlyr_col(:)

    real(r8)  , pointer                  :: vsfm_soilp_col_1d(:)
    real(r8)  , pointer                  :: vsfm_mass_col_1d(:)
    real(r8)  , pointer                  :: vsfm_smpl_col_1d(:)
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp_l(:,:)
    real(r8)  , pointer                  :: e2l_zwt(:)
    real(r8)  , pointer                  :: e2l_mflx_snowlyr_col(:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    integer                              :: jwt
    integer                              :: idx
    integer                              :: soe_auxvar_id
    real(r8)                             :: z_up, z_dn
    integer :: bounds_proc_begc, bounds_proc_endc
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_flux_mflx_snowlyr_col , l2e_mflx_snowlyr_col )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_state_soilp           , l2e_soilp            )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi                , l2e_col_zi           )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_liq            , l2e_h2osoi_liq       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_smp             , e2l_smp_l            )

    ! PreSolve: Allows saturation value to be computed based on ICs and stored
    !           in GE auxvar
    call this%vsfm_mpp%soe%SetDtime(1.d0)
    call this%vsfm_mpp%soe%PreSolve()

    ! PostSolve: Allows saturation value stored in GE auxvar to be copied into
    !            SoE auxvar
    call this%vsfm_mpp%soe%PostSolve()

    allocate(vsfm_soilp_col_1d((bounds_proc_endc-bounds_proc_begc+1)*nlevgrnd))
    allocate(vsfm_mass_col_1d ((bounds_proc_endc-bounds_proc_begc+1)*nlevgrnd))
    allocate(vsfm_smpl_col_1d ((bounds_proc_endc-bounds_proc_begc+1)*nlevgrnd))

    if (restart_vsfm) then

       ! Save 1D array for VSFM (vsfm_soilp_col_1d) and
       ! set initial value of mflx_snowlyr_col for ALM
       do c = bounds_proc_begc, bounds_proc_endc
          do j = 1, nlevgrnd
             idx = (c - bounds_proc_begc)*nlevgrnd + j
             vsfm_soilp_col_1d(idx) = l2e_soilp(c,j)
          end do
          idx = c-bounds_proc_begc+1
          e2l_mflx_snowlyr_col(c) = l2e_mflx_snowlyr_col(c)
       end do

       ! Set the initial conditions
       call this%vsfm_mpp%Restart(vsfm_soilp_col_1d)

       ! PreSolve: Allows saturation value to be computed based on ICs and stored
       !           in GE auxvar
       call this%vsfm_mpp%soe%SetDtime(1.d0)
       call this%vsfm_mpp%soe%PreSolve()

       ! PostSolve: Allows saturation value stored in GE auxvar to be copied into
       !            SoE auxvar
       call this%vsfm_mpp%soe%PostSolve()

    else
       ! Set initial value of mflx_snowlyr_col for ALM
       e2l_mflx_snowlyr_col(:) = 0._r8
    end if

    ! Get total mass
    soe_auxvar_id = 1;
    call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,   &
                                          VAR_MASS,          &
                                          soe_auxvar_id,     &
                                          vsfm_mass_col_1d)

    ! Get liquid soil matrix potential
    soe_auxvar_id = 1;
    call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,       &
                                          VAR_SOIL_MATRIX_POT,   &
                                          soe_auxvar_id,         &
                                          vsfm_smpl_col_1d)

    do c = bounds_proc_begc, bounds_proc_endc
       ! initialization
       jwt = -1

       ! Loops in decreasing j so WTD can be computed in the same loop
       do j = nlevgrnd, 1, -1
          idx = (c-bounds_proc_begc)*nlevgrnd + j

          if (restart_vsfm) then
             e2l_h2osoi_liq(c,j) = l2e_h2osoi_liq(c,j)
             e2l_h2osoi_ice(c,j) = l2e_h2osoi_ice(c,j)
          else
             e2l_h2osoi_liq(c,j) = vsfm_mass_col_1d(idx)
             e2l_h2osoi_ice(c,j) = 0.d0
          end if
          e2l_smp_l(c,j)      = vsfm_smpl_col_1d(idx)*1000._r8      ! [m] --> [mm]

          if (jwt == -1) then
             ! Find the first soil that is unsaturated
             if (e2l_smp_l(c,j) < 0._r8) jwt = j
          end if

       end do

       if (jwt == -1 .or. jwt == nlevgrnd) then
          ! Water table below or in the last layer
          e2l_zwt(c) = l2e_col_zi(c,nlevgrnd)
       else
          z_dn = (l2e_col_zi(c,jwt-1) + l2e_col_zi(c,jwt  ))/2._r8
          z_up = (l2e_col_zi(c,jwt ) + l2e_col_zi(c,jwt+1))/2._r8
          e2l_zwt(c) = (0._r8 - e2l_smp_l(c,jwt))/(e2l_smp_l(c,jwt) - e2l_smp_l(c,jwt+1))*(z_dn - z_up) + z_dn
        endif
     enddo

   end subroutine extract_data_for_alm

    !------------------------------------------------------------------------
   subroutine EM_VSFM_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The VSFM dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    select case (em_stage)
    case (EM_VSFM_SOIL_HYDRO_STAGE)
       call EM_VSFM_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case default
       write(iulog,*)'EM_FATES_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_VSFM_Solve

    !------------------------------------------------------------------------
  subroutine EM_VSFM_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Solve the Variably Saturated Flow Model (VSFM) in soil columns.
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use mpp_varpar                , only : nlevgrnd
    !use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    use petscsnes
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g                                                       ! do loop indices
    integer                              :: pi                                                               ! pft index
    real(r8)                             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)                             :: dtime

    real(r8)  , pointer                  :: mflx_et_col_1d         (:)
    real(r8)  , pointer                  :: mflx_infl_col_1d       (:)
    real(r8)  , pointer                  :: mflx_dew_col_1d        (:)
    real(r8)  , pointer                  :: mflx_drain_col_1d      (:)
    real(r8)  , pointer                  :: mflx_sub_snow_col_1d   (:)
    real(r8)  , pointer                  :: mflx_snowlyr_col_1d    (:)
    real(r8)  , pointer                  :: t_soil_col_1d          (:)

    real(r8)  , pointer                  :: vsfm_fliq_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_mass_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_smpl_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_soilp_col_1d      (:)
    real(r8)  , pointer                  :: vsfm_sat_col_1d        (:)

    real(r8)  , pointer                  :: frac_ice                    (:,:) ! fraction of ice
    real(r8)  , pointer                  :: total_mass_flux_col         (:)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_et_col      (:)            ! ET sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_infl_col    (:)            ! Infiltration source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_dew_col     (:)            ! Dew source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_drain_col   (:)            ! Drainage sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_snowlyr_col (:)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_sub_col     (:)            ! Sublimation sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_lateral_col (:)            ! Lateral flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_seepage_col (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: qflx_seepage                (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: vsfm_mass_prev_col          (:,:) ! Mass of water before a VSFM solve
    real(r8)  , pointer                  :: vsfm_dmass_col              (:)            ! Change in mass of water after a VSFM solve
    real(r8)  , pointer                  :: mass_beg_col                (:)            ! Total mass before a VSFM solve
    real(r8)  , pointer                  :: mass_end_col                (:)            ! Total mass after a VSFM solve
    integer                              :: ier                                                              ! error status

    integer                              :: begc, endc
    integer                              :: idx
    real(r8)                             :: area

    PetscInt                             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscErrorCode                       :: ierr                                                             ! PETSc return error code

    PetscBool                            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscInt                             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal                            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal                            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal                            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt                             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt                             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal                            :: stol                                                             ! solution convergance tolerance
    PetscReal                            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter                  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal                            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal                            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal                            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal                            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt                             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter                  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt                             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt                             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal                            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter                 :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool                            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal , pointer                  :: vsfm_soilp_col_ghosted_1d(:)
    PetscReal , pointer                  :: vsfm_fliq_col_ghosted_1d(:)
    PetscReal , pointer                  :: mflx_lateral_col_1d(:)
    PetscReal , pointer                  :: lat_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_press_1d(:)

    integer                              :: jwt
    real(r8)                             :: z_dn, z_up

    real(r8)  , pointer                  :: l2e_mflux_infil(:)
    real(r8)  , pointer                  :: l2e_mflux_dew(:)
    real(r8)  , pointer                  :: l2e_mflux_sub_snow(:)
    real(r8)  , pointer                  :: l2e_mflux_snowlyr(:)
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_mflux_drain(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_zi(:,:)
    integer   , pointer                  :: l2e_filter_hydrologyc(:)
    integer                              :: l2e_num_hydrologyc

    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp(:,:)
    real(r8)  , pointer                  :: e2l_wtd(:)
    real(r8)  , pointer                  :: e2l_soilp(:,:)
    real(r8)  , pointer                  :: e2l_qrecharge(:)

    integer :: bounds_proc_begc, bounds_proc_endc
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

      ! Get time step

      dtime = dt

      call l2e_list%GetPointerToReal1D(this%index_l2e_flux_infil       , l2e_mflux_infil       )
      call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dew         , l2e_mflux_dew         )
      call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snow_sub    , l2e_mflux_sub_snow    )
      call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snowlyr     , l2e_mflux_snowlyr     )

      call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et          , l2e_mflux_et          )
      call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drainage    , l2e_mflux_drain       )
      call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq , l2e_h2osoi_liq        )
      call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice , l2e_h2osoi_ice        )

      call l2e_list%GetPointerToInt1D(this%index_l2e_filter_hydrologyc , l2e_filter_hydrologyc )
      call l2e_list%GetIntValue(this%index_l2e_filter_num_hydrologyc   , l2e_num_hydrologyc    )

      call l2e_list%GetPointerToReal2D(this%index_l2e_column_zi        , l2e_zi                )

      call e2l_list%GetPointerToReal1D(this%index_e2l_state_wtd        , e2l_wtd               )
      call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq        )
      call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice        )
      call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp               )
      call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp             )

      call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrecharge   , e2l_qrecharge         )

      begc = bounds_proc_begc
      endc = bounds_proc_endc

      allocate(frac_ice                    (begc:endc,1:nlevgrnd))
      allocate(total_mass_flux_col         (begc:endc))
      allocate(total_mass_flux_et_col      (begc:endc))
      allocate(total_mass_flux_infl_col    (begc:endc))
      allocate(total_mass_flux_dew_col     (begc:endc))
      allocate(total_mass_flux_drain_col   (begc:endc))
      allocate(total_mass_flux_snowlyr_col (begc:endc))
      allocate(total_mass_flux_sub_col     (begc:endc))
      allocate(total_mass_flux_lateral_col (begc:endc))
      allocate(total_mass_flux_seepage_col (begc:endc))
      allocate(qflx_seepage                (begc:endc))
      allocate(vsfm_mass_prev_col          (begc:endc,1:nlevgrnd))
      allocate(vsfm_dmass_col              (begc:endc))
      allocate(mass_beg_col                (begc:endc))
      allocate(mass_end_col                (begc:endc))

      allocate(mflx_et_col_1d              ((endc-begc+1)*nlevgrnd))
      allocate(mflx_drain_col_1d           ((endc-begc+1)*nlevgrnd))
      allocate(mflx_infl_col_1d            (endc-begc+1))
      allocate(mflx_dew_col_1d             (endc-begc+1))
      allocate(mflx_sub_snow_col_1d        (endc-begc+1))
      allocate(mflx_snowlyr_col_1d         (endc-begc+1))
      allocate(t_soil_col_1d               ((endc-begc+1)*nlevgrnd))

      allocate(vsfm_mass_col_1d            ((endc-begc+1)*nlevgrnd))
      allocate(vsfm_fliq_col_1d            ((endc-begc+1)*nlevgrnd))
      allocate(vsfm_smpl_col_1d            ((endc-begc+1)*nlevgrnd))
      allocate(vsfm_soilp_col_1d           ((endc-begc+1)*nlevgrnd))
      allocate(vsfm_sat_col_1d             ((endc-begc+1)*nlevgrnd))

      ! initialize
      mflx_et_col_1d(:)                = 0.d0
      mflx_infl_col_1d(:)              = 0.d0
      mflx_dew_col_1d(:)               = 0.d0
      mflx_drain_col_1d(:)             = 0.d0
      mflx_sub_snow_col_1d(:)          = 0.d0
      mflx_snowlyr_col_1d(:)           = 0.d0
      t_soil_col_1d(:)                 = 298.15d0

      mass_beg                         = 0.d0
      mass_end                         = 0.d0
      total_mass_flux                  = 0.d0
      total_mass_flux_et               = 0.d0
      total_mass_flux_infl             = 0.d0
      total_mass_flux_dew              = 0.d0
      total_mass_flux_drain            = 0.d0
      total_mass_flux_snowlyr          = 0.d0
      total_mass_flux_sub              = 0.d0
      total_mass_flux_lateral          = 0.d0

      mass_beg_col(:)                  = 0.d0
      mass_end_col(:)                  = 0.d0
      total_mass_flux_col(:)           = 0.d0
      total_mass_flux_et_col(:)        = 0.d0
      total_mass_flux_infl_col(:)      = 0.d0
      total_mass_flux_dew_col(:)       = 0.d0
      total_mass_flux_drain_col(:)     = 0.d0
      total_mass_flux_snowlyr_col(:)   = 0.d0
      total_mass_flux_sub_col(:)       = 0.d0
      total_mass_flux_lateral_col(:)   = 0.d0

      vsfm_mass_prev_col(:,:)          = 0.d0
      vsfm_dmass_col(:)                = 0.d0

      ! Get total mass
      soe_auxvar_id = 1;
      call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL ,       &
                                            VAR_MASS        ,       &
                                            soe_auxvar_id   ,       &
                                            vsfm_mass_col_1d        &
                                           )

      do fc = 1, l2e_num_hydrologyc
         c = l2e_filter_hydrologyc(fc)

         do j = 1, nlevgrnd

            idx = (c - begc)*nlevgrnd + j
            mflx_et_col_1d(idx)          = l2e_mflux_et(c,j)
            mflx_drain_col_1d(idx)       = l2e_mflux_drain(c,j)

            total_mass_flux_et           = total_mass_flux_et           + mflx_et_col_1d(idx)
            total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + mflx_et_col_1d(idx)

            total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(idx)
            total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(idx)

            mass_beg                     = mass_beg                     + vsfm_mass_col_1d(idx)
            mass_beg_col(c)              = mass_beg_col(c)              + vsfm_mass_col_1d(idx)
            vsfm_mass_prev_col(c,j)      = vsfm_mass_col_1d(idx)
         end do

         idx = c - begc+1

         mflx_dew_col_1d(idx)           = l2e_mflux_dew(c)
         mflx_infl_col_1d(idx)          = l2e_mflux_infil(c)
         mflx_snowlyr_col_1d(idx)       = l2e_mflux_snowlyr(c)
         mflx_sub_snow_col_1d(idx)      = l2e_mflux_sub_snow(c)

         total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(idx)
         total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(idx)

         total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(idx)
         total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(idx)

         total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(idx)
         total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(idx)

         total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(idx)
         total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(idx)

         total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                  total_mass_flux_infl_col(c)    + &
                                  total_mass_flux_dew_col(c)     + &
                                  total_mass_flux_drain_col(c)   + &
                                  total_mass_flux_snowlyr_col(c) + &
                                  total_mass_flux_sub_col(c)     + &
                                  total_mass_flux_lateral_col(c)
      end do
      total_mass_flux        = total_mass_flux_et        + &
                               total_mass_flux_infl      + &
                               total_mass_flux_dew       + &
                               total_mass_flux_drain     + &
                               total_mass_flux_snowlyr   + &
                               total_mass_flux_sub       + &
                               total_mass_flux_lateral

      ! Set temperature
      soe_auxvar_id = 1;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL ,      &
                                             VAR_TEMPERATURE ,      &
                                             soe_auxvar_id   ,      &
                                             t_soil_col_1d          &
                                            )
      ! Set Infiltration
      soe_auxvar_id = this%vsfm_cond_id_for_infil;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_infl_col_1d       &
                                            )
      ! Set ET
      soe_auxvar_id = this%vsfm_cond_id_for_et;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_et_col_1d         &
                                            )
      ! Set Dew
      soe_auxvar_id = this%vsfm_cond_id_for_dew;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_dew_col_1d        &
                                            )
      ! Set Drainage sink
      soe_auxvar_id = this%vsfm_cond_id_for_drainage;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_drain_col_1d      &
                                            )
      ! Set mass flux associated with disappearance of snow layer
      ! from last time step
      soe_auxvar_id = this%vsfm_cond_id_for_snow;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_snowlyr_col_1d    &
                                            )
      ! Set mass flux associated with sublimation of snow
      soe_auxvar_id = this%vsfm_cond_id_for_sublimation;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS            , &
                                             VAR_BC_SS_CONDITION  , &
                                             soe_auxvar_id        , &
                                             mflx_sub_snow_col_1d   &
                                            )

      frac_ice(:,:)       = 0.d0
      vsfm_fliq_col_1d(:) = 1.d0
      do fc = 1, l2e_num_hydrologyc
         c = l2e_filter_hydrologyc(fc)
         do j = 1, nlevgrnd

            frac_ice(c,j) = l2e_h2osoi_ice(c,j)/(l2e_h2osoi_liq(c,j) + l2e_h2osoi_ice(c,j))

            idx = (c - begc)*nlevgrnd + j
            vsfm_fliq_col_1d(idx) = 1._r8 - frac_ice(c,j)
         end do
      end do

      ! Set frac_liq
      soe_auxvar_id = 1;
      call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL  , &
                                             VAR_FRAC_LIQ_SAT , &
                                             soe_auxvar_id    , &
                                             vsfm_fliq_col_1d   &
                                            )


#if 0
      if (vsfm_lateral_model_type == 'source_sink') then

         call get_proc_bounds(bounds_proc)

         allocate(vsfm_soilp_col_ghosted_1d((bounds_proc%endc_all - bounds_proc%begc_all+1)*nlevgrnd))
         allocate(vsfm_fliq_col_ghosted_1d( (bounds_proc%endc_all - bounds_proc%begc_all+1)*nlevgrnd))
         allocate(mflx_lateral_col_1d( (bounds_proc%endc - bounds_proc%begc+1)*nlevgrnd))

         soe_auxvar_id = 1;
         call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL   , &
                                               VAR_PRESSURE      , &
                                               soe_auxvar_id     , &
                                               vsfm_soilp_col_1d   &
                                               )

         call ExchangeColumnLevelGhostData(bounds, nlevgrnd, vsfm_soilp_col_1d, vsfm_soilp_col_ghosted_1d)
         call ExchangeColumnLevelGhostData(bounds, nlevgrnd, vsfm_fliq_col_1d,  vsfm_fliq_col_ghosted_1d )

         soe_auxvar_id = 1;
         call this%vsfm_mpp%soe%SetDataFromCLMForGhost(AUXVAR_INTERNAL           , &
                                                        VAR_PRESSURE              , &
                                                        soe_auxvar_id             , &
                                                        vsfm_soilp_col_ghosted_1d   &
                                                       )

         soe_auxvar_id = 1;
         call this%vsfm_mpp%soe%SetDataFromCLMForGhost(AUXVAR_INTERNAL          , &
                                                        VAR_FRAC_LIQ_SAT         , &
                                                        soe_auxvar_id            , &
                                                        vsfm_fliq_col_ghosted_1d   &
                                                       )

         call this%vsfm_mpp%soe%ComputeLateralFlux(dtime)

         soe_auxvar_id = this%vsfm_cond_id_for_lateral_flux;
         call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS   , &
                                               VAR_BC_SS_CONDITION      , &
                                               soe_auxvar_id     , &
                                               mflx_lateral_col_1d   &
                                               )

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)

            g    = col%gridCell(c)
            area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

            ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
            flux_unit_conversion     = area * denh2o * 1.0d-3

            qflx_lateral(c) = 0._r8
            do j = 1, nlevgrnd
               idx = (c-bounds%begc)*nlevgrnd + j

               total_mass_flux_lateral_col(c) =      &
                    total_mass_flux_lateral_col(c) + &
                    mflx_lateral_col_1d(idx)

               qflx_lateral(c) = qflx_lateral(c) - &
                    mflx_lateral_col_1d(idx)/flux_unit_conversion
            enddo

            total_mass_flux_lateral = total_mass_flux_lateral + &
                 total_mass_flux_lateral_col(c)
         enddo

         deallocate(vsfm_soilp_col_ghosted_1d )
         deallocate(vsfm_fliq_col_ghosted_1d  )
         deallocate(mflx_lateral_col_1d       )

         if (vsfm_include_seepage_bc) then
            allocate(seepage_press_1d( (bounds_proc%endc - bounds_proc%begc+1)))
            seepage_press_1d(:) = 101325.d0
            soe_auxvar_id = 1
            call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
                 VAR_BC_SS_CONDITION, soe_auxvar_id, seepage_press_1d)
            deallocate(seepage_press_1d)
         endif

      else if (vsfm_lateral_model_type == 'three_dimensional') then

         call get_proc_bounds(bounds_proc)

         if (vsfm_include_seepage_bc) then
            allocate(seepage_press_1d( (bounds_proc%endc - bounds_proc%begc+1)))
            seepage_press_1d(:) = 101325.d0
            soe_auxvar_id = 1
            call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
                 VAR_BC_SS_CONDITION, soe_auxvar_id, seepage_press_1d)
            deallocate(seepage_press_1d)
         endif

      endif
#endif

      ! Preform Pre-StepDT operations
      call this%vsfm_mpp%soe%PreStepDT()

      ! Get default SNES settings
      call SNESGetTolerances(this%vsfm_mpp%soe%solver%snes , &
                             atol_default            , &
                             rtol_default            , &
                             stol_default            , &
                             max_it_default          , &
                             max_f_default           , &
                             ierr                      &
                            )
      CHKERRQ(ierr)

      stol = stol_default
      rtol = rtol_default

      !
      ! Solve the VSFM.
      !
      iter_count           = 0
      diverged_count       = 0
      mass_bal_err_count   = 0
      abs_mass_error_col   = 0.d0
      successful_step      = PETSC_FALSE

      do

         iter_count = iter_count + 1

         call SNESSetTolerances(this%vsfm_mpp%soe%solver%snes , &
                                atol_default            , &
                                rtol                    , &
                                stol                    , &
                                max_it_default          , &
                                max_f_default           , &
                                ierr                      &
                               );
         CHKERRQ(ierr)

         call this%vsfm_mpp%soe%StepDT(dtime, nstep, &
              converged, converged_reason, ierr); CHKERRQ(ierr)

         if (.not. converged) then

            ! VSFM solver did not converge, so let's try again with different
            ! solver settings.

            stol             = stol_alternate
            diverged_count   = diverged_count + 1
            successful_step  = PETSC_FALSE

            ! Reduce total run length time by the amount VSFM ran successfully
            ! with previous solver settings
            dtime = dtime - this%vsfm_mpp%soe%time

            if (diverged_count > 1) then
               ! Set frac_liq
               vsfm_fliq_col_1d(:) = 1.d0
               soe_auxvar_id = 1;

               call this%vsfm_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL  , &
                                                      VAR_FRAC_LIQ_SAT , &
                                                      soe_auxvar_id    , &
                                                      vsfm_fliq_col_1d   &
                                                     )
            end if
         else

            ! Solver converged, so let's copy data from VSFM model to
            ! CLM's data structure.

            ! Get Liquid saturation
            soe_auxvar_id = 1;
            call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL , &
                                                  VAR_LIQ_SAT     , &
                                                  soe_auxvar_id   , &
                                                  vsfm_sat_col_1d   &
                                                 )

            ! Get total mass
            soe_auxvar_id = 1;
            call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL  , &
                                                  VAR_MASS         , &
                                                  soe_auxvar_id    , &
                                                  vsfm_mass_col_1d   &
                                                 )

            ! Get liquid soil matrix potential
            soe_auxvar_id = 1;
            call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL     , &
                                                  VAR_SOIL_MATRIX_POT , &
                                                  soe_auxvar_id       , &
                                                  vsfm_smpl_col_1d      &
                                                 )

            ! Get soil liquid pressure. This is the prognostic state of VSFM
            ! and needs to be saved in the restart file.
            soe_auxvar_id = 1;
            call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL   , &
                                                  VAR_PRESSURE      , &
                                                  soe_auxvar_id     , &
                                                  vsfm_soilp_col_1d   &
                                                 )

            qflx_seepage(:) = 0._r8

#if 0
            if (vsfm_lateral_model_type == 'source_sink') then

               ! Get following fluxes from VSFM:
               ! (i) seepage mass exchanged.

               call get_proc_bounds(bounds_proc)

               allocate(seepage_mass_exc_col_1d( (bounds_proc%endc - bounds_proc%begc+1)         ))
               seepage_mass_exc_col_1d = 0.d0

               if (vsfm_include_seepage_bc) then
                  soe_auxvar_id = 1;
                  call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_BC              ,  &
                                                        VAR_BC_MASS_EXCHANGED  ,  &
                                                        soe_auxvar_id          ,  &
                                                        seepage_mass_exc_col_1d   &
                                                        )
               endif

               do fc = 1, num_hydrologyc
                  c = filter_hydrologyc(fc)

                  g    = col%gridCell(c)
                  area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

                  ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
                  flux_unit_conversion     = area * denh2o * 1.0d-3

                  idx = (c-bounds%begc) + 1

                  qflx_seepage(c)                = seepage_mass_exc_col_1d(idx)/flux_unit_conversion/dtime
                  total_mass_flux_seepage_col(c) = -seepage_mass_exc_col_1d(idx)/dtime

                  total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                           total_mass_flux_infl_col(c)    + &
                                           total_mass_flux_dew_col(c)     + &
                                           total_mass_flux_drain_col(c)   + &
                                           total_mass_flux_snowlyr_col(c) + &
                                           total_mass_flux_sub_col(c)     + &
                                           total_mass_flux_lateral_col(c) + &
                                           total_mass_flux_seepage_col(c)
               enddo

            else if (vsfm_lateral_model_type == 'three_dimensional') then

               ! Get following fluxes from VSFM:
               ! (i) lateral mass exchanged, and
               ! (ii) seepage mass exchanged.

               call get_proc_bounds(bounds_proc)

               allocate(lat_mass_exc_col_1d(     (bounds_proc%endc - bounds_proc%begc+1)*nlevgrnd))
               allocate(seepage_mass_exc_col_1d( (bounds_proc%endc - bounds_proc%begc+1)         ))

               lat_mass_exc_col_1d(:) = 0.d0
               seepage_mass_exc_col_1d(:) = 0.d0

               soe_auxvar_id = 1
               call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL            , &
                                                     VAR_LATERAL_MASS_EXCHANGED , &
                                                     soe_auxvar_id              , &
                                                     lat_mass_exc_col_1d          &
                                                     )

               if (vsfm_include_seepage_bc) then
                  soe_auxvar_id = 1;
                  call this%vsfm_mpp%soe%GetDataForCLM(AUXVAR_BC              ,  &
                                                        VAR_BC_MASS_EXCHANGED  ,  &
                                                        soe_auxvar_id          ,  &
                                                        seepage_mass_exc_col_1d   &
                                                        )
               endif

               total_mass_flux_lateral_col(:)   = 0.d0
               total_mass_flux_lateral          = 0.d0

               do fc = 1, num_hydrologyc
                  c = filter_hydrologyc(fc)

                  g    = col%gridCell(c)
                  area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

                  ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
                  flux_unit_conversion     = area * denh2o * 1.0d-3

                  qflx_lateral(c) = 0._r8
                  do j = 1, nlevgrnd
                     idx = (c-bounds%begc)*nlevgrnd + j

                     total_mass_flux_lateral_col(c) = &
                          total_mass_flux_lateral_col(c) + &
                          lat_mass_exc_col_1d(idx)/dtime

                     qflx_lateral(c)  = qflx_lateral(c) - &
                          lat_mass_exc_col_1d(idx)/flux_unit_conversion/dtime

                  enddo

                  idx = (c-bounds%begc) + 1
                  qflx_seepage(c)                = seepage_mass_exc_col_1d(idx)/flux_unit_conversion/dtime
                  total_mass_flux_seepage_col(c) = -seepage_mass_exc_col_1d(idx)/dtime

                  total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                           total_mass_flux_infl_col(c)    + &
                                           total_mass_flux_dew_col(c)     + &
                                           total_mass_flux_drain_col(c)   + &
                                           total_mass_flux_snowlyr_col(c) + &
                                           total_mass_flux_sub_col(c)     + &
                                           total_mass_flux_lateral_col(c) + &
                                           total_mass_flux_seepage_col(c)

                  total_mass_flux_lateral = total_mass_flux_lateral + &
                       total_mass_flux_lateral_col(c)
               enddo

               deallocate(lat_mass_exc_col_1d)
            endif
#endif

            ! Put the data in CLM's data structure
            mass_end        = 0.d0
            area            = 1.d0 ! [m^2]

            do fc = 1, l2e_num_hydrologyc
               c = l2e_filter_hydrologyc(fc)

               !if (lateral_connectivity) then
               !  g    = col%gridCell(c)
               !   area = ldomain_lateral%ugrid%areaGrid_ghosted(g)
               !endif

               ! initialization
               jwt = -1

               ! Loops in decreasing j so WTD can be computed in the same loop
               do j = nlevgrnd, 1, -1
                  idx = (c-begc)*nlevgrnd + j

                  e2l_h2osoi_liq(c,j) = (1.d0 - frac_ice(c,j))*vsfm_mass_col_1d(idx)/area
                  e2l_h2osoi_ice(c,j) = frac_ice(c,j)         *vsfm_mass_col_1d(idx)/area

                  mass_end        = mass_end        + vsfm_mass_col_1d(idx)
                  mass_end_col(c) = mass_end_col(c) + vsfm_mass_col_1d(idx)

                  vsfm_dmass_col(c) = vsfm_dmass_col(c) + &
                                      (vsfm_mass_col_1d(idx)-vsfm_mass_prev_col(c,j))

                  e2l_smp(c,j)    = vsfm_smpl_col_1d(idx)*1000.0_r8      ! [m] --> [mm]

                  if (jwt == -1) then
                     ! Find the first soil that is unsaturated
                     if (e2l_smp(c,j) < 0._r8) jwt = j
                  end if

               end do

               ! Find maximum water balance error over the column
               abs_mass_error_col = max(abs_mass_error_col,                     &
                                        abs(mass_beg_col(c) - mass_end_col(c) + &
                                            total_mass_flux_col(c)*dt))
               e2l_qrecharge(c) = 0._r8

               if (jwt == -1 .or. jwt == nlevgrnd) then
                  ! Water table below or in the last layer
                  e2l_wtd(c) = l2e_zi(c,nlevgrnd)
               else
                  z_dn = (l2e_zi(c,jwt-1) + l2e_zi(c,jwt  ))/2._r8
                  z_up = (l2e_zi(c,jwt ) + l2e_zi(c,jwt+1))/2._r8
                  e2l_wtd(c) = (0._r8 - e2l_smp(c,jwt))/(e2l_smp(c,jwt) - e2l_smp(c,jwt+1))*(z_dn - z_up) + z_dn
               endif
            end do

            ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
            ! soilp_col is used for restarting VSFM.
            do c = begc, endc
               do j = 1, nlevgrnd
                  idx = (c - begc)*nlevgrnd + j
                  e2l_soilp(c,j) = vsfm_soilp_col_1d(idx)
               end do
            end do

            ! For the solution that did converge, is the mass error acceptable?
            if (abs_mass_error_col >= max_abs_mass_error_col) then

               ! For the solution that converged, the mass error
               ! is unacceptable. So let's try again with tighter
               ! solution tolerance (stol) for SNES.

               mass_bal_err_count  = mass_bal_err_count + 1

               if (converged_reason == SNES_CONVERGED_FNORM_RELATIVE) then
                  rtol = rtol/10._r8
               else if (converged_reason == SNES_CONVERGED_SNORM_RELATIVE) then
                  stol = stol/10._r8
               endif

               dtime               = dt
               successful_step     = PETSC_FALSE
               abs_mass_error_col  = 0._r8
               mass_end_col(:)     = 0._r8

               ! Perform Pre-StepDT operations
               call this%vsfm_mpp%soe%PreStepDT()

            else

               successful_step  = PETSC_TRUE

            endif

         endif

         if (successful_step) exit

         if (iter_count >= max_iter_count) then
            write(iulog,*)'In soilwater_vsfm: VSFM failed to converge after multiple attempts.'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

      end do

#if 0
      ! Add seepage flux from VSFM to surface runoff
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         qflx_surf(c) = qflx_surf(c) + qflx_seepage(c)
      enddo
#endif

      call SNESSetTolerances(this%vsfm_mpp%soe%solver%snes, atol_default, rtol_default, stol_default, &
                             max_it_default, max_f_default, ierr); CHKERRQ(ierr)

      call this%vsfm_mpp%soe%PostStepDT()

#if VSFM_DEBUG
      write(iulog,*)'VSFM-DEBUG: nstep                      = ',get_nstep()
      write(iulog,*)'VSFM-DEBUG: dtime                      = ',dt
      write(iulog,*)'VSFM-DEBUG: change in mass between dt  = ',-(mass_beg - mass_end)
      write(iulog,*)'VSFM-DEBUG: change in mass due to flux = ',total_mass_flux*dt
      write(iulog,*)'VSFM-DEBUG: Error in mass conservation = ',mass_beg - mass_end + total_mass_flux*dt
      write(iulog,*)'VSFM-DEBUG: et_flux    * dtime         = ',total_mass_flux_et*dt
      write(iulog,*)'VSFM-DEBUG: infil_flux * dtime         = ',total_mass_flux_infl*dt
      write(iulog,*)'VSFM-DEBUG: dew_flux   * dtime         = ',total_mass_flux_dew*dt
      write(iulog,*)'VSFM-DEBUG: drain_flux * dtime         = ',total_mass_flux_drain*dt
      write(iulog,*)'VSFM-DEBUG: snow_flux  * dtime         = ',total_mass_flux_snowlyr*dt
      write(iulog,*)'VSFM-DEBUG: sub_flux   * dtime         = ',total_mass_flux_sub*dt
      write(iulog,*)'VSFM-DEBUG: lat_flux   * dtime         = ',total_mass_flux_lateral*dt
      write(iulog,*)'VSFM-DEBUG: total_mass_flux            = ',total_mass_flux/flux_unit_conversion
      write(iulog,*)'VSFM-DEBUG: et_flux                    = ',total_mass_flux_et
      write(iulog,*)'VSFM-DEBUG: infil_flux                 = ',total_mass_flux_infl
      write(iulog,*)'VSFM-DEBUG: dew_flux                   = ',total_mass_flux_dew
      write(iulog,*)'VSFM-DEBUG: drain_flux                 = ',total_mass_flux_drain
      write(iulog,*)'VSFM-DEBUG: snow_flux                  = ',total_mass_flux_snowlyr
      write(iulog,*)'VSFM-DEBUG: sub_flux                   = ',total_mass_flux_sub
      write(iulog,*)''
#endif

      deallocate(frac_ice                    )
      deallocate(total_mass_flux_col         )
      deallocate(total_mass_flux_et_col      )
      deallocate(total_mass_flux_infl_col    )
      deallocate(total_mass_flux_dew_col     )
      deallocate(total_mass_flux_drain_col   )
      deallocate(total_mass_flux_snowlyr_col )
      deallocate(total_mass_flux_sub_col     )
      deallocate(total_mass_flux_lateral_col )
      deallocate(total_mass_flux_seepage_col )
      deallocate(qflx_seepage                )
      deallocate(vsfm_mass_prev_col          )
      deallocate(vsfm_dmass_col              )
      deallocate(mass_beg_col                )
      deallocate(mass_end_col                )

      deallocate(mflx_et_col_1d              )
      deallocate(mflx_drain_col_1d           )
      deallocate(mflx_infl_col_1d            )
      deallocate(mflx_dew_col_1d             )
      deallocate(mflx_sub_snow_col_1d        )
      deallocate(mflx_snowlyr_col_1d         )
      deallocate(t_soil_col_1d               )

      deallocate(vsfm_mass_col_1d            )
      deallocate(vsfm_fliq_col_1d            )
      deallocate(vsfm_smpl_col_1d            )
      deallocate(vsfm_soilp_col_1d           )
      deallocate(vsfm_sat_col_1d             )

  end subroutine EM_VSFM_Solve_Soil_Hydro

#endif

end module ExternalModelVSFMMod
