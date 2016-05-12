
module SystemOfEquationsThermalType

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for Thermal system-of-equations
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl      , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType
  use GoveqnThermalKSPTemperatureSnowType
  use GoveqnThermalKSPTemperatureSSWType
  use GoveqnThermalKSPTemperatureSoilType
  use SystemOfEquationsBaseType
  use SystemOfEquationsThermalAuxType
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"

  type, public, extends(sysofeqns_base_type) :: sysofeqns_thermal_type

     type (sysofeqns_thermal_auxvar_type), pointer :: aux_vars_in(:)            !!< Internal state.
     type (sysofeqns_thermal_auxvar_type), pointer :: aux_vars_bc(:)            !!< Boundary conditions.
     type (sysofeqns_thermal_auxvar_type), pointer :: aux_vars_ss(:)            !!< Source-sink.

     PetscInt                                      :: num_calls_to_ifunction
     PetscInt                                      :: num_calls_to_ijacobian

     PetscInt, pointer                             :: soe_auxvars_bc_offset (:) ! Cummulative sum of number of control volumes associated with each boundary condition.
     PetscInt, pointer                             :: soe_auxvars_ss_offset (:) ! Cummulative sum of number of control volumes associated with each source-sink condition.
     PetscInt, pointer                             :: soe_auxvars_bc_ncells (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer                             :: soe_auxvars_ss_ncells (:) ! Number of control volumes associated with each source-sink condition.
     PetscInt                                      :: num_auxvars_in            ! Number of auxvars associated with internal state.
     PetscInt                                      :: num_auxvars_in_local      ! Number of auxvars associated with internal state.
     PetscInt                                      :: num_auxvars_bc            ! Number of auxvars associated with boundary condition.
     PetscInt                                      :: num_auxvars_ss            ! Number of auxvars associated with source-sink condition.

   contains
     procedure, public :: Init            => ThermalSOEInit
     procedure, public :: Setup           => ThermalSOESetup
     procedure, public :: SetSolnPrevCLM  => ThermalSOESetSolnPrevCLM
     procedure, public :: GetSoln         => ThermalSOEGetSoln
     procedure, public :: SetRDataFromCLM => ThermalSOESetRDataFromCLM
     procedure, public :: SetIDataFromCLM => ThermalSOESetIDataFromCLM
     procedure, public :: SetBDataFromCLM => ThermalSOESetBDataFromCLM
     procedure, public :: PreStepDT       => ThermalSOEPreStepDT
     procedure, public :: PreSolve        => ThermalSOEPreSolve
     procedure, public :: ComputeRHS      => ThermalSOEComputeRHS
     procedure, public :: ComputeOperators=> ThermalSOEComputeOperators
   end type sysofeqns_thermal_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalSOEInit(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this

    call SOEBaseInit(this)

    this%num_auxvars_in         = 0
    this%num_auxvars_in_local   = 0
    this%num_auxvars_bc         = 0
    this%num_auxvars_ss         = 0

    nullify(this%aux_vars_in           )
    nullify(this%aux_vars_bc           )
    nullify(this%aux_vars_ss           )

    nullify(this%soe_auxvars_bc_offset )
    nullify(this%soe_auxvars_ss_offset )
    nullify(this%soe_auxvars_bc_ncells )
    nullify(this%soe_auxvars_ss_ncells )

  end subroutine ThermalSOEInit

  !------------------------------------------------------------------------
  subroutine ThermalSOESetup(this, mpi_rank, mpp_type, soe_type, meshes, nmesh, begc, endc, z, zi)
    !
    ! !DESCRIPTION:
    ! Sets up SoE for the Thermal
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_THERMAL_TBASED
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    use MeshType, only                  : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    PetscInt                      :: mpi_rank
    PetscInt                      :: mpp_type
    PetscInt                      :: soe_type
    class(mesh_type) , pointer    :: meshes(:)
    PetscInt                      :: nmesh
    integer          , intent(in) :: begc,endc
    PetscReal        , pointer    :: z(:,:)
    PetscReal        , pointer    :: zi(:,:)

    call this%Init()

    this%mpi_rank = mpi_rank

    select case (soe_type)

    case (SOE_THERMAL_TBASED)

       select case(mpp_type)
       case (MPP_THERMAL_TBASED_KSP_CLM)
          call ThermalSOEThermalKSPSetup(this, meshes, nmesh, begc, endc, z, zi)

       case default
          write(iulog,*) 'ThermalSOESetup: Unknown mpp_type'
          call endrun(msg=errMsg(__FILE__, __LINE__))

       end select

    case default
       write(iulog,*) 'ThermalSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalSOESetup


  !------------------------------------------------------------------------
  subroutine ThermalSOEThermalKSPSetup(this, meshes, nmesh, begc, endc, z, zi)
    !
    ! !DESCRIPTION:
    ! Sets up SoE for the Thermal
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_THERMAL_TBASED
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    use MultiPhysicsProbConstants , only : PETSC_KSP
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MeshType                  , only : mesh_type
    use MeshType                  , only : MeshCreateConnectionSet
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionNew
    use ConditionType             , only : ConditionListAddCondition
    use ConnectionSetType         , only : connection_set_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use GoveqnThermalKSPTemperatureSnowType
    use GoveqnThermalKSPTemperatureSSWType
    use GoveqnThermalKSPTemperatureSoilType
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type)                       :: this
    class(mesh_type)                          , pointer :: meshes(:)
    PetscInt                                            :: nmesh
    integer, intent(in)                                 :: begc,endc
    PetscReal, pointer                                  :: z(:,:)
    PetscReal, pointer                                  :: zi(:,:)
    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_ksp_temp_snow_type) , pointer :: goveq_snow
    class (goveqn_thermal_ksp_temp_ssw_type)  , pointer :: goveq_sh2o
    class (goveqn_thermal_ksp_temp_soil_type) , pointer :: goveq_soil
    type(condition_type)                      , pointer :: bc
    type(condition_type)                      , pointer :: ss
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type(connection_set_type)                 , pointer :: cur_conn_set
    type(condition_type)                      , pointer :: cur_cond
    type(condition_type)                      , pointer :: soil_bc_for_snow
    type(condition_type)                      , pointer :: soil_bc_for_sh2o
    type(condition_type)                      , pointer :: snow_bc_for_soil
    type(condition_type)                      , pointer :: sh2o_bc_for_soil
    
    PetscInt                                            :: snow_bc_auxvar_offset_1
    PetscInt                                            :: sh2o_bc_auxvar_offset_1
    PetscInt                                            :: soil_bc_auxvar_offset_1
    PetscInt                                            :: soil_bc_auxvar_offset_2
    PetscInt                                            :: snow_bc_auxvar_ncells_1
    PetscInt                                            :: sh2o_bc_auxvar_ncells_1
    PetscInt                                            :: soil_bc_auxvar_ncells_1
    PetscInt                                            :: soil_bc_auxvar_ncells_2
    PetscInt                                            :: snow_bc_auxvar_idx_1
    PetscInt                                            :: sh2o_bc_auxvar_idx_1
    PetscInt                                            :: soil_bc_auxvar_idx_1
    PetscInt                                            :: soil_bc_auxvar_idx_2
    PetscInt                                            :: snow_bc_auxvar_idx_of_other_goveqn_1
    PetscInt                                            :: sh2o_bc_auxvar_idx_of_other_goveqn_1
    PetscInt                                            :: soil_bc_auxvar_idx_of_other_goveqn_1
    PetscInt                                            :: soil_bc_auxvar_idx_of_other_goveqn_2
    PetscInt                                            :: icond
    PetscInt                                            :: snow_num_bc
    PetscInt                                            :: sh2o_num_bc
    PetscInt                                            :: soil_num_bc
    PetscInt                                            :: snow_num_ss
    PetscInt                                            :: sh2o_num_ss
    PetscInt                                            :: soil_num_ss
    PetscInt                                            :: iauxvar
    PetscInt                                            :: nvars
    PetscInt                                            :: count
    PetscInt                                            :: total_ncells_for_bc
    PetscInt                                            :: total_ncells_for_ss
    PetscInt                                  , pointer :: snow_ncells_for_bc(:)
    PetscInt                                  , pointer :: sh2o_ncells_for_bc(:)
    PetscInt                                  , pointer :: soil_ncells_for_bc(:)
    PetscInt                                  , pointer :: snow_ncells_for_ss(:)
    PetscInt                                  , pointer :: sh2o_ncells_for_ss(:)
    PetscInt                                  , pointer :: soil_ncells_for_ss(:)
    PetscInt                                  , pointer :: num_cells(:)
    PetscInt                                            :: igoveqn
    PetscInt                                            :: iauxvar_beg
    PetscInt                                            :: iauxvar_end
    PetscInt                                            :: iconn
    PetscBool                                           :: snow_bc_for_soil_found
    PetscBool                                           :: sh2o_bc_for_soil_found
    PetscBool                                           :: soil_bc_for_snow_found
    PetscBool                                           :: soil_bc_for_sh2o_found
    PetscInt                                  , pointer :: snow_ids_for_soil_bc(:)
    PetscInt                                  , pointer :: sh2o_ids_for_soil_bc(:)
    PetscInt                                  , pointer :: soil_ids_for_snow_bc(:)
    PetscInt                                  , pointer :: soil_ids_for_sh2o_bc(:)
    PetscReal                                 , pointer :: soil_dist_for_snow_bc(:)
    PetscReal                                 , pointer :: soil_dist_for_sh2o_bc(:)
    PetscInt                                  , pointer :: snow_soe_auxvars_bc_offset(:)
    PetscInt                                  , pointer :: sh2o_soe_auxvars_bc_offset(:)
    PetscInt                                  , pointer :: soil_soe_auxvars_bc_offset(:)
    PetscInt                                  , pointer :: snow_soe_auxvars_ss_offset(:)
    PetscInt                                  , pointer :: sh2o_soe_auxvars_ss_offset(:)
    PetscInt                                  , pointer :: soil_soe_auxvars_ss_offset(:)

    this%name         = "SOE for Thermal equation based on temperature using KSP"
    this%ngoveqns     = 3
    this%itype        = SOE_THERMAL_TBASED
    this%solver_type  = PETSC_KSP
    
    ! Add governing-equations to the system-of-equations
    allocate(goveq_snow)
    allocate(goveq_sh2o)
    allocate(goveq_soil)

    call goveq_snow%Setup()
    call goveq_sh2o%Setup()
    call goveq_soil%Setup()

    goveq_snow%name = "Thermal equation using temprature formulation in snow (KSP formulation)"
    goveq_sh2o%name = "Thermal equation using temprature formulation in standing surface water (KSP formulation)"
    goveq_soil%name = "Thermal equation using temprature formulation in soil (KSP formulation)"

    goveq_snow%id_in_list  = 1
    goveq_snow%mesh_itype  = MESH_CLM_SNOW_COL

    goveq_sh2o%id_in_list  = 2
    goveq_sh2o%mesh_itype  = MESH_CLM_SSW_COL

    goveq_soil%id_in_list  = 3
    goveq_soil%mesh_itype  = MESH_CLM_THERMAL_SOIL_COL

    this%goveqns      => goveq_snow
    goveq_snow%next   => goveq_sh2o
    goveq_sh2o%next   => goveq_soil

    ! Assign mesh to each governing equation
    call SOESetMeshesOfGoveqns(this, meshes, nmesh)

    ! Set BCs for each governing equation
    !                                GOVERNING EQUATION
    !---------------------------------------------------------------|
    !        |      |     snow      |     sh2o      |     soil      |
    !---------------------------------------------------------------|
    !   B    |      |               |               |               |
    !   O    |  1   | Flux at top   | Flux at top   | Flux at top   |
    !   U    |      | snow layer    | of sh2o       | soil layer    |
    !   N    |      |               |               |               |
    !   D    |------------------------------------------------------|
    !   A    |      |               |               |               |
    !   R    |      | For coupling  | For coupling  | For coupling  |
    !   Y    |  2   | to soil GE.   | to soil GE.   | to snow GE.   |
    !        |      |               |               |               |
    !   C    |------------------------------------------------------|
    !   O    |      |               |               |               |
    !   N    |      |               |               | For coupling  |
    !   D    |   3  |               |               | to sh2o GE.   |
    !   S.   |      |               |               |               |
    !        |      |               |               |               |
    !---------------------------------------------------------------|


    ! Set BCs for each governing equation

    ! Snow: BC-1
    bc               => ConditionNew()
    bc%name          = 'Heat_flux_BC_at_top_of_snow'
    bc%units         = 'W/m^2'
    bc%itype         = COND_HEAT_FLUX
    bc%region_itype  = SNOW_TOP_CELLS

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_snow%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells)

    snow_bc_auxvar_offset_1 = bc%ncells

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_snow%boundary_conditions, bc)
    nullify(bc)

    ! Snow: BC-2
    bc                        => ConditionNew()
    bc%name                   = 'BC_from_soil_governing_equation'
    bc%units                  = '[K]'
    bc%itype                  = COND_DIRICHLET_FRM_OTR_GOVEQ
    bc%region_itype           = SNOW_BOTTOM_CELLS
    bc%list_id_of_other_goveq = goveq_soil%id_in_list

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_snow%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells)

    snow_bc_auxvar_ncells_1               = bc%ncells
    snow_bc_auxvar_idx_1                  = 2
    snow_bc_auxvar_idx_of_other_goveqn_1  = 2

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_snow%boundary_conditions, bc)
    nullify(bc)

    ! Standing surface water: BC-1
    bc               => ConditionNew()
    bc%name          = 'Heat_flux_BC_at_top_of_standing_surface_water'
    bc%units         = 'W/m^2'
    bc%itype         = COND_HEAT_FLUX
    bc%region_itype  = SSW_TOP_CELLS

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_sh2o%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells)

    sh2o_bc_auxvar_offset_1 = bc%ncells

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_sh2o%boundary_conditions, bc)
    nullify(bc)

    ! Standing surface water: BC-2
    bc               => ConditionNew()
    bc%name          = 'BC_from_soil_governing_equation'
    bc%units         = '[K]'
    bc%itype         = COND_DIRICHLET_FRM_OTR_GOVEQ
    bc%region_itype  = SSW_TOP_CELLS
    bc%list_id_of_other_goveq  = goveq_soil%id_in_list

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_sh2o%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells)
    sh2o_bc_auxvar_ncells_1               = bc%ncells
    sh2o_bc_auxvar_idx_1                  = 2
    sh2o_bc_auxvar_idx_of_other_goveqn_1  = 3

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_sh2o%boundary_conditions, bc)
    nullify(bc)

    ! Soil: BC-1
    bc               => ConditionNew()
    bc%name          = 'Heat flux BC at top of soil'
    bc%units         = 'W/m^2'
    bc%itype         = COND_HEAT_FLUX
    bc%region_itype  = SOIL_TOP_CELLS

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_soil%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells)

    soil_bc_auxvar_offset_1 = bc%ncells

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_soil%boundary_conditions, bc)
    nullify(bc)

    ! Soil: BC-2
    bc               => ConditionNew()
    bc%name          = 'BC_from_snow_governing_equation'
    bc%units         = '[K]'
    bc%itype         = COND_DIRICHLET_FRM_OTR_GOVEQ
    bc%region_itype  = SOIL_TOP_CELLS
    bc%list_id_of_other_goveq  = goveq_snow%id_in_list

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_soil%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells,        &
                                 use_clm_dist_to_interface = PETSC_TRUE, &
                                 begc=begc, endc=endc, z=z, zi=zi)

    soil_bc_auxvar_ncells_1               = bc%ncells
    soil_bc_auxvar_idx_1                  = 2
    soil_bc_auxvar_idx_of_other_goveqn_1  = 2
    soil_bc_auxvar_offset_2               = soil_bc_auxvar_offset_1 + bc%ncells

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_soil%boundary_conditions, bc)
    nullify(bc)

    ! Soil: BC-3
    bc               => ConditionNew()
    bc%name          = 'BC_from_standing_surface_water_governing_equation'
    bc%units         = '[K]'
    bc%itype         = COND_DIRICHLET_FRM_OTR_GOVEQ
    bc%region_itype  = SOIL_TOP_CELLS
    bc%list_id_of_other_goveq  = goveq_sh2o%id_in_list

    allocate(bc%conn_set)
    call MeshCreateConnectionSet(goveq_soil%mesh,  &
                                 bc%region_itype,  &
                                 bc%conn_set,      &
                                 bc%ncells,        &
                                 use_clm_dist_to_interface = PETSC_TRUE, &
                                 begc=begc, endc=endc, z=z, zi=zi)

    soil_bc_auxvar_ncells_2               = bc%ncells
    soil_bc_auxvar_idx_2                  = 3
    soil_bc_auxvar_idx_of_other_goveqn_2  = 2

    allocate(bc%value(bc%ncells))
    bc%value(:) = 0.d0
    call ConditionListAddCondition(goveq_soil%boundary_conditions, bc)
    nullify(bc)

    ! Set SSs for each governing equation

    ss               => ConditionNew()
    ss%name          = 'Absorbed_solar_radiation'
    ss%units         = 'W/m^2'
    ss%itype         = COND_HEAT_RATE
    ss%region_itype  = ALL_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(goveq_snow%mesh, &
                                 ss%region_itype,    &
                                 ss%conn_set,        &
                                 ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_snow%source_sinks, ss)
    nullify(ss)

    ss               => ConditionNew()
    ss%name          = 'Absorbed_solar_radiation'
    ss%units         = 'W/m^2'
    ss%itype         = COND_HEAT_RATE
    ss%region_itype  = ALL_CELLS
    allocate(ss%conn_set)
    call MeshCreateConnectionSet(goveq_soil%mesh, &
                                 ss%region_itype,    &
                                 ss%conn_set,        &
                                 ss%ncells)
    allocate(ss%value(ss%ncells))
    ss%value(:) = 0.d0
    call ConditionListAddCondition(goveq_soil%source_sinks, ss)
    nullify(ss)

    ! Allocate memory for aux vars
    call goveq_snow%AllocateAuxVars()
    call goveq_sh2o%AllocateAuxVars()
    call goveq_soil%AllocateAuxVars()


    ! Determine number of SoE auxvars for internal connections
    allocate(num_cells(this%ngoveqns))
    count = 0
    this%num_auxvars_in = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       count = count + 1
       num_cells(count) = cur_goveq%mesh%ncells_all
       this%num_auxvars_in = this%num_auxvars_in + &
                              cur_goveq%mesh%ncells_all
       cur_goveq => cur_goveq%next
    enddo
    
    ! Determine number of SoE auxvars for boundary condition
    ! (excluding bc_typ COND_DIRICHLET_FRM_OTR_GOVEQ)
    this%num_auxvars_bc = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)

          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, snow_num_bc, snow_ncells_for_bc)

          do icond = 1, snow_num_bc
             this%num_auxvars_bc = this%num_auxvars_bc + &
                                       snow_ncells_for_bc(icond)
          enddo

       class is (goveqn_thermal_ksp_temp_ssw_type)

          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, sh2o_num_bc, sh2o_ncells_for_bc)

          do icond = 1, sh2o_num_bc
             this%num_auxvars_bc = this%num_auxvars_bc + &
                                       sh2o_ncells_for_bc(icond)
          enddo

       class is (goveqn_thermal_ksp_temp_soil_type)

          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, soil_num_bc, soil_ncells_for_bc)

          do icond = 1, soil_num_bc
             this%num_auxvars_bc = this%num_auxvars_bc + &
                                       soil_ncells_for_bc(icond)
          enddo

       class default
           write(iulog,*)'ThermalSOEThermalKSPSetup: Unknown class'
           call endrun(msg=errMsg(__FILE__, __LINE__))
        end select

       cur_goveq => cur_goveq%next
    enddo


    ! Determine number of SoE auxvars for source-sinks
    this%num_auxvars_ss = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, snow_num_ss, snow_ncells_for_ss)
          do icond = 1, snow_num_ss
             this%num_auxvars_ss = this%num_auxvars_ss + &
                                       snow_ncells_for_ss(icond)
          enddo

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, sh2o_num_ss, sh2o_ncells_for_ss)
          do icond = 1, sh2o_num_ss
             this%num_auxvars_ss = this%num_auxvars_ss + &
                                       sh2o_ncells_for_ss(icond)
          enddo

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, soil_num_ss, soil_ncells_for_ss)

          do icond = 1, soil_num_ss
             this%num_auxvars_ss = this%num_auxvars_ss + &
                                       soil_ncells_for_ss(icond)
          enddo

        class default
           write(iulog,*)'ThermalSOEThermalKSPSetup: Unknown class'
           call endrun(msg=errMsg(__FILE__, __LINE__))
        end select

       cur_goveq => cur_goveq%next
    enddo

    
    ! Allocate memory
    allocate(this%aux_vars_in(this%num_auxvars_in))
    allocate(this%aux_vars_bc(this%num_auxvars_bc))
    allocate(this%aux_vars_ss(this%num_auxvars_ss))

    ! Initialize SoE auxvars
    iauxvar_beg = 0
    iauxvar_end = 0
    do igoveqn = 1, this%ngoveqns
       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = num_cells(igoveqn) + iauxvar_end

       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_in(iauxvar)%is_in     = PETSC_TRUE
          this%aux_vars_in(iauxvar)%goveqn_id = igoveqn
       enddo

    enddo

    allocate(this%soe_auxvars_bc_offset(this%num_auxvars_bc))
    allocate(this%soe_auxvars_ss_offset(this%num_auxvars_bc))
    allocate(this%soe_auxvars_bc_ncells(this%num_auxvars_bc))
    allocate(this%soe_auxvars_ss_ncells(this%num_auxvars_ss))

    allocate(snow_soe_auxvars_bc_offset(snow_num_bc))
    allocate(sh2o_soe_auxvars_bc_offset(sh2o_num_bc))
    allocate(soil_soe_auxvars_bc_offset(soil_num_bc))
    allocate(snow_soe_auxvars_ss_offset(snow_num_ss))
    allocate(sh2o_soe_auxvars_ss_offset(sh2o_num_ss))
    allocate(soil_soe_auxvars_ss_offset(soil_num_ss))

    count               = 0
    total_ncells_for_bc = 0
    iauxvar_beg         = 0
    iauxvar_end         = 0

    igoveqn = 1
    do icond = 1, snow_num_bc
       count = count + 1
       this%soe_auxvars_bc_offset(count) = total_ncells_for_bc
       snow_soe_auxvars_bc_offset(icond) = total_ncells_for_bc
       this%soe_auxvars_bc_ncells(count) = snow_ncells_for_bc(icond)
       total_ncells_for_bc = total_ncells_for_bc + snow_ncells_for_bc(icond)

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + snow_ncells_for_bc(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_bc(iauxvar)%is_bc     = PETSC_TRUE
          this%aux_vars_bc(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_bc(iauxvar)%condition_id = icond
       enddo
    enddo

    igoveqn = 2
    do icond = 1, sh2o_num_bc
       count = count + 1
       this%soe_auxvars_bc_offset(count) = total_ncells_for_bc
       sh2o_soe_auxvars_bc_offset(icond) = total_ncells_for_bc
       this%soe_auxvars_bc_ncells(count) = sh2o_ncells_for_bc(icond)
       total_ncells_for_bc = total_ncells_for_bc + sh2o_ncells_for_bc(icond)

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + sh2o_ncells_for_bc(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_bc(iauxvar)%is_bc     = PETSC_TRUE
          this%aux_vars_bc(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_bc(iauxvar)%condition_id = icond
       enddo
    enddo

    igoveqn = 3
    do icond = 1, soil_num_bc
       count = count + 1
       this%soe_auxvars_bc_offset(count) = total_ncells_for_bc
       soil_soe_auxvars_bc_offset(icond) = total_ncells_for_bc
       this%soe_auxvars_bc_ncells(count) = soil_ncells_for_bc(icond)
       total_ncells_for_bc = total_ncells_for_bc + soil_ncells_for_bc(icond)
       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + soil_ncells_for_bc(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_bc(iauxvar)%is_bc     = PETSC_TRUE
          this%aux_vars_bc(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_bc(iauxvar)%condition_id = icond
       enddo
    enddo


    count               = 0
    total_ncells_for_ss = 0
    iauxvar_beg         = 0
    iauxvar_end         = 0

    igoveqn = 1
    do icond = 1, snow_num_ss
       count = count + 1
       this%soe_auxvars_ss_offset(count) = total_ncells_for_ss
       snow_soe_auxvars_ss_offset(icond) = total_ncells_for_ss
       this%soe_auxvars_ss_ncells(count) = snow_ncells_for_ss(icond)
       total_ncells_for_ss = total_ncells_for_ss + snow_ncells_for_ss(icond)

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + snow_ncells_for_ss(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_ss(iauxvar)%is_ss     = PETSC_TRUE
          this%aux_vars_ss(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_ss(iauxvar)%condition_id = icond
       enddo
    enddo

    igoveqn = 2
    do icond = 1, sh2o_num_ss
       count = count + 1
       this%soe_auxvars_ss_offset(count) = total_ncells_for_ss
       sh2o_soe_auxvars_ss_offset(icond) = total_ncells_for_ss
       this%soe_auxvars_ss_ncells(count) = sh2o_ncells_for_ss(icond)
       total_ncells_for_ss = total_ncells_for_ss + sh2o_ncells_for_ss(icond)

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + sh2o_ncells_for_ss(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_ss(iauxvar)%is_ss     = PETSC_TRUE
          this%aux_vars_ss(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_ss(iauxvar)%condition_id = icond
       enddo
    enddo

    igoveqn = 3
    do icond = 1, soil_num_ss
       count = count + 1
       this%soe_auxvars_ss_offset(count) = total_ncells_for_ss
       soil_soe_auxvars_ss_offset(icond) = total_ncells_for_ss
       this%soe_auxvars_ss_ncells(count) = soil_ncells_for_ss(icond)
       total_ncells_for_ss = total_ncells_for_ss + soil_ncells_for_ss(icond)

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + soil_ncells_for_ss(icond)
       do iauxvar = iauxvar_beg, iauxvar_end          
          call this%aux_vars_in(iauxvar)%Init()
          this%aux_vars_ss(iauxvar)%is_ss     = PETSC_TRUE
          this%aux_vars_ss(iauxvar)%goveqn_id = igoveqn
          this%aux_vars_ss(iauxvar)%condition_id = icond
       enddo
    enddo

    call goveq_snow%SetSOEAuxVarOffsets(snow_num_bc, snow_soe_auxvars_bc_offset, &
         snow_num_ss, snow_soe_auxvars_ss_offset)
    call goveq_sh2o%SetSOEAuxVarOffsets(sh2o_num_bc, sh2o_soe_auxvars_bc_offset, &
         sh2o_num_ss, sh2o_soe_auxvars_ss_offset)
    call goveq_soil%SetSOEAuxVarOffsets(soil_num_bc, soil_soe_auxvars_bc_offset, &
         soil_num_ss, soil_soe_auxvars_ss_offset)

    if (associated(snow_ncells_for_bc)) deallocate(snow_ncells_for_bc)
    if (associated(sh2o_ncells_for_bc)) deallocate(sh2o_ncells_for_bc)
    if (associated(soil_ncells_for_bc)) deallocate(soil_ncells_for_bc)
    if (associated(snow_ncells_for_ss)) deallocate(snow_ncells_for_ss)
    if (associated(sh2o_ncells_for_ss)) deallocate(sh2o_ncells_for_ss)
    if (associated(soil_ncells_for_ss)) deallocate(soil_ncells_for_ss)

    
    ! Set the variables needed by a given governing equation from
    ! other governing equation

    nvars = 2
    call goveq_snow%AllocVarsFromOtherGEs(nvars)
    goveq_snow%var_ids_needed_from_other_goveqns (1) = VAR_TEMPERATURE
    goveq_snow%ids_of_other_goveqns              (1) = goveq_soil%id_in_list
    goveq_snow%is_bc_auxvar_type                 (1) = PETSC_TRUE
    goveq_snow%bc_auxvar_offset                  (1) = snow_bc_auxvar_offset_1
    goveq_snow%bc_auxvar_ncells                  (1) = snow_bc_auxvar_ncells_1
    goveq_snow%bc_auxvar_idx                     (1) = snow_bc_auxvar_idx_1
    goveq_snow%bc_auxvar_idx_of_other_goveqn     (1) = snow_bc_auxvar_idx_of_other_goveqn_1

    goveq_snow%var_ids_needed_from_other_goveqns (2) = VAR_THERMAL_COND
    goveq_snow%ids_of_other_goveqns              (2) = goveq_soil%id_in_list
    goveq_snow%is_bc_auxvar_type                 (2) = PETSC_TRUE
    goveq_snow%bc_auxvar_offset                  (2) = snow_bc_auxvar_offset_1
    goveq_snow%bc_auxvar_ncells                  (2) = snow_bc_auxvar_ncells_1
    goveq_snow%bc_auxvar_idx                     (2) = snow_bc_auxvar_idx_1
    goveq_snow%bc_auxvar_idx_of_other_goveqn     (2) = snow_bc_auxvar_idx_of_other_goveqn_1

    nvars = 2
    call goveq_sh2o%AllocVarsFromOtherGEs(nvars)
    goveq_sh2o%var_ids_needed_from_other_goveqns (1) = VAR_TEMPERATURE
    goveq_sh2o%ids_of_other_goveqns              (1) = goveq_soil%id_in_list
    goveq_sh2o%is_bc_auxvar_type                 (1) = PETSC_TRUE
    goveq_sh2o%bc_auxvar_offset                  (1) = sh2o_bc_auxvar_offset_1
    goveq_sh2o%bc_auxvar_ncells                  (1) = sh2o_bc_auxvar_ncells_1
    goveq_sh2o%bc_auxvar_idx                     (1) = sh2o_bc_auxvar_idx_1
    goveq_sh2o%bc_auxvar_idx_of_other_goveqn     (1) = sh2o_bc_auxvar_idx_of_other_goveqn_1

    goveq_sh2o%var_ids_needed_from_other_goveqns (2) = VAR_THERMAL_COND
    goveq_sh2o%ids_of_other_goveqns              (2) = goveq_soil%id_in_list
    goveq_sh2o%is_bc_auxvar_type                 (2) = PETSC_TRUE
    goveq_sh2o%bc_auxvar_offset                  (2) = snow_bc_auxvar_offset_1
    goveq_sh2o%bc_auxvar_ncells                  (2) = snow_bc_auxvar_ncells_1
    goveq_sh2o%bc_auxvar_idx                     (2) = snow_bc_auxvar_idx_1
    goveq_sh2o%bc_auxvar_idx_of_other_goveqn     (2) = snow_bc_auxvar_idx_of_other_goveqn_1

    nvars = 10
    call goveq_soil%AllocVarsFromOtherGEs(nvars)
    goveq_soil%var_ids_needed_from_other_goveqns (1) = VAR_TEMPERATURE
    goveq_soil%ids_of_other_goveqns              (1) = goveq_snow%id_in_list
    goveq_soil%is_bc_auxvar_type                 (1) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (1) = soil_bc_auxvar_offset_1
    goveq_soil%bc_auxvar_ncells                  (1) = soil_bc_auxvar_ncells_1
    goveq_soil%bc_auxvar_idx                     (1) = soil_bc_auxvar_idx_1
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (1) = soil_bc_auxvar_idx_of_other_goveqn_1

    goveq_soil%var_ids_needed_from_other_goveqns (2) = VAR_THERMAL_COND
    goveq_soil%ids_of_other_goveqns              (2) = goveq_snow%id_in_list
    goveq_soil%is_bc_auxvar_type                 (2) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (2) = soil_bc_auxvar_offset_1
    goveq_soil%bc_auxvar_ncells                  (2) = soil_bc_auxvar_ncells_1
    goveq_soil%bc_auxvar_idx                     (2) = soil_bc_auxvar_idx_1
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (2) = soil_bc_auxvar_idx_of_other_goveqn_1

    goveq_soil%var_ids_needed_from_other_goveqns (3) = VAR_FRAC
    goveq_soil%ids_of_other_goveqns              (3) = goveq_snow%id_in_list
    goveq_soil%is_bc_auxvar_type                 (3) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (3) = soil_bc_auxvar_offset_1
    goveq_soil%bc_auxvar_ncells                  (3) = soil_bc_auxvar_ncells_1
    goveq_soil%bc_auxvar_idx                     (3) = soil_bc_auxvar_idx_1
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (3) = soil_bc_auxvar_idx_of_other_goveqn_1

    goveq_soil%var_ids_needed_from_other_goveqns (4) = VAR_ACTIVE
    goveq_soil%ids_of_other_goveqns              (4) = goveq_snow%id_in_list
    goveq_soil%is_bc_auxvar_type                 (4) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (4) = soil_bc_auxvar_offset_1
    goveq_soil%bc_auxvar_ncells                  (4) = soil_bc_auxvar_ncells_1
    goveq_soil%bc_auxvar_idx                     (4) = soil_bc_auxvar_idx_1
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (4) = soil_bc_auxvar_idx_of_other_goveqn_1

    goveq_soil%var_ids_needed_from_other_goveqns (5) = VAR_DIST_UP
    goveq_soil%ids_of_other_goveqns              (5) = goveq_snow%id_in_list
    goveq_soil%is_bc_auxvar_type                 (5) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (5) = soil_bc_auxvar_offset_1
    goveq_soil%bc_auxvar_ncells                  (5) = soil_bc_auxvar_ncells_1
    goveq_soil%bc_auxvar_idx                     (5) = soil_bc_auxvar_idx_1
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (5) = soil_bc_auxvar_idx_of_other_goveqn_1

    goveq_soil%var_ids_needed_from_other_goveqns (6) = VAR_TEMPERATURE
    goveq_soil%ids_of_other_goveqns              (6) = goveq_sh2o%id_in_list
    goveq_soil%is_bc_auxvar_type                 (6) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (6) = soil_bc_auxvar_offset_2
    goveq_soil%bc_auxvar_ncells                  (6) = soil_bc_auxvar_ncells_2
    goveq_soil%bc_auxvar_idx                     (6) = soil_bc_auxvar_idx_2
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (6) = soil_bc_auxvar_idx_of_other_goveqn_2

    goveq_soil%var_ids_needed_from_other_goveqns (7) = VAR_THERMAL_COND
    goveq_soil%ids_of_other_goveqns              (7) = goveq_sh2o%id_in_list
    goveq_soil%is_bc_auxvar_type                 (7) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (7) = soil_bc_auxvar_offset_2
    goveq_soil%bc_auxvar_ncells                  (7) = soil_bc_auxvar_ncells_2
    goveq_soil%bc_auxvar_idx                     (7) = soil_bc_auxvar_idx_2
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (7) = soil_bc_auxvar_idx_of_other_goveqn_2

    goveq_soil%var_ids_needed_from_other_goveqns (8) = VAR_FRAC
    goveq_soil%ids_of_other_goveqns              (8) = goveq_sh2o%id_in_list
    goveq_soil%is_bc_auxvar_type                 (8) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (8) = soil_bc_auxvar_offset_2
    goveq_soil%bc_auxvar_ncells                  (8) = soil_bc_auxvar_ncells_2
    goveq_soil%bc_auxvar_idx                     (8) = soil_bc_auxvar_idx_2
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (8) = soil_bc_auxvar_idx_of_other_goveqn_2

    goveq_soil%var_ids_needed_from_other_goveqns (9) = VAR_ACTIVE
    goveq_soil%ids_of_other_goveqns              (9) = goveq_sh2o%id_in_list
    goveq_soil%is_bc_auxvar_type                 (9) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (9) = soil_bc_auxvar_offset_2
    goveq_soil%bc_auxvar_ncells                  (9) = soil_bc_auxvar_ncells_2
    goveq_soil%bc_auxvar_idx                     (9) = soil_bc_auxvar_idx_2
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (9) = soil_bc_auxvar_idx_of_other_goveqn_2

    goveq_soil%var_ids_needed_from_other_goveqns (10) = VAR_DZ
    goveq_soil%ids_of_other_goveqns              (10) = goveq_sh2o%id_in_list
    goveq_soil%is_bc_auxvar_type                 (10) = PETSC_TRUE
    goveq_soil%bc_auxvar_offset                  (10) = soil_bc_auxvar_offset_2
    goveq_soil%bc_auxvar_ncells                  (10) = soil_bc_auxvar_ncells_2
    goveq_soil%bc_auxvar_idx                     (10) = soil_bc_auxvar_idx_2
    goveq_soil%bc_auxvar_idx_of_other_goveqn     (10) = soil_bc_auxvar_idx_of_other_goveqn_2

    !
    ! Update boundary connections coupling various equations
    !
    snow_bc_for_soil_found = PETSC_FALSE
    sh2o_bc_for_soil_found = PETSC_FALSE
    soil_bc_for_snow_found = PETSC_FALSE
    soil_bc_for_sh2o_found = PETSC_FALSE

    cur_cond => goveq_snow%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       if (trim(cur_cond%name) == 'BC_from_soil_governing_equation') then
          snow_bc_for_soil => cur_cond
          snow_bc_for_soil_found = PETSC_TRUE
       endif
       cur_cond => cur_cond%next
    enddo

    cur_cond => goveq_sh2o%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       if (trim(cur_cond%name) == 'BC_from_soil_governing_equation') then
          sh2o_bc_for_soil => cur_cond
          sh2o_bc_for_soil_found = PETSC_TRUE
       endif
       cur_cond => cur_cond%next
    enddo
    
    cur_cond => goveq_soil%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       if (trim(cur_cond%name) == 'BC_from_snow_governing_equation') then
          soil_bc_for_snow => cur_cond
          soil_bc_for_snow_found = PETSC_TRUE
       endif

       if (trim(cur_cond%name) == 'BC_from_standing_surface_water_governing_equation') then
          soil_bc_for_sh2o => cur_cond
          soil_bc_for_sh2o_found = PETSC_TRUE
       endif
       cur_cond => cur_cond%next
    enddo

    !
    ! Get data
    !
    
    cur_conn_set => soil_bc_for_snow%conn_set
    allocate(soil_ids_for_snow_bc(cur_conn_set%num_connections))
    allocate(soil_dist_for_snow_bc(cur_conn_set%num_connections))
    do iconn = 1, cur_conn_set%num_connections
       soil_ids_for_snow_bc(iconn)  = cur_conn_set%id_dn(iconn)
       soil_dist_for_snow_bc(iconn) = cur_conn_set%dist_dn(iconn)
    enddo
    
    cur_conn_set => soil_bc_for_sh2o%conn_set
    allocate(soil_ids_for_sh2o_bc(cur_conn_set%num_connections))
    allocate(soil_dist_for_sh2o_bc(cur_conn_set%num_connections))
    do iconn = 1, cur_conn_set%num_connections
       soil_ids_for_sh2o_bc(iconn)  = cur_conn_set%id_dn(iconn)
       soil_dist_for_sh2o_bc(iconn) = cur_conn_set%dist_dn(iconn)
    enddo

    cur_conn_set => snow_bc_for_soil%conn_set
    allocate(snow_ids_for_soil_bc(cur_conn_set%num_connections))
    do iconn = 1, cur_conn_set%num_connections
       snow_ids_for_soil_bc(iconn)  = cur_conn_set%id_dn(iconn)
    enddo

    cur_conn_set => sh2o_bc_for_soil%conn_set
    allocate(sh2o_ids_for_soil_bc(cur_conn_set%num_connections))
    do iconn = 1, cur_conn_set%num_connections
       sh2o_ids_for_soil_bc(iconn)  = cur_conn_set%id_dn(iconn)
    enddo

    !
    ! Set data
    !

    cur_conn_set => soil_bc_for_snow%conn_set
    do iconn = 1, cur_conn_set%num_connections
       cur_conn_set%id_up(iconn) = snow_ids_for_soil_bc(iconn)
    enddo

    cur_conn_set => soil_bc_for_sh2o%conn_set
    do iconn = 1, cur_conn_set%num_connections
       cur_conn_set%id_up(iconn) = sh2o_ids_for_soil_bc(iconn)
    enddo

    cur_conn_set => snow_bc_for_soil%conn_set
    do iconn = 1, cur_conn_set%num_connections
       cur_conn_set%id_up(iconn)   = soil_ids_for_snow_bc(iconn)
       cur_conn_set%dist_up(iconn) = soil_dist_for_snow_bc(iconn)
    enddo

    cur_conn_set => sh2o_bc_for_soil%conn_set
    do iconn = 1, cur_conn_set%num_connections
       cur_conn_set%id_up(iconn)   = soil_ids_for_sh2o_bc(iconn)
       cur_conn_set%dist_up(iconn) = soil_dist_for_sh2o_bc(iconn)
    enddo

    deallocate(soil_ids_for_snow_bc)
    deallocate(soil_ids_for_sh2o_bc)
    deallocate(snow_ids_for_soil_bc)
    deallocate(sh2o_ids_for_soil_bc)
    deallocate(soil_dist_for_snow_bc)
    deallocate(soil_dist_for_sh2o_bc) 

    !call this%PrintInfo()

  end subroutine ThermalSOEThermalKSPSetup

  !------------------------------------------------------------------------
  subroutine ThermalSOESetSolnPrevCLM(this, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    PetscReal, intent (in)        :: data_1d(:)
    !
    PetscInt                      :: vsize
    PetscErrorCode                :: ierr
    PetscScalar, pointer          :: x_p(:)

    call VecGetSize(this%soln_prev_clm, vsize, ierr); CHKERRQ(ierr)

    if ( vsize /= size(data_1d)) then
       write(iulog,*) 'ThermalSOESetSolnPrevCLM: vsize /= size(data_1d)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayF90(this%soln_prev_clm, x_p, ierr); CHKERRQ(ierr)
    x_p(:) = data_1d(:)
    call VecRestoreArrayF90(this%soln_prev_clm, x_p, ierr); CHKERRQ(ierr)

  end subroutine ThermalSOESetSolnPrevCLM

  !------------------------------------------------------------------------
  subroutine ThermalSOEGetSoln(this, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    PetscReal, intent (out)       :: data_1d(:)
    !
    PetscInt                      :: vsize
    PetscErrorCode                :: ierr
    PetscScalar, pointer          :: x_p(:)

    call VecGetSize(this%soln, vsize, ierr); CHKERRQ(ierr)

    if ( vsize /= size(data_1d)) then
       write(iulog,*) 'ThermalSOESetSolnPrevCLM: vsize /= size(data_1d)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayF90(this%soln, x_p, ierr); CHKERRQ(ierr)
    data_1d(:) = x_p(:)
    call VecRestoreArrayF90(this%soln, x_p, ierr); CHKERRQ(ierr)

  end subroutine ThermalSOEGetSoln
  !------------------------------------------------------------------------

  subroutine ThermalSOESetRDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    use SystemOfEquationsThermalAuxMod, only : SOEThermalAuxSetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type)              :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal, pointer                         :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in

       call SOEThermalAuxSetRData(this%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case(AUXVAR_BC)
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)

       call SOEThermalAuxSetRData(this%aux_vars_bc, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case(AUXVAR_SS)
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
 
       call SOEThermalAuxSetRData(this%aux_vars_ss, var_type, &
            nauxvar, iauxvar_off, data_1d)

   case default
       write(iulog,*) 'ThermalSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalSOESetRDataFromCLM
  
  !------------------------------------------------------------------------
  subroutine ThermalSOESetIDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    use SystemOfEquationsThermalAuxMod, only : SOEThermalAuxSetIData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type)              :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscInt, pointer                          :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in

       call SOEThermalAuxSetIData(this%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case default
       write(iulog,*) 'ThermalSOESetIDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalSOESetIDataFromCLM
  
  !------------------------------------------------------------------------
  subroutine ThermalSOESetBDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants      , only : SOE_RE_ODE
    use MultiPhysicsProbConstants      , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants      , only : AUXVAR_BC
    use MultiPhysicsProbConstants      , only : AUXVAR_SS
    use SystemOfEquationsThermalAuxMod , only : SOEThermalAuxSetBData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    PetscInt, intent(in)          :: var_type
    PetscInt                      :: soe_auxvar_type
    PetscInt                      :: soe_auxvar_id
    PetscBool, pointer            :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                      :: iauxvar
    PetscInt                      :: iauxvar_off
    PetscInt                      :: nauxvar

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in

       call SOEThermalAuxSetBData(this%aux_vars_in, var_type, nauxvar, &
            iauxvar_off, data_1d)

    case default
       write(iulog,*) 'ThermalSOESetBDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalSOESetBDataFromCLM
  
  !------------------------------------------------------------------------
  subroutine ThermalSOEPreStepDT(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    !
    PetscErrorCode                :: ierr

    call VecCopy(this%soln_prev_clm, this%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(this%soln_prev_clm, this%soln     , ierr); CHKERRQ(ierr)

  end subroutine ThermalSOEPreStepDT

  !------------------------------------------------------------------------
  subroutine ThermalSOEPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use MultiPhysicsProbConstants           , only : SOE_THERMAL_TBASED
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSnowType , only : goveqn_thermal_ksp_temp_snow_type
    use GoveqnThermalKSPTemperatureSSWType  , only : goveqn_thermal_ksp_temp_ssw_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt :: offset

    offset = 0

    select case (this%itype)
    case(SOE_THERMAL_TBASED)

       ! 1) {soln_prev}  ---> sim_aux()
       call ThermalSOEUpdateAuxVarsTBased(this, this%soln_prev)

       ! 2) GE ---> GetFromSimAux()
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)

          class is (goveqn_thermal_ksp_temp_snow_type)
             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             offset = offset + cur_goveq%mesh%ncells_local

             call cur_goveq%UpdateInternalConn()
             call cur_goveq%UpdateBoundaryConn()

             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

          class is (goveqn_thermal_ksp_temp_ssw_type)
             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             offset = offset + cur_goveq%mesh%ncells_local

             call cur_goveq%UpdateInternalConn()
             call cur_goveq%UpdateBoundaryConn()

             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

          class is (goveqn_thermal_ksp_temp_soil_type)
             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             offset = offset + cur_goveq%mesh%ncells_local
             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

             call cur_goveq%UpdateBoundaryConn()
          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'ThermalSOEPreSolve: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalSOEPreSolve

  !------------------------------------------------------------------------
  subroutine ThermalSOEUpdateAuxVarsTBased(therm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE based on the input
    ! vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: therm_soe
    Vec                           :: X
    !
    ! !LOCAL VARIABLES:
    PetscInt                   :: dm_id
    PetscInt                   :: nDM
    DM, pointer                :: dms(:)
    Vec, pointer               :: X_subvecs(:)
    PetscInt                   :: size
    PetscInt                   :: offset
    PetscErrorCode             :: ierr

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(therm_soe%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(therm_soe%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(therm_soe%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Update the SoE auxvars
    offset = 0
    do dm_id = 1, nDM
       call ThermalSOESetAuxVars(therm_soe, AUXVAR_INTERNAL, VAR_TEMPERATURE, &
            X_subvecs(dm_id), offset)
       call VecGetSize(X_subvecs(dm_id), size, ierr); CHKERRQ(ierr)
       offset = offset + size
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(therm_soe%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)


  end subroutine ThermalSOEUpdateAuxVarsTBased

  !------------------------------------------------------------------------
  subroutine ThermalSOESetAuxVars(therm_soe, auxvar_type, var_type, &
       var_vec, offset)
    !
    ! !DESCRIPTION:
    ! Set values in SoE auxvars.
    !
    ! !USES:
    use MultiPhysicsProbConstants      , only : AUXVAR_INTERNAL
    use SystemOfEquationsThermalAuxMod , only : SOEThermalAuxSetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: therm_soe
    PetscInt                      :: auxvar_type
    PetscInt, intent(in)          :: var_type
    Vec                           :: var_vec
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer            :: var_p(:)
    PetscInt                      :: nauxvar
    PetscInt                      :: nvar
    PetscInt, optional            :: offset
    PetscInt                      :: iauxvar
    PetscInt                      :: iauxvar_off
    PetscErrorCode                :: ierr

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    select case(auxvar_type)
    case (AUXVAR_INTERNAL)

       nauxvar = size(therm_soe%aux_vars_in)

       call VecGetLocalSize(var_vec, nvar, ierr); CHKERRQ(ierr)

       if (nvar+iauxvar_off > nauxvar) then
          write(iulog,*) 'ThermalSOESetAuxVars: nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call VecGetArrayF90(var_vec, var_p, ierr); CHKERRQ(ierr)

       call SOEThermalAuxSetRData(therm_soe%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, var_p)

       call VecRestoreArrayF90(var_vec, var_p, ierr); CHKERRQ(ierr)

    case default
       write(iulog,*) 'ThermalSOESetAuxVars: auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select


  end subroutine ThermalSOESetAuxVars

  !------------------------------------------------------------------------

  subroutine ThermalSOEComputeRHS(this, ksp, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use GoverningEquationBaseType, only : goveqn_base_type
    use MultiPhysicsProbConstants, only : SOE_THERMAL_TBASED
    use GoverningEquationBaseType     , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    KSP                           :: ksp
    Vec                           :: B
    PetscErrorCode                :: ierr
    !
    !
    class(goveqn_base_type),pointer :: cur_goveq
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt                      :: dm_id
    PetscInt                      :: row
    PetscInt                      :: col
    PetscInt                      :: nDM
    DM, pointer                   :: dms(:)
    Vec, pointer                  :: B_subvecs(:)

    
    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(B_subvecs(nDM))

    ! Get vectors for individual GEs
    call DMCompositeGetAccessArray(this%dm, B, nDM, PETSC_NULL_INTEGER, &
         B_subvecs, ierr); CHKERRQ(ierr)

    ! 1) GE: UpdateAuxVars
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       call cur_goveq%UpdateAuxVarsIntrn()

       cur_goveq => cur_goveq%next
    enddo
    
    ! 2) GE_1 <---> GE_2 exchange AuxVars()
    do row = 1,nDM
       do col = row+1,nDM
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call ThermalSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
          call ThermalSOEGovEqnExchangeAuxVars(cur_goveq_2, cur_goveq_1)
       enddo
    enddo

    ! 3) Call ComputeRHS
    dm_id = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       dm_id = dm_id + 1

       call VecZeroEntries(B_subvecs(dm_id), ierr); CHKERRQ(ierr)       

       select type (cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          !call cur_goveq%UpdateBoundaryConn()
       class is (goveqn_thermal_ksp_temp_ssw_type) 
       class is (goveqn_thermal_ksp_temp_soil_type)          
          call cur_goveq%UpdateBoundaryConn()
       end select

       call cur_goveq%ComputeRHS(B_subvecs(dm_id),    &
            ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! Restore vectors for individual GEs
    call DMCompositeRestoreAccessArray(this%dm, B, nDM, PETSC_NULL_INTEGER, &
         B_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(B_subvecs)

  end subroutine ThermalSOEComputeRHS

  !------------------------------------------------------------------------

  subroutine ThermalSOEComputeOperators(this, ksp, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc KSP.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_type) :: this
    KSP                        :: ksp
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                      :: row
    PetscInt                      :: col
    PetscInt                      :: nDM

    IS,pointer                    :: is(:)
    DM, pointer                   :: dms(:)
    Mat, pointer                  :: B_submats(:,:)
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%dm, dms, ierr); CHKERRQ(ierr)

    ! Initialize the matrix
    call MatZeroEntries(B, ierr); CHKERRQ(ierr)

    ! Get submatrices
    allocate(is(nDM))
    allocate(B_submats(nDM,nDM))
    call DMCompositeGetLocalISs(this%dm, is, ierr); CHKERRQ(ierr)
    do row = 1,nDM
      do col = 1,nDM
        call MatGetLocalSubMatrix(B, is(row), is(col), B_submats(row,col), ierr); CHKERRQ(ierr)
      enddo
    enddo

    ! Operators and OperatorsOffDiag
    row = 0
    cur_goveq_1 => this%goveqns
    do
       if (.not.associated(cur_goveq_1)) exit

       row = row + 1

       call cur_goveq_1%ComputeOperatorsDiag(B_submats(row,row),    &
                                             B_submats(row,row),    &
                                             ierr); CHKERRQ(ierr)

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          call cur_goveq_1%ComputeOperatorsOffDiag(B_submats(row,col),    &
                                                   B_submats(row,col),    &
                                                   cur_goveq_2%id,        &
                                                   cur_goveq_2%id_in_list,&
                                                   ierr); CHKERRQ(ierr)

          call cur_goveq_2%ComputeOperatorsOffDiag(B_submats(col,row),    &
                                                   B_submats(col,row),    &
                                                   cur_goveq_1%id,        &
                                                   cur_goveq_1%id_in_list,&
                                                   ierr); CHKERRQ(ierr)

          cur_goveq_2 => cur_goveq_2%next
       enddo

       cur_goveq_1 => cur_goveq_1%next
    enddo

    ! Restore submatrices
    do row = 1,nDM
      do col = 1,nDM
        call MatRestoreLocalSubMatrix(B, is(row), is(col), B_submats(row,col), ierr); CHKERRQ(ierr)
      enddo
    enddo

    ! Destroy IS
    do row = 1,nDM
      call ISDestroy(is(row), ierr); CHKERRQ(ierr)
    enddo

    ! Assemble matrix
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

    ! Free memory
    deallocate(dms         )
    deallocate(is          )
    deallocate(B_submats   )
  end subroutine ThermalSOEComputeOperators

  !------------------------------------------------------------------------

  subroutine ThermalSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSnowType , only : goveqn_thermal_ksp_temp_snow_type
    use GoveqnThermalKSPTemperatureSSWType  , only : goveqn_thermal_ksp_temp_ssw_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    use ConnectionSetType                   , only : connection_set_type
    use ConditionType                       , only : condition_type
    use MultiPhysicsProbConstants           , only : VAR_DZ
    use MultiPhysicsProbConstants           , only : VAR_THERMAL_COND
    use ThermalKSPTemperatureSnowAuxMod     , only : ThermKSPTempSnowAuxVarSetRValues
    use ThermalKSPTemperatureSnowAuxMod     , only : ThermKSPTempSnowAuxVarGetRValues
    use ThermalKSPTemperatureSSWAuxMod      , only : ThermKSPTempSSWAuxVarSetRValues
    use ThermalKSPTemperatureSSWAuxMod      , only : ThermKSPTempSSWAuxVarGetRValues
    use ThermalKSPTemperatureSoilAuxMod     , only : ThermKSPTempSoilAuxVarSetRValues
    use ThermalKSPTemperatureSoilAuxMod     , only : ThermKSPTempSoilAuxVarGetRValues
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)   , pointer :: cur_goveq_1
    class(goveqn_base_type)   , pointer :: cur_goveq_2
    !
    type(connection_set_type) , pointer :: cur_conn_set_2
    type(condition_type)      , pointer :: cur_cond_2
    PetscInt                            :: idx
    PetscInt, pointer                   :: ids(:)
    PetscInt                            :: iauxvar
    PetscInt                            :: ivar
    PetscInt                            :: var_type
    PetscInt                            :: bc_idx
    PetscInt                            :: bc_offset
    PetscInt                            :: bc_auxvar_idx_of_other_goveqn
    PetscReal                           :: var_value
    PetscReal, pointer                  :: var_values(:)
    PetscBool                           :: bc_found
    PetscBool                           :: bc_type

    do ivar = 1,cur_goveq_1%nvars_needed_from_other_goveqns

       ! Does cur_goveq_1 needs ivar-th variable from cur_goveq_2?
       if (cur_goveq_1%ids_of_other_goveqns(ivar) == &
            cur_goveq_2%id_in_list) then

          var_type                      = cur_goveq_1%var_ids_needed_from_other_goveqns(ivar)
          bc_type                       = cur_goveq_1%is_bc_auxvar_type(ivar)
          bc_offset                     = cur_goveq_1%bc_auxvar_offset(ivar)
          bc_auxvar_idx_of_other_goveqn = cur_goveq_1%bc_auxvar_idx_of_other_goveqn(ivar)

          if (.not.bc_type) then
             
             write(iulog,*) 'ThermalSOEGovEqnExchangeAuxVars: Extend code to ' // &
                  'exchange non-boundary condition data'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          else

             ! Get the appropriate pointer to the BC from cur_goveq_2
             bc_idx = 1
             bc_found = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit
                cur_conn_set_2 => cur_cond_2%conn_set
                
                ! Is this the appropriate BC?
                if (bc_idx == bc_auxvar_idx_of_other_goveqn) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'ThermalSOEGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_conn_set_2%num_connections /= &
                  cur_goveq_1%bc_auxvar_ncells(ivar) ) then
                write(iulog,*) 'conn_set_2%num_connections        = ', cur_conn_set_2%num_connections
                write(iulog,*) 'cur_goveq_1%bc_auxvar_ncells(ivar)= ', cur_goveq_1%bc_auxvar_ncells(ivar)
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             allocate(ids       (cur_goveq_1%bc_auxvar_ncells(ivar)))
             allocate(var_values(cur_goveq_1%bc_auxvar_ncells(ivar)))

             ! Save the IDs to get the data from
             do iauxvar = 1, cur_goveq_1%bc_auxvar_ncells(ivar)
                ids(iauxvar) = cur_conn_set_2%id_dn(iauxvar)
             enddo

             ! Get the data
             select type(cur_goveq_2)
             class is (goveqn_thermal_ksp_temp_snow_type)
                call ThermKSPTempSnowAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     size(cur_goveq_2%aux_vars_in), ids, var_values)
                
             class is (goveqn_thermal_ksp_temp_ssw_type)
                call ThermKSPTempSSWAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     size(cur_goveq_2%aux_vars_in), ids, var_values)

             class is (goveqn_thermal_ksp_temp_soil_type)
                call ThermKSPTempSoilAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     size(cur_goveq_2%aux_vars_in), ids, var_values)

             class default
                write(iulog,*)'ThermalSOEGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             ! Save the IDs to set the data to
             do iauxvar = 1, cur_goveq_1%bc_auxvar_ncells(ivar)
                ids(iauxvar) = iauxvar + bc_offset
             enddo

             ! Set the data
             select type(cur_goveq_1)
             class is (goveqn_thermal_ksp_temp_snow_type)
                call ThermKSPTempSnowAuxVarSetRValues(cur_goveq_1%aux_vars_bc, var_type, &
                     size(cur_goveq_1%aux_vars_bc), ids, var_values)
                
             class is (goveqn_thermal_ksp_temp_ssw_type)
                call ThermKSPTempSSWAuxVarSetRValues(cur_goveq_1%aux_vars_bc, var_type, &
                     size(cur_goveq_1%aux_vars_bc), ids, var_values)

             class is (goveqn_thermal_ksp_temp_soil_type)
                call ThermKSPTempSoilAuxVarSetRValues(cur_goveq_1%aux_vars_bc, var_type, &
                     size(cur_goveq_1%aux_vars_bc), ids, var_values)

             class default
                write(iulog,*)'ThermalSOEGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             deallocate(ids       )
             deallocate(var_values)

          endif

       endif
    enddo

  end subroutine ThermalSOEGovEqnExchangeAuxVars

#endif

end module SystemOfEquationsThermalType
