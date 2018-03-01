module GoveqnThermalKSPTemperatureSoilType

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                       , only : iulog
  use mpp_abortutils                   , only : endrun
  use mpp_shr_log_mod                  , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType        , only : goveqn_base_type
  use ThermalKSPTemperatureSoilAuxType , only : therm_ksp_temp_soil_auxvar_type
  use SystemOfEquationsThermalAuxType  , only : sysofeqns_thermal_auxvar_type
  use petscsys
  use petscvec
  use petscmat
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_thermal_ksp_temp_soil_type
     type (therm_ksp_temp_soil_auxvar_type), pointer :: aux_vars_in(:)  ! Internal state.
     type (therm_ksp_temp_soil_auxvar_type), pointer :: aux_vars_bc(:)  ! Boundary conditions.
     type (therm_ksp_temp_soil_auxvar_type), pointer :: aux_vars_ss(:)  ! Source-sink.

     PetscInt, pointer                             :: soe_auxvars_bc_offset (:) ! SoE auxvar offset corresponding to BCs
     PetscInt, pointer                             :: soe_auxvars_ss_offset (:) ! SoE auxvar offset corresponding to SSs

   contains

     procedure, public :: Setup                   => ThermKSPTempSoilSetup
     procedure, public :: AllocateAuxVars         => ThermKSPTempSoilAllocateAuxVars

     procedure, public :: GetFromSOEAuxVarsIntrn  => ThermKSPTempSoilGetFromSOEAuxVarsIntrn
     procedure, public :: GetFromSOEAuxVarsBC     => ThermKSPTempSoilGetFromSOEAuxVarsBC
     procedure, public :: GetFromSOEAuxVarsSS     => ThermKSPTempSoilGetFromSOEAuxVarsSS
     procedure, public :: GetDataFromSOEAuxVar    => ThermKSPTempSoilGetDataFromSOEAuxVar

     procedure, public :: SetDataInSOEAuxVar      => ThermKSPTempSoilSetDataInSOEAuxVar
     procedure, public :: SetSOEAuxVarOffsets     => ThermKSPTempSoilSetSOEAuxVarOffsets
     procedure, public :: UpdateBoundaryConn      => ThermKSPTempSoilUpdateBoundaryConn
     procedure, public :: UpdateAuxVarsIntrn      => ThermKSPTempSoilUpdateAuxVarsIntrn

     procedure, public :: ComputeRHS              => ThermKSPTempSoilComputeRHS
     procedure, public :: ComputeOperatorsDiag    => ThermKSPTempSoilComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag => ThermKSPTempSoilComputeOperatorsOffDiag

  end type goveqn_thermal_ksp_temp_soil_type

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation for Thermal equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_TBASED
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this

    call this%Create()

    this%name       = "Soil thermal equation based on temperature"
    this%id         = GE_THERM_SOIL_TBASED
    this%mesh_itype = MESH_CLM_SOIL_COL

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_ss)

    nullify(this%soe_auxvars_bc_offset)
    nullify(this%soe_auxvars_ss_offset)

  end subroutine ThermKSPTempSoilSetup

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !   + Boundary condtions,
    !   + Source-sink condition.
    !
    ! !USES:
    use ConditionType, only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: ncells_cond
    PetscInt                                 :: ncond
    PetscInt                                 :: icond

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    do icond = 1,this%mesh%ncells_all
       call this%aux_vars_in(icond)%Init()
    enddo

    ! Allocate memory and initialize aux vars: For boundary connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_bc(ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_bc(icond)%Init()
    enddo
    allocate(this%soe_auxvars_bc_offset(ncond))
    this%soe_auxvars_bc_offset(:) = -1

    ! Allocate memory and initialize aux vars: For source sink connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_ss(ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_ss(icond)%Init()
    enddo
    allocate(this%soe_auxvars_ss_offset(ncond))
    this%soe_auxvars_ss_offset(:) = -1

  end subroutine ThermKSPTempSoilAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilGetDataFromSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    !
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type), intent(inout)        :: this
    PetscInt, intent(in)                                           :: soe_avar_type
    type (sysofeqns_thermal_auxvar_type), dimension(:), intent(in) :: soe_avars
    PetscInt, intent(in), optional                                 :: offset
    !
    ! !LOCAL VARIABLES
    PetscInt                                                       :: iauxvar_off

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       if (present(offset)) then
          iauxvar_off = offset
       else
          iauxvar_off = 0
       endif
       call ThermKSPTempSoilGetFromSOEAuxVarsIntrn(this, soe_avars, iauxvar_off)
    case (AUXVAR_BC)
       call ThermKSPTempSoilGetFromSOEAuxVarsBC(this, soe_avars)
    case (AUXVAR_SS)
       call ThermKSPTempSoilGetFromSOEAuxVarsSS(this, soe_avars)
    case default
       write(iulog,*) 'ThermKSPTempSoilGetDataFromSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSoilGetDataFromSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilGetFromSOEAuxVarsIntrn(this, soe_avars, offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for internal nodes
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) , intent(inout)            :: this
    type(sysofeqns_thermal_auxvar_type)      , dimension(:), intent(in) :: soe_avars
    PetscInt                                 , intent(in)               :: offset
    !
    ! LOCAL VARIABLES
    PetscInt :: iauxvar
    PetscInt :: nauxvar

    nauxvar = size(this%aux_vars_in)
    if( nauxvar > size(soe_avars) ) then
       write(iulog,*) 'size(this%aux_vars_in) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, nauxvar

       this%aux_vars_in(iauxvar)%temperature    = soe_avars(iauxvar+offset)%temperature
       this%aux_vars_in(iauxvar)%liq_areal_den  = soe_avars(iauxvar+offset)%liq_areal_den
       this%aux_vars_in(iauxvar)%ice_areal_den  = soe_avars(iauxvar+offset)%ice_areal_den
       this%aux_vars_in(iauxvar)%snow_water     = soe_avars(iauxvar+offset)%snow_water
       this%aux_vars_in(iauxvar)%num_snow_layer = soe_avars(iauxvar+offset)%num_snow_layer
       this%aux_vars_in(iauxvar)%tuning_factor  = soe_avars(iauxvar+offset)%tuning_factor
       this%aux_vars_in(iauxvar)%frac           = soe_avars(iauxvar+offset)%frac
       this%aux_vars_in(iauxvar)%dz             = soe_avars(iauxvar+offset)%dz

    enddo

  end subroutine ThermKSPTempSoilGetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilGetFromSOEAuxVarsBC(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) , intent(inout)            :: this
    type(sysofeqns_thermal_auxvar_type)      , dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                             :: iauxvar
    integer                             :: iauxvar_off 
    integer                             :: iconn
    integer                             :: nauxVar_ge
    integer                             :: nauxVar_soe
    integer                             :: condition_id
    integer                             :: sum_conn
    integer                             :: cell_id
    type(condition_type)      , pointer :: cur_cond
    type(connection_set_type) , pointer :: cur_conn_set
    character(len=256)                  :: string

    nauxVar_ge = size(this%aux_vars_bc)

    nauxVar_soe = size(soe_avars)
    if( nauxVar_ge > nauxVar_soe ) then
       write(iulog,*) 'size(this%aux_vars_bc) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_cond => cur_cond%next
          cycle
       endif

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = this%soe_auxvars_bc_offset(condition_id)

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          select case(cur_cond%itype)
          case (COND_HEAT_FLUX)
             this%aux_vars_bc(sum_conn)%condition_value =  &
                  soe_avars(iconn + iauxvar_off)%condition_value

             ! H - dH/dT * T
             cur_cond%value(iconn) =                               &
                  soe_avars(iconn + iauxvar_off)%condition_value - &
                  soe_avars(iconn + iauxvar_off)%dhsdT *           &
                  this%aux_vars_in(cell_id)%temperature

             this%aux_vars_bc(sum_conn)%dhsdT = &
                  soe_avars(iconn + iauxvar_off)%dhsdT

             this%aux_vars_bc(sum_conn)%frac = &
                  soe_avars(iconn + iauxvar_off)%frac

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing

          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilGetFromSOEAuxVarsBC


  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilGetFromSOEAuxVarsSS(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type), intent(inout)       :: this
    type(sysofeqns_thermal_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                             :: iauxvar
    integer                             :: iauxvar_off
    integer                             :: iconn
    integer                             :: nauxVar_ge
    integer                             :: nauxVar_soe
    integer                             :: condition_id
    integer                             :: sum_conn
    PetscReal                           :: var_value
    type(condition_type)      , pointer :: cur_cond
    type(connection_set_type) , pointer :: cur_conn_set
    character(len=256)                  :: string

    nauxVar_ge  = size(this%aux_vars_ss)
    nauxVar_soe = size(soe_avars)

    if( nauxVar_ge > nauxVar_soe ) then
       write(iulog,*) 'size(this%aux_vars_ss) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = this%soe_auxvars_ss_offset(condition_id)

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          select case(cur_cond%itype)
          case (COND_HEAT_RATE)

             var_value = soe_avars(iconn + iauxvar_off)%condition_value

             this%aux_vars_ss(sum_conn)%condition_value = var_value
             cur_cond%value(iconn)                      = var_value
          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilGetFromSOEAuxVarsSS

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilSetDataInSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies data GE auxiliary variable into SoE auxiliary variable.
    ! This is done as a part of post solve.
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type)                        :: this
    PetscInt                                                        :: soe_avar_type
    type (sysofeqns_thermal_auxvar_type), dimension(:), intent(out) :: soe_avars
    PetscInt, optional                                              :: offset
    !
    ! !LOCAL VARIABLES
    PetscInt                                                       :: iauxvar
    PetscInt                                                       :: iauxvar_off

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       this%aux_vars_in => this%aux_vars_in

       if ( size(this%aux_vars_in) > size(soe_avars)) then
          write(iulog,*) 'size(this%aux_vars_in) > size(soe_avars)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       do iauxvar = 1, size(this%aux_vars_in)
          if (this%mesh%is_active(iauxvar)) then
             soe_avars(iauxvar+iauxvar_off)%temperature =  &
                  this%aux_vars_in(iauxvar)%temperature
          endif
       enddo

    case default
       write(iulog,*) 'RichardsODESetDataInSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSoilSetDataInSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilSetSOEAuxVarOffsets(this, bc_offset_count, bc_offsets, &
    ss_offset_count, ss_offsets)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    PetscInt                                 :: bc_offset_count
    PetscInt, pointer                        :: bc_offsets(:)
    PetscInt                                 :: ss_offset_count
    PetscInt, pointer                        :: ss_offsets(:)
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: cond_count

    call this%GetNConditionsExcptCondItype(COND_BC, -1, cond_count)
    if (bc_offset_count > cond_count) then
       write(iulog,*) 'ERROR: bc_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call this%GetNConditionsExcptCondItype(COND_SS, -1, cond_count)
    if (ss_offset_count > cond_count) then
       write(iulog,*) 'ERROR: ss_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cond_count = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          this%soe_auxvars_bc_offset(cond_count) = bc_offsets(cond_count)
       endif
       cur_cond => cur_cond%next
    enddo

    cond_count = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       this%soe_auxvars_ss_offset(cond_count) = ss_offsets(cond_count)
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilSetSOEAuxVarOffsets

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilUpdateAuxVarsIntrn(this)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt :: ghosted_id

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_local
       call this%aux_vars_in(ghosted_id)%AuxVarCompute( &
            this%mesh%dz(ghosted_id), this%mesh%vol(ghosted_id))
    enddo

  end subroutine ThermKSPTempSoilUpdateAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    Vec                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscReal, dimension(:), pointer         :: b_p
    
    call VecGetArrayF90(B, b_p, ierr); CHKERRQ(ierr);

    call ThermalKSPTempSoilAccum(this, b_p)
    call ThermalKSPTempSoilDivergence(this, b_p)

    call VecRestoreArrayF90(B, b_p, ierr); CHKERRQ(ierr)

  end subroutine ThermKSPTempSoilComputeRHS

  !------------------------------------------------------------------------

  subroutine ThermalKSPTempSoilAccum(geq_soil, b_p)
    !
    ! !DESCRIPTION:
    !
    ! \frac{\partial \xi}{\partial dt} dV
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: geq_soil
    PetscReal, dimension(:), intent(out)     :: b_p
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: cell_id
    PetscReal                                :: T
    PetscReal                                :: heat_cap
    PetscReal                                :: tfactor
    PetscReal                                :: vol
    PetscReal                                :: dt

    dt = geq_soil%dtime

    ! Interior cells
    do cell_id = 1, geq_soil%mesh%ncells_local

       if (geq_soil%aux_vars_in(cell_id)%is_active) then

          T        = geq_soil%aux_vars_in(cell_id)%temperature
          heat_cap = geq_soil%aux_vars_in(cell_id)%heat_cap_pva
          tfactor  = geq_soil%aux_vars_in(cell_id)%tuning_factor
          vol      = geq_soil%mesh%vol(cell_id)

#ifdef MATCH_CLM_FORMULATION
          b_p(cell_id) = T
#else
          b_p(cell_id) = heat_cap*vol/(dt*tfactor)*T
#endif

       endif
    enddo

  end subroutine ThermalKSPTempSoilAccum

  !------------------------------------------------------------------------

  subroutine ThermalKSPTempSoilDivergence(geq_soil, b_p)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: geq_soil
    PetscReal, dimension(:), intent(out)     :: b_p
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscReal                                :: flux
    PetscReal                                :: area
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set
    PetscReal                                :: factor
    PetscReal                                :: T
    PetscReal                                :: dt
    PetscReal                                :: heat_cap
    PetscReal                                :: tfactor
    PetscReal                                :: vol
    PetscBool                                :: is_bc_sh2o
    PetscReal                                :: dist, dist_up, dist_dn
    PetscReal                                :: therm_cond_aveg, therm_cond_up, therm_cond_dn
    PetscReal                                :: T_up, T_dn

    dt = geq_soil%dtime

    ! Interior cells
    cur_conn_set => geq_soil%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          if ((.not.geq_soil%aux_vars_in(cell_id_up)%is_active) .or. &
              (.not.geq_soil%aux_vars_in(cell_id_dn)%is_active)) cycle

          call DiffHeatFlux(geq_soil%aux_vars_in(cell_id_up)%temperature,  &
                            geq_soil%aux_vars_in(cell_id_up)%therm_cond,   &
                            geq_soil%aux_vars_in(cell_id_dn)%temperature,  &
                            geq_soil%aux_vars_in(cell_id_dn)%therm_cond,   &
                            cur_conn_set%conn(iconn)%GetDistUp(),                   &
                            cur_conn_set%conn(iconn)%GetDistDn(),                   &
                            flux                                           &
                            )

          area = cur_conn_set%conn(iconn)%GetArea()

          T        = geq_soil%aux_vars_in(cell_id_up)%temperature
          heat_cap = geq_soil%aux_vars_in(cell_id_up)%heat_cap_pva
          tfactor  = geq_soil%aux_vars_in(cell_id_up)%tuning_factor
          vol      = geq_soil%mesh%vol(cell_id_up)
#ifdef MATCH_CLM_FORMULATION
          factor =  (dt*tfactor)/(heat_cap*vol)
#else
          factor = 1.d0
#endif
          b_p(cell_id_up) = b_p(cell_id_up) + cnfac*flux*area*factor


          T        = geq_soil%aux_vars_in(cell_id_dn)%temperature
          heat_cap = geq_soil%aux_vars_in(cell_id_dn)%heat_cap_pva
          tfactor  = geq_soil%aux_vars_in(cell_id_dn)%tuning_factor
          vol      = geq_soil%mesh%vol(cell_id_dn)
#ifdef MATCH_CLM_FORMULATION
          factor =  (dt*tfactor)/(heat_cap*vol)
#else
          factor = 1.d0
#endif
          b_p(cell_id_dn) = b_p(cell_id_dn) - cnfac*flux*area*factor

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => geq_soil%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       is_bc_sh2o = PETSC_FALSE
       do ieqn = 1, cur_cond%num_other_goveqs
          if (cur_cond%itype_of_other_goveqs(ieqn) == GE_THERM_SSW_TBASED) then
             is_bc_sh2o = PETSC_TRUE
          endif
       enddo

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
          sum_conn = sum_conn + 1

          if (.not.geq_soil%aux_vars_in(cell_id )%is_active) cycle
          
          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ)

             if (.not.geq_soil%aux_vars_bc(sum_conn)%is_active) cycle

                if (is_bc_sh2o) then
                   dist_up       = cur_conn_set%conn(iconn)%GetIDUp()
                   dist_dn       = cur_conn_set%conn(iconn)%GetIDDn()
                   dist          = dist_up + dist_dn

                   therm_cond_up = geq_soil%aux_vars_bc(sum_conn)%therm_cond
                   therm_cond_dn = geq_soil%aux_vars_in(cell_id )%therm_cond

                   T_up = geq_soil%aux_vars_bc(sum_conn)%temperature
                   T_dn = geq_soil%aux_vars_in(cell_id )%temperature

                   dist_dn = geq_soil%aux_vars_in(cell_id)%dz/2.d0
                   therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                        (therm_cond_up*dist_dn + therm_cond_dn*dist_up)
                   dist = dist_dn + max(1.d-6, dist_up*2.d0)/2.d0

                   flux = -therm_cond_aveg*(T_up - T_dn)/(dist)

                else

                   call DiffHeatFlux(geq_soil%aux_vars_bc(sum_conn)%temperature, &
                                     geq_soil%aux_vars_bc(sum_conn)%therm_cond,  &
                                     geq_soil%aux_vars_in(cell_id )%temperature, &
                                     geq_soil%aux_vars_in(cell_id )%therm_cond,  &
                                     cur_conn_set%conn(iconn)%GetDistUp(),                &
                                     cur_conn_set%conn(iconn)%GetDistDn(),                &
                                     flux                                        &
                                    )
                endif

                T        = geq_soil%aux_vars_in(cell_id)%temperature
                heat_cap = geq_soil%aux_vars_in(cell_id)%heat_cap_pva
                tfactor  = geq_soil%aux_vars_in(cell_id)%tuning_factor
                vol      = geq_soil%mesh%vol(cell_id)
#ifdef MATCH_CLM_FORMULATION
                factor =  (dt*tfactor)/(heat_cap*vol)
#else
                factor = 1.d0
#endif

                b_p(cell_id) = b_p(cell_id) - &
                     geq_soil%aux_vars_bc(sum_conn)%frac* &
                     cnfac * flux * area * factor

          case (COND_HEAT_FLUX)             
             area = cur_conn_set%conn(iconn)%GetArea()

             T        = geq_soil%aux_vars_in(cell_id)%temperature
             heat_cap = geq_soil%aux_vars_in(cell_id)%heat_cap_pva
             tfactor  = geq_soil%aux_vars_in(cell_id)%tuning_factor
             vol      = geq_soil%mesh%vol(cell_id)
#ifdef MATCH_CLM_FORMULATION
             factor =  (dt*tfactor)/(heat_cap*vol)
#else
             factor = 1.d0
#endif
             b_p(cell_id) = b_p(cell_id) + &
                  factor*cur_cond%value(iconn)               * &
                  geq_soil%aux_vars_bc(sum_conn)%frac * &
                  area

          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown boundary condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

       enddo
       cur_cond => cur_cond%next
    enddo    

    ! Source-sink cells
    cur_cond => geq_soil%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          if ((.not.geq_soil%aux_vars_in(cell_id)%is_active)) cycle

          select case(cur_cond%itype)
          case(COND_HEAT_RATE)
             T        = geq_soil%aux_vars_in(cell_id)%temperature
             heat_cap = geq_soil%aux_vars_in(cell_id)%heat_cap_pva
             tfactor  = geq_soil%aux_vars_in(cell_id)%tuning_factor
             vol      = geq_soil%mesh%vol(cell_id)
#ifdef MATCH_CLM_FORMULATION
             factor =  (dt*tfactor)/(heat_cap*vol)
#else
             factor = 1.d0
#endif

             b_p(cell_id) = b_p(cell_id) + cur_cond%value(iconn)*factor
          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown source-sink condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

     enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermalKSPTempSoilDivergence

  !------------------------------------------------------------------------

  subroutine DiffHeatFlux(T_up, therm_cond_up, T_dn, therm_cond_dn, &
                          dist_up, dist_dn, flux)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: T_up
    PetscReal, intent(in)  :: therm_cond_up
    PetscReal, intent(in)  :: T_dn
    PetscReal, intent(in)  :: therm_cond_dn
    PetscReal, intent(in)  :: dist_up
    PetscReal, intent(in)  :: dist_dn
    PetscReal, intent(out) :: flux
    !
    ! !LOCAL VARIABLES
    PetscReal :: therm_cond

    ! Distance weighted harmonic average
    therm_cond = therm_cond_up*therm_cond_dn*(dist_up+dist_dn)/ &
                 (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

    flux = -therm_cond*(T_up - T_dn)/(dist_up + dist_dn)

  end subroutine DiffHeatFlux

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    Mat                                      :: A
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id_dn
    PetscReal                                :: dist
    PetscReal                                :: dist_up
    PetscReal                                :: dist_dn
    PetscReal                                :: area
    PetscReal                                :: therm_cond_aveg
    PetscReal                                :: therm_cond_up
    PetscReal                                :: therm_cond_dn
    PetscReal                                :: heat_cap
    PetscReal                                :: tfactor
    PetscReal                                :: vol
    PetscReal                                :: dhsdT
    PetscReal                                :: frac
    PetscReal                                :: dt
    PetscReal                                :: value
    type(connection_set_type), pointer       :: cur_conn_set
    type(condition_type),pointer             :: cur_cond
    PetscReal :: factor
    PetscReal                                :: T
    PetscBool :: is_bc_sh2o

    dt = this%dtime

    ! Diagonal term
    do cell_id = 1, this%mesh%ncells_local

       heat_cap   = this%aux_vars_in(cell_id)%heat_cap_pva
       tfactor    = this%aux_vars_in(cell_id)%tuning_factor
       vol        = this%mesh%vol(cell_id)

       if (this%aux_vars_in(cell_id)%is_active) then
#ifdef MATCH_CLM_FORMULATION
          value = 1.d0
#else
          value = heat_cap*vol/(dt*tfactor)
#endif
       else
          value = 1.d0
       endif

       call MatSetValuesLocal(B, 1, cell_id-1, 1, cell_id-1, value, &
            ADD_VALUES, ierr); CHKERRQ(ierr)
    enddo
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          if ((.not.this%aux_vars_in(cell_id_up)%is_active) .or. &
              (.not.this%aux_vars_in(cell_id_dn)%is_active)) cycle

          area          = cur_conn_set%conn(iconn)%GetArea()
          dist_up       = cur_conn_set%conn(iconn)%GetDistUp()
          dist_dn       = cur_conn_set%conn(iconn)%GetDistDn()
          dist          = dist_up + dist_dn

          therm_cond_up = this%aux_vars_in(cell_id_up)%therm_cond
          therm_cond_dn = this%aux_vars_in(cell_id_dn)%therm_cond

          ! Distance weighted harmonic average
          therm_cond_aveg = therm_cond_up*therm_cond_dn*dist/ &
                 (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

          value = (1.d0 - cnfac)*therm_cond_aveg/dist*area

          T        = this%aux_vars_in(cell_id_up)%temperature
          heat_cap = this%aux_vars_in(cell_id_up)%heat_cap_pva
          tfactor  = this%aux_vars_in(cell_id_up)%tuning_factor
          vol      = this%mesh%vol(cell_id_up)
#ifdef MATCH_CLM_FORMULATION
          factor =  (dt*tfactor)/(heat_cap*vol)
#else
          factor = 1.d0
#endif

          call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_up-1,  value*factor, ADD_VALUES, ierr); CHKERRQ(ierr)
          call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_dn-1, -value*factor, ADD_VALUES, ierr); CHKERRQ(ierr)

          T        = this%aux_vars_in(cell_id_dn)%temperature
          heat_cap = this%aux_vars_in(cell_id_dn)%heat_cap_pva
          tfactor  = this%aux_vars_in(cell_id_dn)%tuning_factor
          vol      = this%mesh%vol(cell_id_dn)
#ifdef MATCH_CLM_FORMULATION
          factor =  (dt*tfactor)/(heat_cap*vol)
#else
          factor = 1.d0
#endif

          call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -value*factor, ADD_VALUES, ierr); CHKERRQ(ierr)
          call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_dn-1,  value*factor, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       is_bc_sh2o = PETSC_FALSE
       do ieqn = 1, cur_cond%num_other_goveqs
          if (cur_cond%itype_of_other_goveqs(ieqn) == GE_THERM_SSW_TBASED) then
             is_bc_sh2o = PETSC_TRUE
          endif
       enddo

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()
          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          sum_conn = sum_conn + 1

          if ((.not.this%aux_vars_in(cell_id_dn)%is_active)) cycle

          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ)
             if (.not. this%aux_vars_bc(sum_conn)%is_active) cycle

             frac          = this%aux_vars_bc(sum_conn)%frac
             area          = cur_conn_set%conn(iconn)%GetArea()
             dist_up       = cur_conn_set%conn(iconn)%GetDistUp()
             dist_dn       = cur_conn_set%conn(iconn)%GetDistDn()
             dist          = dist_up + dist_dn

             therm_cond_up = this%aux_vars_bc(sum_conn)%therm_cond
             therm_cond_dn = this%aux_vars_in(cell_id_dn)%therm_cond

             ! Distance weighted harmonic average
             if (is_bc_sh2o) then
                therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                     (therm_cond_up*dist_dn + therm_cond_dn*dist_up)
                dist = dist_dn + max(1.d-6, dist_up*2.d0)/2.d0
             else
                therm_cond_aveg = therm_cond_up*therm_cond_dn*dist/ &
                     (therm_cond_up*dist_dn + therm_cond_dn*dist_up)
             endif

             T        = this%aux_vars_in(cell_id_dn)%temperature
             heat_cap = this%aux_vars_in(cell_id_dn)%heat_cap_pva
             tfactor  = this%aux_vars_in(cell_id_dn)%tuning_factor
             vol      = this%mesh%vol(cell_id_dn)
#ifdef MATCH_CLM_FORMULATION
             factor =  (dt*tfactor)/(heat_cap*vol)
#else
             factor = 1.d0
#endif

             value = frac*(1.d0 - cnfac)*therm_cond_aveg/dist*area*factor

             call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_dn-1, value, &
                  ADD_VALUES, ierr); CHKERRQ(ierr)
             
          case (COND_HEAT_FLUX)

             dhsdT = this%aux_vars_bc(sum_conn)%dhsdT
             frac  = this%aux_vars_bc(sum_conn)%frac
             area  = cur_conn_set%conn(iconn)%GetArea()

             T        = this%aux_vars_in(cell_id_dn)%temperature
             heat_cap = this%aux_vars_in(cell_id_dn)%heat_cap_pva
             tfactor  = this%aux_vars_in(cell_id_dn)%tuning_factor
             vol      = this%mesh%vol(cell_id_dn)
#ifdef MATCH_CLM_FORMULATION
             factor =  (dt*tfactor)/(heat_cap*vol)
#else
             factor = 1.d0
#endif

             value = -frac*dhsdT**area*factor

             call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_dn-1, value, &
                  ADD_VALUES, ierr); CHKERRQ(ierr)

          case default
             write(iulog,*) 'ThermKSPTempSoilComputeOperatorsDiag: Unknown cond%itype'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilComputeOperatorsDiag

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, list_id_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    Mat                                      :: A
    Mat                                      :: B
    PetscInt                                 :: itype_of_other_goveq
    PetscInt                                 :: list_id_of_other_goveq
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id_dn
    PetscReal                                :: dist
    PetscReal                                :: dist_up
    PetscReal                                :: dist_dn
    PetscReal                                :: area
    PetscReal                                :: therm_cond_aveg
    PetscReal                                :: therm_cond_up
    PetscReal                                :: therm_cond_dn
    PetscReal                                :: frac
    PetscReal                                :: value
    PetscReal                                :: factor
    PetscReal                                :: T
    PetscReal                                :: heat_cap
    PetscReal                                :: tfactor
    PetscReal                                :: vol
    PetscReal                                :: dt
    PetscBool                                :: is_bc_sh2o
    type(connection_set_type), pointer       :: cur_conn_set
    type(condition_type)     , pointer       :: cur_cond

    dt = this%dtime

    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       is_bc_sh2o = PETSC_FALSE
       do ieqn = 1, cur_cond%num_other_goveqs
          if (cur_cond%itype_of_other_goveqs(ieqn) == GE_THERM_SSW_TBASED) then
             is_bc_sh2o = PETSC_TRUE
          endif
       enddo
       cur_conn_set => cur_cond%conn_set

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          do ieqn = 1, cur_cond%num_other_goveqs
             if (cur_cond%list_id_of_other_goveqs(ieqn) == list_id_of_other_goveq) then

                do iconn = 1, cur_conn_set%num_connections

                   cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()
                   cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()

                   sum_conn = sum_conn + 1

                   if ((.not.this%aux_vars_in(cell_id_dn)%is_active)) cycle

                   if (.not. this%aux_vars_bc(sum_conn)%is_active) cycle

                   frac          = this%aux_vars_bc(sum_conn)%frac
                   area          = cur_conn_set%conn(iconn)%GetArea()
                   dist_up       = cur_conn_set%conn(iconn)%GetDistUp()
                   dist_dn       = cur_conn_set%conn(iconn)%GetDistDn()
                   dist          = dist_up + dist_dn

                   therm_cond_up = this%aux_vars_bc(sum_conn)%therm_cond
                   therm_cond_dn = this%aux_vars_in(cell_id_dn)%therm_cond

                   ! Distance weighted harmonic average
                   if (is_bc_sh2o) then
                      therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                           (therm_cond_up*dist_dn + therm_cond_dn*dist_up)
                      dist = dist_dn + max(1.d-6, dist_up*2.d0)/2.d0
                   else
                      therm_cond_aveg = therm_cond_up*therm_cond_dn*dist/ &
                           (therm_cond_up*dist_dn + therm_cond_dn*dist_up)
                   endif

                   T        = this%aux_vars_in(cell_id_dn)%temperature
                   heat_cap = this%aux_vars_in(cell_id_dn)%heat_cap_pva
                   tfactor  = this%aux_vars_in(cell_id_dn)%tuning_factor
                   vol      = this%mesh%vol(cell_id_dn)
#ifdef MATCH_CLM_FORMULATION
                   factor =  (dt*tfactor)/(heat_cap*vol)
#else
                   factor = 1.d0
#endif

                   value = frac*(1.d0 - cnfac)*therm_cond_aveg/dist*area*factor

                   call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -value, &
                        ADD_VALUES, ierr); CHKERRQ(ierr)

                enddo
             else
                sum_conn = sum_conn + cur_conn_set%num_connections
             endif
          enddo
       else
          sum_conn = sum_conn + cur_conn_set%num_connections
       endif

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilComputeOperatorsOffDiag


  !------------------------------------------------------------------------

  subroutine ThermKSPTempSoilUpdateBoundaryConn(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ConnectionSetType, only          : connection_set_type
    use ConditionType, only              : condition_type
    use mpp_varpar, only                 : nlevsno
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SNOW_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_soil_type) :: this
    !
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: cell_id
    PetscInt                                 :: sum_conn
    PetscBool                                :: is_bc_snow
    PetscBool                                :: is_bc_sh2o
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set
      
    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       is_bc_snow = PETSC_FALSE
       is_bc_sh2o = PETSC_FALSE

       do ieqn = 1, cur_cond%num_other_goveqs
       if (cur_cond%itype_of_other_goveqs(ieqn) == GE_THERM_SSW_TBASED) then
          is_bc_sh2o = PETSC_TRUE
       endif
       enddo
       
       do ieqn = 1, cur_cond%num_other_goveqs
       if (cur_cond%itype_of_other_goveqs(ieqn) == GE_THERM_SNOW_TBASED) then
          is_bc_snow = PETSC_TRUE
       endif
       enddo

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
          sum_conn = sum_conn + 1

          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ)
             if (.not. this%aux_vars_bc(sum_conn)%is_active) cycle
             
             if (is_bc_snow) then
                call cur_conn_set%conn(iconn)%SetDistUp(this%aux_vars_bc(sum_conn)%dist_up)
             elseif (is_bc_sh2o) then
                call cur_conn_set%conn(iconn)%SetDistUp(this%aux_vars_bc(sum_conn)%dz/2.d0)
             endif

          case default

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSoilUpdateBoundaryConn
#endif
  
end module GoveqnThermalKSPTemperatureSoilType
