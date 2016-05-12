module GoveqnThermalKSPTemperatureSSWType

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl                       , only : iulog
  use mpp_abortutils                       , only : endrun
  use mpp_shr_log_mod                      , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType        , only : goveqn_base_type
  use ThermalKSPTemperatureSSWAuxType  , only : therm_ksp_temp_ssw_auxvar_type
  use SystemOfEquationsThermalAuxType  , only : sysofeqns_thermal_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscviewer.h"

  type, public, extends(goveqn_base_type) :: goveqn_thermal_ksp_temp_ssw_type
     type (therm_ksp_temp_ssw_auxvar_type), pointer :: aux_vars_in(:)  ! Internal state.
     type (therm_ksp_temp_ssw_auxvar_type), pointer :: aux_vars_bc(:)  ! Boundary conditions.
     type (therm_ksp_temp_ssw_auxvar_type), pointer :: aux_vars_ss(:)  ! Source-sink.

     PetscInt, pointer                             :: soe_auxvars_bc_offset (:) ! SoE auxvar offset corresponding to BCs
     PetscInt, pointer                             :: soe_auxvars_ss_offset (:) ! SoE auxvar offset corresponding to SSs

  contains
     procedure, public :: Setup                     => ThermKSPTempSSWSetup
     procedure, public :: AllocateAuxVars           => ThermKSPTempSSWAllocateAuxVars

     procedure, public :: GetFromSOEAuxVarsIntrn    => ThermKSPTempSSWGetFromSOEAuxVarsIntrn
     procedure, public :: GetFromSOEAuxVarsBC       => ThermKSPTempSSWGetFromSOEAuxVarsBC
     procedure, public :: GetFromSOEAuxVarsSS       => ThermKSPTempSSWGetFromSOEAuxVarsSS
     procedure, public :: GetDataFromSOEAuxVar      => ThermKSPTempSSWGetDataFromSOEAuxVar
     procedure, public :: GetConditionNames         => ThermKSPTempSSWGetConditionNames
     procedure, public :: GetNumConditions          => ThermKSPTempSSWGetNumConditions
     procedure, public :: GetNumCellsInConditions   => ThermKSPTempSSWGetNumCellsInConditions

     procedure, public :: SetDataInSOEAuxVar        => ThermKSPTempSSWSetDataInSOEAuxVar
     procedure, public :: SetSOEAuxVarOffsets       => ThermKSPTempSSWSetSOEAuxVarOffsets
     procedure, public :: UpdateInternalConn        => ThermKSPTempSSWUpdateInternalConn
     procedure, public :: UpdateBoundaryConn        => ThermKSPTempSSWUpdateBoundaryConn
     procedure, public :: UpdateAuxVarsIntrn        => ThermKSPTempSSWUpdateAuxVarsIntrn

     procedure, public :: ComputeRHS                => ThermKSPTempSSWComputeRHS
     procedure, public :: ComputeOperatorsDiag      => ThermKSPTempSSWComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag   => ThermKSPTempSSWComputeOperatorsOffDiag

  end type goveqn_thermal_ksp_temp_ssw_type

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation for Thermal equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants, only : MESH_CLM_SSW_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this

    call this%Init()

    this%name         = "Standing surface water thermal equation based on temperature"
    this%id           = GE_THERM_SSW_TBASED
    this%mesh_itype   = MESH_CLM_SSW_COL

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_ss)

    nullify(this%soe_auxvars_bc_offset)
    nullify(this%soe_auxvars_ss_offset)

  end subroutine ThermKSPTempSSWSetup

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWAllocateAuxVars(this)
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
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
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

  end subroutine ThermKSPTempSSWAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetDataFromSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    !
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : AUXVAR_BC
    use MultiPhysicsProbConstants       , only : AUXVAR_SS
    use ConditionType                   , only : condition_type
    use ConnectionSetType               , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type), intent(inout)        :: this
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
       call ThermKSPTempSSWGetFromSOEAuxVarsIntrn(this, soe_avars, iauxvar_off)
    case (AUXVAR_BC)
       call ThermKSPTempSSWGetFromSOEAuxVarsBC(this, soe_avars)
    case (AUXVAR_SS)
       call ThermKSPTempSSWGetFromSOEAuxVarsSS(this, soe_avars)
    case default
       write(iulog,*) 'ERROR: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSSWGetDataFromSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetFromSOEAuxVarsIntrn(this, soe_avars, offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for internal nodes
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) , intent(inout)            :: this
    type(sysofeqns_thermal_auxvar_type)     , dimension(:), intent(in) :: soe_avars
    PetscInt                                , intent(in)               :: offset
    !
    ! LOCAL VARIABLES
    PetscInt :: iauxvar
    PetscInt :: nauxvar

    nauxvar = size(this%aux_vars_in)
    if( nauxvar > size(soe_avars) ) then
       write(iulog,*) 'ERROR: size(ge_avars) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, nauxvar
       this%aux_vars_in(iauxvar)%temperature = soe_avars(iauxvar+offset)%temperature
       this%aux_vars_in(iauxvar)%is_active   = soe_avars(iauxvar+offset)%is_active
       this%aux_vars_in(iauxvar)%dz          = soe_avars(iauxvar+offset)%dz
       this%aux_vars_in(iauxvar)%frac        = soe_avars(iauxvar+offset)%frac
    enddo

  end subroutine ThermKSPTempSSWGetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWGetFromSOEAuxVarsBC(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType                   , only : condition_type
    use ConnectionSetType               , only : connection_set_type
    use MultiPhysicsProbConstants       , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants       , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants       , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type), intent(inout)        :: this
    type(sysofeqns_thermal_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                                                       :: iauxvar, iauxvar_off, iconn
    integer                                                       :: auxVarCt_ge, auxVarCt_soe
    integer                                                       :: condition_id, sum_conn
    integer                                                       :: cell_id
    type(condition_type), pointer                                 :: cur_cond
    type(connection_set_type), pointer                            :: cur_conn_set
    character(len=256)                                            :: string

    auxVarCt_ge = size(this%aux_vars_bc)

    auxVarCt_soe = size(soe_avars)
    if( auxVarCt_ge > auxVarCt_soe ) then
       write(iulog,*) 'ERROR: size(this%aux_vars_bc) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          sum_conn = sum_conn + cur_cond%conn_set%num_connections
          cur_cond => cur_cond%next
          cycle
       endif

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = this%soe_auxvars_bc_offset(condition_id)

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1
          cell_id = cur_conn_set%id_dn(iconn)

          select case(cur_cond%itype)
          case (COND_HEAT_FLUX)
             ! H
             this%aux_vars_bc(sum_conn)%condition_value =  &
                  soe_avars(iconn + iauxvar_off)%condition_value

             ! H - dH/dT * T
             cur_cond%value(iconn) = &
                  soe_avars(iconn + iauxvar_off)%condition_value - &
                  soe_avars(iconn + iauxvar_off)%dhsdT * &
                  this%aux_vars_in(cell_id)%temperature

             this%aux_vars_bc(sum_conn)%dhsdT = &
                  soe_avars(iconn + iauxvar_off)%dhsdT

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing
          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'ERROR: Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWGetFromSOEAuxVarsBC


  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetFromSOEAuxVarsSS(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType                   , only : condition_type
    use ConnectionSetType               , only : connection_set_type
    use MultiPhysicsProbConstants       , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants       , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type), intent(inout)       :: this
    type(sysofeqns_thermal_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                                                       :: iauxvar
    integer                                                       :: iauxvar_off
    integer                                                       :: iconn
    integer                                                       :: auxVarCt_ge, auxVarCt_soe
    integer                                                       :: condition_id, sum_conn
    PetscReal                                                     :: var_value
    type(condition_type), pointer                                 :: cur_cond
    type(connection_set_type), pointer                            :: cur_conn_set
    character(len=256)                                            :: string

    auxVarCt_ge = size(this%aux_vars_ss)

    auxVarCt_soe = size(soe_avars)
    if( auxVarCt_ge > auxVarCt_soe ) then
       write(iulog,*) 'ERROR: size(this%aux_vars_ss) > size(soe_avars)'
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

       if (trim(cur_cond%name) /= 'Lateral_flux') then
          !
          ! Do not update values associated with lateral flux because:
          !  - SOE auxvars associated with lateral flux source-sink
          !    have zero values, AND
          !  - GE auxvars already have pre-computed values of lateral flux.
          !
          cur_conn_set => cur_cond%conn_set
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             select case(cur_cond%itype)
             case (COND_HEAT_RATE)
                var_value = soe_avars(iconn + iauxvar_off)%condition_value
                this%aux_vars_ss(sum_conn)%condition_value = var_value
                cur_cond%value(iconn) = var_value
             case default
                write(string,*) cur_cond%itype
                write(iulog,*) 'ERROR: Unknown cur_cond%itype = ' // trim(string)
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select
          enddo
       endif

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWGetFromSOEAuxVarsSS

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetConditionNames(this, cond_type, &
                cond_type_to_exclude, num_conds, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the total number and names of conditions (eg. boundary condition
    ! or source-sink) present.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    PetscInt, intent(in)                     :: cond_type
    PetscInt, intent(in)                     :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
    character (len=256), pointer             :: cond_names(:)
    !
    type(condition_type),pointer             :: cur_cond
    character(len=256)                       :: string

    ! Find number of BCs
    call this%GetNumConditions(cond_type, COND_NULL, num_conds)

    if (num_conds == 0) then
       nullify(cond_names)
       return
    endif

    allocate(cond_names(num_conds))

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       cur_cond => this%boundary_conditions%first
    case (COND_SS)
      cur_cond => this%source_sinks%first
    case default
       write(string,*) cond_type
       write(iulog,*) 'ERROR: Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    num_conds = 0
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_type_to_exclude) then
          num_conds = num_conds + 1
          cond_names(num_conds) = cur_cond%name
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWGetConditionNames

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetNumConditions(this, cond_type, &
              cond_type_to_exclude, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    PetscInt                                 :: cond_type
    PetscInt                                 :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
    !
    type(condition_type),pointer             :: cur_cond
    character(len=256)                       :: string

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       cur_cond => this%boundary_conditions%first
    case (COND_SS)
      cur_cond => this%source_sinks%first
    case default
       write(string,*) cond_type
       write(iulog,*) 'ERROR: Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    num_conds = 0
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_type_to_exclude) then
          num_conds = num_conds + 1
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWGetNumConditions

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWGetNumCellsInConditions(this, cond_type, &
                cond_type_to_exclude, num_conds, ncells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions (eg. boundary condition or
    ! source-sink) and number of control volumes associated with each condition
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    PetscInt, intent(in)                     :: cond_type
    PetscInt, intent(in)                     :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
    PetscInt, intent(out), pointer           :: ncells_for_conds(:)
    !
    type(condition_type),pointer             :: cur_cond
    character(len=256)                       :: string

    ! Find number of BCs
    call this%GetNumConditions(cond_type, COND_NULL, num_conds)

    if (num_conds == 0) then
       nullify(ncells_for_conds)
       return
    endif

    allocate(ncells_for_conds(num_conds))

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       cur_cond => this%boundary_conditions%first
    case (COND_SS)
      cur_cond => this%source_sinks%first
    case default
       write(string,*) cond_type
       write(iulog,*) 'ERROR: Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    num_conds = 0
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_type_to_exclude) then
          num_conds = num_conds + 1
          ncells_for_conds(num_conds) = cur_cond%ncells
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWGetNumCellsInConditions

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWSetDataInSOEAuxVar(this, soe_avar_type, soe_avars, &
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
    class(goveqn_thermal_ksp_temp_ssw_type)                          :: this
    PetscInt                                                         :: soe_avar_type
    type (sysofeqns_thermal_auxvar_type) , dimension(:), intent(out) :: soe_avars
    PetscInt                             , optional                  :: offset
    !
    ! !LOCAL VARIABLES
    type (therm_ksp_temp_ssw_auxvar_type), pointer :: ge_avars(:)
    PetscInt                                       :: iauxvar
    PetscInt                                       :: iauxvar_off

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       ge_avars => this%aux_vars_in

       if ( size(ge_avars) > size(soe_avars)) then
          write(iulog,*) 'ERROR: size(ge_avars) > size(soe_avars)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       do iauxvar = 1, size(ge_avars)
          if (this%mesh%is_active(iauxvar)) then
             soe_avars(iauxvar+iauxvar_off)%temperature =  &
                  ge_avars(iauxvar)%temperature
          endif
       enddo

    case default
       write(iulog,*) 'ERROR: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSSWSetDataInSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWSetSOEAuxVarOffsets(this, bc_offset_count, bc_offsets, &
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
    class(goveqn_thermal_ksp_temp_ssw_type)  :: this
    PetscInt                                 :: bc_offset_count
    PetscInt, pointer                        :: bc_offsets(:)
    PetscInt                                 :: ss_offset_count
    PetscInt, pointer                        :: ss_offsets(:)
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: cond_count

    call this%GetNumConditions(COND_BC, -1, cond_count)
    if (bc_offset_count > cond_count) then
       write(iulog,*) 'ERROR: bc_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call this%GetNumConditions(COND_SS, -1, cond_count)
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
       this%soe_auxvars_ss_offset(cond_count) = bc_offsets(cond_count)
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWSetSOEAuxVarOffsets

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWUpdateInternalConn(this)

    !
    ! !DESCRIPTION:
    !
    use ConnectionSetType, only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type)        :: this
    !
    ! !LOCAL VARIABLES
    type (therm_ksp_temp_ssw_auxvar_type), pointer :: aux_vars_in(:)
    PetscInt                                       :: icell

    aux_vars_in   => this%aux_vars_in

    do icell = 1, this%mesh%ncells_all
       if (aux_vars_in(icell)%is_active) then
          this%mesh%dz(icell)  = max(1.0d-6, aux_vars_in(icell)%dz)
          this%mesh%vol(icell) = this%mesh%dx(icell)* &
                                 this%mesh%dy(icell)* &
                                 this%mesh%dz(icell)          
       else
          this%mesh%vol(icell) = 0.d0
       endif
    enddo

  end subroutine ThermKSPTempSSWUpdateInternalConn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWUpdateBoundaryConn(this)
    !
    ! !DESCRIPTION:
    !
    use ConnectionSetType, only          : connection_set_type
    use ConditionType, only              : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type)        :: this
    !
    ! !LOCAL VARIABLES
    type(connection_set_type), pointer             :: conn_set
    PetscInt                                       :: cell_id
    PetscInt                                       :: iconn
    type(condition_type),pointer                   :: cur_cond
    type(condition_type),pointer                   :: soil_temp_cond
    PetscBool                                      :: soil_temp_cond_found

    ! Boundary cells
    soil_temp_cond_found = PETSC_FALSE

    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       if (trim(cur_cond%name) == 'BC_from_soil_governing_equation') then
          soil_temp_cond => cur_cond
          soil_temp_cond_found = PETSC_TRUE
       endif

       cur_cond => cur_cond%next
    enddo

    if (.not.soil_temp_cond_found) then
       write(iulog,*)'ERROR: BC with soil not found.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Update 'dist' on which the top heat flux BC will be applied
    conn_set => soil_temp_cond%conn_set
    do iconn = 1, conn_set%num_connections

       cell_id = conn_set%id_dn(iconn)

       if (this%aux_vars_in(cell_id)%is_active ) then
          conn_set%dist_dn(iconn) = this%mesh%dz(cell_id)/2.d0
       endif       
    enddo

  end subroutine ThermKSPTempSSWUpdateBoundaryConn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWUpdateAuxVarsIntrn(this)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt :: ghosted_id

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_local
       call this%aux_vars_in(ghosted_id)%AuxVarCompute( &
            this%mesh%dz(ghosted_id), this%mesh%vol(ghosted_id))
    enddo

  end subroutine ThermKSPTempSSWUpdateAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    Vec                                     :: B
    PetscErrorCode                          :: ierr
    !
    ! !LOCAL VARIABLES
    PetscReal, dimension(:), pointer         :: b_p

    call VecGetArrayF90(B, b_p, ierr); CHKERRQ(ierr);

    call ThermalKSPTempSSWAccum(this, b_p)
    call ThermalKSPTempSSWDivergence(this, b_p)
    
    call VecRestoreArrayF90(B, b_p, ierr); CHKERRQ(ierr)

  end subroutine ThermKSPTempSSWComputeRHS

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    Mat                                     :: A
    Mat                                     :: B
    PetscErrorCode                          :: ierr
    !
    !
    ! !LOCAL VARIABLES
    PetscInt                                   :: iconn, sum_conn
    PetscInt                                   :: cell_id
    PetscReal                                  :: dist, dist_up, dist_dn
    PetscReal                                  :: area, vol
    PetscReal                                  :: therm_cond_aveg, therm_cond_up, therm_cond_dn
    PetscReal                                  :: heat_cap
    PetscReal                                  :: dhsdT
    PetscReal                                  :: dt
    PetscReal                                  :: coeff
    PetscReal                                  :: value
    type(connection_set_type), pointer         :: cur_conn_set
    type(condition_type),pointer               :: cur_cond

    dt = this%dtime

    ! Diagonal term
    do cell_id = 1, this%mesh%ncells_local

       heat_cap   = this%aux_vars_in(cell_id)%heat_cap_pva
       vol        = this%mesh%vol(cell_id)

       if (this%aux_vars_in(cell_id)%is_active) then
          value = heat_cap*vol/dt
          !value = 1.d0
       else
          value = 1.d0
       endif

       call MatSetValuesLocal(B, 1, cell_id-1, 1, cell_id-1, value, &
            ADD_VALUES, ierr); CHKERRQ(ierr)
    enddo
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id = cur_conn_set%id_dn(iconn)
          sum_conn = sum_conn + 1

          if ((.not.this%aux_vars_in(cell_id)%is_active)) cycle

          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ)
             !if (.not. this%aux_vars_bc(sum_conn)%is_active) cycle

             area          = cur_conn_set%area(iconn)
             dist_up       = cur_conn_set%dist_up(iconn)
             dist_dn       = cur_conn_set%dist_dn(iconn)
             dist          = dist_up + dist_dn

             therm_cond_up = this%aux_vars_bc(sum_conn)%therm_cond
             therm_cond_dn = this%aux_vars_in(cell_id )%therm_cond

             ! Distance weighted harmonic average
             dist_dn = this%aux_vars_in(cell_id)%dz/2.d0
             therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                  (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

             coeff = (1.d0 - cnfac)*therm_cond_aveg/dist*area

             call MatSetValuesLocal(B, 1, cell_id-1, 1, cell_id-1, coeff, &
                  ADD_VALUES, ierr); CHKERRQ(ierr)
             
          case (COND_HEAT_FLUX)

             dhsdT = this%aux_vars_bc(sum_conn)%dhsdT
             area  = cur_conn_set%area(iconn)

             value = -dhsdT*area

             value = value

             call MatSetValuesLocal(B, 1, cell_id-1, 1, cell_id-1, value, &
                  ADD_VALUES, ierr); CHKERRQ(ierr)

          case default
             write(iulog,*) 'ERROR: Unknown cond%itype'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWComputeOperatorsDiag

  !------------------------------------------------------------------------

  subroutine ThermKSPTempSSWComputeOperatorsOffDiag(this, A, B, &
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: this
    Mat                                     :: A
    Mat                                     :: B
    PetscInt                                :: itype_of_other_goveq
    PetscInt                                :: list_id_of_other_goveq
    PetscErrorCode                          :: ierr

    !
    ! !LOCAL VARIABLES
    PetscInt                                   :: iconn, sum_conn
    PetscInt                                   :: cell_id_up, cell_id_dn
    PetscReal                                  :: dist, dist_up, dist_dn
    PetscReal                                  :: area
    PetscReal                                  :: therm_cond_aveg, therm_cond_up, therm_cond_dn
    PetscReal                                  :: coeff
    type(connection_set_type), pointer         :: cur_conn_set
    type(condition_type),pointer               :: cur_cond

    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ .and. &
           cur_cond%list_id_of_other_goveq == list_id_of_other_goveq) then

          do iconn = 1, cur_conn_set%num_connections

             cell_id_dn = cur_conn_set%id_dn(iconn)
             cell_id_up = cur_conn_set%id_up(iconn)
             sum_conn   = sum_conn + 1

             if ((.not.this%aux_vars_in(cell_id_dn)%is_active)) cycle

             area          = cur_conn_set%area(iconn)
             dist_up       = cur_conn_set%dist_up(iconn)
             dist_dn       = cur_conn_set%dist_dn(iconn)
             dist          = dist_up + dist_dn

             therm_cond_up = this%aux_vars_bc(sum_conn)%therm_cond
             therm_cond_dn = this%aux_vars_in(cell_id_dn)%therm_cond

             ! Distance weighted harmonic average
             dist_dn = this%aux_vars_in(cell_id_dn)%dz/2.d0
             therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                  (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

             coeff = (1.d0 - cnfac)*therm_cond_aveg/dist*area

             call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -coeff, &
                  ADD_VALUES, ierr); CHKERRQ(ierr)
             
          enddo
       else
          sum_conn = sum_conn + cur_conn_set%num_connections
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermKSPTempSSWComputeOperatorsOffDiag

  !------------------------------------------------------------------------

  subroutine ThermalKSPTempSSWAccum(geq_ssw, b_p)
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
    class(goveqn_thermal_ksp_temp_ssw_type) :: geq_ssw
    PetscReal, dimension(:), intent(out)     :: b_p
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: cell_id
    PetscReal                                :: T
    PetscReal                                :: heat_cap
    PetscReal                                :: vol
    PetscReal                                :: dt

    dt = geq_ssw%dtime

    ! Interior cells
    do cell_id = 1, geq_ssw%mesh%ncells_local

       if (geq_ssw%aux_vars_in(cell_id)%is_active) then

          T        = geq_ssw%aux_vars_in(cell_id)%temperature
          heat_cap = geq_ssw%aux_vars_in(cell_id)%heat_cap_pva
          vol      = geq_ssw%mesh%vol(cell_id)

          b_p(cell_id) = heat_cap*vol/(dt)*T
          !b_p(cell_id) = geq_ssw%aux_vars_in(cell_id)%temperature
       endif
    enddo

  end subroutine ThermalKSPTempSSWAccum

  !------------------------------------------------------------------------

  subroutine ThermalKSPTempSSWDivergence(geq_ssw, b_p)
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_ksp_temp_ssw_type) :: geq_ssw
    PetscReal, dimension(:), intent(out)     :: b_p
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id
    PetscReal                                :: flux
    PetscReal                                :: area
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set
    PetscReal                                :: dist, dist_up, dist_dn
    PetscReal                                :: therm_cond_aveg, therm_cond_up, therm_cond_dn
    PetscReal                                :: T_up, T_dn

    ! Interior cells
    cur_conn_set => geq_ssw%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       write(iulog,*) 'ERROR: No internal connections allowed in standing surface water'
       call endrun(msg=errMsg(__FILE__, __LINE__))

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => geq_ssw%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id  = cur_conn_set%id_dn(iconn)
          sum_conn = sum_conn + 1

          if (.not.geq_ssw%aux_vars_in(cell_id )%is_active) cycle
          
          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ)
             area = cur_conn_set%area(iconn)

             dist_up       = cur_conn_set%dist_up(iconn)
             dist_dn       = cur_conn_set%dist_dn(iconn)
             dist          = dist_up + dist_dn

             therm_cond_up = geq_ssw%aux_vars_bc(sum_conn)%therm_cond
             therm_cond_dn = geq_ssw%aux_vars_in(cell_id )%therm_cond

             T_up = geq_ssw%aux_vars_bc(sum_conn)%temperature
             T_dn = geq_ssw%aux_vars_in(cell_id )%temperature

             dist_dn = geq_ssw%aux_vars_in(cell_id)%dz/2.d0
             therm_cond_aveg = therm_cond_up*therm_cond_dn*(dist_up + dist_dn)/ &
                  (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

             flux = -therm_cond_aveg*(T_up - T_dn)/(dist)

             b_p(cell_id) = b_p(cell_id) - cnfac*flux*area


          case (COND_HEAT_FLUX)             
             area = cur_conn_set%area(iconn)

             b_p(cell_id) = b_p(cell_id) + &
                  cur_cond%value(iconn) * area

          case default
           write(iulog,*)'ERROR: Unknown boundary condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

       enddo
       cur_cond => cur_cond%next
    enddo

    ! Source-sink cells
    cur_cond => geq_ssw%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          cell_id = cur_conn_set%id_dn(iconn)

          if ((.not.geq_ssw%aux_vars_in(cell_id)%is_active)) cycle


          select case(cur_cond%itype)
          case(COND_HEAT_RATE)
             b_p(cell_id) = b_p(cell_id) + cur_cond%value(iconn)

          case default
           write(iulog,*)'ERROR: Unknown source-sink condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

     enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermalKSPTempSSWDivergence

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

#endif
  
end module GoveqnThermalKSPTemperatureSSWType
