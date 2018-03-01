
module SystemOfEquationsTHType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                              , only : iulog
  use mpp_abortutils                          , only : endrun
  use mpp_shr_log_mod                         , only : errMsg => shr_log_errMsg
  use SystemOfEquationsThermalEnthalpyAuxType , only : sysofeqns_thermal_enthalpy_auxvar_type
  use SystemOfEquationsBaseType               , only : sysofeqns_base_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(sysofeqns_base_type) :: sysofeqns_th_type

     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_in(:)            !!< Internal state.
     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_bc(:)            !!< Boundary conditions.
     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_ss(:)            !!< Source-sink.

     PetscInt                                      :: num_calls_to_ifunction
     PetscInt                                      :: num_calls_to_ijacobian

     PetscInt, pointer                             :: soe_auxvars_in_offset (:) ! Cummulative sum of number of control volumes associated with internal condition.
     PetscInt, pointer                             :: soe_auxvars_bc_offset (:) ! Cummulative sum of number of control volumes associated with each boundary condition.
     PetscInt, pointer                             :: soe_auxvars_ss_offset (:) ! Cummulative sum of number of control volumes associated with each source-sink condition.

     PetscInt, pointer                             :: soe_auxvars_in_ncells (:) ! Number of control volumes associated with each internal condition.
     PetscInt, pointer                             :: soe_auxvars_bc_ncells (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer                             :: soe_auxvars_ss_ncells (:) ! Number of control volumes associated with each source-sink condition.
     PetscInt, pointer                             :: soe_auxvars_bc_ncells_per_goveqn (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer                             :: soe_auxvars_ss_ncells_per_goveqn (:) ! Number of control volumes associated with each source-sink condition.

     PetscInt                                      :: num_auxvars_in            ! Number of auxvars associated with internal state.
     PetscInt                                      :: num_auxvars_in_local      ! Number of auxvars associated with internal state.
     PetscInt                                      :: num_auxvars_bc            ! Number of auxvars associated with boundary condition.
     PetscInt                                      :: num_auxvars_ss            ! Number of auxvars associated with source-sink condition.

   contains
     procedure, public :: Init                   => SOETHInit
     procedure, public :: AddGovEqn              => SOETHAddGovEqn

     procedure, public :: CreateVectorsForGovEqn => SOETHCreateVectorsForGovEqn

     procedure, public :: PreSolve               => SOETHPreSolve

     procedure, public :: Residual               => SOETHResidual
     procedure, public :: Jacobian               => SOETHJacobian

     procedure, public :: SetDataFromCLM         => SOETHSetDataFromCLM
     procedure, public :: GetDataForCLM          => SOETHGetDataForCLM

  end type sysofeqns_th_type

  public :: SOETHSetAuxVars, &
            SOETHUpdateConnections

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine SOETHInit(this)
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
    class(sysofeqns_th_type) :: this

    call SOEBaseInit(this)

    this%num_auxvars_in         = 0
    this%num_auxvars_in_local   = 0
    this%num_auxvars_bc         = 0
    this%num_auxvars_ss         = 0

    nullify(this%aux_vars_in           )
    nullify(this%aux_vars_bc           )
    nullify(this%aux_vars_ss           )

    nullify(this%soe_auxvars_in_offset )
    nullify(this%soe_auxvars_bc_offset )
    nullify(this%soe_auxvars_ss_offset )
    nullify(this%soe_auxvars_in_ncells )
    nullify(this%soe_auxvars_bc_ncells )
    nullify(this%soe_auxvars_ss_ncells )
    nullify(this%soe_auxvars_bc_ncells_per_goveqn)
    nullify(this%soe_auxvars_ss_ncells_per_goveqn)

  end subroutine SOETHInit

  !------------------------------------------------------------------------
  subroutine SOETHAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : SOEBaseInit
    use MultiPhysicsProbConstants     , only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants     , only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants     , only : GE_RE
    use MultiPhysicsProbConstants     , only : MESH_CLM_SOIL_COL
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)                            :: this
    PetscInt                                            :: geq_type
    character(len=*)                                    :: name
    PetscInt                                            :: mesh_itype
    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_enthalpy_soil_type) , pointer :: goveq_soil
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    class(goveqn_base_type)                   , pointer :: cur_goveqn
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
       case (GE_THERM_SOIL_EBASED)

          allocate(goveq_soil)
          call goveq_soil%Setup()

          goveq_soil%name              = trim(name)
          goveq_soil%rank_in_soe_list  = this%ngoveqns
          goveq_soil%mesh_itype        = mesh_itype

          if (this%ngoveqns == 1) then
             this%goveqns => goveq_soil
          else
             cur_goveqn%next => goveq_soil
          endif

       case (GE_RE)

          allocate(goveq_richards)
          call goveq_richards%Setup()

          goveq_richards%name              = trim(name)
          goveq_richards%rank_in_soe_list  = this%ngoveqns
          goveq_richards%mesh_itype        = mesh_itype

          if (this%ngoveqns == 1) then
             this%goveqns => goveq_richards
          else
             cur_goveqn%next => goveq_richards
          endif

       case default
          write(iulog,*) 'Unknown governing equation type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

     end subroutine SOETHAddGovEqn

  !------------------------------------------------------------------------
  subroutine SOETHCreateVectorsForGovEqn(this)
    !
    ! !DESCRIPTION:
    ! Creates vectors required by each governing equation
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type) :: this
    !
    class(goveqn_base_type), pointer       :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%CreateVectors()

       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%CreateVectors()

       class default
          write(iulog,*) 'Unsupported cur_goveq type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine SOETHCreateVectorsForGovEqn

  !------------------------------------------------------------------------
  subroutine SOETHPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_TH
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use MultiPhysicsProbConstants     , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)          :: this
    class(goveqn_base_type) , pointer :: cur_goveq
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    PetscInt                          :: offset_in
    PetscInt                          :: offset_bc
    PetscInt                          :: offset_ss
    PetscInt                          :: ndata
    PetscInt                          :: var_type
    PetscInt                          :: row
    PetscInt                          :: col
    PetscInt                          :: igoveqn
    PetscReal               , pointer :: data(:)

    igoveqn   = 0
    offset_in = 0
    offset_bc = 0
    offset_ss = 0

    select case (this%itype)
    case(SOE_TH)

       ! 1) {soln_prev}  ---> sim_aux()
       call SOETHUpdateAuxVars(this, this%solver%soln_prev)

       ! 2) GE ---> GetFromSimAux()
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit

          igoveqn = igoveqn + 1
          select type(cur_goveq)

          class is (goveqn_richards_ode_pressure_type)

             ! Internal auxvars
             ndata = cur_goveq%mesh%ncells_local
             allocate(data(ndata))

             var_type = VAR_PRESSURE
             call SOETHGetAuxVars(this, AUXVAR_INTERNAL, var_type, &
                  offset_in, ndata, data)
             call cur_goveq%SetFromSOEAuxVarsIntrn(var_type, ndata, data)

             offset_in = offset_in + ndata
             deallocate(data)

             ! Boundary auxvars
             ndata = this%soe_auxvars_bc_ncells_per_goveqn(igoveqn)
             if (ndata > 0) then
                write(*,*)'SOETHPreSolve: Add code to support richards BC'
                stop
             endif

             ! Source/sink auxvars
             ndata = this%soe_auxvars_ss_ncells_per_goveqn(igoveqn)
             if (ndata > 0) then
                write(*,*)'SOETHPreSolve: Add code to support richards SS'
                stop
             endif

          class is (goveqn_thermal_enthalpy_soil_type)

             ! Internal auxvars
             ndata = cur_goveq%mesh%ncells_local
             allocate(data(ndata))

             var_type = VAR_TEMPERATURE
             call SOETHGetAuxVars(this, AUXVAR_INTERNAL, var_type, &
                  offset_in, ndata, data)
             call cur_goveq%SetFromSOEAuxVarsIntrn(var_type, ndata, data)

             offset_in = offset_in + ndata
             deallocate(data)

             ! Boundary auxvars
             ndata = this%soe_auxvars_bc_ncells_per_goveqn(igoveqn)
             if (ndata > 0) then
                allocate(data(ndata))

                var_type = VAR_BC_SS_CONDITION
                call SOETHGetAuxVars(this, AUXVAR_BC, var_type, &
                     offset_bc, ndata, data)
                call cur_goveq%SetFromSOEAuxVarsBC(var_type, ndata, data)

                offset_bc = offset_bc + ndata
                deallocate(data)

             endif

             ! Source/sink auxvars
             ndata = this%soe_auxvars_ss_ncells_per_goveqn(igoveqn)
             if (ndata > 0) then
                write(*,*)'SOETHPreSolve: Add code to support energy SS'
                stop
             endif

          class default
             write(iulog,*) 'Unknown goveqn_type'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select

          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'SOEThermalEnthalpyPreSolve: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! 3) GE_1 <---> GE_2 exchange AuxVars()
    do row = 1, this%ngoveqns
       do col = row+1, this%ngoveqns
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call SOETHGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
          call SOETHGovEqnExchangeAuxVars(cur_goveq_2, cur_goveq_1)
       enddo
    enddo

    ! 4) UpdateAuxVars
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%UpdateAuxVars()

       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%UpdateAuxVars()
       end select
          
       call cur_goveq%PreSolve()

       cur_goveq => cur_goveq%next
    enddo

  end subroutine SOETHPreSolve

  !------------------------------------------------------------------------
  subroutine SOETHUpdateAuxVars(therm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE based on the input
    ! vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)          :: therm_soe
    Vec                               :: X
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: dm_id
    PetscInt                          :: var_type
    PetscInt                          :: nDM
    DM, pointer                       :: dms(:)
    Vec, pointer                      :: X_subvecs(:)
    PetscInt                          :: size
    PetscInt                          :: offset
    PetscErrorCode                    :: ierr
    class(goveqn_base_type) , pointer :: cur_goveq

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(therm_soe%solver%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(therm_soe%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(therm_soe%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Update the SoE auxvars
    offset = 0
    dm_id = 0
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)

       class is (goveqn_thermal_enthalpy_soil_type)
          var_type = VAR_TEMPERATURE

       class is (goveqn_richards_ode_pressure_type)
          var_type = VAR_PRESSURE

       class default
          write(iulog,*) 'Unknown goveqn_type'
          call endrun(msg=errMsg(__FILE__, __LINE__))

       end select

       dm_id = dm_id + 1
       call SOETHSetAuxVars(therm_soe, AUXVAR_INTERNAL, var_type, &
            X_subvecs(dm_id), offset)
       call VecGetSize(X_subvecs(dm_id), size, ierr); CHKERRQ(ierr)
       offset = offset + size

       cur_goveq => cur_goveq%next
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(therm_soe%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)


  end subroutine SOETHUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine SOETHSetAuxVars(therm_soe, auxvar_type, var_type, &
       var_vec, offset)
    !
    ! !DESCRIPTION:
    ! Set values in SoE auxvars.
    !
    ! !USES:
    use MultiPhysicsProbConstants              , only : AUXVAR_INTERNAL
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxSetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type) :: therm_soe
    PetscInt                 :: auxvar_type
    PetscInt, intent(in)     :: var_type
    Vec                      :: var_vec
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer       :: var_p(:)
    PetscInt                 :: nauxvar
    PetscInt                 :: nvar
    PetscInt, optional       :: offset
    PetscInt                 :: iauxvar
    PetscInt                 :: iauxvar_off
    PetscErrorCode           :: ierr

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
          write(iulog,*) 'nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call VecGetArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)

       call SOEThermalEnthalpyAuxSetRData(therm_soe%aux_vars_in, var_type, &
            nvar, iauxvar_off, var_p)

       call VecRestoreArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)

    case default
       write(iulog,*) 'auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select


  end subroutine SOETHSetAuxVars

  !------------------------------------------------------------------------
  subroutine SOETHGetAuxVars(therm_soe, auxvar_type, var_type, &
       offset, nvar, data_1d)
    !
    ! !DESCRIPTION:
    ! Set values in SoE auxvars.
    !
    ! !USES:
    use MultiPhysicsProbConstants              , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants              , only : AUXVAR_BC
    use MultiPhysicsProbConstants              , only : AUXVAR_SS
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxGetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type) :: therm_soe
    PetscInt                 :: auxvar_type
    PetscInt, intent(in)     :: var_type
    PetscInt                 :: nvar
    PetscReal, pointer       :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer       :: var_p(:)
    PetscInt                 :: nauxvar
    PetscInt, optional       :: offset
    PetscInt                 :: iauxvar
    PetscInt                 :: iauxvar_off
    PetscErrorCode           :: ierr

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    select case(auxvar_type)
    case (AUXVAR_INTERNAL)

       nauxvar = size(therm_soe%aux_vars_in)

       if (nvar+iauxvar_off > nauxvar) then
          write(iulog,*) 'nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call SOEThermalEnthalpyAuxGetRData(therm_soe%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case (AUXVAR_BC)

       nauxvar = size(therm_soe%aux_vars_bc)

       if (nvar+iauxvar_off > nauxvar) then
          write(iulog,*) 'nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call SOEThermalEnthalpyAuxGetRData(therm_soe%aux_vars_bc, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case (AUXVAR_SS)

       nauxvar = size(therm_soe%aux_vars_ss)

       if (nvar+iauxvar_off > nauxvar) then
          write(iulog,*) 'nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call SOEThermalEnthalpyAuxGetRData(therm_soe%aux_vars_ss, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case default
       write(iulog,*) 'auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select


  end subroutine SOETHGetAuxVars

  !------------------------------------------------------------------------
  subroutine SOETHGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType              , only : goveqn_base_type
    use ConnectionSetType                      , only : connection_set_type
    use ConditionType                          , only : condition_type
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxSetRData
    use ThermalEnthalpySoilAuxMod              , only : ThermalEnthalpySoilAuxVarSetRValues
    use ThermalEnthalpySoilAuxMod              , only : ThermalEnthalpySoilAuxVarGetRValues
    use GoveqnThermalEnthalpySoilType          , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType          , only : goveqn_richards_ode_pressure_type
    use RichardsODEPressureAuxMod              , only : RichardsODEPressureAuxVarSetRValues
    use RichardsODEPressureAuxMod              , only : RichardsODEPressureAuxVarGetRValues
    use CouplingVariableType                   , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)   , pointer :: cur_goveq_1
    class(goveqn_base_type)   , pointer :: cur_goveq_2
    !
    type(connection_set_type) , pointer :: cur_conn_set_2
    type(condition_type)      , pointer :: cur_cond_2
    type(coupling_variable_type), pointer :: cpl_var_1
    PetscInt                            :: idx
    PetscInt, pointer                   :: ids(:)
    PetscInt                            :: nauxvar
    PetscInt                            :: iauxvar
    PetscInt                            :: var_type
    PetscInt                            :: bc_idx
    PetscInt                            :: bc_offset
    PetscInt                            :: bc_rank_in_cpl_eqn
    PetscReal                           :: var_value
    PetscReal, pointer                  :: var_values(:)
    PetscBool                           :: bc_found
    PetscBool                           :: is_bc

    cpl_var_1 => cur_goveq_1%coupling_vars%first
    do
       if (.not.associated(cpl_var_1)) exit

       ! Does cur_goveq_1 needs ivar-th variable from cur_goveq_2?
       if (cpl_var_1%rank_of_coupling_goveqn == &
            cur_goveq_2%rank_in_soe_list) then

          var_type           = cpl_var_1%variable_type
          is_bc              = cpl_var_1%variable_is_bc_in_coupling_goveqn
          bc_offset          = cpl_var_1%offset_of_bc_in_current_goveqn
          bc_rank_in_cpl_eqn = cpl_var_1%rank_of_bc_in_coupling_goveqn

          if (.not.is_bc) then
             
             ! Get the data
             select type(cur_goveq_2)
             class is (goveqn_thermal_enthalpy_soil_type)

                nauxvar = size(cur_goveq_2%aux_vars_in)

                allocate(var_values(nauxvar))
                allocate(ids(nauxvar))

                do iauxvar = 1, nauxvar
                   ids(iauxvar) = iauxvar
                enddo

                call ThermalEnthalpySoilAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     nauxvar, ids, var_values)

             class is (goveqn_richards_ode_pressure_type)
                nauxvar = size(cur_goveq_2%aux_vars_in)

                allocate(var_values(nauxvar))
                allocate(ids(nauxvar))

                do iauxvar = 1, nauxvar
                   ids(iauxvar) = iauxvar
                enddo

                call RichardsODEPressureAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     nauxvar, ids, var_values)
                
             class default
                write(iulog,*)'Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))

             end select

             ! Set the data
             select type(cur_goveq_1)
             class is (goveqn_thermal_enthalpy_soil_type)
                call ThermalEnthalpySoilAuxVarSetRValues(cur_goveq_1%aux_vars_in, var_type, &
                     nauxvar, ids, var_values)

             class is (goveqn_richards_ode_pressure_type)
                call RichardsODEPressureAuxVarSetRValues(cur_goveq_1%aux_vars_in, var_type, &
                     nauxvar, ids, var_values)
                
                
             class default
                write(iulog,*)'SOEThermalEnthalpyGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             deallocate(ids       )
             deallocate(var_values)

          else

             ! Get the appropriate pointer to the BC from cur_goveq_2
             bc_idx   = 1
             bc_found = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit
                cur_conn_set_2 => cur_cond_2%conn_set
                
                ! Is this the appropriate BC?
                if (bc_idx == bc_rank_in_cpl_eqn) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'SOEThermalEnthalpyGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_conn_set_2%num_connections /= &
                  cpl_var_1%num_cells ) then
                write(iulog,*) 'conn_set_2%num_connections = ', cur_conn_set_2%num_connections
                write(iulog,*) 'cpl_var_1%num_cells        = ',cpl_var_1%num_cells
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             allocate(ids       (cpl_var_1%num_cells))
             allocate(var_values(cpl_var_1%num_cells))

             ! Get the data
             select type(cur_goveq_2)
             class is (goveqn_richards_ode_pressure_type)
                do iauxvar = 1, cpl_var_1%num_cells
                   idx = cur_conn_set_2%conn(iauxvar)%GetIDDn()
                   call cur_goveq_2%aux_vars_in(idx)%GetValue(var_type, var_value)
                   var_values(iauxvar) = var_value
                enddo

             class is (goveqn_thermal_enthalpy_soil_type)

                ! Save the IDs to get the data from
                do iauxvar = 1, cpl_var_1%num_cells
                   ids(iauxvar) = cur_conn_set_2%conn(iauxvar)%GetIDDn()
                enddo

                call ThermalEnthalpySoilAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     size(cur_goveq_2%aux_vars_in), ids, var_values)
                
             class default
                write(iulog,*)'Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             ! Set the data
             select type(cur_goveq_1)
             class is (goveqn_richards_ode_pressure_type)
                do iauxvar = 1, cpl_var_1%num_cells
                   idx = iauxvar + bc_offset
                   var_value = var_values(iauxvar)
                   call cur_goveq_1%aux_vars_bc(idx)%SetValue(var_type, var_value)
                enddo

             class is (goveqn_thermal_enthalpy_soil_type)

                ! Save the IDs to set the data to
                do iauxvar = 1, cpl_var_1%num_cells
                   ids(iauxvar) = iauxvar + bc_offset
                enddo

                call ThermalEnthalpySoilAuxVarSetRValues(cur_goveq_1%aux_vars_bc, var_type, &
                     size(cur_goveq_1%aux_vars_bc), ids, var_values)
                
             class default
                write(iulog,*)'SOEThermalEnthalpyGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             deallocate(ids       )
             deallocate(var_values)

          endif

       endif

       cpl_var_1 => cpl_var_1%next
    enddo

  end subroutine SOETHGovEqnExchangeAuxVars

  !------------------------------------------------------------------------
  subroutine SOETHResidual(this, snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Performs residual function evaluation for the VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_TH
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)          :: this
    SNES                              :: snes
    Vec                               :: X
    Vec                               :: F
    PetscErrorCode                    :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: dm_id
    PetscInt                          :: nDM
    DM                      , pointer :: dms(:)
    Vec                     , pointer :: X_subvecs(:)
    Vec                     , pointer :: F_subvecs(:)
    class(goveqn_base_type) , pointer :: cur_goveq
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    PetscInt                          :: row , col
    PetscViewer                       :: viewer
    character(len=256)                :: string
    PetscReal               , pointer :: data(:)
    PetscInt                          :: ndata
    PetscInt                          :: var_type
    PetscInt                          :: offset

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))
    allocate(F_subvecs(    nDM))

    ! Get vectors (X,F) for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, X_subvecs, &
         ierr); CHKERRQ(ierr)
    call DMCompositeGetAccessArray(this%solver%dm, F, nDM, PETSC_NULL_INTEGER, F_subvecs, &
         ierr); CHKERRQ(ierr)


    offset = 0

    select case (this%itype)
    case(SOE_TH)

       ! 1) {soln_prev}  ---> sim_aux()
       call SOETHUpdateAuxVars(this, X)

       ! 2) GE ---> GetFromSimAux()
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit

          select type(cur_goveq)

          class is (goveqn_thermal_enthalpy_soil_type)

             ndata = cur_goveq%mesh%ncells_local
             allocate(data(ndata))

             var_type = VAR_TEMPERATURE
             call SOETHGetAuxVars(this, AUXVAR_INTERNAL, var_type, &
                  offset, ndata, data)
             call cur_goveq%SetFromSOEAuxVarsIntrn(var_type, ndata, data)

             offset = offset + ndata

          class is (goveqn_richards_ode_pressure_type)

             ndata = cur_goveq%mesh%ncells_local
             allocate(data(ndata))

             var_type = VAR_PRESSURE
             call SOETHGetAuxVars(this, AUXVAR_INTERNAL, var_type, &
                  offset, ndata, data)
             call cur_goveq%SetFromSOEAuxVarsIntrn(var_type, ndata, data)

             offset = offset + ndata

             deallocate(data)

          class default
             write(iulog,*) 'Unknown goveqn_type'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select

          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'SOEThermalEnthalpyPreSolve: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! 3) GE_1 <---> GE_2 exchange AuxVars()
    do row = 1, this%ngoveqns
       do col = row+1, this%ngoveqns
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call SOETHGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
          call SOETHGovEqnExchangeAuxVars(cur_goveq_2, cur_goveq_1)
       enddo
    enddo

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)

       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%UpdateAuxVars()

       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%UpdateAuxVars()

       end select
          cur_goveq => cur_goveq%next
    enddo

    ! Call Residual
    dm_id = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       dm_id = dm_id + 1

       call VecZeroEntries(F_subvecs(dm_id), ierr); CHKERRQ(ierr)

       call cur_goveq%ComputeResidual( &
            X_subvecs(dm_id),          &
            F_subvecs(dm_id),          &
            ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)
    call DMCompositeRestoreAccessArray(this%solver%dm, F, nDM, PETSC_NULL_INTEGER, &
         F_subvecs, ierr); CHKERRQ(ierr)
    
    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)
    deallocate(F_subvecs)

  end subroutine SOETHResidual

  !------------------------------------------------------------------------
  subroutine SOETHJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes jacobian for the
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type) :: this
    SNES                            :: snes
    Vec                             :: X
    Mat                             :: A
    Mat                             :: B
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: row
    PetscInt                        :: col
    PetscInt                        :: nDM
    PetscInt                        :: icell

    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    IS,pointer                      :: is(:)
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Mat, pointer                    :: B_submats(:,:)
    PetscViewer                     :: viewer
    character(len=256)              :: string


    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Initialize the matrix
    call MatZeroEntries(B, ierr); CHKERRQ(ierr)

    ! Get submatrices
    allocate(is(nDM))
    allocate(B_submats(nDM,nDM))
    call DMCompositeGetLocalISs(this%solver%dm, is, ierr); CHKERRQ(ierr)
    do row = 1,nDM
       do col = 1,nDM
          call MatGetLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
       enddo
    enddo

    ! Jacobian and JacobianOffDiag
    row = 0
    cur_goveq_1 => this%goveqns
    do
       if (.not.associated(cur_goveq_1)) exit

       row = row + 1

       call cur_goveq_1%ComputeJacobian( &
            X_subvecs(row),              &
            B_submats(row,row),          &
            B_submats(row,row),          &
            ierr); CHKERRQ(ierr)

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          ! J = dF_1/dx_2
          call cur_goveq_1%ComputeOffDiagJacobian( &
               X_subvecs(row),                     &
               X_subvecs(col),                     &
               B_submats(row,col),                 &
               B_submats(row,col),                 &
               cur_goveq_2%id,                     &
               cur_goveq_2%rank_in_soe_list,       &
               ierr); CHKERRQ(ierr)

          ! J = dF_2/dx_1
          call cur_goveq_2%ComputeOffDiagJacobian( &
               X_subvecs(col),                     &
               X_subvecs(row),                     &
               B_submats(col,row),                 &
               B_submats(col,row),                 &
               cur_goveq_1%id,                     &
               cur_goveq_1%rank_in_soe_list,       &
               ierr); CHKERRQ(ierr)

          cur_goveq_2 => cur_goveq_2%next
       enddo

       cur_goveq_1 => cur_goveq_1%next
    enddo

    ! Restore vectors (X) for individual GEs
    call DMCompositeRestoreAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Restore submatrices
    do row = 1,nDM
       do col = 1,nDM
          call MatRestoreLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
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
    deallocate(dms       )
    deallocate(X_subvecs )
    deallocate(is        )
    deallocate(B_submats )

  end subroutine SOETHJacobian

  !------------------------------------------------------------------------
  subroutine SOETHSetDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the VSFM solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    use SystemOfEquationsThermalEnthalpyAuxMod, only : SOEThermalEnthalpyAuxSetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)          :: this
    PetscInt, intent(in)                                   :: var_type
    PetscInt                                               :: soe_auxvar_type
    PetscInt                                               :: soe_auxvar_id
    PetscReal                                              :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                               :: iauxvar
    PetscInt                                               :: iauxvar_off
    PetscInt                                               :: nauxvar
    type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
    case default
       write(iulog,*) 'VSFMSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOESetDataFromCLM: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call SOEThermalEnthalpyAuxSetRData(auxvars, var_type, &
         nauxvar, iauxvar_off, data_1d)

  end subroutine SOETHSetDataFromCLM

  !------------------------------------------------------------------------
  subroutine SOETHGetDataForCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the VSFM solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    use SystemOfEquationsThermalEnthalpyAuxMod, only : SOEThermalEnthalpyAuxGetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)          :: this
    PetscInt, intent(in)                                   :: var_type
    PetscInt                                               :: soe_auxvar_type
    PetscInt                                               :: soe_auxvar_id
    PetscReal                                              :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                               :: iauxvar
    PetscInt                                               :: iauxvar_off
    PetscInt                                               :: nauxvar
    type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
    case default
       write(iulog,*) 'VSFMSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOESetDataFromCLM: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call SOEThermalEnthalpyAuxGetRData(auxvars, var_type, &
         nauxvar, iauxvar_off, data_1d)

  end subroutine SOETHGetDataForCLM

  !------------------------------------------------------------------------
  subroutine SOETHUpdateConnections(this, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use SystemOfEquationsVSFMType     , only : VSFMSOEUpdateBCConnections
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_th_type)      :: this
    PetscInt, intent(in)          :: mpp_id
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type), pointer :: cur_goveq_1
    class(goveqn_base_type), pointer :: cur_goveq_2
    PetscInt                         :: row, col
    PetscBool                        :: goveq_1_is_rich
    PetscBool                        :: goveq_1_is_ther
    PetscBool                        :: goveq_2_is_rich
    PetscBool                        :: goveq_2_is_ther

    do row = 1, this%ngoveqns
       do col = row+1, this%ngoveqns
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)

          goveq_1_is_rich = PETSC_FALSE
          goveq_1_is_ther = PETSC_FALSE
          goveq_2_is_rich = PETSC_FALSE
          goveq_2_is_ther = PETSC_FALSE

          select type(cur_goveq_1)
             class is (goveqn_richards_ode_pressure_type)
                goveq_1_is_rich = PETSC_TRUE
             class is (goveqn_thermal_enthalpy_soil_type)
                goveq_1_is_ther = PETSC_TRUE
             class default
                write(iulog,*)'Unknown goveqn class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
           end select

          select type(cur_goveq_2)
             class is (goveqn_richards_ode_pressure_type)
                goveq_2_is_rich = PETSC_TRUE
             class is (goveqn_thermal_enthalpy_soil_type)
                goveq_2_is_ther = PETSC_TRUE
             class default
                write(iulog,*)'Unknown goveqn class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
           end select

          if (goveq_1_is_rich .and. goveq_2_is_rich) then
             call VSFMSOEUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
          else if (goveq_1_is_ther .and. goveq_2_is_ther) then
             call SOETHUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
          endif

       enddo
    enddo

  end subroutine SOETHUpdateConnections

  !------------------------------------------------------------------------

  subroutine SOETHUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use ConnectionSetType             , only : connection_set_type
    use ConditionType                 , only : condition_type
    use MultiPhysicsProbConstants     , only : MPP_TH_SNES_CLM
    use ThermalEnthalpySoilAuxType    , only : therm_enthalpy_soil_auxvar_type
    use CouplingVariableType          , only : coupling_variable_type
    use RichardsODEPressureAuxType    , only : RichODEPressureAuxVarCopy
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)               , pointer    :: cur_goveq_1
    class(goveqn_base_type)               , pointer    :: cur_goveq_2
    PetscInt                              , intent(in) :: mpp_id
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type)             , pointer    :: cur_conn_set_1
    type(connection_set_type)             , pointer    :: cur_conn_set_2
    type(condition_type)                  , pointer    :: cur_cond_1
    type(condition_type)                  , pointer    :: cur_cond_2
    type(therm_enthalpy_soil_auxvar_type) , pointer    :: aux_vars_bc_1(:)
    type(therm_enthalpy_soil_auxvar_type) , pointer    :: aux_vars_bc_2(:)
    type(therm_enthalpy_soil_auxvar_type)              :: tmp_aux_var_bc_1
    type(therm_enthalpy_soil_auxvar_type)              :: tmp_aux_var_bc_2
    type(coupling_variable_type)          , pointer    :: cpl_var_1
    PetscInt                                           :: iauxvar
    PetscInt                                           :: sum_conn_1
    PetscInt                                           :: sum_conn_2
    PetscInt                                           :: bc_idx
    PetscInt                                           :: rank_of_bc_in_cur_eqn
    PetscInt                                           :: bc_rank_in_cpl_eqn
    PetscInt                                           :: id_dn_1
    PetscInt                                           :: id_dn_2
    PetscReal                                          :: x_up   , y_up, z_up
    PetscReal                                          :: x_dn   , y_dn, z_dn
    PetscReal                                          :: x_half , y_half, z_half
    PetscReal                                          :: dx     , dy, dz
    PetscReal                                          :: dist   , dist_up, dist_dn
    PetscReal                                          :: area_1 , area_2
    PetscBool                                          :: bc_found
    PetscBool                                          :: is_bc

    PetscReal, dimension(3) :: vec3Swap
    PetscReal :: valSwap

    select type(cur_goveq_1)
    class is (goveqn_thermal_enthalpy_soil_type)
        aux_vars_bc_1 => cur_goveq_1%aux_vars_bc
    class default
        write(iulog,*)'VSFMSOEUpdateBCConnections: Unknown class'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    select type(cur_goveq_2)
    class is (goveqn_thermal_enthalpy_soil_type)
        aux_vars_bc_2 => cur_goveq_2%aux_vars_bc
    class default
        write(iulog,*)'VSFMSOEUpdateBCConnections: Unknown class'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    cpl_var_1 => cur_goveq_1%coupling_vars%first

    do
       if (.not.associated(cpl_var_1)) exit
    
       if (cpl_var_1%rank_of_coupling_goveqn == &
           cur_goveq_2%rank_in_soe_list) then

          is_bc                 = cpl_var_1%variable_is_bc_in_coupling_goveqn
          rank_of_bc_in_cur_eqn = cpl_var_1%rank_of_bc_in_current_goveqn
          bc_rank_in_cpl_eqn    = cpl_var_1%rank_of_bc_in_coupling_goveqn

          if (is_bc) then

             bc_idx = 1
             sum_conn_1 = 0
             bc_found = PETSC_FALSE
             cur_cond_1 => cur_goveq_1%boundary_conditions%first
             do
                if (.not.associated(cur_cond_1)) exit
                cur_conn_set_1 => cur_cond_1%conn_set
                if (bc_idx == rank_of_bc_in_cur_eqn) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                sum_conn_1 = sum_conn_1 + cur_conn_set_1%num_connections
                cur_cond_1 => cur_cond_1%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             bc_idx = 1
             sum_conn_2 = 0
             bc_found = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit
                cur_conn_set_2 => cur_cond_2%conn_set
                if (bc_idx == bc_rank_in_cpl_eqn) then
                   bc_found = PETSC_TRUE
                   exit
                endif

                bc_idx = bc_idx + 1
                sum_conn_2 = sum_conn_2 + cur_conn_set_2%num_connections
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_conn_set_1%num_connections /= &
                 cur_conn_set_2%num_connections ) then
                write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: Number of cells in BC of eq-1 and eq-2 are not same.'
                write(iulog,*) 'Eq-1'
                write(iulog,*) '   Name      : ',trim(cur_goveq_1%name)
                write(iulog,*) '   BC name   : ',trim(cur_cond_1%name)
                write(iulog,*) '   BC ncells : ',cur_conn_set_1%num_connections
                write(iulog,*) 'Eq-2'
                write(iulog,*) '   Name      : ',trim(cur_goveq_2%name)
                write(iulog,*) '   BC name   : ',trim(cur_cond_2%name)
                write(iulog,*) '   BC ncells : ',cur_conn_set_2%num_connections
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (cur_goveq_2%rank_in_soe_list > cur_goveq_1%rank_in_soe_list) then
                cur_cond_2%swap_order = PETSC_TRUE
             else
                cur_cond_1%swap_order = PETSC_TRUE
             endif

             do iauxvar = 1, cur_conn_set_1%num_connections


               !
               ! Eq-1
               !              unit_vec
               !            <----------
               !        "dn"           "up"
               !    _____________ _____________
               !   |             |             |
               !   |             |             |
               !   |      x      o      x      |
               !   |             |             |
               !   |   mesh-eq1  |   mesh-eq2  |
               !   |_____________|_____________|
               !
               ! Eq-2
               !        "up"           "dn"
               !            ---------->
               !              unit_vec
               !

                id_dn_1 = cur_conn_set_1%conn(iauxvar)%GetIDDn()
                id_dn_2 = cur_conn_set_2%conn(iauxvar)%GetIDDn()

                x_dn = cur_goveq_1%mesh%x(id_dn_1)
                y_dn = cur_goveq_1%mesh%y(id_dn_1)
                z_dn = cur_goveq_1%mesh%z(id_dn_1)

                x_up = cur_goveq_2%mesh%x(id_dn_2)
                y_up = cur_goveq_2%mesh%y(id_dn_2)
                z_up = cur_goveq_2%mesh%z(id_dn_2)

                dx = -(x_up - x_dn)
                dy = -(y_up - y_dn)
                dz = -(z_up - z_dn)

                dist = (dx**2.d0 + dy**2.d0 + dz**2.d0)**0.5d0

                call cur_conn_set_1%conn(iauxvar)%SetDistUnitVec(dx/dist, dy/dist, dz/dist)

                x_half = (x_dn + x_up)/2.d0
                y_half = (y_dn + y_up)/2.d0
                z_half = (z_dn + z_up)/2.d0

                call cur_conn_set_1%conn(iauxvar)%SetIDUp(id_dn_2)

                ! unit_vector for eq1 = -unit_vector for eq2
                call cur_conn_set_2%conn(iauxvar)%SetDistUnitVec(    &
                    -cur_conn_set_1%conn(iauxvar)%GetDistUnitVecX(), &
                    -cur_conn_set_1%conn(iauxvar)%GetDistUnitVecY(), &
                    -cur_conn_set_1%conn(iauxvar)%GetDistUnitVecZ() )

                ! dist_up eq2 = dist_dn eq1
                ! dist_up eq1 = dist_dn eq2
                call cur_conn_set_2%conn(iauxvar)%SetDistUp(cur_conn_set_1%conn(iauxvar)%GetDistDn())
                call cur_conn_set_1%conn(iauxvar)%SetDistUp(cur_conn_set_2%conn(iauxvar)%GetDistDn())

                call cur_conn_set_2%conn(iauxvar)%SetIDUp(id_dn_1)

                sum_conn_1 = sum_conn_1 + 1
                sum_conn_2 = sum_conn_2 + 1

                call RichODEPressureAuxVarCopy(tmp_aux_var_bc_1, aux_vars_bc_1(sum_conn_1))
                call RichODEPressureAuxVarCopy(tmp_aux_var_bc_2, aux_vars_bc_2(sum_conn_2))

                call RichODEPressureAuxVarCopy(aux_vars_bc_1(sum_conn_1), tmp_aux_var_bc_2)
                call RichODEPressureAuxVarCopy(aux_vars_bc_2(sum_conn_2), tmp_aux_var_bc_1)

             enddo

          endif
       endif

       cpl_var_1 => cpl_var_1%next
    enddo

  end subroutine SOETHUpdateBCConnections

#endif

end module SystemOfEquationsTHType
