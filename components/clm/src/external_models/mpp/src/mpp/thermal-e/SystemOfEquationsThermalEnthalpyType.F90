
module SystemOfEquationsThermalEnthalpyType

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for Thermal system-of-equations
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                              , only : iulog
  use mpp_abortutils                          , only : endrun
  use mpp_shr_log_mod                         , only : errMsg => shr_log_errMsg
  use SystemOfEquationsThermalEnthalpyAuxType , only : sysofeqns_thermal_enthalpy_auxvar_type
  use SystemOfEquationsBaseType               , only : sysofeqns_base_type
  use GoveqnThermalEnthalpySoilType           , only : goveqn_thermal_enthalpy_soil_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(sysofeqns_base_type) :: sysofeqns_thermal_enthalpy_type

     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_in(:)            !!< Internal state.
     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_bc(:)            !!< Boundary conditions.
     type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: aux_vars_ss(:)            !!< Source-sink.

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
     procedure, public :: Init                   => SOEThermalEnthalpyInit
     procedure, public :: SetSolnPrevCLM         => SOEThermalEnthalpySetSolnPrevCLM
     procedure, public :: GetSoln                => SOEThermalEnthalpyGetSoln
     procedure, public :: SetRDataFromCLM        => SOEThermalEnthalpySetRDataFromCLM
     procedure, public :: SetBDataFromCLM        => SOEThermalEnthalpySetBDataFromCLM
     procedure, public :: PreStepDT              => SOEThermalEnthalpyPreStepDT
     procedure, public :: PreSolve               => SOEThermalEnthalpyPreSolve
     procedure, public :: AddGovEqn              => SOEThermalEnthalpyAddGovEqn
     procedure, public :: SetDataFromCLM         => SOEThermalEnthalpySetDataFromCLM
     procedure, public :: GetDataForCLM          => SOEThermalEnthalpyGetDataForCLM
     procedure, public :: CreateVectorsForGovEqn => SOEThermalEnthalpyCreateVectorsForGovEqn
     procedure, public :: Residual               => SOEThermalEnthalpyResidual
     procedure, public :: Jacobian               => SOEThermalEnthalpyJacobian
     procedure, public :: PostSolve              => SOEThermalEnthalpyPostSolve
     procedure, public :: PostStepDT             => SOEThermalEnthalpyPostStepDT

  end type sysofeqns_thermal_enthalpy_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyInit(this)
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
    class(sysofeqns_thermal_enthalpy_type) :: this

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

  end subroutine SOEThermalEnthalpyInit

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : SOEBaseInit
    use MultiPhysicsProbConstants     , only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants     , only : GE_RE
    use MultiPhysicsProbConstants     , only : MESH_CLM_SOIL_COL
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type)              :: this
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

          goveq_soil%name                 = trim(name)
          goveq_soil%rank_in_soe_list     = this%ngoveqns
          goveq_soil%mesh_itype           = mesh_itype

          if (this%ngoveqns               == 1) then
             this%goveqns                 => goveq_soil
          else
             cur_goveqn%next              => goveq_soil
          endif

       case (GE_RE)

          allocate(goveq_richards)
          call goveq_richards%Setup()

          goveq_richards%name             = trim(name)
          goveq_richards%rank_in_soe_list = this%ngoveqns
          goveq_richards%mesh_itype       = mesh_itype

          if (this%ngoveqns == 1) then
             this%goveqns => goveq_richards
          else
             cur_goveqn%next => goveq_richards
          endif

       case default
          write(iulog,*) 'Unknown governing equation type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

  end subroutine SOEThermalEnthalpyAddGovEqn

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpySetSolnPrevCLM(this, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    PetscReal, intent (in)                 :: data_1d(:)
    !
    PetscInt                               :: vsize
    PetscErrorCode                         :: ierr
    PetscScalar, pointer                   :: x_p(:)

    call VecGetSize(this%solver%soln_prev_clm, vsize, ierr); CHKERRQ(ierr)

    if ( vsize /= size(data_1d)) then
       write(iulog,*) 'SOEThermalEnthalpySetSolnPrevCLM: vsize /= size(data_1d)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayF90(this%solver%soln_prev_clm, x_p, ierr); CHKERRQ(ierr)
    x_p(:) = data_1d(:)
    call VecRestoreArrayF90(this%solver%soln_prev_clm, x_p, ierr); CHKERRQ(ierr)

  end subroutine SOEThermalEnthalpySetSolnPrevCLM

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyGetSoln(this, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    PetscReal, intent (out)                :: data_1d(:)
    !
    PetscInt                               :: vsize
    PetscErrorCode                         :: ierr
    PetscScalar, pointer                   :: x_p(:)

    call VecGetSize(this%solver%soln, vsize, ierr); CHKERRQ(ierr)

    if ( vsize /= size(data_1d)) then
       write(iulog,*) 'SOEThermalEnthalpySetSolnPrevCLM: vsize /= size(data_1d)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayF90(this%solver%soln, x_p, ierr); CHKERRQ(ierr)
    data_1d(:) = x_p(:)
    call VecRestoreArrayF90(this%solver%soln, x_p, ierr); CHKERRQ(ierr)

  end subroutine SOEThermalEnthalpyGetSoln

  !------------------------------------------------------------------------

  subroutine SOEThermalEnthalpySetRDataFromCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants              , only : SOE_RE_ODE
    use MultiPhysicsProbConstants              , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants              , only : AUXVAR_BC
    use MultiPhysicsProbConstants              , only : AUXVAR_SS
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxSetRData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    PetscInt, intent(in)                   :: var_type
    PetscInt                               :: soe_auxvar_type
    PetscInt                               :: soe_auxvar_id
    PetscReal, pointer                     :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_off
    PetscInt                               :: nauxvar

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in

       call SOEThermalEnthalpyAuxSetRData (this%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case(AUXVAR_BC)
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)

       call SOEThermalEnthalpyAuxSetRData(this%aux_vars_bc, var_type, &
            nauxvar, iauxvar_off, data_1d)

    case(AUXVAR_SS)
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
 
       call SOEThermalEnthalpyAuxSetRData(this%aux_vars_ss, var_type, &
            nauxvar, iauxvar_off, data_1d)

   case default
       write(iulog,*) 'SOEThermalEnthalpySetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SOEThermalEnthalpySetRDataFromCLM
  
  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpySetBDataFromCLM(this, soe_auxvar_type, var_type, &
       data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants              , only : SOE_RE_ODE
    use MultiPhysicsProbConstants              , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants              , only : AUXVAR_BC
    use MultiPhysicsProbConstants              , only : AUXVAR_SS
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxSetBData
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    PetscInt, intent(in)                   :: var_type
    PetscInt                               :: soe_auxvar_type
    PetscBool, pointer                     :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_off
    PetscInt                               :: nauxvar

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in

       call SOEThermalEnthalpyAuxSetBData(this%aux_vars_in, var_type, nauxvar, &
            iauxvar_off, data_1d)

    case default
       write(iulog,*) 'SOEThermalEnthalpySetBDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SOEThermalEnthalpySetBDataFromCLM
  
  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyPreStepDT(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    !
    PetscErrorCode                         :: ierr
    class(goveqn_base_type),pointer        :: cur_goveq

    call VecCopy(this%solver%soln_prev_clm, this%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(this%solver%soln_prev_clm, this%solver%soln     , ierr); CHKERRQ(ierr)

  end subroutine SOEThermalEnthalpyPreStepDT

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_THERMAL_EBASED
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    class(goveqn_base_type),pointer        :: cur_goveq
    PetscInt                               :: offset

    offset = 0

    select case (this%itype)
    case(SOE_THERMAL_EBASED)

       ! 1) {soln_prev}  ---> sim_aux()
       call SOEThermalEnthalpyUpdateAuxVars(this, this%solver%soln_prev)

       ! 2) GE ---> GetFromSimAux()
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)

          class is (goveqn_thermal_enthalpy_soil_type)
             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             offset = offset + cur_goveq%mesh%ncells_local

             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

             call cur_goveq%UpdateAuxVars()
             call cur_goveq%PreSolve()

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

  end subroutine SOEThermalEnthalpyPreSolve

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyUpdateAuxVars(therm_soe, X)
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
    class(sysofeqns_thermal_enthalpy_type) :: therm_soe
    Vec                                    :: X
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: dm_id
    PetscInt                               :: nDM
    DM, pointer                            :: dms(:)
    Vec, pointer                           :: X_subvecs(:)
    PetscInt                               :: size
    PetscInt                               :: offset
    PetscErrorCode                         :: ierr

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
    do dm_id = 1, nDM
       call SOEThermalEnthalpySetAuxVars(therm_soe, AUXVAR_INTERNAL, VAR_TEMPERATURE, &
            X_subvecs(dm_id), offset)
       call VecGetSize(X_subvecs(dm_id), size, ierr); CHKERRQ(ierr)
       offset = offset + size
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(therm_soe%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)


  end subroutine SOEThermalEnthalpyUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpySetAuxVars(therm_soe, auxvar_type, var_type, &
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
    class(sysofeqns_thermal_enthalpy_type) :: therm_soe
    PetscInt                               :: auxvar_type
    PetscInt, intent(in)                   :: var_type
    Vec                                    :: var_vec
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer                     :: var_p(:)
    PetscInt                               :: nauxvar
    PetscInt                               :: nvar
    PetscInt, optional                     :: offset
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_off
    PetscErrorCode                         :: ierr

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
          write(iulog,*) 'SOEThermalEnthalpySetAuxVars: nvar+iauxvar_off > nauxvar.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call VecGetArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)

       call SOEThermalEnthalpyAuxSetRData(therm_soe%aux_vars_in, var_type, &
            nauxvar, iauxvar_off, var_p)

       call VecRestoreArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)

    case default
       write(iulog,*) 'SOEThermalEnthalpySetAuxVars: auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select


  end subroutine SOEThermalEnthalpySetAuxVars

  !------------------------------------------------------------------------

  subroutine SOEThermalEnthalpyResidual(this, snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use GoverningEquationBaseType, only : goveqn_base_type
    use MultiPhysicsProbConstants, only : SOE_THERMAL_EBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    SNES                                   :: snes
    Vec                                    :: X
    Vec                                    :: F
    PetscErrorCode                         :: ierr
    !
    !
    class(goveqn_base_type),pointer        :: cur_goveq
    class(goveqn_base_type),pointer        :: cur_goveq_1
    class(goveqn_base_type),pointer        :: cur_goveq_2
    PetscInt                               :: dm_id
    PetscInt                               :: row
    PetscInt                               :: col
    PetscInt                               :: nDM
    DM, pointer                            :: dms(:)
    Vec, pointer                           :: X_subvecs(:)
    Vec, pointer                           :: F_subvecs(:)
    PetscInt                               :: offset
    PetscViewer :: viewer
    
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

    ! 1) {X}  ---> sim_aux()
    call SOEThermalEnthalpyUpdateAuxVars(this, X)

    ! 2.1) GE ---> GetFromSimAux()
    ! Get pointers to governing-equations
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)

          call cur_goveq%UpdateAuxVarsIntrn()
          offset = offset + cur_goveq%mesh%ncells_local

       end select

       cur_goveq => cur_goveq%next
    enddo

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
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

    ! Restore vectors for individual GEs
    call DMCompositeRestoreAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)
    call DMCompositeRestoreAccessArray(this%solver%dm, F, nDM, PETSC_NULL_INTEGER, &
         F_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)
    deallocate(F_subvecs)

  end subroutine SOEThermalEnthalpyResidual

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyJacobian(this, snes, X, A, B, ierr)
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
    class(sysofeqns_thermal_enthalpy_type) :: this
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

  end subroutine SOEThermalEnthalpyJacobian

  !------------------------------------------------------------------------

  subroutine SOEThermalEnthalpyGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType              , only : goveqn_base_type
    use ConnectionSetType                      , only : connection_set_type
    use ConditionType                          , only : condition_type
    use SystemOfEquationsThermalEnthalpyAuxMod , only : SOEThermalEnthalpyAuxSetRData
    use ThermalEnthalpySoilAuxMod              , only : ThermalEnthalpySoilAuxVarSetRValues
    use ThermalEnthalpySoilAuxMod              , only : ThermalEnthalpySoilAuxVarGetRValues
    use CouplingVariableType                   , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)       , pointer :: cur_goveq_1
    class(goveqn_base_type)       , pointer :: cur_goveq_2
    !
    type(connection_set_type)     , pointer :: cur_conn_set_2
    type(condition_type)          , pointer :: cur_cond_2
    type (coupling_variable_type) , pointer :: cpl_var_1
    PetscInt                                :: idx
    PetscInt                      , pointer :: ids(:)
    PetscInt                                :: iauxvar
    PetscInt                                :: var_type
    PetscInt                                :: bc_idx
    PetscInt                                :: bc_offset
    PetscInt                                :: bc_rank_in_cpl_eqn
    PetscReal                               :: var_value
    PetscReal                     , pointer :: var_values(:)
    PetscBool                               :: bc_found
    PetscBool                               :: is_bc

    cpl_var_1 => cur_goveq_1%coupling_vars%first
    do
       if (.not.associated(cpl_var_1)) exit

       ! Does cur_goveq_1 needs ivar-th variable from cur_goveq_2?
       if (cpl_var_1%rank_of_coupling_goveqn == &
            cur_goveq_2%rank_in_soe_list) then

          var_type                      = cpl_var_1%variable_type
          is_bc                       = cpl_var_1%variable_is_bc_in_coupling_goveqn
          bc_offset                     = cpl_var_1%offset_of_bc_in_current_goveqn
          bc_rank_in_cpl_eqn = cpl_var_1%rank_of_bc_in_coupling_goveqn

          if (.not.is_bc) then
             
             write(iulog,*) 'SOEThermalEnthalpyGovEqnExchangeAuxVars: Extend code to ' // &
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
                write(iulog,*) 'conn_set_2%num_connections        = ', cur_conn_set_2%num_connections
                write(iulog,*) 'cpl_var_1%num_cells               = ', cpl_var_1%num_cells
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             allocate(ids       (cpl_var_1%num_cells))
             allocate(var_values(cpl_var_1%num_cells))

             ! Save the IDs to get the data from
             do iauxvar = 1, cpl_var_1%num_cells
                ids(iauxvar) = cur_conn_set_2%conn(iauxvar)%GetIDDn()
             enddo

             ! Get the data
             select type(cur_goveq_2)
             class is (goveqn_thermal_enthalpy_soil_type)
                call ThermalEnthalpySoilAuxVarGetRValues(cur_goveq_2%aux_vars_in, var_type, &
                     size(cur_goveq_2%aux_vars_in), ids, var_values)
                
             class default
                write(iulog,*)'SOEThermalEnthalpyGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             ! Save the IDs to set the data to
             do iauxvar = 1, cpl_var_1%num_cells
                ids(iauxvar) = iauxvar + bc_offset
             enddo

             ! Set the data
             select type(cur_goveq_1)
             class is (goveqn_thermal_enthalpy_soil_type)
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

  end subroutine SOEThermalEnthalpyGovEqnExchangeAuxVars

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpySetDataFromCLM(this, soe_auxvar_type, var_type, &
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
    class(sysofeqns_thermal_enthalpy_type)                 :: this
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

  end subroutine SOEThermalEnthalpySetDataFromCLM

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyGetDataForCLM(this, soe_auxvar_type, var_type, &
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
    class(sysofeqns_thermal_enthalpy_type)                 :: this
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
       write(iulog,*) 'size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call SOEThermalEnthalpyAuxGetRData(auxvars, var_type, &
         nauxvar, iauxvar_off, data_1d)

  end subroutine SOEThermalEnthalpyGetDataForCLM

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyCreateVectorsForGovEqn(this)
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
    class(sysofeqns_thermal_enthalpy_type) :: this
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

  end subroutine SOEThermalEnthalpyCreateVectorsForGovEqn


  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyPostSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations after a successful call to the PETSc solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_THERMAL_EBASED
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode                  :: ierr

    call VecCopy(this%solver%soln, this%solver%soln_prev,ierr); CHKERRQ(ierr)

    select case (this%itype)
    case(SOE_THERMAL_EBASED)

       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
          class is (goveqn_thermal_enthalpy_soil_type)

             call cur_goveq%SetDataInSOEAuxVar(AUXVAR_INTERNAL, this%aux_vars_in)
             !call cur_goveq%SetDataInSOEAuxVar(AUXVAR_BC      , this%aux_vars_bc)

          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SOEThermalEnthalpyPostSolve

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyPostStepDT(this)
    !
    ! !DESCRIPTION:
    ! This subroutines make a copy of solution vector post StepDT is
    ! called.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_type) :: this
    PetscErrorCode             :: ierr

    call VecCopy(this%solver%soln_prev, this%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine SOEThermalEnthalpyPostStepDT


#endif

end module SystemOfEquationsThermalEnthalpyType
