
module SystemOfEquationsVSFMType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for the Variable Saturate Flow Model (VSFM) system-of-equations
  !
#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use SystemOfEquationsVSFMAuxType   , only : sysofeqns_vsfm_auxvar_type
  use SystemOfEquationsBaseType      , only : sysofeqns_base_type
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

  type, public, extends(sysofeqns_base_type) :: sysofeqns_vsfm_type
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_in(:)            ! Internal state.
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_bc(:)            ! Boundary conditions.
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_ss(:)            ! Source-sink.
     type (sysofeqns_vsfm_auxvar_type), pointer :: aux_vars_conn_in(:)       ! Internal connections

     
     PetscInt, pointer                          :: soe_auxvars_bc_offset (:) ! Cummulative sum of number of control volumes associated with each boundary condition.
     PetscInt, pointer                          :: soe_auxvars_ss_offset (:) ! Cummulative sum of number of control volumes associated with each source-sink condition.
     PetscInt, pointer                          :: soe_auxvars_bc_ncells (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer                          :: soe_auxvars_ss_ncells (:) ! Number of control volumes associated with each source-sink condition.
     PetscInt                                   :: num_auxvars_in            ! Number of auxvars associated with internal state.
     PetscInt                                   :: num_auxvars_in_local      ! Number of auxvars associated with internal state.
     PetscInt                                   :: num_auxvars_bc            ! Number of auxvars associated with boundary condition.
     PetscInt                                   :: num_auxvars_ss            ! Number of auxvars associated with source-sink condition.
     PetscInt                                   :: num_auxvars_conn_in       ! Number of auxvars associated with internal connections
   contains
     procedure, public :: Init                   => VSFMSOEInit
     procedure, public :: Residual               => VSFMSOEResidual
     procedure, public :: Jacobian               => VSFMJacobian
     procedure, public :: PreSolve               => VSFMSOEPreSolve
     procedure, public :: SetDataFromCLM         => VSFMSOESetDataFromCLM
     procedure, public :: GetDataForCLM          => VSFMSOEGetDataForCLM
     procedure, public :: PostSolve              => VSFMSOEPostSolve
     procedure, public :: PostStepDT             => VSFMSPostStepDT
     procedure, public :: PreStepDT              => VSFMSPreStepDT
     procedure, public :: GetConditionNames      => VSFMSGetConditionNames
     procedure, public :: SetDataFromCLMForGhost => VSFMSOESetDataFromCLMForGhost
     procedure, public :: ComputeLateralFlux     => VSFMComputeLateralFlux
     procedure, public :: AddGovEqn              => VSFMAddGovEqn
     procedure, public :: AddGovEqnWithMeshRank  => VSFMAddGovEqnWithMeshRank
     procedure, public :: CreateVectorsForGovEqn => VSFMCreateVectorsForGovEqn
  end type sysofeqns_vsfm_type

  public :: VSFMSOESetAuxVars,         &
            VSFMSOEUpdateConnections,  &
            VSFMSOEUpdateBCConnections

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine VSFMSOEInit(this)
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
    class(sysofeqns_vsfm_type) :: this

    call SOEBaseInit(this)

    this%num_auxvars_in         = 0
    this%num_auxvars_in_local   = 0
    this%num_auxvars_bc         = 0
    this%num_auxvars_ss         = 0

    nullify(this%aux_vars_in           )
    nullify(this%aux_vars_bc           )
    nullify(this%aux_vars_ss           )
    nullify(this%aux_vars_conn_in      )

    nullify(this%soe_auxvars_bc_offset )
    nullify(this%soe_auxvars_ss_offset )
    nullify(this%soe_auxvars_bc_ncells )
    nullify(this%soe_auxvars_ss_ncells )

  end subroutine VSFMSOEInit

  !------------------------------------------------------------------------
  subroutine VSFMSOEResidual(this, snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Performs residual function evaluation for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    SNES                            :: snes
    Vec                             :: X
    Vec                             :: F
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: dm_id
    PetscInt                        :: nDM
    PetscInt                        :: offset
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Vec, pointer                    :: F_subvecs(:)
    class(goveqn_base_type),pointer :: cur_goveq
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt                        :: row, col
    PetscViewer                     :: viewer
    character(len=256)              :: string

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
    call VSFMSOEUpdateAuxVarsODE(this, X)

    ! 2.1) GE ---> GetFromSimAux()
    ! Get pointers to governing-equations
    offset = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)

          call cur_goveq%UpdateAuxVarsIntrn()
          call cur_goveq%UpdateAuxVarsBC()

          offset = offset + cur_goveq%mesh%ncells_local
       end select

       cur_goveq => cur_goveq%next
    enddo

    ! 3  ) GE_1 <---> GE_2 exchange AuxVars()
    do row = 1,nDM
       do col = row+1,nDM
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call VSFMSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
          call VSFMSOEGovEqnExchangeAuxVars(cur_goveq_2, cur_goveq_1)
       enddo
    enddo

    if (nDM > 1) then
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
             class is (goveqn_richards_ode_pressure_type)
             call cur_goveq%UpdateAuxVars()
          end select
          cur_goveq => cur_goveq%next
       enddo
    end if

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

  end subroutine VSFMSOEResidual

  !------------------------------------------------------------------------
  subroutine VSFMJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes jacobian for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
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

    IS,pointer                      :: is(:)
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Mat, pointer                    :: B_submats(:,:)
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
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

  end subroutine VSFMJacobian

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateAuxVarsODE(vsfm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE based on the input
    ! vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)    :: vsfm_soe
    Vec                           :: X
    !
    ! !LOCAL VARIABLES:

    select case (vsfm_soe%itype)
    case (SOE_RE_ODE)
       call VSFMSOEUpdateAuxVarsRichEqnODE(vsfm_soe, X)
    case default
       write(iulog,*) 'VSFMSOEUpdateAuxVars: unknown vsfm_soe%itype'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEUpdateAuxVarsODE

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateAuxVarsRichEqnODE(vsfm_soe, X)
    !
    ! !DESCRIPTION:
    ! Updates the SoE vars for the discretized ODE of Richards equation
    ! based on the input vector X
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: vsfm_soe
    Vec                        :: X
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
    call DMCompositeGetNumberDM(vsfm_soe%solver%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(vsfm_soe%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(vsfm_soe%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Update the SoE auxvars
    offset = 0
    do dm_id = 1, nDM
       call VSFMSOESetAuxVars(vsfm_soe, AUXVAR_INTERNAL, VAR_PRESSURE, &
            X_subvecs(dm_id), offset)
       call VecGetSize(X_subvecs(dm_id), size, ierr); CHKERRQ(ierr)
       offset = offset + size
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(vsfm_soe%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)

  end subroutine VSFMSOEUpdateAuxVarsRichEqnODE

  !------------------------------------------------------------------------
  subroutine VSFMSOESetAuxVars(vsfm_soe, auxvar_type, var_type, &
       var_vec, offset)
    !
    ! !DESCRIPTION:
    ! Set values in SoE auxvars.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: vsfm_soe
    PetscInt                                   :: auxvar_type
    PetscInt, intent(in)                       :: var_type
    Vec                                        :: var_vec
    !
    ! !LOCAL VARIABLES:
    PetscReal, pointer                         :: var_p(:)
    type (sysofeqns_vsfm_auxvar_type), pointer :: avars(:)
    PetscInt                                   :: nauxvar
    PetscInt                                   :: nvar
    PetscInt, optional                         :: offset
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscErrorCode                             :: ierr

    select case(auxvar_type)
    case (AUXVAR_INTERNAL)
       avars => vsfm_soe%aux_vars_in
    case default
       write(iulog,*) 'VSFMSOESetAuxVars: auxvar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    nauxvar = size(avars)

    call VecGetLocalSize(var_vec, nvar, ierr); CHKERRQ(ierr)

    if (nvar+iauxvar_off > nauxvar) then
       write(iulog,*) 'VSFMSOESetAuxVars: nvar+iauxvar_off > nauxvar.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call VecGetArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)
    do iauxvar = 1, nvar
       if (avars(iauxvar + iauxvar_off)%is_active) then
          call avars(iauxvar + iauxvar_off)%SetValue(var_type, var_p(iauxvar))
       end if
    enddo

    call VecRestoreArrayReadF90(var_vec, var_p, ierr); CHKERRQ(ierr)

  end subroutine VSFMSOESetAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMSOEBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Display information about the VSFM problem
    !
    ! !USES:
    use GoverningEquationBaseType       , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveqn

    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(iulog,*)'  SystemOfEqns_name  : ',trim(this%name)
    write(iulog,*)'  SystemOfEqns_itype : ',this%itype
    write(iulog,*)''

    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       call cur_goveqn%PrintInfo()
       cur_goveqn => cur_goveqn%next
    enddo

    write(iulog,*)'  No. of aux_vars_in : ',size(this%aux_vars_in)
    write(iulog,*)''
    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  end subroutine VSFMSOEBasePrintInfo

  !------------------------------------------------------------------------
  subroutine VSFMSOEPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations before PETSc solver are called.
    ! - For ODE based Richards equation VSFM solver, this subroutine calls
    !   indiviaul governing equations to compute the accumulation term
    !   that depends on the solution from previous time step.
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_RE_ODE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt :: offset

    select case (this%itype)
    case(SOE_RE_ODE)

       ! 1) {soln_prev}  ---> sim_aux()
       call VSFMSOEUpdateAuxVarsODE(this, this%solver%soln_prev)

       ! 2) GE ---> GetFromSimAux()
       offset = 0
       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
             class is (goveqn_richards_ode_pressure_type)

             call cur_goveq%GetFromSOEAuxVarsIntrn(this%aux_vars_in, offset)
             call cur_goveq%GetFromSOEAuxVarsBC(this%aux_vars_bc)
             call cur_goveq%GetFromSOEAuxVarsSS(this%aux_vars_ss)

             call cur_goveq%UpdateAuxVars()

             call cur_goveq%PreSolve()

             offset = offset + cur_goveq%mesh%ncells_local

          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEPreSolve

  !------------------------------------------------------------------------
  subroutine VSFMSOEPostSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations after a successful call to the PETSc solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : SOE_RE_ODE
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    use MultiPhysicsProbConstants     , only : AUXVAR_CONN_INTERNAL
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: iauxvar_in_off
    PetscInt                        :: iauxvar_bc_off
    PetscInt                        :: iauxvar_ss_off
    PetscInt                        :: iauxvar_conn_in_off
    PetscInt                        :: num_auxvars_filled
    PetscErrorCode                  :: ierr

    call VecCopy(this%solver%soln, this%solver%soln_prev,ierr); CHKERRQ(ierr)

    iauxvar_in_off       = 0
    iauxvar_bc_off       = 0
    iauxvar_ss_off       = 0
    iauxvar_conn_in_off  = 0

    select case (this%itype)
    case(SOE_RE_ODE)

       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
             class is (goveqn_richards_ode_pressure_type)

                call cur_goveq%SetDataInSOEAuxVar(AUXVAR_INTERNAL, this%aux_vars_in, &
                     iauxvar_in_off, num_auxvars_filled)
                iauxvar_in_off = iauxvar_in_off + num_auxvars_filled

                call cur_goveq%SetDataInSOEAuxVar(AUXVAR_BC      , this%aux_vars_bc, &
                     iauxvar_bc_off, num_auxvars_filled)
                iauxvar_bc_off = iauxvar_bc_off + num_auxvars_filled

                call cur_goveq%SetDataInSOEAuxVar(AUXVAR_SS      , this%aux_vars_ss, &
                     iauxvar_ss_off, num_auxvars_filled)
                iauxvar_ss_off = iauxvar_ss_off + num_auxvars_filled

                call cur_goveq%SetDataInSOEAuxVar(AUXVAR_CONN_INTERNAL, this%aux_vars_conn_in, &
                     iauxvar_conn_in_off, num_auxvars_filled)
                iauxvar_conn_in_off = iauxvar_conn_in_off + num_auxvars_filled


             end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'VSFMSOESetup: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEPostSolve

  !------------------------------------------------------------------------
  subroutine VSFMSOESetDataFromCLM(this, soe_auxvar_type, var_type, &
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
       nauxvar      = min(size(data_1d), this%soe_auxvars_bc_ncells(soe_auxvar_id))
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
       nauxvar      = min(size(data_1d), this%soe_auxvars_ss_ncells(soe_auxvar_id))
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

    do iauxvar = 1, nauxvar
       if (auxvars(iauxvar + iauxvar_off)%is_active) then
          call auxvars(iauxvar + iauxvar_off)%SetValue(var_type, data_1d(iauxvar))
       end if
    enddo

  end subroutine VSFMSOESetDataFromCLM

  !------------------------------------------------------------------------
  subroutine VSFMSOESetDataFromCLMForGhost(this, soe_auxvar_type, var_type, &
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
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_beg
    PetscInt                                   :: iauxvar_end
    PetscInt                                   :: nauxvar
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       nauxvar      = this%num_auxvars_in
       iauxvar_beg  = this%num_auxvars_in_local + 1
       iauxvar_end  = nauxvar
    case default
       write(iulog,*) 'VSFMSOESetDataFromCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOESetDataFromCLMForGhost: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = iauxvar_beg, iauxvar_end
       call auxvars(iauxvar)%SetValue(var_type, data_1d(iauxvar))
    enddo

  end subroutine VSFMSOESetDataFromCLMForGhost

  !------------------------------------------------------------------------
  subroutine VSFMSOEGetDataForCLM(this, soe_auxvar_type, var_type, &
       soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to extracted values from the VSFM solver
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : SOE_RE_ODE
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    use MultiPhysicsProbConstants, only : AUXVAR_CONN_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)                 :: this
    PetscInt, intent(in)                       :: var_type
    PetscInt                                   :: soe_auxvar_type
    PetscInt                                   :: soe_auxvar_id
    PetscReal                                  :: data_1d(:)
    PetscInt                                   :: nsize
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: iauxvar
    PetscInt                                   :: iauxvar_off
    PetscInt                                   :: nauxvar
    PetscReal                                  :: var_value
    type (sysofeqns_vsfm_auxvar_type), pointer :: auxvars(:)

    select case(soe_auxvar_type)
    case(AUXVAR_INTERNAL)
       auxvars      => this%aux_vars_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_in
    case(AUXVAR_SS)
       auxvars      => this%aux_vars_ss
       iauxvar_off  = this%soe_auxvars_ss_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_ss_ncells(soe_auxvar_id)
    case(AUXVAR_BC)
       auxvars      => this%aux_vars_bc
       iauxvar_off  = this%soe_auxvars_bc_offset(soe_auxvar_id)
       nauxvar      = this%soe_auxvars_bc_ncells(soe_auxvar_id)
    case(AUXVAR_CONN_INTERNAL)
       auxvars      => this%aux_vars_conn_in
       iauxvar_off  = 0
       nauxvar      = this%num_auxvars_conn_in
    case default
       write(iulog,*) 'VSFMSOEGetDataForCLM: Unknown soe_auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'VSFMSOEGetDataForCLM: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, size(data_1d)
       call auxvars(iauxvar + iauxvar_off)%GetValue(var_type, var_value)
       data_1d(iauxvar) = var_value
    enddo

  end subroutine VSFMSOEGetDataForCLM

  !------------------------------------------------------------------------
  subroutine VSFMComputeLateralFlux(this, dt)
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use MultiPhysicsProbConstants     , only : AUXVAR_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)    :: this
    PetscReal                     :: dt
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: iauxvar_ss_off
    PetscInt                        :: num_auxvars_filled

    ! Set timestep
    call this%SetDtime(dt)

    ! Do any pre-solve operations
    call this%PreSolve()

    iauxvar_ss_off = 0

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)

             call cur_goveq%ComputeLateralFlux()
             call cur_goveq%SetDataInSOEAuxVar(AUXVAR_SS, this%aux_vars_ss, &
                  iauxvar_ss_off, num_auxvars_filled)

             iauxvar_ss_off = iauxvar_ss_off + num_auxvars_filled

          end select

       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMComputeLateralFlux

  !------------------------------------------------------------------------
  subroutine VSFMSPreStepDT(this)
    !
    ! !DESCRIPTION:
    ! This subroutines copies solution vector at previous CLM time step
    ! before StepDT is called.
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode                  :: ierr
    PetscInt                        :: iauxvar

    call VecCopy(this%solver%soln_prev_clm, this%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(this%solver%soln_prev_clm, this%solver%soln     , ierr); CHKERRQ(ierr)

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%PreStepDT()
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMSPreStepDT

  !------------------------------------------------------------------------
  subroutine VSFMSPostStepDT(this)
    !
    ! !DESCRIPTION:
    ! This subroutines make a copy of solution vector post StepDT is
    ! called.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscErrorCode             :: ierr

    call VecCopy(this%solver%soln_prev, this%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine VSFMSPostStepDT

  !------------------------------------------------------------------------
  subroutine VSFMSGetConditionNames(this, cond_type, cond_type_to_exclude, &
       num_conds, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the total number and names of conditions (eg. boundary condition
    ! or source-sink) present
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    PetscInt, intent(in)            :: cond_type
    PetscInt, intent(in)            :: cond_type_to_exclude
    character (len=256), pointer    :: cond_names(:)
    PetscInt                        :: num_conds
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: nn
    PetscErrorCode                  :: ierr
    PetscInt                        :: num_conds_tmp
    character (len=256), pointer    :: cond_names_tmp(:)

    num_conds = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
             call cur_goveq%NumConditions(cond_type, cond_type_to_exclude, num_conds_tmp)
          class default
             write(iulog,*) 'VSFMSGetConditionNames: Supported for only goveqn_richards_ode_pressure_type'
             call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
       num_conds = num_conds + num_conds_tmp
       cur_goveq => cur_goveq%next
    enddo

    allocate(cond_names(num_conds))

    num_conds = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       call cur_goveq%GetCondNamesExcptCondItype( &
            cond_type, cond_type_to_exclude, num_conds_tmp, cond_names_tmp)
       if (num_conds_tmp > 0) then
          do nn = 1, num_conds_tmp
             num_conds = num_conds + 1
             cond_names(num_conds) = cond_names_tmp(nn)
          enddo
          deallocate(cond_names_tmp)
       endif

       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMSGetConditionNames

  !------------------------------------------------------------------------
  subroutine VSFMAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : SOEBaseInit
    use GoverningEquationBaseType     , only : goveqn_base_type
    use MultiPhysicsProbConstants     , only : GE_RE
    use MultiPhysicsProbConstants     , only : MESH_CLM_SOIL_COL
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscInt                   :: geq_type
    character(len=*)           :: name
    PetscInt                   :: mesh_itype
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    class(goveqn_base_type),pointer                     :: cur_goveqn
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
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

  end subroutine VSFMAddGovEqn

  !------------------------------------------------------------------------
  subroutine VSFMAddGovEqnWithMeshRank(this, geq_type, name, mesh_rank)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : SOEBaseInit
    use GoverningEquationBaseType     , only : goveqn_base_type
    use MultiPhysicsProbConstants     , only : GE_RE
    use MultiPhysicsProbConstants     , only : MESH_CLM_SOIL_COL
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type) :: this
    PetscInt                   :: geq_type
    character(len=*)           :: name
    PetscInt                   :: mesh_rank
    !
    ! !LOCAL VARIABLES:
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    class(goveqn_base_type), pointer                    :: cur_goveqn
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
    case (GE_RE)

       allocate(goveq_richards)
       call goveq_richards%Setup()

       goveq_richards%name              = trim(name)
       goveq_richards%rank_in_soe_list  = this%ngoveqns
       goveq_richards%mesh_rank         = mesh_rank

       if (this%ngoveqns == 1) then
          this%goveqns => goveq_richards
       else
          cur_goveqn%next => goveq_richards
       endif

    case default
       write(iulog,*) 'Unknown governing equation type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMAddGovEqnWithMeshRank

  !------------------------------------------------------------------------
  subroutine VSFMCreateVectorsForGovEqn(this)
    !
    ! !DESCRIPTION:
    ! Creates vectors required by each governing equation
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_type)      :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call cur_goveq%CreateVectors()
       class default
          write(iulog,*) 'Unsupported cur_goveq type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
       cur_goveq => cur_goveq%next
    enddo

  end subroutine VSFMCreateVectorsForGovEqn

  !------------------------------------------------------------------------

  subroutine VSFMSOEGovEqnExchangeAuxVars(cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use ConnectionSetType             , only : connection_set_type
    use ConditionType                 , only : condition_type
    use CouplingVariableType          , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type), pointer                  :: cur_conn_set_2
    type(condition_type),pointer                        :: cur_cond_2
    type (coupling_variable_type), pointer              :: cpl_var_1
    PetscInt                                            :: idx
    PetscInt                                            :: iauxvar
    PetscInt                                            :: var_type
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: bc_idx
    PetscInt                                            :: bc_offset
    PetscInt                                            :: bc_rank_in_cpl_eqn
    PetscReal                                           :: var_value
    PetscBool                                           :: bc_found
    PetscBool                                           :: is_bc

    cpl_var_1 => cur_goveq_1%coupling_vars%first
    do
       if (.not.associated(cpl_var_1)) exit

       if (cpl_var_1%rank_of_coupling_goveqn == &
           cur_goveq_2%rank_in_soe_list) then

          var_type           = cpl_var_1%variable_type
          is_bc              = cpl_var_1%variable_is_bc_in_coupling_goveqn
          bc_offset          = cpl_var_1%offset_of_bc_in_current_goveqn
          bc_rank_in_cpl_eqn = cpl_var_1%rank_of_bc_in_coupling_goveqn

          !GB: Presently assuming that auxvars for ALL cells is needed from
          !    the other goveqn. Additionally, it is assumed that the retrieved
          !    data is needed for 'aux_vars_in' and not for boundary-conditions or
          !    source-sinks
          if (.not.is_bc) then
             do ghosted_id = 1,cur_goveq_1%mesh%ncells_local

                select type(cur_goveq_2)
                class is (goveqn_richards_ode_pressure_type)
                   call cur_goveq_2%aux_vars_in(ghosted_id)%GetValue(var_type, var_value)
                class default
                   write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                select type(cur_goveq_1)
                class is (goveqn_richards_ode_pressure_type)
                   if (.not.is_bc) then
                      call cur_goveq_1%aux_vars_in(ghosted_id)%SetValue(var_type, var_value)
                   else
                      call cur_goveq_1%aux_vars_bc(ghosted_id+bc_offset)%SetValue(var_type, var_value)
                   endif

                class default
                   write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

             enddo

           else

              bc_idx = 1
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
                 cur_cond_2 => cur_cond_2%next
              enddo

              if (.not.bc_found) then
                 write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: BC not found'
                 call endrun(msg=errMsg(__FILE__, __LINE__))
              endif

              if (cur_conn_set_2%num_connections /= &
                  cpl_var_1%num_cells ) then
                 write(iulog,*) 'VSFMSOEGovEqnExchangeAuxVars: Number of '
                 call endrun(msg=errMsg(__FILE__, __LINE__))
              endif

              do iauxvar = 1, cpl_var_1%num_cells

                 ! Get value from cur_goveq_2
                 select type(cur_goveq_2)
                 class is (goveqn_richards_ode_pressure_type)
                    idx = cur_conn_set_2%conn(iauxvar)%GetIDDn()
                    call cur_goveq_2%aux_vars_in(idx)%GetValue(var_type, var_value)
                 class default
                    write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                    call endrun(msg=errMsg(__FILE__, __LINE__))
                 end select

                 idx = iauxvar + bc_offset

                 ! Set value from cur_goveq_1
                 select type(cur_goveq_1)
                 class is (goveqn_richards_ode_pressure_type)
                    call cur_goveq_1%aux_vars_bc(idx)%SetValue(var_type, var_value)
                 class default
                    write(iulog,*)'VSFMSOEGovEqnExchangeAuxVars: Unknown class'
                    call endrun(msg=errMsg(__FILE__, __LINE__))
                 end select
              enddo
           endif ! if (.not.is_bc)

        endif ! (

        cpl_var_1 => cpl_var_1%next
        
     enddo

  end subroutine VSFMSOEGovEqnExchangeAuxVars

  !------------------------------------------------------------------------
  subroutine VSFMSOEUpdateConnections(this, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type)    :: this
    PetscInt, intent(in)          :: mpp_id
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt :: row, col
    PetscBool :: check

    check = PETSC_TRUE

    do row = 1, this%ngoveqns
       do col = row+1, this%ngoveqns
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call VSFMSOEUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
       enddo
    enddo

  end subroutine VSFMSOEUpdateConnections

  !------------------------------------------------------------------------

  subroutine VSFMSOEUpdateBCConnections(cur_goveq_1, cur_goveq_2, mpp_id)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use ConnectionSetType             , only : connection_set_type
    use ConditionType                 , only : condition_type
    use RichardsODEPressureAuxType
    use SaturationFunction
    use MultiPhysicsProbConstants     , only : MPP_VSFM_SNES_CLM
    use CouplingVariableType          , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscInt, intent(in)            :: mpp_id
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type), pointer                  :: cur_conn_set_1
    type(connection_set_type), pointer                  :: cur_conn_set_2
    type(condition_type),pointer                        :: cur_cond_1
    type(condition_type),pointer                        :: cur_cond_2
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc_1(:)
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc_2(:)
    type (rich_ode_pres_auxvar_type)                    :: tmp_aux_var_bc_1
    type (rich_ode_pres_auxvar_type)                    :: tmp_aux_var_bc_2
    type (coupling_variable_type), pointer              :: cpl_var_1
    PetscInt                                            :: iauxvar
    PetscInt                                            :: sum_conn_1
    PetscInt                                            :: sum_conn_2
    PetscInt                                            :: bc_idx
    PetscInt                                            :: rank_of_bc_in_cur_eqn
    PetscInt                                            :: bc_rank_in_cpl_eqn
    PetscInt                                            :: id_dn_1
    PetscInt                                            :: id_dn_2
    PetscReal                                           :: x_up, y_up, z_up
    PetscReal                                           :: x_dn, y_dn, z_dn
    PetscReal                                           :: x_half, y_half, z_half
    PetscReal                                           :: dx, dy, dz
    PetscReal                                           :: dist, dist_up, dist_dn
    PetscReal                                           :: area_1, area_2
    PetscBool                                           :: bc_found
    PetscBool                                           :: is_bc

    PetscReal, dimension(3) :: vec3Swap
    PetscReal :: valSwap
    type(saturation_params_type) :: satParamsSwap

    select type(cur_goveq_1)
    class is (goveqn_richards_ode_pressure_type)
        aux_vars_bc_1 => cur_goveq_1%aux_vars_bc
    class default
        write(iulog,*)'VSFMSOEUpdateBCConnections: Unknown class'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    select type(cur_goveq_2)
    class is (goveqn_richards_ode_pressure_type)
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

             enddo ! iauxvar

          endif ! if bc
       endif ! ids are equal

       cpl_var_1 => cpl_var_1%next
       
    enddo

  end subroutine VSFMSOEUpdateBCConnections

#endif

end module SystemOfEquationsVSFMType
