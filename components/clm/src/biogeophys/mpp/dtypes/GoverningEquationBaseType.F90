module GoverningEquationBaseType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils     , only : endrun
  use mpp_shr_log_mod    , only : errMsg => shr_log_errMsg
  use MeshType           , only : mesh_type
  use ConditionType      , only : condition_list_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type, public :: goveqn_base_type

     character (len=256)             :: name                                 ! name of governing equation (GE)
     PetscInt                        :: id                                   ! identifier
     PetscInt                        :: id_in_list                           ! order in the list
     class(mesh_type),pointer        :: mesh                                 ! pointer to the mesh
     PetscInt                        :: mesh_itype                           ! type of mesh
     type(condition_list_type)       :: boundary_conditions                  ! boundary conditions to the GE
     type(condition_list_type)       :: source_sinks                         ! source/sinks in the GE

     PetscReal                       :: dtime                                ! time step [sec]

                                                                             ! Track variables supplied by other governing equations.
     PetscInt                        :: nvars_needed_from_other_goveqns      !
     PetscInt, pointer               :: var_ids_needed_from_other_goveqns(:) !

     PetscInt, pointer               :: ids_of_other_goveqns(:)              ! index of the other governing equation in the list
     PetscBool, pointer              :: is_bc_auxvar_type(:)                 !
     PetscInt, pointer               :: bc_auxvar_offset(:)                  !
     PetscInt, pointer               :: bc_auxvar_idx(:)                     !
     PetscInt, pointer               :: bc_auxvar_idx_of_other_goveqn(:)     !
     PetscInt, pointer               :: bc_auxvar_ncells(:)                  !

     class(goveqn_base_type),pointer :: next
   contains
     procedure, public :: Init                    => GoveqnBaseInit
     procedure, public :: Clean                   => GoveqnBaseClean
     procedure, public :: AllocVarsFromOtherGEs   => GoveqnBaseAllocVarsFromOtherGEs
     procedure, public :: DeallocVarsFromOtherGEs => GoveqnBaseDeallocVarsFromOtherGEs
     procedure, public :: PrintInfo               => GoveqnBasePrintInfo
     procedure, public :: UpdateAuxVars           => GoveqnBaseUpdateAuxVars
     procedure, public :: UpdateAuxVarsIntrn      => GoveqnBaseUpdateAuxVarsIntrn
     procedure, public :: UpdateAuxVarsBC         => GoveqnBaseUpdateAuxVarsBC
     procedure, public :: UpdateAuxVarsSS         => GoveqnBaseUpdateAuxVarsSS
     procedure, public :: PreSolve                => GoveqnBasePreSolve
     procedure, public :: IFunction               => GoveqnBaseIFunction
     procedure, public :: IJacobian               => GoveqnBaseIJacobian
     procedure, public :: IJacobianOffDiag        => GoveqnBaseIJacobianOffDiag
     procedure, public :: JacobianOffDiag         => GoveqnBaseJacobianOffDiag
     procedure, public :: Jacobian                => GoveqnBaseJacobian
     procedure, public :: Residual                => GoveqnBaseResidual
     procedure, public :: ComputeRHS              => GoveqnBaseComputeRHS
     procedure, public :: ComputeOperatorsDiag    => GoveqnBaseComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag => GoveqnBaseComputeOperatorsOffDiag
     procedure, public :: SetDtime                => GoveqnBaseSetDtime
     procedure, public :: GetNumConditions        => GoveqnBaseGetNumConditions
     procedure, public :: AddCondition            => GoveqnBaseAddCondition
     procedure, public :: UpdateConditionConnSet  => GoveqnBaseUpdateConditionConnSet
  end type goveqn_base_type
  !------------------------------------------------------------------------

  public :: GoveqnBasePrintInfo

contains

  !------------------------------------------------------------------------
  subroutine GoveqnBaseInit(this)
    !
    ! !DESCRIPTION:
    ! Initialze a GE object
    !
    ! !USES:
    use ConditionType ,only : ConditionListInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    this%name                            = ""
    this%id                              = -1
    this%id_in_list                      = -1
    this%dtime                           = 0.d0
    this%nvars_needed_from_other_goveqns = 0


    nullify(this%mesh)
    call ConditionListInit(this%boundary_conditions )
    call ConditionListInit(this%source_sinks        )

    nullify(this%var_ids_needed_from_other_goveqns  )
    nullify(this%ids_of_other_goveqns               )
    nullify(this%is_bc_auxvar_type                  )
    nullify(this%bc_auxvar_offset                   )
    nullify(this%bc_auxvar_idx                      )
    nullify(this%bc_auxvar_idx_of_other_goveqn      )
    nullify(this%bc_auxvar_ncells                   )
    nullify(this%next                               )

  end subroutine GoveqnBaseInit

  !------------------------------------------------------------------------
  subroutine GoveqnBaseClean(this)
    !
    ! !DESCRIPTION:
    ! Release allocated memory
    !
    ! !USES:
    use ConditionType ,only : ConditionListClean
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    call ConditionListClean(this%boundary_conditions )
    call ConditionListClean(this%source_sinks        )

  end subroutine GoveqnBaseClean


  !------------------------------------------------------------------------
  subroutine GoveqnBaseAllocVarsFromOtherGEs(this, nvars)
    !
    ! !DESCRIPTION:
    ! Allocate memory for tracking variables provided by other governing equations
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type), intent(inout):: this
    PetscInt, intent(in):: nvars

    if (this%nvars_needed_from_other_goveqns /= 0 ) then
       write(iulog,*) 'GoveqnBaseAllocVarsFromOtherGEs: Bad initialization'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    this%nvars_needed_from_other_goveqns = nvars

    allocate(this%var_ids_needed_from_other_goveqns(nvars)); this%var_ids_needed_from_other_goveqns(:) = 0
    allocate(this%ids_of_other_goveqns             (nvars)); this%ids_of_other_goveqns             (:) = 0
    allocate(this%is_bc_auxvar_type                (nvars)); this%is_bc_auxvar_type                (:) = PETSC_FALSE
    allocate(this%bc_auxvar_offset                 (nvars)); this%bc_auxvar_offset                 (:) = 0
    allocate(this%bc_auxvar_idx                    (nvars)); this%bc_auxvar_idx                    (:) = 0
    allocate(this%bc_auxvar_idx_of_other_goveqn    (nvars)); this%bc_auxvar_idx_of_other_goveqn    (:) = 0
    allocate(this%bc_auxvar_ncells                 (nvars)); this%bc_auxvar_ncells                 (:) = 0

  end subroutine GoveqnBaseAllocVarsFromOtherGEs


  !------------------------------------------------------------------------
  subroutine GoveqnBaseDeallocVarsFromOtherGEs(this)
    !
    ! !DESCRIPTION:
    ! Release allocated memory
    !
    implicit none

    ! !ARGUMENTS
    class(goveqn_base_type), intent(inout):: this

    if (this%nvars_needed_from_other_goveqns <= 0 ) then
       write(iulog,*) 'GoveqnBaseAllocVarsFromOtherGEs: Not allocated'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    this%nvars_needed_from_other_goveqns = 0

    if (associated(this%var_ids_needed_from_other_goveqns)) deallocate(this%var_ids_needed_from_other_goveqns )
    if (associated(this%ids_of_other_goveqns             )) deallocate(this%ids_of_other_goveqns              )
    if (associated(this%is_bc_auxvar_type                )) deallocate(this%is_bc_auxvar_type                 )
    if (associated(this%bc_auxvar_offset                 )) deallocate(this%bc_auxvar_offset                  )
    if (associated(this%bc_auxvar_idx                    )) deallocate(this%bc_auxvar_idx                     )
    if (associated(this%bc_auxvar_idx_of_other_goveqn    )) nullify(this%bc_auxvar_idx_of_other_goveqn        )
    if (associated(this%bc_auxvar_ncells                 )) nullify(this%bc_auxvar_ncells                     )

    nullify(this%var_ids_needed_from_other_goveqns )
    nullify(this%ids_of_other_goveqns              )
    nullify(this%is_bc_auxvar_type                 )
    nullify(this%bc_auxvar_offset                  )
    nullify(this%bc_auxvar_idx                     )
    nullify(this%bc_auxvar_idx_of_other_goveqn     )
    nullify(this%bc_auxvar_ncells                  )

  end subroutine GoveqnBaseDeallocVarsFromOtherGEs


  !------------------------------------------------------------------------
  subroutine GoveqnBaseIFunction(this, U, Udot, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS IFunction
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: U
    Vec                     :: Udot
    Vec                     :: F
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseIFunction must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseIFunction

  !------------------------------------------------------------------------
  subroutine GoveqnBaseIJacobian(this, U, Udot, shift, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS IJacobian
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: U
    Vec                     :: Udot
    PetscReal               :: shift
    Mat                     :: A
    Mat                     :: B
    PetscErrorCode                :: ierr

    write(iulog,*)'GoveqnBaseJFunction must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseIJacobian

  !------------------------------------------------------------------------
  subroutine GoveqnBaseIJacobianOffDiag(this, U_1, Udot_1, U_2, Udot_2, &
       shift, A, B, id_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS IJacobian corresponding to off-diagonal
    ! matrix coupling between two GEs
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: U_1
    Vec                     :: Udot_1
    Vec                     :: U_2
    Vec                     :: Udot_2
    PetscReal               :: shift
    Mat                     :: A
    Mat                     :: B
    PetscInt                :: id_of_other_goveq
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseIJacobianOffDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseIJacobianOffDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseJacobianOffDiag(this, X_1, X_2, A, B, &
       id_of_other_goveq, &
       list_id_of_other_goveq, &
       ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Jacobian corresponding to off-diagonal
    ! matrix coupling between two GEs
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X_1
    Vec                     :: X_2
    Mat                     :: A
    Mat                     :: B
    PetscInt                :: id_of_other_goveq
    PetscInt                :: list_id_of_other_goveq
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseJacobianOffDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseJacobianOffDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseResidual(this, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Function evaluation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X
    Vec                     :: F
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseResidual must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseResidual

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeRHS must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeRHS

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc KSP Operator matrix
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Mat                     :: A
    Mat                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeOperatorsDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOperatorsDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, list_id_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc KSP Operator matrix
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Mat                     :: A
    Mat                     :: B
    PetscInt                :: itype_of_other_goveq
    PetscInt                :: list_id_of_other_goveq
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeOperatorsOffDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOperatorsOffDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseSetDtime(this, dtime)
    !
    ! !DESCRIPTION:
    ! Sets time step
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    PetscReal               :: dtime

    this%dtime = dtime

  end subroutine GoveqnBaseSetDtime


  !------------------------------------------------------------------------
  subroutine GoveqnBaseJacobian(this, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Jacobian
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X
    Mat                     :: A
    Mat                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseJacobian must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseJacobian

  !------------------------------------------------------------------------
  subroutine GoveqnBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType, only : ConditionListPrintInfo
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*)'    ---------------------------------------------------------'
    write(iulog,*)'    Goveqn_name       : ',trim(this%name)
    write(iulog,*)'    Goveqn_id         : ',this%id
    write(iulog,*)'    Goveqn_mesh_itype : ',this%mesh_itype
    write(iulog,*)'    Num_vars_needed   : ',this%nvars_needed_from_other_goveqns
    write(iulog,*)' '

    write(iulog,*)'    BC'
    call ConditionListPrintInfo(this%boundary_conditions)
    write(iulog,*)'    SS'
    call ConditionListPrintInfo(this%source_sinks)

  end subroutine GoveqnBasePrintInfo

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine to update auxiliary variables associated with the GE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBaseUpdateAuxVars: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateAuxVarsIntrn(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine to update auxiliary variables for internal cells
    ! associated with the GE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBaseUpdateAuxVarsIntrn: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseUpdateAuxVarsIntrn

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateAuxVarsBC(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine to update auxiliary variables for boundary condition
    ! associated with the GE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBaseUpdateAuxVarsBC: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseUpdateAuxVarsBC

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateAuxVarsSS(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine to update auxiliary variables for source-sink
    ! associated with the GE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBaseUpdateAuxVarsSS: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseUpdateAuxVarsSS

  !------------------------------------------------------------------------
  subroutine GoveqnBasePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PreSolve
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBasePreSolve: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBasePreSolve

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetNumConditions(this, cond_type, &
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
    class(goveqn_base_type) :: this
    PetscInt                                 :: cond_type
    PetscInt                                 :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
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
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
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

  end subroutine GoveqnBaseGetNumConditions

  !------------------------------------------------------------------------
  subroutine GoveqnBaseAddCondition(this, ss_or_bc_type, name, unit,    &
       cond_type, region_type, id_of_other_goveq, itype_of_other_goveq, &
       num_other_goveqs, id_of_other_goveqs, itype_of_other_goveqs,     &
       icoupling_of_other_goveqns, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds boundary/source-sink condition to governing equation
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionListAddCondition
    use ConditionType             , only : ConditionNew
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)                     :: this
    PetscInt                                    :: ss_or_bc_type
    character(len =*)                           :: name
    character(len =*)                           :: unit
    PetscInt                                    :: cond_type
    PetscInt                                    :: region_type
    PetscInt, optional                          :: id_of_other_goveq
    PetscInt, optional                          :: itype_of_other_goveq
    PetscInt, optional                          :: num_other_goveqs
    PetscInt, optional, pointer                 :: id_of_other_goveqs(:)
    PetscInt, optional, pointer                 :: itype_of_other_goveqs(:)
    PetscBool, optional, pointer                :: icoupling_of_other_goveqns(:)
    type(connection_set_type),pointer, optional :: conn_set
    !
    type(condition_type), pointer :: cond

    cond => ConditionNew()

    cond%name         = trim(name)
    cond%units        = trim(unit)
    cond%itype        = cond_type
    cond%region_itype = region_type

    allocate(cond%conn_set)
    if (.not. present(conn_set)) then
       call MeshCreateConnectionSet(this%mesh, cond%region_itype, cond%conn_set, cond%ncells)
    else
       cond%conn_set => conn_set
       cond%ncells = conn_set%num_connections
       nullify(conn_set)
    endif

    allocate(cond%value(cond%ncells))
    cond%value(:) = 0.d0

    select case (ss_or_bc_type)
    case (COND_BC)
       if (cond_type == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          if (.not. present(num_other_goveqs)) then

             if (.not.present(id_of_other_goveq) ) then
                write(iulog,*) 'BC type = COND_DIRICHLET_FRM_OTR_GOVEQ ' // &
                     ' but id_of_other_goveq is absent'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             if (.not.present(itype_of_other_goveq) ) then
                write(iulog,*) 'BC type = COND_DIRICHLET_FRM_OTR_GOVEQ ' // &
                     ' but itype_of_other_goveq is absent'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             cond%num_other_goveqs = 1

             allocate(cond%list_id_of_other_goveqs(                 cond%num_other_goveqs))
             allocate(cond%itype_of_other_goveqs(                   cond%num_other_goveqs))
             allocate(cond%swap_order_of_other_goveqs(              cond%num_other_goveqs))
             allocate(cond%coupled_via_intauxvar_with_other_goveqns(cond%num_other_goveqs))

             cond%list_id_of_other_goveqs(1)                  = id_of_other_goveq
             cond%itype_of_other_goveqs(1)                    = itype_of_other_goveq
             cond%coupled_via_intauxvar_with_other_goveqns(1) = PETSC_FALSE
          else

             if (.not.present(id_of_other_goveqs) ) then
                write(iulog,*) 'BC type = COND_DIRICHLET_FRM_OTR_GOVEQ ' // &
                     ' but id_of_other_goveqs is absent'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             else
                if (size(id_of_other_goveqs) /= num_other_goveqs) then
                   write(iulog,*) 'size(id_of_other_goveqs) /= num_other_goveqs'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
             endif

             if (.not.present(itype_of_other_goveqs) ) then
                write(iulog,*) 'BC type = COND_DIRICHLET_FRM_OTR_GOVEQ ' // &
                     ' but itype_of_other_goveqs is absent'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             else
                if (size(itype_of_other_goveqs) /= num_other_goveqs) then
                   write(iulog,*) 'size(itype_of_other_goveqs) /= num_other_goveqs'
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
             endif

             cond%num_other_goveqs = num_other_goveqs

             allocate(cond%list_id_of_other_goveqs(                 cond%num_other_goveqs))
             allocate(cond%itype_of_other_goveqs(                   cond%num_other_goveqs))
             allocate(cond%swap_order_of_other_goveqs(              cond%num_other_goveqs))
             allocate(cond%coupled_via_intauxvar_with_other_goveqns(cond%num_other_goveqs))

             cond%list_id_of_other_goveqs(:)                  = id_of_other_goveqs(:)
             cond%itype_of_other_goveqs(:)                    = itype_of_other_goveqs(:)
             cond%coupled_via_intauxvar_with_other_goveqns(:) = icoupling_of_other_goveqns(:)
          endif
       endif
       call ConditionListAddCondition(this%boundary_conditions, cond)

    case (COND_SS)
       call ConditionListAddCondition(this%source_sinks, cond)

    case default
       write(iulog,*) 'Unknown condition type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine GoveqnBaseAddCondition

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateConditionConnSet(this, icond, &
       ss_or_bc_type, nconn,  conn_id_up, conn_id_dn, &
       conn_dist_up, conn_dist_dn,  conn_unitvec, conn_area)
    !
    ! !DESCRIPTION:
    ! Updates connection of a source-sink/boundary condition
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionListAddCondition
    use ConditionType             , only : ConditionNew
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use ConnectionSetType         , only : ConnectionSetNew
    use ConnectionSetType         , only : ConnectionSetDestroy
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)       :: this
    PetscInt                      :: icond
    PetscInt                      :: ss_or_bc_type
    PetscInt                      :: nconn
    PetscInt, pointer             :: conn_id_up(:)   !
    PetscInt, pointer             :: conn_id_dn(:)   !
    PetscReal, pointer            :: conn_dist_up(:) !
    PetscReal, pointer            :: conn_dist_dn(:) !
    PetscReal, pointer            :: conn_area(:)    !
    PetscReal, pointer            :: conn_type(:)    !
    PetscReal, pointer            :: conn_unitvec(:,:)
    !
    type(condition_type), pointer :: cond
    PetscInt                      :: iconn, ii
    PetscInt                      :: num_conditions


    num_conditions = 0
    select case (ss_or_bc_type)
    case (COND_BC)

       cond => this%boundary_conditions%first
       do
          if (.not.associated(cond)) exit
          num_conditions = num_conditions + 1
          cond => cond%next
       enddo

       if (icond > num_conditions) then
          write(iulog,*),'Could not find the BC whose ' // &
               'connection set needs to be updated.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       cond => this%boundary_conditions%first
       do ii = 1,icond-1
          cond => cond%next
       enddo

    case (COND_SS)
       cond => this%source_sinks%first
       do
          if (.not.associated(cond)) exit
          num_conditions = num_conditions + 1
          cond => cond%next
       enddo

       if (icond > num_conditions) then
          write(iulog,*),'Could not find the SS whose ' // &
               'connection set needs to be updated.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       cond => this%source_sinks%first
       do ii = 1,icond-1
          cond => cond%next
       enddo

    case default
       write(iulog,*) 'Unknown condition type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

    call ConnectionSetDestroy(cond%conn_set)

    allocate(cond%conn_set)
    cond%conn_set => ConnectionSetNew(nconn)

    do iconn = 1, nconn
       cond%conn_set%id_up(iconn)               = conn_id_up(iconn)
       cond%conn_set%id_dn(iconn)               = conn_id_dn(iconn)
       cond%conn_set%area(iconn)                = conn_area(iconn)
       cond%conn_set%dist_up(iconn)             = conn_dist_up(iconn)
       cond%conn_set%dist_dn(iconn)             = conn_dist_dn(iconn)
       cond%conn_set%dist_unitvec(iconn)%arr(1) = conn_unitvec(iconn,1)
       cond%conn_set%dist_unitvec(iconn)%arr(2) = conn_unitvec(iconn,2)
       cond%conn_set%dist_unitvec(iconn)%arr(3) = conn_unitvec(iconn,3)
    enddo

  end subroutine GoveqnBaseUpdateConditionConnSet

#endif

end module GoverningEquationBaseType
