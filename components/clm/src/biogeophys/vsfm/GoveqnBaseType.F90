#ifdef USE_PETSC_LIB


module GoverningEquationBaseType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

  ! !USES:
  use clm_varctl         , only : iulog
  use abortutils         , only : endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg
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

     PetscInt, pointer               :: ids_of_other_goveqns(:)              !
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
     procedure, public :: ComputeOperators        => GoveqnBaseComputeOperators
     procedure, public :: SetDtime                => GoveqnBaseSetDtime
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
  subroutine GoveqnBaseComputeRHS(this, R, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: R
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeRHS must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeRHS

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOperators(this, A, B, ierr)
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

    write(iulog,*)'GoveqnBaseComputeOperators must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOperators

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

end module GoverningEquationBaseType
#endif
