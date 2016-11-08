#ifdef USE_PETSC_LIB


module SystemOfEquationsBaseType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Base object for a system-of-equations
  !-----------------------------------------------------------------------

  ! !USES:
  use abortutils                , only : endrun
  use clm_varctl                , only : iulog
  use shr_log_mod               , only : errMsg => shr_log_errMsg
  use MeshType                  , only : mesh_type
  use GoverningEquationBaseType , only : goveqn_base_type
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
#include "finclude/petscsnes.h"
#include "finclude/petscsnes.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscksp.h"
#include "finclude/petscksp.h90"

  type, public :: sysofeqns_base_type
     character(len =256)             :: name                         ! name for system-of-equations (SoE)
     PetscInt                        :: itype                        ! identifier for SoE

     class(goveqn_base_type),pointer :: goveqns                      ! pointer to governing equations within SoE
     PetscInt                        :: ngoveqns                     ! number of governing equations within SoE


     PetscReal                       :: time                         ! [sec]
     PetscReal                       :: dtime                        ! [sec]

     PetscInt                        :: cumulative_newton_iterations ! Total number of Newton iterations
     PetscInt                        :: cumulative_linear_iterations ! Total number of Linear iterations

     PetscInt                        :: solver_type                  ! type of PETSc equation being solved (KSP, SNES, TS)
     DM                              :: dm                           ! PETSc DM
     TS                              :: ts                           ! PETSc TS
     SNES                            :: snes                         ! PETSc SNES
     KSP                             :: ksp                          ! PETSc KSP

     Vec                             :: soln                         ! solution at current iteration + time step
     Vec                             :: soln_prev                    ! solution vector at previous time step
     Vec                             :: soln_prev_clm                ! solution vector at previous CLM time step
     Vec                             :: rhs                          ! used if SoE is a PETSc TS
     Mat                             :: jac                          ! used if SoE is a PETSc TS/SNES
     Mat                             :: Amat                         ! used if SoE is a PETSc KSP

   contains
     procedure, public :: Init                  => SOEBaseInit
     procedure, public :: Clean                 => SOEBaseClean
     procedure, public :: IFunction             => SOEBaseIFunction
     procedure, public :: IJacobian             => SOEBaseIJacobian
     procedure, public :: Residual              => SOEBaseResidual
     procedure, public :: Jacobian              => SOEBaseJacobian
     procedure, public :: ComputeRHS            => SOEComputeRHS
     procedure, public :: ComputeOperators      => SOEComputeOperators
     procedure, public :: StepDT                => SOEBaseStepDT
     procedure, public :: PreSolve              => SOEBasePreSolve
     procedure, public :: PostSolve             => SOEBasePostSolve
     procedure, public :: PreStepDT             => SOEBasePreStepDT
     procedure, public :: PostStepDT            => SOEBasePostStepDT
     procedure, public :: PrintInfo             => SOEBasePrintInfo
     procedure, public :: SetPointerToIthGovEqn => SOEBaseSetPointerToIthGovEqn
     procedure, public :: SetDtime              => SOEBaseSetDtime
     procedure, public :: SetDataFromCLM        => SOEBaseSetDataFromCLM
     procedure, public :: GetDataForCLM         => SOEBaseGetDataForCLM
  end type sysofeqns_base_type

  public :: SOEBaseInit
  public :: SOESetMeshesOfGoveqns

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine SOEBaseInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize a SoE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    this%name                         = ""
    this%itype                        = 0
    this%ngoveqns                     = 0

    this%time                         = 0.d0
    this%dtime                        = 0.d0

    this%cumulative_newton_iterations = 0
    this%cumulative_linear_iterations = 0

    this%solver_type                  = 0
    this%dm                           = 0
    this%ts                           = 0
    this%snes                         = 0
    this%ksp                          = 0

    this%soln                         = 0
    this%soln_prev                    = 0
    this%soln_prev_clm                = 0
    this%rhs                          = 0
    this%jac                          = 0
    this%Amat                         = 0

    nullify(this%goveqns)

  end subroutine SOEBaseInit

  !------------------------------------------------------------------------
  subroutine SOEBaseIFunction(this, ts, t, U, Udot, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc TS.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    TS                         :: ts
    PetscReal                  :: t
    Vec                        :: U
    Vec                        :: Udot
    Vec                        :: F
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseIFunction must be extended')

  end subroutine SOEBaseIFunction

  !------------------------------------------------------------------------
  subroutine SOEBaseIJacobian(this, ts, t, U, Udot, shift, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc TS.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    TS                         :: ts
    PetscReal                  :: t
    Vec                        :: U
    Vec                        :: Udot
    PetscReal                  :: shift
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseIJacobian must be extended')

  end subroutine SOEBaseIJacobian

  !------------------------------------------------------------------------
  subroutine SOEBaseResidual(this,snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc SNES.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    SNES                       :: snes
    PetscReal                  :: t
    Vec                        :: X
    Vec                        :: F
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseResidual must be extended')

  end subroutine SOEBaseResidual

  !------------------------------------------------------------------------
  subroutine SOEBaseJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc SNES.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    SNES                       :: snes
    Vec                        :: X
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseJacobian must be extended')

  end subroutine SOEBaseJacobian

  !------------------------------------------------------------------------
  subroutine SOEComputeRHS(this, ksp, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc TS.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    KSP                        :: ksp
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SOEComputeRHS: '//&
         'SOEComputeRHS must be extended')

  end subroutine SOEComputeRHS

  !------------------------------------------------------------------------
  subroutine SOEComputeOperators(this, ksp, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc KSP.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    KSP                        :: ksp
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SOEComputeOperators: '//&
         'SOEComputeOperators must be extended')

  end subroutine SOEComputeOperators

  !------------------------------------------------------------------------
  subroutine SOEBasePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations before calling
    ! PETSc solover.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBasePreSolve must be extended')

  end subroutine SOEBasePreSolve

  !------------------------------------------------------------------------
  subroutine SOEBasePreStepDT(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations before calling
    ! StepDT.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SOEBasePreStepDT: '//&
         'SOEBasePreStepDT must be extended')

  end subroutine SOEBasePreStepDT

  !------------------------------------------------------------------------
  subroutine SOEBasePostStepDT(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations post StepDT.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SOEBasePostStepDT: '//&
         'SOEBasePostStepDT must be extended')

  end subroutine SOEBasePostStepDT

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT(this, dt, converged, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE by calling appropriate subroutine dependning on the choice
    ! of PETSc solver set for the SoE.
    !
    ! !USES
    use MultiPhysicsProbConstants, only : PETSC_TS
    use MultiPhysicsProbConstants, only : PETSC_SNES
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscReal                  :: dt
    PetscBool,intent(out)      :: converged
    PetscErrorCode             :: ierr

    select case(this%solver_type)
    case (PETSC_TS)
       call SOEBaseStepDT_TS(this, dt, ierr)
    case (PETSC_SNES)
       call SOEBaseStepDT_SNES(this, dt, converged, ierr)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SOEBaseStepDT

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT_TS(soe, dt, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE via PETSc TS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: soe
    PetscReal                  :: dt
    PetscErrorCode             :: ierr
    !
    ! !LOCAL VARIABLES:
    SNESConvergedReason        :: snes_reason

    call TSSetTime(soe%ts, 0.d0, ierr); CHKERRQ(ierr)
    call TSSetDuration(soe%ts, 100000, dt, ierr); CHKERRQ(ierr)
    call TSsetInitialTimeStep(soe%ts, 0.0d0, 3600.0d0, ierr); CHKERRQ(ierr)

    call TSSetFromOptions(soe%ts, ierr); CHKERRQ(ierr)

    call TSSolve(soe%ts, soe%soln, ierr); CHKERRQ(ierr)

  end subroutine SOEBaseStepDT_TS

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT_SNES(soe, dt, converged, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE via PETSc SNES
    !
    ! !USES
    use spmdMod          , only : masterproc, iam
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : vsfm_use_dynamic_linesearch
    !
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscsnes.h"
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: soe
    PetscErrorCode             :: ierr
    PetscBool,intent(out)      :: converged
    PetscReal                  :: dt
    !
    ! !LOCAL VARIABLES:
    SNESConvergedReason        :: snes_reason
    PetscInt                   :: num_newton_iterations
    PetscInt                   :: num_linear_iterations
    PetscInt                   :: num_time_cuts
    PetscInt, parameter        :: max_num_time_cuts = 20
    PetscReal                  :: target_time
    PetscReal                  :: dt_iter
    PetscInt                   :: linesearch_iter
    SNESLineSearch             :: linesearch
    PetscInt                   :: max_linesearch_iter
    PetscInt, pointer          :: linesearch_iter_types(:)
    PetscInt, parameter        :: LS_BASIC = 1
    PetscInt, parameter        :: LS_BT    = 2
    PetscInt, parameter        :: LS_L2    = 3
    PetscInt, parameter        :: LS_CP    = 4
    PetscBool                  :: is_default_linesearch_basic
    PetscBool                  :: is_default_linesearch_bt
    PetscBool                  :: is_default_linesearch_l2
    PetscBool                  :: is_default_linesearch_cp
    character(len=32)          :: linesearch_name

    ! initialize
    linesearch_iter = 0
    num_time_cuts  = 0
    soe%time       = 0.d0
    target_time    = dt
    dt_iter        = dt

    ! Determine the default linesearch option
    call SNESGetLineSearch(soe%snes, linesearch, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHBASIC, is_default_linesearch_basic, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHBT, is_default_linesearch_bt, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHL2, is_default_linesearch_l2, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHCP, is_default_linesearch_cp, ierr); CHKERRQ(ierr)

    ! Check if the default linesearch option was a known option
    if ( (.not. is_default_linesearch_basic) .and. &
         (.not. is_default_linesearch_bt   ) .and. &
         (.not. is_default_linesearch_l2   ) .and. &
         (.not. is_default_linesearch_cp   ) ) then
       write(iulog,*) 'Unknown default linesearch option'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Create a list of linesearch options that will be tried before
    ! cutting the timestep
    if (is_default_linesearch_bt .or. is_default_linesearch_l2) then
       max_linesearch_iter = 2
       allocate(linesearch_iter_types(max_linesearch_iter))

       if (is_default_linesearch_bt) then
          linesearch_iter_types(1) = LS_BT
          linesearch_iter_types(2) = LS_L2
       else
          linesearch_iter_types(1) = LS_L2
          linesearch_iter_types(2) = LS_BT
       endif
    else
       max_linesearch_iter = 3
       allocate(linesearch_iter_types(max_linesearch_iter))

       if (is_default_linesearch_basic) linesearch_iter_types(1) = LS_BASIC
       if (is_default_linesearch_cp   ) linesearch_iter_types(1) = LS_CP
       linesearch_iter_types(2) = LS_L2
       linesearch_iter_types(3) = LS_BT
    endif

    do

       ! Set timestep
       call soe%SetDtime(dt_iter)

       ! Do any pre-solve operations
       call soe%PreSolve()

       select case (linesearch_iter_types(linesearch_iter+1))
       case (LS_BASIC)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_basic'
       case (LS_BT)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHBT, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_bt'
       case (LS_L2)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHL2, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_l2'
       case (LS_CP)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHCP, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_cp'
       end select

       ! Solve the nonlinear equation
       call SNESSolve(soe%snes, PETSC_NULL_OBJECT, soe%soln, ierr); CHKERRQ(ierr)

       ! Get reason why SNES iteration stopped
       call SNESGetConvergedReason(soe%snes, snes_reason, ierr); CHKERRQ(ierr)

       ! Did SNES converge?
       if (snes_reason < 0) then

          linesearch_iter = linesearch_iter + 1

          if (vsfm_use_dynamic_linesearch .and. linesearch_iter < max_linesearch_iter) then
             ! Let's try another linesearch
             write(iulog,*),'On proc ', iam, ' time_step = ', get_nstep(), &
                  linesearch_name // ' unsuccessful. Trying another one.'
             call VecCopy(soe%soln_prev, soe%soln, ierr); CHKERRQ(ierr)
          else
             ! All linesearch options failed, reset the linesearch iteration
             ! counter
             linesearch_iter = 0

             ! SNES diverged, so let's cut the timestep and try again.
             num_time_cuts = num_time_cuts + 1
             dt_iter = 0.5d0*dt_iter
             write(iulog,*),'On proc ', iam, ' time_step = ', get_nstep(), &
                  'snes_reason = ',snes_reason,' cutting dt to ',dt_iter
          endif

          call VecCopy(soe%soln_prev, soe%soln, ierr); CHKERRQ(ierr)
       else
          ! SNES converged.
          converged = PETSC_TRUE
          soe%time = soe%time + dt_iter

          call SNESGetIterationNumber(soe%snes,       &
               num_newton_iterations, ierr); CHKERRQ(ierr)
          call SNESGetLinearSolveIterations(soe%snes, &
               num_linear_iterations, ierr); CHKERRQ(ierr)

          soe%cumulative_newton_iterations = soe%cumulative_newton_iterations + &
               num_newton_iterations
          soe%cumulative_linear_iterations = soe%cumulative_linear_iterations + &
               num_linear_iterations

          ! Do any post-solve operations
          call soe%PostSolve()
       endif

       ! Do number of time cuts exceed maximum allowable number of cuts?
       if (num_time_cuts > max_num_time_cuts) then
          converged = PETSC_FALSE
          return
       endif

       if (soe%time >= target_time) exit
    enddo

    ! Set the linsearch type to be the default setting
    select case (linesearch_iter_types(1))
    case (LS_BASIC)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr); CHKERRQ(ierr)
    case (LS_BT)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHBT, ierr); CHKERRQ(ierr)
    case (LS_L2)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHL2, ierr); CHKERRQ(ierr)
    case (LS_CP)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHCP, ierr); CHKERRQ(ierr)
    end select
    deallocate(linesearch_iter_types)


  end subroutine SOEBaseStepDT_SNES

  !------------------------------------------------------------------------
  subroutine SOEBasePostSolve(this)
    !
    ! !DESCRIPTION:
    ! Subroutine that makes the copy of current solution after a successful
    ! PETSc solove.
    ! This subroutines may be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    !
    ! !LOCAL VARIABLES:
    PetscErrorCode             :: ierr

    call VecCopy(this%soln, this%soln_prev,ierr); CHKERRQ(ierr)

  end subroutine SOEBasePostSolve

  !------------------------------------------------------------------------
  subroutine SOEBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Displays information about SoE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
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
    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  end subroutine SOEBasePrintInfo

  !------------------------------------------------------------------------
  subroutine SOESetMeshesOfGoveqns(soe, meshes, nmesh)
    !
    ! !DESCRIPTION:
    ! Match the meshes in `meshes` with the governing equations in `soe`.
    !
    ! !USES
    use MeshType                     , only : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type), intent(inout)  :: soe
    class(mesh_type), pointer, intent(in)      :: meshes(:)
    PetscInt, intent(in)                       :: nmesh
    !
    ! !LOCAL VARIABLES:
    PetscInt                                   :: imesh
    PetscInt                                   :: mesh_itype
    PetscBool                                  :: mesh_found
    class(mesh_type), pointer                  :: cur_mesh
    class(goveqn_base_type),pointer            :: cur_goveqn

    cur_goveqn => soe%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       mesh_itype = cur_goveqn%mesh_itype

       mesh_found = PETSC_FALSE
       do imesh = 1, nmesh
          cur_mesh => meshes(imesh)
          if (mesh_itype == cur_mesh%itype) then
             cur_goveqn%mesh => cur_mesh
             mesh_found = PETSC_TRUE
             exit
          endif
       enddo

       if (.not.mesh_found) then
          call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
               'Mesh not found for Governing equation within the list')
       endif

       cur_goveqn => cur_goveqn%next
    enddo

  end subroutine SOESetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine SOEBaseSetPointerToIthGovEqn(this, goveqn_id, goveqn_ptr)
    !
    ! !DESCRIPTION:
    ! Returns pointer to the i-th governing equation present with SoE
    !
    ! !USES
    use MeshType , only : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type)      :: this
    PetscInt                        :: goveqn_id
    class(goveqn_base_type),pointer :: goveqn_ptr
    !
    ! !LOCAL VARIABLES
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: sum_goveqn
    PetscBool                       :: found

    sum_goveqn = 0
    found      = PETSC_FALSE

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       sum_goveqn = sum_goveqn + 1
       if (sum_goveqn == goveqn_id) then
          found = PETSC_TRUE
          goveqn_ptr => cur_goveq
          exit
       endif

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*) 'In SOEBaseSetPointerToIthGovEqn: ith-goveqn not found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine SOEBaseSetPointerToIthGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseSetDtime(this, dtime)
    !
    ! !DESCRIPTION:
    ! Sets timestep for SoE and all governing equations present within SoE.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscReal                  :: dtime
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer   :: cur_goveqn

    this%dtime = dtime

    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       call cur_goveqn%SetDtime(dtime)
       cur_goveqn => cur_goveqn%next
    enddo

  end subroutine SOEBaseSetDtime

  !------------------------------------------------------------------------
  subroutine SOEBaseSetDataFromCLM(this, soe_auxvar_type, var_type, soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that sets data from CLM before proceeding to solve the SoE.
    ! E.g. Getting infilitration source from CLM for VSFM.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt, intent(in)       :: var_type
    PetscInt                   :: soe_auxvar_type
    PetscInt                   :: soe_auxvar_id
    PetscReal                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseSetDataFromCLM must be extended')

  end subroutine SOEBaseSetDataFromCLM

  !------------------------------------------------------------------------
  subroutine SOEBaseGetDataForCLM(this, soe_auxvar_type, var_type, soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that extracts data from SoE's data structure for CLM after
    ! a successfull SoE solve. E.g. Getting updated soil moisture value from VSFM
    ! for CLM.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt, intent(in)       :: var_type
    PetscInt                   :: soe_auxvar_type
    PetscInt                   :: soe_auxvar_id
    PetscReal                  :: data_1d(:)
    PetscInt                   :: nsize
    !
    ! !LOCAL VARIABLES:

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseGetDataForCLM must be extended')

  end subroutine SOEBaseGetDataForCLM

  !------------------------------------------------------------------------
  subroutine SOEBaseClean(this)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    !call this%goveqns%Clean

    nullify(this%goveqns)

  end subroutine SOEBaseClean

end module SystemOfEquationsBaseType

#endif
