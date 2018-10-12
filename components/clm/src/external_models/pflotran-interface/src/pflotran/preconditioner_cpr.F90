module CPR_Preconditioner_module
! Implements a CPR preconditioner using the PCSHELL
! funcitonality of PETSC.
! Daniel Stone and Sebastien Loisel
#include <petsc/finclude/petscsys.h>
#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscviewer.h"
  use petscmat
  use petscpc
  use petscvec
  use PFLOTRAN_Constants_module
  use Option_module
  implicit none

  private

  type, public :: cpr_pc_type 
    Mat :: A, Ap
    ! The two stage CPR preconditioner calls two other preconditioners:
    PC :: T1_PC        ! T1 will be a SHELL pc that extracts the pressure system residual from the
                    ! input residual, then approximates the inverse of Ap acting on that residual
                    ! using AMG.
    KSP :: T1_KSP   ! Usually PCONLY, the PC for this KSP is hypre/boomeramg. This is called
                    ! by T1 (above) as the AMG part.
    PC :: T2_PC        ! A regular PC, usually BJACOBI
    Vec :: T1r,r2, s, z, factors1vec, factors2vec
    PetscInt ::  t2fillin, timescalled, asmoverlap, exrow_offset
    PetscBool :: firstT1Call, firstT2call, asmfactorinplace, &
                 t2shiftinblocks, zeroing, useGAMG, mv_output, t2_zeroing, &
                 skip_T1, amg_report, amg_manual
    character(len=MAXWORDLENGTH) :: T1_type, T2_type, extract_type
    ! following are buffers/workers for the pressure system extraction.
    ! 1d arrays:
    PetscReal, dimension(:), allocatable :: vals, insert_vals
    PetscInt, dimension(:), allocatable :: colIdx, colIdx_keep, insert_colIdx
    ! 2d arrays:
    PetscReal, dimension(:,:), allocatable :: all_vals
    ! point at the option object, needed for error outputs
    type(option_type), pointer :: option
  end type cpr_pc_type 

  ! interfaces need to make the set and get context routines
  ! work
  interface
    subroutine PCShellSetContext (P_in, ctx_in, ierr_in)
      use petscksp
      Import :: cpr_pc_type 
      PC :: P_in
      type(cpr_pc_type) :: ctx_in
      PetscErrorCode :: ierr_in
    end subroutine PCShellSetContext
  end interface

  interface
    subroutine PCShellGetContext (P_in, ctx_in, ierr_in)
      use petscksp
      Import :: cpr_pc_type 
      PC :: P_in
      type(cpr_pc_type), pointer :: ctx_in
      PetscErrorCode :: ierr_in
    end subroutine PCShellGetContext
  end interface

public :: DeallocateWorkersInCPRStash, &
          CPRMake

contains

! ************************************************************************** !
!  CPR Apply routines

subroutine CPRApply(p, r, y,ierr)
  ! To be used as a PCSHELL apply routine

  ! 
  ! Applies the CPR preconditioner to r
  ! (output y)
  !
  ! Author: Sebastien Loisel, Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  ! fixed signature: will approximate
  ! y = inv(A)r using shell preconditioner P:
  PC :: P
  Vec :: r
  Vec :: y
  PetscErrorCode :: ierr

  ! worker vectors:
  Vec :: t1r, r2
  ! PC context holds workers:
  type(cpr_pc_type), pointer :: ctx
  ! misc variables and parms:
  PetscReal:: one
  Mat :: a 

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  one = 1.0
  r2 = ctx%r2
  t1r = ctx%T1r
  a = ctx%a

  if (ctx%skip_T1) then
    call VecZeroEntries(t1r,ierr); CHKERRQ(ierr)
  else
    call PCApply(ctx%T1_PC,r,t1r,ierr); CHKERRQ(ierr)
  endif

  call MatMult(a,t1r,r2,ierr); CHKERRQ(ierr)
  call VecAYPX(r2,-one,r,ierr); CHKERRQ(ierr) ! r1 <- r-r2

  call PCApply(ctx%T2_PC,r2,y,ierr); CHKERRQ(ierr)

  call VecAYPX(y,one,t1r,ierr); CHKERRQ(ierr)

end subroutine CPRApply

! ************************************************************************** !

subroutine CPRT1Apply(p, x, y,ierr)
  ! To be used as a PCSHELL apply routine

  ! 
  ! Applies the T1 part of the CPR preconditioner
  ! to r (ourput y)
  !
  ! Author: Sebastien Loisel, Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 
  
  implicit none

  ! fixed signature: will approximate
  ! y = inv(A)r using shell preconditioner P:
  PC :: p
  Vec :: x
  Vec :: y
  PetscErrorCode :: ierr

  ! PC context holds workers:
  type(cpr_pc_type), pointer :: ctx
  ! some other KSPs and PCs we'll use:
  KSP :: ksp
  PC :: amg_pc
  ! misc workers, etc:
  PetscInt :: b,start,end,k,its
  Vec :: s, z

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  s = ctx%s
  z = ctx%z

  call QIRHS(ctx%factors1Vec, ctx%factors2vec, x, s, ierr)

  ksp = ctx%T1_KSP

  if (ctx%T1_type /= "NONE") then
    call KSPSolve(ksp,s,z,ierr); CHKERRQ(ierr)
    call KSPGetIterationNumber(ksp, its, ierr); CHKERRQ(ierr)
  else
    call KSPGetPC(ksp, amg_pc, ierr);CHKERRQ(ierr)
    call PCApply(amg_pc,s,z,ierr); CHKERRQ(ierr)

    if (ctx%amg_report) then
      call PCView(amg_pc, PETSC_VIEWER_STDOUT_SELF, ierr);CHKERRQ(ierr)
    endif
  endif

  call VecZeroEntries(y,ierr); CHKERRQ(ierr)
  call VecStrideScatter(z,0,y,INSERT_VALUES,ierr); CHKERRQ(ierr)

end subroutine CPRT1Apply

!  end of CPR apply routines
! ************************************************************************** !

! ************************************************************************** !
!  CPR Setup Routines (every time Jaocobian updates) 

subroutine CPRSetup(p,ierr)
  ! to be used as a PCSHELL setup routine

  ! 
  ! sets up the CPR preconditioner, is called
  ! every time the matrix updates
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  ! fixed signature:
  PC :: p
  PetscErrorCode :: ierr

  ! PC context holds workers:
  type(cpr_pc_type), pointer ::ctx

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  call CPRSetupT1(ctx, ierr)
  call CPRSetupT2(ctx, ierr)

end subroutine CPRSetup

! ************************************************************************** !

subroutine CPRSetupT1(ctx,  ierr)
  ! 
  ! set up the T1 part of the CPR preconditioner,
  ! mostly by extracting the new pressure system
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PetscInt :: b, mx
  Mat :: A
  KSP :: ksp
  
  A = ctx%A

  call MatGetBlockSize(A,b,ierr); CHKERRQ(ierr)

  if (ctx%firstT1Call) then ! some final initializations that need knowlege
                            ! of the fully assembled Jacobian
    ! let T1 and its inner ksp solver know what their inner operators are
    ksp = ctx%T1_KSP
    call PCSetOperators(ctx%T1_PC,A,A,ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp,ctx%Ap,ctx%Ap,ierr); CHKERRQ(ierr)


    ! now sparsity pattern of the Jacobian is defined, 
    ! allocate worker arrays that are big enough to 
    ! be used in the extraction routines coming up.
    call MatGetMaxRowCount(A, mx, ierr)
    call AllocateWorkersInCPRStash(ctx, mx, b)

    ctx%firstT1Call = .FALSE.
  end if

  call MatZeroEntries(ctx%Ap, ierr); CHKERRQ(ierr)

  select case(ctx%extract_type)
    case('QIMPES_VARIABLE')
       if (b == 3) then
        ! we have a more efficient version for 3x3 blocks so do that if we can instead
        call MatGetSubQIMPES(A, ctx%Ap, ctx%factors1vec,  ierr, ctx)
      else
        call MatGetSubQIMPES_var(A, ctx%Ap, ctx%factors1vec,  ierr,  b, ctx)
      endif
    case('QIMPES_VARIABLE_FORCE')
      ! force to use the variables block size implementation even for 3x3 blocks
      call MatGetSubQIMPES_var(A, ctx%Ap, ctx%factors1vec,  ierr,   b,  ctx)
    case default
      ctx%option%io_buffer = 'CPRSetupT1, extraction type not defined'
      call printErrMsg(ctx%option)
  end select

end subroutine CPRSetupT1

! ************************************************************************** !

subroutine CPRSetupT2(ctx, ierr)
  ! 
  ! set up the T2 part of the CPR preconditioner,
  ! if anything needs to be done
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  use String_module

  implicit none

  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PC :: T2, pc_inner
  Mat :: A
  PetscInt :: nsub_ksp, first_sub_ksp, i
  KSP, pointer :: sub_ksps(:)

  ! set operator on first call.

  ! also if first call then time is right
  ! to adjust things like fill in, etc. 

  ! if BJACOBI, this is the best place to try
  ! modifying sub-ksps
  if (ctx%firstT2Call) then

  T2 = ctx%T2_PC

    if (StringCompare(ctx%T2_type, 'PCASM')) then
      call PCASMSetOverlap(T2, ctx%asmoverlap, ierr);CHKERRQ(ierr)
    endif

    A = ctx%A
    call PCSetOperators(T2,A,A,ierr); CHKERRQ(ierr)
    call PCSetFromOptions(T2, ierr);CHKERRQ(ierr)
    call PCSetUp(T2, ierr); CHKERRQ(ierr)

    if (StringCompare(ctx%T2_type, 'PCASM')) then

      ! default should be preonly
      call PCASMGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
      !                         ksp array
      ! allocate ksp array now number known:
      allocate(sub_ksps(nsub_ksp))
      ! call again:
      call PCASMGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               sub_ksps, ierr); CHKERRQ(ierr)

      do i = 1,nsub_ksp
        call KSPGetPC(sub_ksps(i), pc_inner);CHKERRQ(ierr)
        if (ctx%t2shiftinblocks) then
          call PCFactorSetShiftType(pc_inner,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
        endif
        call PCFactorSetLevels(pc_inner, ctx%t2fillin, ierr); CHKERRQ(ierr)
        if (ctx%asmfactorinplace) then
          call PCFactorSetUseInPlace(pc_inner, PETSC_TRUE,  ierr);CHKERRQ(ierr)
        endif
      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    ! specifically for block jacobi:
    elseif (StringCompare(ctx%T2_type, 'PCBJACOBI')) then

      ! default should be preonly
      call PCBJacobiGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
      !                         ksp array
      ! allocate ksp array now number known:
      allocate(sub_ksps(nsub_ksp))
      ! call again:
      call PCBJacobiGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               sub_ksps, ierr); CHKERRQ(ierr)

      do i = 1,nsub_ksp
        ! default should be ilu with 0 fill
        call KSPGetPC(sub_ksps(i), pc_inner, ierr); CHKERRQ(ierr)

        ! MUST DO THIS: the pc is otherwise not setup by this
        ! point and will not accept changes like factor levels.
        call PCSetType(pc_inner, PCILU, ierr); CHKERRQ(ierr)
        call PCSetUp(PC_inner, ierr); CHKERRQ(ierr)

        call PCFactorSetLevels(pc_inner, ctx%t2fillin, ierr); CHKERRQ(ierr)

      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    endif

    ctx%firstT2call = PETSC_FALSE
  endif

end subroutine CPRSetupT2

!  end of CPR Setup Routines 
! ************************************************************************** !

! ************************************************************************** !
!  CPR Creation Routines  

subroutine CPRMake(p, ctx, c, ierr, option)
  ! 
  ! make the CPR preconditioner
  ! create all necessary sub PC/KSP objects
  ! and set up as much as can be done at this
  ! point 
  !
  ! Author: Sebastien Loisel,  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 
  use String_module

  implicit none

  PC :: p
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx
  MPI_Comm :: c
  type(option_type), target :: option

  ctx%option => option

  call PCSetType(p,PCSHELL,ierr); CHKERRQ(ierr)
  call PCShellSetApply(p,CPRapply,ierr); CHKERRQ(ierr)
  call PCShellSetSetUp(p,CPRSetup,ierr); CHKERRQ(ierr)

  call PCShellSetContext(p, ctx, ierr); CHKERRQ(ierr)

#if !defined(PETSC_HAVE_LIBHYPRE)
  if (.not. ctx%useGAMG .or. &
      (StringCompare(ctx%T2_type,'SAILS') .or. &
       StringCompare(ctx%T2_type,'PILUT') .or. &
       StringCompare(ctx%T2_type,'EUCLID'))) then
    option%io_buffer = 'CPR solver settings require that PETSc be &
      &configured with hypre (--download-hypre=yes).'
    call printErrMsg(option)
  endif
#endif

  call CPRCreateT1(c, ctx,  ierr); CHKERRQ(ierr)
  call CPRCreateT2(c, ctx,  ierr); CHKERRQ(ierr)

end subroutine CPRMake

! ************************************************************************** !

subroutine CPRCreateT1(c,  ctx,   ierr)
  ! 
  ! This will perform bare minimum
  ! creation of objects for T1 that do not
  ! need an existing system (i.e. matrix A).
  ! The system dependent setup must be repeated
  ! every step and therefore is seperated out.
  !
  ! Author: Sebastien Loisel,  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  MPI_Comm :: c
  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  KSP :: ksp
  PC :: prec

  ! nice default options for boomeramg:
  if (.NOT. ctx%amg_manual) then 
    call SetCPRDefaults(ierr)
  endif

  call KSPCreate(c,ksp,ierr); CHKERRQ(ierr)

  select case(ctx%T1_type)
    case('RICHARDSON')
      call KSPSetType(ksp,KSPRICHARDSON,ierr); CHKERRQ(ierr)
    case('FGMRES')
      call KSPSetType(ksp,KSPFGMRES,ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0d-3, 1.d0-3, PETSC_DEFAULT_REAL, &
                            PETSC_DEFAULT_INTEGER, ierr);CHKERRQ(ierr) 
    case('GMRES')
      call KSPSetType(ksp,KSPGMRES,ierr); CHKERRQ(ierr)
    case default
      ! really we will just skip over the ksp entirely
      ! in this case, but for completeness..
      call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
  end select

  call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)

  if (ctx%useGAMG) then
    call PCSetType(prec,PCGAMG,ierr); CHKERRQ(ierr)
  else
    call PCSetType(prec,PCHYPRE,ierr); CHKERRQ(ierr)
#if defined(PETSC_HAVE_LIBHYPRE)
    call PCHYPRESetType(prec,"boomeramg",ierr); CHKERRQ(ierr)
#endif
  endif

  call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
  call PetscObjectSetName(ksp,"T1",ierr); CHKERRQ(ierr)

  ctx%T1_KSP = ksp
  call PCCreate(C,ctx%T1_PC,ierr); CHKERRQ(ierr)
  call PCSetType(ctx%T1_PC,PCSHELL,ierr); CHKERRQ(ierr)
  call PCShellSetApply(ctx%T1_PC,CPRT1apply,ierr); CHKERRQ(ierr)

  call PCShellSetContext(ctx%T1_PC, ctx, ierr); CHKERRQ(ierr)

end subroutine CPRCreateT1

! ************************************************************************** !

subroutine CPRCreateT2(c, ctx, ierr)
  ! 
  ! create the T2 preconditioner for the CPR PC
  !
  ! Author: Sebastien Loisel,  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  MPI_Comm :: c
  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PC :: t2

  call PCCreate(c,t2,ierr); CHKERRQ(ierr)

  !select case(ctx%T2_type)
  select case(trim(ctx%T2_type))
    case('SAILS')
      call PCSetType(T2,PCHYPRE,ierr); CHKERRQ(ierr)
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"parasails",ierr); CHKERRQ(ierr)
#endif
    case('PBJ')
      call PCSetType(t2,PCPBJACOBI,ierr); CHKERRQ(ierr)
    case('NONE')
      call PCSetType(t2,PCNONE,ierr); CHKERRQ(ierr)
    case('PCASM')
      call PCSetType(t2,PCASM,ierr); CHKERRQ(ierr)
    case('PCGASM')
      call PCSetType(t2,PCGASM,ierr); CHKERRQ(ierr)
    case('PILUT')
      call PCSetType(t2,PCHYPRE,ierr); CHKERRQ(ierr)
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"pilut",ierr); CHKERRQ(ierr)
#endif
    case('EUCLID')
      call PCSetType(t2,PCHYPRE,ierr); CHKERRQ(ierr)
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"euclid",ierr); CHKERRQ(ierr)
#endif
    case('ILU')
      call PCSetType(t2,PCILU,ierr); CHKERRQ(ierr)
    case default
      call PCSetType(t2,PCBJACOBI,ierr); CHKERRQ(ierr)
  end select

  ctx%T2_PC = t2

end subroutine CPRCreateT2

!  end of CPR Creation Routines  
! ************************************************************************** !

! ************************************************************************** !
! suplementary setup/init/deinit/routines  

subroutine SetCPRDefaults(ierr)
  ! 
  ! set in Petsc's options tables some good
  ! default options for th HYPRE BOOMERAMG
  ! AMG preconditioner.
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string, word

  ! set sensible defualts for the boomeramg pc here,
  ! the defaults it comes with are rarely preferable to
  ! us.

  ! strong threshold, 0.5 is ok for 3D
  string =  '-pc_hypre_boomeramg_strong_threshold'
  word   =  '0.5'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! coarsen type, PMIS is generally more efficient 
  string =  '-pc_hypre_boomeramg_coarsen_type'
  word   =  'PMIS'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! interpolation type, ext+i is reccomended 
  string =  '-pc_hypre_boomeramg_interp_type'
  word   =  'ext+i'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! relaxer, Jacobi should be super weak but cheap 
  string =  '-pc_hypre_boomeramg_relax_type_all'
  word   =  'Jacobi'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

end subroutine SetCPRDefaults

! ************************************************************************** !

subroutine AllocateWorkersInCPRStash(ctx, n, b)
  ! 
  ! allocate the various arrays that are used to hold data 
  ! during the pressure extraction phase of the CPR pc.
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none


  type(cpr_pc_type) :: ctx
  PetscInt :: entriesInARow, entriesInAReducedRow

  PetscInt :: n, b

  ! n: the maxiumum number of nonzero columns in any row
  !    of the matrix we will work with
  ! b: the block size of the matrix we will work with

  ! entriesInARow: the maximum number of entries in a row in
  !               the matrix we will work with, plus a
  !               safety buffer in case it was wrong
  entriesInARow = n + 10

  ! entriesInAReducedRow: the maximum number of entries in a row
  !                       of the extracted matrix we will work with.
  !
  entriesInAReducedRow = n/b
  !                       plus a safety buffer in case it was wrong
  entriesInAReducedRow = entriesInAReducedRow + 10

  ! vals: a real buffer to hold the matrix values output from
  !       matgetrows
  allocate(ctx%vals (0:entriesInARow))
  ! insert_vals: a real buffer to hold the values we will
  !               insert to the reduced matrix
  allocate(ctx%insert_vals (0:entriesInAReducedRow))
  ! colIdx: an integer buffer to hold the column indexes
  !         output from matgetvals
  allocate(ctx%colIdx (0:entriesInARow))
  ! colIdx_keep: will copy colIdx into here, since matrestorerows
  !              resets colIdx
  allocate(ctx%colIdx_keep (0:entriesInARow))
  ! insert_colIdx: an integer buffer of the column indexes
  !                to be input to the reduced matrix we will
  !                work with
  allocate(ctx%insert_colIdx (0:entriesInAReducedRow))
  ! all_vals: will just copy every value from the current row block
  !           (i.e. set of b rows) into here to work on more
  !           easilly. Can very much be improved.
  allocate(ctx%all_vals (0:b-1, 0:entriesInARow))

end subroutine AllocateWorkersInCPRStash

! ************************************************************************** !

subroutine DeallocateWorkersInCPRStash(ctx)
  ! 
  ! DEallocate the various arrays that are used to hold data 
  ! during the pressure extraction phase of the CPR pc.
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none


  type(cpr_pc_type) :: ctx

  deallocate(ctx%vals)
  deallocate(ctx%insert_vals)
  deallocate(ctx%colIdx)
  deallocate(ctx%colIdx_keep)
  deallocate(ctx%insert_colIdx)
  deallocate(ctx%all_vals)

  !! also the pointers:
  nullify(ctx%option)

end subroutine DeallocateWorkersInCPRStash

! end of suplementary setup/init/deinit/routines  
! ************************************************************************** !


! ************************************************************************** !
!        pressure system extraction routines  

subroutine MatGetSubQIMPES(a, ap, factors1Vec,  ierr, ctx)
  ! 
  ! extraction of the pressure system matrix of for the 
  ! CPR preconditioner, and store the pivoting factors
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  ! 3x3 blocks ONLY

  implicit none


  Mat :: a, ap
  Vec :: factors1Vec
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  ! the elements of the diagonal block:
  PetscReal :: aa,bb,cc,dd,ee,ff,gg,hh,ii
  ! misc workers:
  PetscReal :: det, fac0, fac1, fac2, sm
  PetscInt :: b, rws, cls, nblks, nblks_l, firstRow, cur_coldex, ncolblks, &
              firstrowdex, loopdex, i, j, numcols, numcols_keep
  PetscMPIInt :: rnk, r_st, r_nd



  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rnk, ierr)
  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)
  call MatGetBlockSize(a,b,ierr); CHKERRQ(ierr)
  call MatGetSize(A, rws, cls, ierr); CHKERRQ(ierr)
  if (rws /= cls) then
    ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call printErrMsg(ctx%option)
  endif
  nblks = rws/b
  nblks_l = (r_nd-r_st)/b

  sm = 0.d0

  ! loop over the row blocks
  do i = 0,nblks_l-1

    firstRow = i*b + r_st

    ! a) extract first row
    call MatGetRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)
    ! store vals since we have to put them back
    ! store colIdx this time as well
    do j = 0,numcols-1
      ctx%all_vals(0, j) = ctx%vals(j)
      ctx%colIdx_keep(j) = ctx%colIdx(j)
    end do
    ! b) we can get index of diagonal block here
    firstrowdex = -1
    do loopdex = 0,numcols-1,3
      if (ctx%colIdx(loopdex) == firstrow) then
        firstrowdex = loopdex
      endif
    enddo
    if (firstrowdex == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPES, cannot find diagonal entry, check matrix'
      call printErrMsg(ctx%option)
    endif
    aa = ctx%vals(firstrowdex)
    bb = ctx%vals(firstrowdex+1)
    cc = ctx%vals(firstrowdex+2)
    ! restore
    call MatRestoreRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! c) second row
    call MatGetRow(a, firstRow+1, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                  ierr); CHKERRQ(ierr)
    do j = 0,numcols-1
      ctx%all_vals(1, j) = ctx%vals(j)
    end do
    dd = ctx%vals(firstrowdex)
    ee = ctx%vals(firstrowdex+1)
    ff = ctx%vals(firstrowdex+2)
    call MatRestoreRow(a, firstRow+1, numcols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! d) third row
    call MatGetRow(a, firstRow+2, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                 ierr); CHKERRQ(ierr)
    do j = 0,numcols-1
      ctx%all_vals(2, j) =ctx%vals(j)
    end do
    gg = ctx%vals(firstrowdex)
    hh = ctx%vals(firstrowdex+1)
    ii = ctx%vals(firstrowdex+2)
    numcols_keep = numcols
    call MatRestoreRow(a, firstRow+2, numcols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! e) factors
    sm = abs(aa)+abs(dd)+abs(gg)
    det = aa*(ee*ii - ff*hh) - bb*(dd*ii-ff*gg) + cc*(dd*hh-ee*gg)
    fac0 = sm*(ee*ii-ff*hh)/det
    fac1 = sm*(cc*hh-bb*ii)/det
    fac2 = sm*(bb*ff-cc*ee)/det

    ! f) store vectors
    call VecSetValue(factors1Vec, firstRow, fac0, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, firstRow+1, fac1, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, firstRow+2, fac2, INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    ! g) prepare to set values
    insert_rows(0) = i + r_st/b
    ncolblks = numcols_keep/b
    do j = 0,ncolblks-1
      cur_coldex = j*b
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_coldex)/b


      ctx%insert_vals(j) = fac0*ctx%all_vals(0, cur_coldex)
      ctx%insert_vals(j) = ctx%insert_vals(j) + fac1*ctx%all_vals(1,cur_coldex)
      ctx%insert_vals(j) = ctx%insert_vals(j) + fac2*ctx%all_vals(2,cur_coldex)
    end do

    ! h) set values
    call MatSetValues(ap, 1, insert_rows, ncolblks, &
                      ctx%insert_colIdx(0:ncolblks-1), &
                      ctx%insert_vals(0:ncolblks-1), INSERT_VALUES, ierr)
    CHKERRQ(ierr)


  end do
  call MatAssemblyBegin(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

end subroutine MatGetSubQIMPES

! ************************************************************************** !

subroutine MatGetSubQIMPES_var(a, ap, factors1Vec,  ierr, &
                              b, ctx                       )
  ! 
  ! extraction of the pressure system matrix of for the 
  ! CPR preconditioner, and store the pivoting factors
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  ! for arbitrary block size b

  implicit none

  Mat :: a, ap
  Vec :: factors1vec
  PetscErrorCode :: ierr
  PetscInt :: b
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  PetscReal, dimension(0:b-1,0:b-1) :: diag_block
  PetscReal, dimension(0:b-1) :: local_factors
  PetscMPIInt :: rnk, r_st, r_nd
  PetscInt :: rws, cls, nblks, nblks_l, firstRow, cur_coldex, ncolblks, &
              firstrowdex, loopdex, i, j, k, numcols, numcols_keep
  PetscReal :: sm, offdiagsum, diagpart
  ! for lapack inversion routine:
  integer, dimension(1:b) :: ipiv
  PetscReal, dimension(1:b) :: work
  PetscInt :: lwork, invinfo, luinfo

  lwork = b

  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rnk, ierr)
  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)
  call MatGetSize(a, rws, cls, ierr); CHKERRQ(ierr)
  if (rws /= cls) then
    ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call printErrMsg(ctx%option)
  endif

  nblks = rws/b
  nblks_l = (r_nd-r_st)/b

  ! loop over the row blocks
  do i = 0,nblks_l-1

    firstRow = i*b + r_st

    ! first row special treatment
    call MatGetRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)

    ! store both values of row and the col indexs
    do k = 0,numcols-1
      ctx%all_vals(0, k) = ctx% vals(k)
      ctx%colIdx_keep(k) = ctx%colIdx(k)
    end do
    ! get index of diagonal block
    firstrowdex = -1
    do loopdex = 0,numcols-1,b
      if (ctx%colIdx(loopdex) == firstrow) then
        firstrowdex = loopdex
      endif
    enddo
    if (firstrowdex == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPES_var, cannot find diagonal entry, check matrix'
      call printErrMsg(ctx%option)
    endif
    numcols_keep = numcols
    ! restore first row
    call MatRestoreRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! loop over remaining rows
    do j = 1,b-1
      call MatGetRow(a, firstRow+j, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                     ierr); CHKERRQ(ierr)
      ! harvest values
      do k = 0,numcols-1
        ctx%all_vals(j, k) = ctx%vals(k)
      end do
      call MatRestoreRow(a, firstRow+j, numcols, PETSC_NULL_INTEGER, &
                         ctx%vals, ierr); CHKERRQ(ierr)
    enddo

    ! get inverse of block
    diag_block = ctx%all_vals(0:b-1, firstrowdex:firstrowdex+b-1)
    ! invert
    ! NOTE: this can fail and seemingly just return permutation
    call DGETRF(b, b, diag_block, b, IPIV, luinfo ) ! factorize first
    call DGETRI(b,   diag_block,       b,   IPIV,   WORK,     LWORK, invinfo)

    if (invinfo > 0) then
      ctx%option%io_buffer = 'MatGetSubQIMPES_var, singular diagonal block'
      call printErrMsg(ctx%option)
    endif
    ! scaling factor: the sum of abs of the first column of
    ! diagonal
    sm = 0.d0
    do j = 0,b-1
      !sm = sm + abs(ctx%all_vals(j, firstrowdex))
      sm = sm + abs(ctx%all_vals(j, firstrowdex+ctx%exrow_offset))
    end do
    ! factors: take the top row of the inverse block
    ! and scale
    do j = 0,b-1
        local_factors(j) =  sm*diag_block(0, j) ! diag block has been
                                                ! replaced with inverse by this point

    end do

    offdiagsum = 0.d0
    diagpart   = 0.d0

    ! prepare to set values
    insert_rows(0) = i + r_st/b
    ncolblks = numcols_keep/b

    do j = 0,ncolblks-1
      cur_coldex = j*b
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_coldex)/b

      ctx%insert_vals(j) = 0.d0
      do k = 0,b-1
        ctx%insert_vals(j) = ctx%insert_vals(j) + local_factors(k)*ctx%all_vals(k, cur_coldex)
      enddo

      if (ctx%insert_colIdx(j) == insert_rows(0)) then
        diagpart = abs(ctx%insert_vals(j))
      else
        offdiagsum = offdiagsum + abs(ctx%insert_vals(j))
      endif

    end do


    do j = 0,b-1
       call VecSetValue(factors1Vec, firstRow+j, local_factors(j), &
                         INSERT_VALUES, ierr); CHKERRQ(ierr)
    end do


    ! set values
    call MatSetValues(ap, 1, insert_rows, ncolblks, &
                      ctx%insert_colIdx(0:ncolblks-1), &
                      ctx%insert_vals(0:ncolblks-1), INSERT_VALUES, ierr)
     CHKERRQ(ierr)

  end do
  call MatAssemblyBegin(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

end subroutine MatGetSubQIMPES_var

! ************************************************************************** !

subroutine QIRHS(factors, worker, r, rhat, ierr)
  ! 
  ! extract the RHS of the pressure system for the CPR
  ! preconditioner, given pivoting factors (factors)
  ! and full system rhs (r).
  ! output is rhat.
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  Vec :: factors, worker, r, rhat
  PetscErrorCode :: ierr

  PetscInt :: b, k

  call VecPointwiseMult(worker, factors, r, ierr);CHKERRQ(ierr)
  call VecGetBlockSize(worker,b,ierr); CHKERRQ(ierr)
  k = 0
  call VecStrideGather(worker, k, rhat, INSERT_VALUES, ierr);CHKERRQ(ierr)
  do k = 1,b-1
   call VecStrideGather(worker, k, rhat, ADD_VALUES, ierr);CHKERRQ(ierr)
  end do

end subroutine QIRHS

! end of pressure system extraction routines  
! ************************************************************************** !

! ************************************************************************** !
!       misc routines  

subroutine MatGetMaxRowCount(a, mx, ierr)
  ! 
  ! loop over rows of matrix, return the greatest number
  ! of nonzeros in a row
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  Mat :: a
  PetscInt :: mx
  PetscErrorCode :: ierr

  PetscInt :: mx_loc, r_st, r_nd, i, numcols

  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)

  mx = 0
  mx_loc = 0

  do i = r_st,r_nd-1
    call MatGetRow(a, i, numcols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, &
                   ierr); CHKERRQ(ierr)
    if (numcols > mx_loc) then
      mx_loc = numcols
    endif
    call MatRestoreRow(a, i, numcols, PETSC_NULL_INTEGER, &
                       PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
  end do

  call MPI_Allreduce(mx_loc, mx, ONE_INTEGER_MPI, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)


end subroutine MatGetMaxRowCount

! end of misc routines  
! ************************************************************************** !

end module CPR_Preconditioner_module
