!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_compute_nonhydrostatic_pressure
!
!> \brief MPAS ocean module to solve for the nonhydrostatic pressure correction
!> \author Luke Van Roekel
!> \date January 2021
!> \details
!    This module will derive the nonhydrostatic pressure correction
!     uses petsc with the ILU preconditioner for solving
!
!-----------------------------------------------------------------------

#ifdef USE_PETSC
#define PETSC_USE_FORTRAN_MODULES 1 
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscmat.h"
!#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petsc.h"
#endif

module ocn_compute_nonhydrostatic_pressure

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_constants
  use mpas_log

  use ocn_constants
  use ocn_config
  use ocn_mesh
  use ocn_diagnostics_variables

#ifdef USE_PETSC
  !PETSC includes for the matrix solve
  use petscsys
  use petscvec
  use petscmat
  use petscpc
  use petscksp
  use petsc
#endif
  implicit none
  private
  save

  !--------------------------------------------------------------------
  !
  ! Public member functions
  !
  !--------------------------------------------------------------------

  public :: ocn_nonhydrostatic_pressure_tend, &
            ocn_nonhydrostatic_solver_init,   &
            ocn_nonhydrostatic_solver_finalize, &
            ocn_nonhydrostatic_solver_everything, &
            ocn_nonhydrostatic_pressure_update_velocity

  !--------------------------------------------------------------------
  !
  ! Private module variables
  !
  !--------------------------------------------------------------------

  double precision :: norm ! norm of solution error
#ifdef USE_PETSC
  PetscInt ::  i, j, II, JJ, m, n, its, Istart, Iend
  PetscInt ::  row, col, ione, globalM, first, last, firstVal
  PetscInt :: overlap,nlocal
  PetscErrorCode :: petsc_err, ierr2
  PetscMPIInt :: petsc_rank, petsc_size
  PetscBool :: flag
  PetscScalar :: v, one, neg_one
  Vec :: x, b, u, negOneVec ! approximate solution, right hand side vector, exact solution vector
  Mat :: Amat ! Matrix that defines the system
  Mat :: Pmat ! Preconditioner Matrix
  KSP :: ksp ! krylov subspace method context
  KSP,allocatable,dimension(:) :: subksp
  PetscRandom :: rctx ! random number generator
  PC :: pc, subpc ! preconditioner context
  PCType :: ptype ! flag to set type of preconditioner MAYBE make this a NL option?
  PetscReal :: tol 

  !Stuff for the new vecscatter
  IS :: origIS, destIS !PETSc index sets
  VecScatter :: scatter ! PETSc vec scatter context computed once and reused
  Vec :: local_x ! local value of things to grab from global vector
  PetscInt,allocatable,dimension(:) :: destIndex, origIndex
  PetscScalar,pointer :: values(:)
#endif

  real(kind=RKIND) :: surfaceBoundaryConditionTop, surfaceBoundaryConditionBottom

  contains

    subroutine ocn_nonhydrostatic_pressure_tend(normalVelocityNew, &
      vertVelocityNonhydroNew, layerThickness, layerThicknessEdge, &
      nhPressureCorrection, mpi_comm, dt)

    real,dimension(:,:), intent(in) :: normalVelocityNew, vertVelocityNonhydroNew, &
      layerThickness, layerThicknessEdge

    real(kind=RKIND), intent(in) :: dt

    real,dimension(:,:), intent(inout) :: nhPressureCorrection

    integer, intent(in) :: mpi_comm

    logical :: reuse_preconditioner 

    reuse_preconditioner = .false.

    call ocn_nonhydrostatic_solver_everything(normalVelocityNew, &
       vertVelocityNonhydroNew, layerThickness, layerThicknessEdge, mpi_comm, &
       dt, nhPressureCorrection)

  end subroutine ocn_nonhydrostatic_pressure_tend

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_pressure_update_velocity
  !
  !> \brief  Update the normal and vertical velocity 
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  creates KSP solver context, gets preconditioner and solves
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_pressure_update_velocity(normalVelocityNew, &
      vertVelocityNonhydroNew, layerThickness, layerThicknessEdge,           &
      nonhydrostaticPressure, nonhydrostaticPressureOld,                 &
      nhPressureCorrection, dt)

    real(KIND=RKIND),intent(in) :: dt

    real(KIND=RKIND),dimension(:,:), intent(inout) ::  &
      normalVelocityNew,    &
      vertVelocityNonhydroNew

    real(KIND=RKIND),dimension(:,:),intent(in) :: &
      layerThicknessEdge, &
      layerThickness

    real(KIND=RKIND),dimension(:,:),intent(inout) :: &
      nonhydrostaticPressure, &
      nonhydrostaticPressureOld, &
      nhPressureCorrection

    integer :: cell1, cell2, nEdges, iCell, iEdge, nCells, i, k

    real(kind=RKIND), dimension(nVertLevels) :: div_hu

    real(kind=RKIND) :: invAreaCell, r_tmp

    nEdges = nEdgesHalo(2)

    call mpas_timer_start('NH solve - update velocity',.false.)
    ! loop over edges to update normalVelocity
    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(k, cell1, cell2)
    do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       do k = 1,maxLevelEdgeTop(iEdge)
          normalVelocityNew(k,iEdge) = normalVelocityNew(k,iEdge) - dt * &
            (nhPressureCorrection(k,cell2) - nhPressureCorrection(k,cell1)) / &
            dcEdge(iEdge)*edgeMask(k,iEdge) 
       end do
    end do
    !$omp end do
    !$omp end parallel

    nCells = nCellsHalo( 1 )

    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(i, r_tmp, iEdge, invAreaCell, k, div_hu)
    do iCell = 1, nCells
       div_hu(:) = 0.0_RKIND
       invAreaCell = 1.0_RKIND / areaCell(iCell)
       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          do k = 1,maxLevelEdgeTop(iEdge)
             r_tmp = dvEdge(iEdge)*normalVelocityNew(k,iEdge)*invAreaCell
             div_hu(k) = div_hu(k) - layerThicknessEdge(k,iEdge)*edgeSignOnCell(i,iCell)*r_tmp
          end do
       end do

       vertVelocityNonhydroNew(maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
       do k=maxLevelCell(iCell),1,-1
          vertVelocityNonhydroNew(k,iCell) = vertVelocityNonhydroNew(k+1,iCell) - div_hu(k)
       end do
    end do
    !$omp end do
    !$omp end parallel

    !Update the nonhydrostatic pressure to be second order accurate

    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(k)
    do iCell=1,nCells
       do k = 1,maxLevelCell(iCell)
          nonhydrostaticPressure(k,iCell) = nonhydrostaticPressureOld(k,iCell) + &
                  nhPressureCorrection(k,iCell)
       end do 
    end do
    !$omp end do
    !$omp end parallel

    call mpas_timer_stop('NH solve - update velocity')

  end subroutine ocn_nonhydrostatic_pressure_update_velocity

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_solver_everything
  !
  !> \brief   Create the matrix and right-hand-side and then solve the
  !>          linear system
  !> \author  Sara Calandrini, Darren Engwirda, Luke Van Roekel
  !> \date    November 2021
  !> \details
  !>  creates KSP solver context, gets preconditioner and solves
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_everything(normalVelocityNew, &
      vertVelocityNonhydroNew, layerThickness, layerThicknessEdge, mpi_comm, &
      dt, nhPressureCorrection)

     implicit none

     real(kind=RKIND), dimension(:,:), intent(in) :: normalVelocityNew, &
                                                     vertVelocityNonhydroNew
     real(kind=RKIND), dimension(:,:), intent(in) :: layerThickness, &
                                                     layerThicknessEdge
     real(kind=RKIND), intent(in) :: dt
     real(kind=RKIND), dimension(:,:), intent(inout) :: nhPressureCorrection
     integer, intent(in) :: mpi_comm

  !--------------------------------------------------------------------
  !
  ! subroutine variables
  !
  !--------------------------------------------------------------------

     integer :: k, iCell, jCell, iEdge, jEdge, iCell3d, jCell3d
     integer :: nnzRowMax, vecSpot
     real(kind=RKIND) :: xEdge, thicknessUpper, thicknessLower
     real(kind=RKIND), dimension(nVertLevels) :: div   

#ifdef USE_PETSC
     PetscScalar :: insertVal
     PetscInt :: itNum

     call mpas_timer_start('NH solve - matrix',.false.)

     ! form matrix: del^2(q) = 1/dt * div(u*)

     !nnzRowMax = 9 * 2 + 4  ! do this better
     nnzRowMax = 1 + 2 + 7

     call MatZeroEntries(Amat, petsc_err)

     do iCell = 1, nCellsOwned

        ! x-part of del^2 operator

        do jEdge = 1, nEdgesOnCell(iCell)

           iEdge = edgesOnCell(jEdge,iCell)
           jCell = cellsOnCell(jEdge,iCell)

           xEdge = dvEdge(iEdge) / dcEdge(iEdge)

           do k = 1, maxLevelCell(iCell)
              if (boundaryEdge(k,iEdge) .ge. 1) cycle

              ! this map could actually create unreferenced rows, etc
              ! when cells are masked out under bathymetry?
              iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
              jCell3d = (indexToCellId(jCell) - 1) * nVertLevels + k - 1   
     
              insertVal = +1 * layerThicknessEdge(k,iEdge) * xEdge
              call MatSetValues(Amat, 1, iCell3d, 1, jCell3d, insertVal, ADD_VALUES, petsc_err)
              insertVal = -1 * layerThicknessEdge(k,iEdge) * xEdge
              call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)
           end do

        end do

        ! for non-flat bathymetry
        do k = maxLevelCell(iCell)+1, nVertLevels
           iCell3D = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
           insertVal = 1.0_RKIND
           call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)
        end do

        ! z-part of del^2 operator

        ! first layer upper q = 0
        ! first layer lower
        k = 1 !minLevelCell(iCell)

        thicknessUpper = 1.0E-15_RKIND + 0.5_RKIND * ( &
           layerThickness(k,iCell) )
        thicknessLower = 1.0E-15_RKIND + 0.5_RKIND * ( &
           layerThickness(k,iCell) + layerThickness(k+1,iCell))

        iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1

        insertVal = surfaceBoundaryConditionTop*(-1 * areaCell(iCell) / thicknessUpper)
        call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)

        iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
        jCell3d = (indexToCellId(iCell) - 1) * nVertLevels + (k+1) - 1

        insertVal = +1 * areaCell(iCell) / thicknessLower
        call MatSetValues(Amat, 1, iCell3d, 1, jCell3d, insertVal, ADD_VALUES, petsc_err)
        insertVal = -1 * areaCell(iCell) / thicknessLower
        call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)


        do k = 2, maxLevelCell(iCell) - 1

           thicknessUpper = 1.0E-15_RKIND + 0.5_RKIND * ( &
              layerThickness(k,iCell) + layerThickness(k-1,iCell))
           thicknessLower = 1.0E-15_RKIND + 0.5_RKIND * ( &
              layerThickness(k,iCell) + layerThickness(k+1,iCell))

           ! upper
           iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
           jCell3d = (indexToCellId(iCell) - 1) * nVertLevels + (k-1) - 1   
      
           insertVal = +1 * areaCell(iCell) / thicknessUpper
           call MatSetValues(Amat, 1, iCell3d, 1, jCell3d, insertVal, ADD_VALUES, petsc_err)
           insertVal = -1 * areaCell(iCell) / thicknessUpper
           call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)
       
           ! lower
           iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
           jCell3d = (indexToCellId(iCell) - 1) * nVertLevels + (k+1) - 1
 
           insertVal = +1 * areaCell(iCell) / thicknessLower
           call MatSetValues(Amat, 1, iCell3d, 1, jCell3d, insertVal, ADD_VALUES, petsc_err)
           insertVal = -1 * areaCell(iCell) / thicknessLower
           call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)   

        end do

        ! final layer upper
        ! final layer lower dq/dz = 0
        k = maxLevelCell(iCell)

        thicknessUpper = 1.0E-15_RKIND + 0.5_RKIND * ( &
           layerThickness(k,iCell) + layerThickness(k-1,iCell))
        thicknessLower = 1.0E-15_RKIND + 0.5_RKIND * ( &
           layerThickness(k,iCell) )

        iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1

        insertVal = surfaceBoundaryConditionBottom*(-1 * areaCell(iCell) / thicknessLower)
        call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)

        iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
        jCell3d = (indexToCellId(iCell) - 1) * nVertLevels + (k-1) - 1
         
        insertVal = +1 * areaCell(iCell) / thicknessUpper
        call MatSetValues(Amat, 1, iCell3d, 1, jCell3d, insertVal, ADD_VALUES, petsc_err)
        insertVal = -1 * areaCell(iCell) / thicknessUpper
        call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)
        
        ! for non-flat bathymetry
        do k = maxLevelCell(iCell)+1, nVertLevels
           iCell3D = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
           insertVal = 1.0_RKIND
           call MatSetValues(Amat, 1, iCell3d, 1, iCell3d, insertVal, ADD_VALUES, petsc_err)
        end do
 
     end do

     call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,petsc_err)
     call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,petsc_err) 

     call mpas_timer_stop('NH solve - matrix')

     ! form vector: del^2(q) = 1/dt * div(u*) 

     call mpas_timer_start('NH solve - rhs',.false.)

     call vecZeroEntries(b,petsc_err)

     do iCell = 1, nCellsOwned

        div(:) = 0.0_RKIND

        ! x-part of 1/dt * div(u*)

        do jEdge = 1, nEdgesOnCell(iCell)

           iEdge = edgesOnCell(jEdge,iCell)

           do k = 1, maxLevelCell(iCell)

              div(k) = div(k) - &
                 edgeSignOnCell(jEdge,iCell) * dvEdge(iEdge) * &
                    normalVelocityNew(k,iEdge) * layerThicknessEdge(k,iEdge)

           end do

        end do

        ! z-part of 1/dt * div(u*)

        do k = 1, maxLevelCell(iCell)

           div(k) = div(k) + areaCell(iCell) * ( &
              vertVelocityNonhydroNew(k,iCell) - vertVelocityNonhydroNew(k+1,iCell))

        end do

       ! copy into full rhs vector

        do k = 1, maxLevelCell(iCell)
           iCell3d = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
           insertVal = div(k) / dt
           call VecSetValues(b, 1, iCell3d, insertVal, INSERT_VALUES, petsc_err)
        end do 

     end do

     call VecAssemblyBegin(b,petsc_err)
     call VecAssemblyEnd(b,petsc_err)

     call mpas_timer_stop('NH solve - rhs')

     ! PCG solver: del^2(q) = 1/dt * div(u*)

     call mpas_timer_start('NH solve - solver',.false.)

     call KSPCreate(mpi_comm, ksp, petsc_err)
     call KSPSetType(ksp, config_nonhydrostatic_solver_type, petsc_err)
     call KSPSetOperators(ksp,Amat,Amat,petsc_err)
     call KSPGetPC(ksp,pc,petsc_err)
     call PCSetType(pc,config_nonhydrostatic_preconditioner,petsc_err)
     call KSPSetFromOptions(ksp,petsc_err)
     call KSPSetTolerances(ksp,config_petsc_rtol,config_petsc_atol,1.0E4,config_petsc_maxit,petsc_err)
     call KSPSolve(ksp, b, x, petsc_err)
     call KSPGetIterationNumber(ksp, itNum, petsc_err)
     print *, 'iteration count = ',itNum
     call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,petsc_err)
     call KSPDestroy(ksp,ierr2)

     !now update the nonhydrostaticPressure array unpack into k,iCell
     call VecScatterBegin(scatter, x, local_x, INSERT_VALUES, SCATTER_FORWARD, petsc_err)
     call VecScatterEnd(scatter, x, local_x, INSERT_VALUES, SCATTER_FORWARD, petsc_err)
     call VecGetArrayReadF90(local_x, values, petsc_err)

     !$omp parallel
     !$omp do schedule(runtime) private(k,vecSpot)
     do iCell=1,nCellsOwned
        do k=1,maxLevelCell(iCell)
          !vecSpot = (indexToCellId(iCell) - 1) * nVertLevels + k - 1
          !nhPressureCorrection(k,iCell) = values(vecSpot+1)
          vecSpot = (iCell - 1) * nVertLevels + k - 1
          nhPressureCorrection(k,iCell) = values(vecSpot+1)
        end do 
     end do
     !print*, nhPressureCorrection(10,100)
     !$omp end do
     !$omp end parallel 

     call mpas_timer_stop('NH solve - solver') 
#endif

  end subroutine ocn_nonhydrostatic_solver_everything


  !--------------------------------------------------------------------
  ! IIIIII
  ! routine ocn_nonhydrostatic_solver_init
  !
  !> \brief   Initializes matricesIIIIII and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  Initializese the matrix for the nonhydrostatic solve, automatically
  !>  determines what is owned by processor.  Does not fill values.
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_init(mpi_comm, ierr)

    integer, intent(in) :: mpi_comm
    integer, intent(inout) :: ierr
    integer :: commsize, globalCells
    real(kind=RKIND) :: decimals, version
#ifdef USE_PETSC
    PetscInt :: major, minor, subminor, release, locRow, locCol, locals(nCellsOwned)
    PetscReal :: val
    ISLocalToGlobalMapping :: mapping
#endif
    integer :: k, iCell, iter, nnzRowMax

    !nnzRowMax = 9 * 2 + 4
    nnzRowMax = 1 + 2 + 7

#ifndef USE_PETSC
    call mpas_log_write("nonhydrostatic solver requires PETSC library compiled and linked", &
                        mpas_log_crit)
#endif

    if(config_nonhydrostatic_solve_surface_boundary_condition == 'noGradient') then
       surfaceBoundaryConditionTop = 0.0_RKIND
       surfaceBoundaryConditionBottom = 0.0_RKIND
    elseif(config_nonhydrostatic_solve_surface_boundary_condition == 'noPressure') then
       surfaceBoundaryConditionTop = 1.0_RKIND
       surfaceBoundaryConditionBottom = 1.0_RKIND
    elseif(config_nonhydrostatic_solve_surface_boundary_condition == 'pressureTopGradientBottom') then
       surfaceBoundaryConditionTop = 1.0_RKIND
       surfaceBoundaryConditionBottom = 0.0_RKIND
    else
       surfaceBoundaryConditionTop = 0.0_RKIND
       surfaceBoundaryConditionBottom = 1.0_RKIND
    endif

    call MPI_comm_size(mpi_comm, commsize, ierr)

#ifdef USE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("error: petsc initialize failed, error code = $i", mpas_log_crit, &
        intargs=(/petsc_err/))
      ierr = petsc_err
    endif

    call PetscGetVersionNumber(major, minor, subminor, release, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc GetVersionNumber failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    end if

    !Solver routines require version to be greater than 3.5
    decimals = floor(log10(float(minor))) + 1
    version = float(major) + float(minor) / 10.0_RKIND**decimals
    if(version < 3.05) then
      call mpas_log_write("ERROR: Petsc Version required to be greater than 3.5, you are using $r", &
        MPAS_LOG_CRIT, realArgs=(/version/))
      ierr = petsc_err
    end if

    m = nCellsOwned * nVertLevels
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flag,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc OptionsGetInt failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatCreate(mpi_comm, Amat, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatCreate failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatSetSizes(Amat, m, m, PETSC_DETERMINE, PETSC_DETERMINE, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetSizes failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatSetType( Amat, MATAIJ, petsc_err )
    if( commsize == 1) then
      call MatSetType( Amat, MATAIJ, petsc_err )
      call MatSetFromOptions(Amat,petsc_err)
      call MatSeqAIJSetPreallocation(Amat, nnzRowMax, PETSC_NULL_INTEGER, petsc_err)
    else
      call MatSetType( Amat, MATMPIAIJ, petsc_err )
      call MatSetFromOptions(Amat,petsc_err)
      call MatMPIAIJSetPreallocation(Amat,2*nnzRowMax, PETSC_NULL_INTEGER, 2*nnzRowMax, PETSC_NULL_INTEGER, petsc_err)
    end if

    !call MatSetFromOptions(Amat,petsc_err)
    !if(petsc_err .ne. 0) then
    !  call mpas_log_write("ERROR: Petsc MatSetFromOptions failed, error_code = $i", MPAS_LOG_CRIT, &
    !    intArgs=(/petsc_err/))
    !  ierr = petsc_err
    !endif

    !call MatSetUp(Amat,petsc_err)
    !if(petsc_err .ne. 0) then
    !  call mpas_log_write("ERROR: Petsc MatSetUp failed, error_code = $i", MPAS_LOG_CRIT, &
    !    intArgs=(/petsc_err/))
    !  ierr = petsc_err
    !endif

    call MatGetOwnershipRange(Amat,Istart,Iend,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetOwnershipRange failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatGetSize(Amat, globalM, globalM, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetSize failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatGetLocalSize(Amat,locRow,locCol,petsc_err)

    print*, m, globalM, locRow, locCol

    call VecCreate(mpi_comm, b, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecSetSizes(b, m, globalM, petsc_err)
    if ( commsize == 1 ) then
       call VecSetType(b, VECSEQ, petsc_err)
    else
       call VecSetType(b, VECMPI, petsc_err)
    end if

    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecDuplicate(b, x, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate Duplicate, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecSetSizes(x, m, globalM, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecGetOwnershipRange(x,first,last,petsc_err)

    call VecDuplicate(b, negOneVec, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate Duplicate, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecSetSizes(negOneVec, m, globalM, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecGetOwnershipRange(negOneVec,first,last,petsc_err)

    call VecCreate(mpi_comm, local_x, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecSetSizes(local_x, m, PETSC_DETERMINE, petsc_err)
    if ( commsize == 1 ) then
       call VecSetType(local_x, VECSEQ, petsc_err)
    else
       call VecSetType(local_x, VECMPI, petsc_err)
    end if

    allocate(destIndex(globalM), origIndex(globalM))
    ! Try a better way with Index Sets and Vec Scatter for some communication
    ! later move this first part to init as we don't have to do it every time FIXME!
    iter = 1
    do iCell = 1,nCellsOwned
      do k=1,nVertLevels
         destIndex(iter) = (iCell-1)*nVertLevels + k - 1 + first
         origIndex(iter) = (indexToCellID(iCell)-1)*nVertLevels + k - 1
         call VecSetValues(negOneVec, 1, origIndex(iter), -1.0_RKIND, INSERT_VALUES, petsc_err)
         iter = iter + 1
      end do
    end do

    call VecAssemblyBegin(negOneVec, petsc_err)
    call VecAssemblyEnd(negOneVec, petsc_err)

    call ISCreateGeneral(mpi_comm, globalM, origIndex, PETSC_COPY_VALUES, origIS, petsc_err)
    call ISCreateGeneral(mpi_comm, globalM, destIndex, PETSC_COPY_VALUES, destIS, petsc_err)
    call VecScatterCreate(x, origIS, local_x, destIS, scatter, petsc_err)
    
#endif

    end subroutine ocn_nonhydrostatic_solver_init

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_solver_finalize
  !
  !> \brief   Destroys matrices and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>   Destroys matrices and frees memory, finalizes Petsc 
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_finalize()

#ifdef USE_PETSC
    !stuff for PETSC_finalize from the scatter routine
    call ISDestroy(origIS, petsc_err)
    call ISDestroy(destIS, petsc_err)
    call VecScatterDestroy(scatter, petsc_err)

    call MatDestroy(Amat,petsc_err)
    call VecDestroy(b,petsc_err)
    call VecDestroy(x,petsc_err)
    call VecDestroy(negOneVec,petsc_err)
    !call KSPDestroy(ksp,ierr2)
    call VecDestroy(local_x,petsc_err)
 
    !LPV FIXME : for some reason this causes error "Corrupt argument"
    !disable for now
    !call PetscFinalize(petsc_err)
#endif

  end subroutine ocn_nonhydrostatic_solver_finalize  

end module ocn_compute_nonhydrostatic_pressure

