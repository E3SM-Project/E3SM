module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_flush
  use RDycoreSpmdMod, only : masterproc, iam, mpicom_rof
  use RDycore_varctl, only : iulog

  implicit none

  private

  type(RDy)                               :: rdy_                      ! RDycore data structure

  PetscInt              , public          :: num_cells_owned           ! number of cells that locally owned
  PetscInt              , public          :: num_cells_global          ! total number of cells in the mesh
  PetscInt              , public, pointer :: natural_id_cells_owned(:) ! natural IDs of cells that are locally owned

  PetscInt                                :: nstep

  PetscReal             , public, pointer :: total_runoff_data(:)      ! the water source to RDycore's SWE


  public :: rdycore_init
  public :: rdycore_run
  public :: rdycore_final

contains

  !-----------------------------------------------------------------------
  subroutine rdycore_init()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !
    use petsc
    implicit none
    !
    ! !LOCAL VARIABLES:
    PetscViewer           :: viewer
    character(len=1024)   :: config_file
    PetscInt              :: g
    PetscMPIInt           :: myrank
    Vec                   :: owner_mpi, owner_seq
    Vec                   :: nid_owned_seq
    PetscScalar, pointer  :: vec_ptr(:), rank(:)
    PetscInt   , pointer  :: int_ptr(:)
    IS                    :: is_from, is_to
    VecScatter            :: scatter
    character(len=100)    :: string
    PetscErrorCode        :: ierr

    config_file = 'rdycore.yaml'

    if (masterproc) then
       write(iulog,*)'RDycore model initialization'
    end if

    ! set PETSc's communicator
    PETSC_COMM_WORLD = mpicom_rof
    PetscCallA(PetscInitialize(ierr))
 
    ! initialize subsystems
    PetscCallA(RDyInit(ierr))
 
    ! create rdycore and set it up with the given file
    PetscCallA(RDyCreate(PETSC_COMM_WORLD, config_file, rdy_, ierr))
    PetscCallA(RDySetup(rdy_, ierr))
 
    ! allocate memory for grid-level runoff dataset
    PetscCallA(RDyGetNumOwnedCells(rdy_, num_cells_owned, ierr))
    allocate(total_runoff_data(num_cells_owned))

    ! get natural ID of owned cells
    allocate(natural_id_cells_owned(num_cells_owned))
    PetscCallA(RDyGetOwnedCellNaturalIDs(rdy_, num_cells_owned, natural_id_cells_owned, ierr))

    ! get number of global cell
    PetscCallA(RDyGetNumGlobalCells(rdy_, num_cells_global, ierr))

    nstep = 0

    if (masterproc) then
       write(iulog,*)'RDycore model initialization completed'
    end if

  end subroutine rdycore_init

  !-----------------------------------------------------------------------
  subroutine rdycore_run(coupling_dt_in_sec, rstwr, rdate)
    !
    use RDycoreRestFile, only : RDycoreRestFileNameBase
    !
    implicit none
    !
    integer          , intent(in) :: coupling_dt_in_sec ! runtime for rdycore before returning
    logical          , intent(in) :: rstwr              ! if .true., then write restart file
    character(len=*) , intent(in) :: rdate              ! date char string for restart file name
    !
    character(len=256)   :: dateTimeString
    character(len=1024)  :: restart_filename_base
    real(r8)             :: dtime
    PetscInt             :: t!, nstep
    integer(RDyTimeUnit) :: time_unit
    PetscReal            :: time_dn, time_up, cur_time, cur_rain
    PetscBool            :: found
    PetscInt             :: g, idx
    PetscErrorCode       :: ierr

    dtime    = coupling_dt_in_sec * 1._r8
    nstep    = nstep + 1
    cur_time = (nstep-1)*dtime

    !call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(iulog,*)'Beginning timestep of RDycore  : '!,trim(dateTimeString)
       call shr_sys_flush(iulog)
    end if

    ! Set the coupling time step
    PetscCallA(RDyGetTimeUnit(rdy_, time_unit, ierr))
    PetscCallA(RDySetCouplingInterval(rdy_, time_unit, dtime, ierr))

    ! Set water source term in RDycore
    PetscCallA(RDySetDomainWaterSource(rdy_, num_cells_owned, total_runoff_data, ierr))

    ! Run the simulation to completion.
    PetscCallA(RDyAdvance(rdy_, ierr))

    if (rstwr) then
       ! determine the name of the restart file
       restart_filename_base = RDycoreRestFileNameBase(rdate)

       PetscCallA(RDyWriteHDF5CheckpointFile(rdy_, restart_filename_base, ierr))
    end if
  end subroutine rdycore_run

  !-----------------------------------------------------------------------
  subroutine rdycore_final()
    !
    implicit none
    !
    !
    ! !LOCAL VARIABLES:
    PetscErrorCode :: ierr

    ! deallocate memory
    deallocate(natural_id_cells_owned)
    deallocate(total_runoff_data)

    ! destroy RDy object
    PetscCallA(RDyDestroy(rdy_, ierr));

    ! finalize
    PetscCallA(RDyFinalize(ierr));

  end subroutine rdycore_final

end module rdycoreMod
