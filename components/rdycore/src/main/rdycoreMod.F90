module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  use rdycore
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_flush
  use RDycoreSpmdMod, only : masterproc, iam, mpicom_rof

  implicit none

  private

  type(RDy)                               :: rdy_                      ! RDycore data structure

  PetscInt              , public          :: num_cells_owned           ! number of cells that locally owned
  PetscInt              , public          :: num_cells_global          ! total number of cells in the mesh
  PetscInt              , public, pointer :: natural_id_cells_owned(:) ! natural IDs of cells that are locally owned

  PetscInt                                :: nstep

  PetscReal             , public, pointer :: total_runoff_data(:)      ! the water source to RDycore's SWE

  integer               , public          :: iulog = 6
  character(len=16)     , public          :: inst_name
  character(len=16)     , public          :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer               , public          :: inst_index          ! number of current instance (ie. 1)


  public :: rdycore_init
  public :: rdycore_run
  public :: rdycore_final

contains

  !-----------------------------------------------------------------------
  subroutine rdycore_init(iulog)
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !
    use petsc
    implicit none
    !
    integer, intent(in)   :: iulog
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
    PetscCallA(RDyGetNumLocalCells(rdy_, num_cells_owned, ierr))
    allocate(total_runoff_data(num_cells_owned))

    ! get natural ID of owned cells
    allocate(natural_id_cells_owned(num_cells_owned))
    PetscCallA(RDyGetLocalCellNaturalIDs(rdy_, num_cells_owned, natural_id_cells_owned, ierr))

    ! get number of global cell
    PetscCallA(RDyGetNumGlobalCells(rdy_, num_cells_global, ierr))

    nstep = 0

    if (masterproc) then
       write(iulog,*)'RDycore model initialization completed'
    end if

  end subroutine rdycore_init

  !-----------------------------------------------------------------------
  subroutine rdycore_run(iulog)
    !
    implicit none
    !
    integer, intent(in)   :: iulog
    !
    character(len=256)   :: dateTimeString
    real(r8)             :: dtime
    PetscInt             :: t!, nstep
    integer(RDyTimeUnit) :: time_unit
    PetscReal            :: time_dn, time_up, cur_time, cur_rain
    PetscBool            :: found
    PetscInt             :: g, idx
    PetscErrorCode       :: ierr

    dtime    = 1800._r8 !get_step_size()
    nstep    = nstep + 1
    cur_time = (nstep-1)*dtime

    !call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(*,*)'Beginning timestep of RDycore  : '!,trim(dateTimeString)
       call shr_sys_flush(iulog)
    end if

    ! Set the coupling time step
    PetscCallA(RDyGetTimeUnit(rdy_, time_unit, ierr))
    PetscCallA(RDySetCouplingInterval(rdy_, time_unit, dtime, ierr))

    ! Run the simulation to completion.
    PetscCallA(RDyAdvance(rdy_, ierr))

    PetscCallA(RDySetDomainWaterSource(rdy_, num_cells_owned, total_runoff_data, ierr))

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
