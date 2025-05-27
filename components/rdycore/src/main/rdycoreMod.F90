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

  integer               , public, pointer :: rdycore_pocn(:)           ! PE rank for each grid cell

  PetscReal             , pointer         :: total_runoff_data(:)      ! the water source to RDycore's SWE

  integer               , public          :: iulog = 6
  character(len=16)     , public          :: inst_name
  character(len=16)     , public          :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer               , public          :: inst_index          ! number of current instance (ie. 1)


  public :: rdycore_init

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
 
    ! allocate memory for grid-level rain data
    PetscCallA(RDyGetNumLocalCells(rdy_, num_cells_owned, ierr))
    allocate(total_runoff_data(num_cells_owned))

    allocate(natural_id_cells_owned(num_cells_owned))
    PetscCallA(RDyGetLocalCellNaturalIDs(rdy_, num_cells_owned, natural_id_cells_owned, ierr))

    ! find the rank for the processor
    PetscCallA(mpi_comm_rank(PETSC_COMM_WORLD, myrank, ierr))

#if 0
    PetscCallA(VecCreateSeq(PETSC_COMM_SELF, num_cells_owned, nid_owned_seq, ierr))

    write(string,*)myrank
    PetscCallA(VecGetArrayF90(nid_owned_seq, vec_ptr,ierr))
    vec_ptr(:) = natural_id_cells_owned(:) * 1.d0
    PetscCallA(VecRestoreArrayF90(nid_owned_seq, vec_ptr,ierr))

    string = 'natural_id_seq_' // trim(adjustl(string)) // '.txt'
    PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_SELF, trim(string), viewer, ierr))
    PetscCallA(VecView(nid_owned_seq, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))
    PetscCallA(VecDestroy(nid_owned_seq, ierr))
#endif

    PetscCallA(RDyGetNumGlobalCells(rdy_, num_cells_global, ierr))

    ! create a MPI and sequential Vec
    PetscCallA(VecCreateMPI(PETSC_COMM_WORLD, num_cells_owned, PETSC_DETERMINE, owner_mpi, ierr))
    PetscCallA(VecCreateSeq(PETSC_COMM_SELF, num_cells_global, owner_seq, ierr))

    allocate(rank(num_cells_owned))
    write(iulog,*)'myrank ',myrank
    rank(:) = myrank
    PetscCallA(VecSetValues(owner_mpi, num_cells_owned, natural_id_cells_owned, rank, INSERT_VALUES, ierr))
    deallocate(rank)
    PetscCallA(VecAssemblyBegin(owner_mpi, ierr))
    PetscCallA(VecAssemblyEnd(owner_mpi, ierr))

#if 0
    string = 'owner_mpi.txt'
    PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(string), viewer, ierr))
    PetscCallA(VecView(owner_mpi, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))
#endif

    ! create index set to scatter from the MPI Vec and to sequential Vec
    allocate(int_ptr(num_cells_global))
    do g = 1, num_cells_global
       int_ptr(g) = g - 1
    end do
    PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD, num_cells_global, int_ptr, PETSC_COPY_VALUES, is_from, ierr))
    PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD, num_cells_global, int_ptr, PETSC_COPY_VALUES, is_to, ierr))
    deallocate(int_ptr)

    ! create the VecScatter
    PetscCallA(VecScatterCreate(owner_mpi, is_from, owner_seq, is_to, scatter, ierr))
    PetscCallA(ISDestroy(is_from, ierr))
    PetscCallA(ISDestroy(is_to, ierr))

    ! scatter th data
    PetscCallA(VecScatterBegin(scatter, owner_mpi, owner_seq, INSERT_VALUES, SCATTER_FORWARD, ierr))
    PetscCallA(VecScatterEnd(scatter, owner_mpi, owner_seq, INSERT_VALUES, SCATTER_FORWARD, ierr))

#if 0
    write(string,*)myrank
    string = 'owner_seq_' // trim(adjustl(string)) // '.txt'
    PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_SELF, trim(string), viewer, ierr))
    PetscCallA(VecView(owner_seq, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))
#endif

    ! save PE number for each grid cell that MOSART will use
    allocate(rdycore_pocn(num_cells_global))

#if 0
    ! determine the begin/end bounds of grid cells
    PetscCallA(VecGetArrayF90(owner_seq, vec_ptr, ierr))
    rdy_bounds%begg = 1
    do g = 1, num_cells_global
       rdycore_pocn(g) = int(vec_ptr(g))
       if (int(vec_ptr(g)) < myrank) then
          rdy_bounds%begg = rdy_bounds%begg + 1
       end if
    end do
    rdy_bounds%endg = rdy_bounds%begg - 1 + num_cells_owned
    PetscCallA(VecRestoreArrayF90(owner_seq, vec_ptr, ierr))

    ! allocate data structure for exchanging data from land to rdycore
    call lnd2rdy_vars%Init(rdy_bounds)
#endif

    ! free up memory
    PetscCallA(VecDestroy(owner_mpi, ierr))
    PetscCallA(VecDestroy(owner_seq, ierr))

    if (masterproc) then
       write(iulog,*)'RDycore model initialization completed'
    end if

  end subroutine rdycore_init

end module rdycoreMod
