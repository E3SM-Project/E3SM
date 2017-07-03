module Factory_Hydrogeophysics_module

  use Simulation_Hydrogeophysics_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: HydrogeophysicsInitialize

contains

! ************************************************************************** !
subroutine HydrogeophysicsInitialize(simulation)
  ! 
  ! Sets up hydrogeophysics simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  use Option_module
  use Wrapper_Hydrogeophysics_module
  use Input_Aux_module
  use Simulation_Base_class 
  use PM_Base_class  
  use Discretization_module
  use String_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscviewer.h"
  
  class(simulation_hydrogeophysics_type) :: simulation

  type(option_type), pointer :: option
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscInt :: i, num_pflotran_processes, offset, num_slaves
  PetscInt :: local_size
  PetscBool :: option_found
  PetscMPIInt :: mpi_int, process_range(3)
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscInt, allocatable :: int_array(:)
  IS :: is_natural
  IS :: is_petsc
  PetscInt :: istart
  PetscInt :: temp_int
  PetscViewer :: viewer
  Vec :: pflotran_tracer_vec_mpi, pflotran_tracer_vec_seq
  Vec :: pflotran_saturation_vec_mpi, pflotran_saturation_vec_seq
  Vec :: pflotran_temperature_vec_mpi, pflotran_temperature_vec_seq
  VecScatter :: pflotran_scatter
  
#ifdef PETSC_HAVE_MPIUNI
  option%io_buffer = 'HYDROGEOPHYSICs not supported with MPIUNI.'
  call printErrMsg(option)
#else

  option => simulation%option
  ! initialize PETSc Vecs to 0
  pflotran_tracer_vec_mpi = 0
  pflotran_tracer_vec_seq = 0
  pflotran_saturation_vec_mpi = 0
  pflotran_saturation_vec_seq = 0
  pflotran_temperature_vec_mpi = 0
  pflotran_temperature_vec_seq = 0

  string = '-num_slaves'
  num_slaves = UNINITIALIZED_INTEGER
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) num_slaves = i

  ! NOTE: PETSc must already have been initialized here!
  if (option%mycommsize < 3) then
    option%io_buffer = 'At least 3 processes must be allocated to ' // &
      'simulation in order to run hydrogeophysics.'
    call printErrMsg(option)
  endif

  ! store original communicator settings to be set back in HydrogeophysicsStrip
  ! in support of multirealization scenarios.
  simulation%mycomm_save = option%mycomm
  simulation%myrank_save = option%myrank
  simulation%mycommsize_save = option%mycommsize
  simulation%mygroup_save = option%mygroup
  simulation%mygroup_id_save = option%mygroup_id
  
  if (Uninitialized(num_slaves)) then
    num_pflotran_processes = simulation%mycommsize_save / 2
    num_slaves = simulation%mycommsize_save - num_pflotran_processes - 1
  else if (num_slaves <= 0) then
    option%io_buffer = 'Number of slaves must be greater than zero. ' // &
      'Currently set to ' // StringFormatInt(num_slaves) // '.'
    call printErrMsg(option)
  else
    if (num_slaves > simulation%mycommsize_save - 2) then
      option%io_buffer = 'Too many slave processes allocated to ' // &
        'simulation: ' // StringFormatInt(num_slaves)
      call printErrMsg(option)
    endif
    num_pflotran_processes = simulation%mycommsize_save - num_slaves - 1
  endif
  
  write(option%io_buffer,*) 'Number of E4D processes: ', &
    StringFormatInt(num_slaves+1)
  call printMsg(option)
  write(option%io_buffer,*) 'Number of PFLOTRAN processes: ', &
    StringFormatInt(num_pflotran_processes)
  call printMsg(option)
  
  ! split the communicator
  option%mygroup_id = 0
  offset = 0
  if (simulation%myrank_save > num_pflotran_processes-1) then
    option%mygroup_id = 1
    offset = num_pflotran_processes
  endif

  if (option%mygroup_id == 0) then
    simulation%pflotran_process = PETSC_TRUE
  endif

  mycolor_mpi = option%mygroup_id
  mykey_mpi = simulation%myrank_save - offset
  call MPI_Comm_split(simulation%mycomm_save,mycolor_mpi,mykey_mpi, &
                      option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)  
  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)
  
  ! create group shared by both master processes
  call MPI_Comm_group(simulation%mycomm_save,simulation%mygroup_save,ierr)
  call MPI_Group_size(option%mygroup,mpi_int,ierr)
!print *, 1, simulation%myrank_save, simulation%mygroup_save, simulation%mycomm_save
!print *, 2, simulation%myrank_save, option%myrank, option%mygroup, option%mycomm, mpi_int
  mpi_int = 1
  process_range(1) = 0
  process_range(2) = num_pflotran_processes ! includes e4d master due to 
  process_range(3) = 1                        ! zero-based indexing
  call MPI_Group_range_incl(simulation%mygroup_save,mpi_int,process_range, &
                            simulation%pf_e4d_scatter_grp,ierr)
!print *, 3, simulation%myrank_save, simulation%pf_e4d_scatter_grp
  call MPI_Comm_create(simulation%mycomm_save,simulation%pf_e4d_scatter_grp, &
                       simulation%pf_e4d_scatter_comm,ierr)
!print *, 4, simulation%myrank_save, simulation%pf_e4d_scatter_comm, simulation%pf_e4d_master_comm
  if (simulation%pf_e4d_scatter_comm /= MPI_COMM_NULL) then
    call MPI_Comm_rank(simulation%pf_e4d_scatter_comm, &
                       simulation%pf_e4d_scatter_rank,ierr)
    call MPI_Comm_size(simulation%pf_e4d_scatter_comm, &
                       simulation%pf_e4d_scatter_size,ierr)
!    PFE4D_COMM = simulation%pf_e4d_scatter_comm
!print *, 5, simulation%myrank_save, simulation%pf_e4d_scatter_rank, simulation%pf_e4d_scatter_size
    ! remove processes between pf_master and e4d_master
    mpi_int = 1
    process_range(1) = 1
    process_range(2) = num_pflotran_processes-1
    process_range(3) = 1
    ! if there are no process ranks to remove, set mpi_int to zero
    if (process_range(2) - process_range(1) < 0) mpi_int = 0
    call MPI_Group_range_excl(simulation%pf_e4d_scatter_grp,mpi_int, &
                              process_range,simulation%pf_e4d_master_grp,ierr)
!print *, 6, simulation%myrank_save, simulation%pf_e4d_master_grp, ierr
    call MPI_Comm_create(simulation%pf_e4d_scatter_comm, &
                         simulation%pf_e4d_master_grp, &
                         simulation%pf_e4d_master_comm,ierr)
!print *, 7, simulation%myrank_save, simulation%pf_e4d_master_comm
    if (simulation%pf_e4d_master_comm /= MPI_COMM_NULL) then
      call MPI_Comm_rank(simulation%pf_e4d_master_comm, &
                         simulation%pf_e4d_master_rank,ierr)
      call MPI_Comm_size(simulation%pf_e4d_master_comm, &
                         simulation%pf_e4d_master_size,ierr)
!      PFE4D_MASTER_COMM = simulation%pf_e4d_master_comm
!print *, 8, simulation%myrank_save, simulation%pf_e4d_master_rank, simulation%pf_e4d_master_size
    endif
  endif

!print *, 9, simulation%myrank_save, simulation%pf_e4d_scatter_comm, simulation%pf_e4d_master_comm, MPI_COMM_NULL
  
  if (simulation%pflotran_process) then
    call HydrogeophysicsInitPostPetsc(simulation)
  else
    option%io_rank = -1 ! turn off I/O from E4D processes.
  endif
#endif

!#define DEBUG

  !   PFLOTRAN subsurface processes      E4D master process
  if (simulation%pflotran_process .or. option%myrank == 0) then 
    if (simulation%pflotran_process) then
      simulation%hydrogeophysics_coupler%pf_to_e4d_master_comm = &
        simulation%pf_e4d_master_comm
    endif

    ! create mpi Vec that includes all PFLOTRAN processes and the E4D master
    call VecCreate(simulation%pf_e4d_scatter_comm,pflotran_tracer_vec_mpi, &
                   ierr);CHKERRQ(ierr)
    if (simulation%pflotran_process) then
      call VecGetLocalSize(simulation%realization%field%work,local_size, &
                           ierr);CHKERRQ(ierr)
    else ! E4D master process
      local_size = 0
    endif
    call VecSetSizes(pflotran_tracer_vec_mpi,local_size,PETSC_DECIDE, &
                     ierr);CHKERRQ(ierr)
    call VecSetFromOptions(pflotran_tracer_vec_mpi,ierr);CHKERRQ(ierr)

    allocate(int_array(local_size))
    if (simulation%pflotran_process) then
      int_array = 0
      do i = 1, local_size
        int_array(i) =  simulation%realization%patch%grid%nG2A( &
                          simulation%realization%patch%grid%nL2G(i))

      enddo
      int_array = int_array - 1
    endif
    call ISCreateGeneral(simulation%pf_e4d_scatter_comm,local_size,int_array, &
!    call ISCreateGeneral(PETSC_COMM_SELF,local_size,int_array, &
                         PETSC_COPY_VALUES,is_natural,ierr);CHKERRQ(ierr)
    deallocate(int_array)
    if (simulation%pflotran_process) then
      call VecGetOwnershipRange(simulation%realization%field%work,istart, &
                                PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
    endif
!    istart = 0
    call ISCreateStride(simulation%pf_e4d_scatter_comm,local_size,istart,1, &
                        is_petsc,ierr);CHKERRQ(ierr)

#ifdef DEBUG
  string = 'is_petsc.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  string = 'is_natural.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
!  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call ISView(is_natural,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif


    ! create seq Vec on each (all PFLOTRAN processes and the E4D master)
    call VecGetSize(pflotran_tracer_vec_mpi,local_size,ierr);CHKERRQ(ierr)
    if (simulation%pflotran_process) local_size = 0
    ! make global size the local size of pflotran_tracer_vec_seq on E4D master
!    call VecCreate(PETSC_COMM_SELF,pflotran_tracer_vec_seq,ierr)
    call VecCreate(simulation%pf_e4d_scatter_comm,pflotran_tracer_vec_seq, &
                   ierr);CHKERRQ(ierr)
    call VecSetSizes(pflotran_tracer_vec_seq,local_size,PETSC_DECIDE, &
                     ierr);CHKERRQ(ierr)
    call VecSetFromOptions(pflotran_tracer_vec_seq,ierr);CHKERRQ(ierr)

#ifdef DEBUG
    string = 'pflotran_tracer_vec_mpi.txt'
    call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(pflotran_tracer_vec_mpi,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    string = 'pflotran_tracer_vec_seq' // trim(StringFormatInt(option%global_rank)) // '.txt'
  !  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
    call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(pflotran_tracer_vec_seq,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
    ! create scatter between mpi Vec and local seq Vecs (only E4D master is 
    ! relevant)
    call VecScatterCreate(pflotran_tracer_vec_mpi,is_petsc, &
                          pflotran_tracer_vec_seq,is_natural, &
                          pflotran_scatter,ierr);CHKERRQ(ierr)
#ifdef DEBUG
    string = 'scatter.txt'
    call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer, &
                              ierr);CHKERRQ(ierr)
    call VecScatterView(pflotran_scatter,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
    call ISDestroy(is_natural,ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)

    ! duplicate tracer vecs into saturation version
    call VecDuplicate(pflotran_tracer_vec_mpi, &
                      pflotran_saturation_vec_mpi,ierr);CHKERRQ(ierr)
    call VecDuplicate(pflotran_tracer_vec_seq, &
                      pflotran_saturation_vec_seq,ierr);CHKERRQ(ierr)
    ! have to broadcast the flow mode to the e4d master
    temp_int = option%iflowmode
    if (simulation%pf_e4d_master_comm /= MPI_COMM_NULL) then
      call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER, &
                     ZERO_INTEGER_MPI,simulation%pf_e4d_master_comm,ierr)
    endif
    if (temp_int /= NULL_MODE .and. temp_int /= RICHARDS_MODE) then
      call VecDuplicate(pflotran_tracer_vec_mpi, &
                        pflotran_temperature_vec_mpi,ierr);CHKERRQ(ierr)
      call VecDuplicate(pflotran_tracer_vec_seq, &
                        pflotran_temperature_vec_seq,ierr);CHKERRQ(ierr)
    endif
  endif
!print *, 'End  -----------'

  if (simulation%pflotran_process) then
    ! tracer
    simulation%hydrogeophysics_coupler%tracer_mpi = pflotran_tracer_vec_mpi
    simulation%hydrogeophysics_coupler%tracer_seq = pflotran_tracer_vec_seq
    simulation%tracer_mpi = pflotran_tracer_vec_mpi
    ! saturation
    simulation%hydrogeophysics_coupler%saturation_mpi = &
      pflotran_saturation_vec_mpi
    simulation%hydrogeophysics_coupler%saturation_seq = &
      pflotran_saturation_vec_seq
    simulation%saturation_mpi = pflotran_saturation_vec_mpi
    ! temperature
    simulation%hydrogeophysics_coupler%temperature_mpi = &
      pflotran_temperature_vec_mpi
    simulation%hydrogeophysics_coupler%temperature_seq = &
      pflotran_temperature_vec_seq
    simulation%temperature_mpi = pflotran_temperature_vec_mpi
    simulation%hydrogeophysics_coupler%pf_to_e4d_scatter = pflotran_scatter
    simulation%hydrogeophysics_coupler%pf_to_e4d_master_comm = &
      simulation%pf_e4d_master_comm
  else
    call HydrogeophysicsWrapperInit(option, &
                                    pflotran_tracer_vec_mpi, &
                                    pflotran_tracer_vec_seq, &
                                    pflotran_saturation_vec_mpi, &
                                    pflotran_saturation_vec_seq, &
                                    pflotran_temperature_vec_mpi, &
                                    pflotran_temperature_vec_seq, &
                                    pflotran_scatter, &
                                    simulation%pf_e4d_master_comm)
  endif
#undef DEBUG

end subroutine HydrogeophysicsInitialize

! ************************************************************************** !

subroutine HydrogeophysicsInitPostPetsc(simulation)
  ! 
  ! HydrogeophysicsInitializePostPetsc: Sets up hydrogeophysics simulation
  ! framework after to PETSc initialization
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  use Factory_Subsurface_module
  use PMC_Hydrogeophysics_class
  use Option_module
  use Logging_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(simulation_hydrogeophysics_type) :: simulation
  
  class(pmc_hydrogeophysics_type), pointer :: hydrogeophysics_coupler
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  ! Init() is called in SubsurfaceInitializePostPetsc
  call SubsurfaceInitializePostPetsc(simulation)
  
  ! add hydrogeophysics coupler to list
  hydrogeophysics_coupler => PMCHydrogeophysicsCreate()
  hydrogeophysics_coupler%option => simulation%option
  hydrogeophysics_coupler%realization => simulation%realization
  simulation%hydrogeophysics_coupler => hydrogeophysics_coupler
  ! set up logging stage
  string = 'Hydrogeophysics'
  call LoggingCreateStage(string,hydrogeophysics_coupler%stage)
  if (associated(simulation%process_model_coupler_list%child)) then
    simulation%process_model_coupler_list%child%child => hydrogeophysics_coupler
  else
    simulation%process_model_coupler_list%child => hydrogeophysics_coupler
  endif

end subroutine HydrogeophysicsInitPostPetsc

! ************************************************************************** !

subroutine HydrogeoInitCommandLineSettings(option)
  ! 
  ! Initializes hydrogeophysics settings
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/17/13
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
!  string = '-dummy'
!  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
!  if (option_found) then
!    option%subsurface_simulation_type = dummy
!  endif
  
end subroutine HydrogeoInitCommandLineSettings

end module Factory_Hydrogeophysics_module
