module Multi_Simulation_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  
  type, public :: multi_simulation_type
    PetscInt :: num_groups
    PetscInt :: num_realizations
    PetscInt :: num_local_realizations
    PetscInt :: cur_realization
    PetscInt, pointer :: realization_ids(:)
  end type multi_simulation_type

  public :: MultiSimulationCreate, &
            MultiSimulationInitialize, &
            MultiSimulationIncrement, &
            MultiSimulationDone, &
            MultiSimulationDestroy
  
contains

! ************************************************************************** !

function MultiSimulationCreate()
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14
  implicit none
  
  type(multi_simulation_type), pointer :: MultiSimulationCreate
  
  type(multi_simulation_type), pointer :: multisimulation
  
  allocate(multisimulation)
  multisimulation%num_realizations = 0
  multisimulation%num_groups = 0
  multisimulation%num_local_realizations = 0
  multisimulation%cur_realization = 0
  MultiSimulationCreate => multisimulation

end function MultiSimulationCreate

! ************************************************************************** !

subroutine MultiSimulationInitialize(multisimulation,option)
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14
  use Option_module
  use Input_Aux_module
  
  implicit none

  type(multi_simulation_type) :: multisimulation
  type(option_type) :: option

  PetscInt :: i
  PetscInt :: offset, delta, remainder

  PetscInt :: realization_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscInt, pointer :: realization_ids_from_file(:)
  character(len=MAXSTRINGLENGTH) :: filename
  type(input_type), pointer :: input
  PetscErrorCode :: ierr
  
  ! query user for number of communicator groups and realizations
  string = '-num_groups'
  call InputGetCommandLineInt(string,multisimulation%num_groups, &
                              option_found,option)

  string = '-num_realizations'
  call InputGetCommandLineInt(string,multisimulation%num_realizations, &
                              option_found,option)

  ! read realization ids from a file - contributed by Xingyuan
  string = '-realization_ids_file'
  call InputGetCommandLineString(string,filename,option_found,option)
  if (option_found) then
    input => InputCreate(IUNIT_TEMP,filename,option)
    if (multisimulation%num_realizations == 0) then
      option%io_buffer = '"-num_realizations <int>" must be specified ' // &
        'with an integer value matching the number of ids in ' // &
        '"-realization_ids_file <string>".'
      call printErrMsg(option)
    endif
    allocate(realization_ids_from_file(multisimulation%num_realizations))
    realization_ids_from_file = 0
    string = &
      '# of realization ids read from file may be too few in StochasticInit()'
    do i = 1, multisimulation%num_realizations
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,string)
      call InputReadInt(input,option,realization_ids_from_file(i))
      call InputErrorMsg(input,option,'realization id', &
                         'MultiSimulationInitialize')
    enddo
    call InputDestroy(input)
  else
    nullify(realization_ids_from_file)
  endif
    
  ! Realization offset contributed by Xingyuan.  This allows one to specify the
  ! smallest/lowest realization id (other than zero) in a stochastic simulation
  string = '-realization_offset'
  call InputGetCommandLineInt(string,offset,option_found,option)
  if (.not.option_found) then
    offset = 0
  endif

  ! error checking
  if (multisimulation%num_groups == 0) then
    option%io_buffer = 'Number of stochastic processor groups not ' // &
                       'initialized. Setting to 1.'
    call printWrnMsg(option)
    multisimulation%num_groups = 1
  endif
  if (multisimulation%num_realizations == 0) then
    option%io_buffer = 'Number of stochastic realizations not ' // &
                       'initialized. Setting to 1.'
    call printWrnMsg(option)
    multisimulation%num_realizations = 1
  endif
  
  call OptionCreateProcessorGroups(option,multisimulation%num_groups)
  
  ! divvy up the realizations
  multisimulation%num_local_realizations = multisimulation%num_realizations / &
                                           multisimulation%num_groups
  remainder = multisimulation%num_realizations - multisimulation%num_groups * &
                                         multisimulation%num_local_realizations
  
  ! offset is initialized above after check for '-realization_offset'
  do i = 1, option%mygroup_id-1
    delta = multisimulation%num_local_realizations
    if (i < remainder) delta = delta + 1
    offset = offset + delta
  enddo
  
  if (option%mygroup_id < remainder) multisimulation%num_local_realizations = &
                                     multisimulation%num_local_realizations + 1
  allocate(multisimulation%realization_ids( &
                                  multisimulation%num_local_realizations))
  multisimulation%realization_ids = 0
  do i = 1, multisimulation%num_local_realizations
    multisimulation%realization_ids(i) = offset + i
  enddo
  
  ! map ids from file - contributed by Xingyuan
  if (associated(realization_ids_from_file)) then
    do i = 1, multisimulation%num_local_realizations
      multisimulation%realization_ids(i) = &
        realization_ids_from_file(multisimulation%realization_ids(i))
    enddo
  endif

end subroutine MultiSimulationInitialize

! ************************************************************************** !

subroutine MultiSimulationIncrement(multisimulation,option)
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14
  use Option_module
  
  implicit none
  
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (.not.associated(multisimulation)) return
  
!#define RESTART_MULTISIMULATION  
#if defined(RESTART_MULTISIMULATION)
  ! code for restarting multisimulation runs; may no longer need this.
  do
#endif

    call OptionInitRealization(option)

    multisimulation%cur_realization = multisimulation%cur_realization + 1
    ! Set group prefix based on id
    option%id = &
      multisimulation%realization_ids(multisimulation%cur_realization)
    write(string,'(i6)') option%id
    option%group_prefix = 'R' // trim(adjustl(string))  

#if defined(RESTART_MULTISIMULATION)
    string = 'restart' // trim(adjustl(option%group_prefix)) // '.chk.info'
    open(unit=86,file=string,status="old",iostat=status)
    ! if file found, cycle
    if (status == 0) then
      close(86)
      cycle
    else
      exit
    endif
  enddo
#endif

end subroutine MultiSimulationIncrement

! ************************************************************************** !

function MultiSimulationDone(multisimulation)
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14
  implicit none
  
  type(multi_simulation_type), pointer :: multisimulation
  
  PetscBool :: MultiSimulationDone

  MultiSimulationDone = PETSC_FALSE
  if (.not.associated(multisimulation)) then
    MultiSimulationDone = PETSC_TRUE
  else if (multisimulation%cur_realization >= &
           multisimulation%num_local_realizations) then
    MultiSimulationDone = PETSC_TRUE
  endif
  
end function MultiSimulationDone

! ************************************************************************** !

subroutine MultiSimulationDestroy(multisimulation)
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14
  implicit none
  
  type(multi_simulation_type), pointer :: multisimulation
  
  if (.not.associated(multisimulation)) return
  
  if (associated(multisimulation%realization_ids)) &
    deallocate(multisimulation%realization_ids)
  nullify(multisimulation%realization_ids)
  
  deallocate(multisimulation)
  nullify(multisimulation)

end subroutine MultiSimulationDestroy

end module Multi_Simulation_module
