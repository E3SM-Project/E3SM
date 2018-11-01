module Time_Storage_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
  
  type, public :: time_storage_type
    PetscReal, pointer :: times(:)
    PetscReal :: cur_time
    PetscReal :: cur_time_fraction
    PetscInt :: cur_time_index
    PetscInt :: max_time_index
    PetscBool :: is_cyclic
    PetscReal :: time_shift    ! shift for cyclic data sets 
    PetscBool :: cur_time_index_changed
    PetscBool :: cur_time_fraction_changed
    PetscInt :: time_interpolation_method
    PetscBool :: force_update
  end type time_storage_type
  
  public :: TimeStorageCreate, &
            TimeStorageGetTimes, &
            TimeStorageVerify, &
            TimeStorageUpdate, &
            TimeStoragePrint, &
            TimeStorageDestroy

contains

! ************************************************************************** !

function TimeStorageCreate()
  ! 
  ! Initializes a time storage
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 
  use petscsys
  implicit none
  
  type(time_storage_type), pointer :: time_storage
  type(time_storage_type), pointer :: TimeStorageCreate

  allocate(time_storage)
  nullify(time_storage%times)
  time_storage%cur_time = 0.d0
  time_storage%cur_time_fraction = 0.d0
  time_storage%cur_time_index = 0
  time_storage%max_time_index = 0
  time_storage%is_cyclic = PETSC_FALSE
  time_storage%time_shift = 0.d0
  time_storage%cur_time_index_changed = PETSC_FALSE
  time_storage%cur_time_fraction_changed = PETSC_FALSE
  time_storage%time_interpolation_method = INTERPOLATION_NULL
  time_storage%force_update = PETSC_FALSE
  
  TimeStorageCreate => time_storage
    
end function TimeStorageCreate

! ************************************************************************** !

subroutine TimeStorageVerify(default_time, time_storage, &
                             default_time_storage, header, option)
  ! 
  ! Verifies the data in a time storage
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 

  use petscsys
  use Option_module

  implicit none
  
  PetscReal :: default_time
  type(time_storage_type), pointer :: time_storage
  type(time_storage_type), pointer :: default_time_storage
  character(len=MAXSTRINGLENGTH) :: header
  type(option_type) :: option
  
  PetscInt :: array_size
  
  if (.not.associated(time_storage)) return
  
  if (associated(default_time_storage)) then
    if (default_time_storage%is_cyclic) time_storage%is_cyclic = PETSC_TRUE
    if (time_storage%time_interpolation_method == INTERPOLATION_NULL) then
      time_storage%time_interpolation_method = &
        default_time_storage%time_interpolation_method
    endif
  endif
  
  if (time_storage%time_interpolation_method == INTERPOLATION_NULL) then
    option%io_buffer = 'Time interpolation method must be specified in ' // &
      trim(header) // '.'
    call printErrMsg(option)
  endif
  
  time_storage%max_time_index = 1
  if (.not.associated(time_storage%times)) then
    if (associated(default_time_storage)) then
      if (.not.associated(default_time_storage%times)) then
        array_size = 1
        allocate(time_storage%times(array_size))
        time_storage%times = default_time
      else
        array_size = size(default_time_storage%times,1)
        allocate(time_storage%times(array_size))
        time_storage%times(1:array_size) = &
          default_time_storage%times(1:array_size)
      endif
    else
      array_size = 1
      allocate(time_storage%times(array_size))
      time_storage%times = default_time
    endif
  endif
  time_storage%max_time_index = size(time_storage%times,1) 
  time_storage%cur_time_index = 1
  
  time_storage%time_shift = time_storage%times(time_storage%max_time_index)

end subroutine TimeStorageVerify

! ************************************************************************** !

subroutine TimeStorageGetTimes(time_storage, option, max_sim_time, time_array)
  ! 
  ! Fills an array of times based on time storage
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 

  use Option_module

  implicit none
  
  type(time_storage_type), pointer :: time_storage
  type(option_type) :: option
  PetscReal :: max_sim_time
  PetscReal, pointer :: time_array(:)
  
  PetscInt :: num_times
  PetscInt :: itime
  PetscReal :: time_shift
  PetscReal, allocatable :: temp_times(:)

  if (.not.associated(time_storage)) then
    nullify(time_array)
    return
  endif
  
  if (.not.time_storage%is_cyclic .or. time_storage%max_time_index == 1) then
    allocate(time_array(time_storage%max_time_index))
    time_array =  time_storage%times
  else ! cyclic
    num_times = (int(max_sim_time / &
                     time_storage%times(time_storage%max_time_index))+1)* &
                time_storage%max_time_index
    allocate(temp_times(num_times))
    temp_times = 0.d0

    num_times = 0
    itime = 0
    time_shift = 0.d0
    do
      num_times = num_times + 1
      itime = itime + 1
      ! exit for non-cyclic - but is will never enter conditional given 
      ! conditional above.
      if (itime > time_storage%max_time_index) exit
      temp_times(num_times) = time_storage%times(itime) + time_shift
      if (mod(itime,time_storage%max_time_index) == 0) then
        itime = 0
        time_shift = time_shift + time_storage%times(time_storage%max_time_index) 
      endif 
      ! exit for cyclic
      if (temp_times(num_times) >= max_sim_time) exit
    enddo

    allocate(time_array(num_times))
    time_array(:) = temp_times(1:num_times)
    deallocate(temp_times)
  endif
 
end subroutine TimeStorageGetTimes

! ************************************************************************** !

subroutine TimeStoragePrint(time_storage,option)
  ! 
  ! Prints time storage info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 

  use Option_module

  implicit none
  
  type(time_storage_type) :: time_storage
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  write(option%fid_out,'(8x,''Time Storage'')')
  if (time_storage%is_cyclic) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(8x,''Is cyclic: '',a)') trim(string)
  if (size(time_storage%times) > 1) then  
    write(option%fid_out,'(8x,''  Number of values: '', i7)') &
      time_storage%max_time_index
    write(option%fid_out,'(8x,''Start value:'',es16.8)') &
      time_storage%times(1)
    write(option%fid_out,'(8x,''End value:'',es16.8)') &
      time_storage%times(time_storage%max_time_index)
  else
    write(option%fid_out,'(8x,''Value:'',es16.8)') time_storage%times(1)
  endif

            
end subroutine TimeStoragePrint

! ************************************************************************** !

subroutine TimeStorageUpdate(time_storage)
  ! 
  ! Updates a time storage
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 

  use petscsys
  use Option_module
  
  implicit none
  
  type(time_storage_type) :: time_storage
  
  PetscInt :: irank
  PetscInt :: cur_time_index
  PetscInt :: next_time_index
  
  ! cycle times if at max_time_index and cyclic
  if (time_storage%cur_time_index == time_storage%max_time_index .and. &
      time_storage%is_cyclic .and. time_storage%max_time_index > 1) then
    do cur_time_index = 1, time_storage%max_time_index
      time_storage%times(cur_time_index) = &
        time_storage%times(cur_time_index) + time_storage%time_shift
    enddo
    time_storage%cur_time_index = 1
  endif
 
  cur_time_index = time_storage%cur_time_index
  next_time_index = min(time_storage%cur_time_index+1, &
                        time_storage%max_time_index)

  ! initialize to no change
  time_storage%cur_time_index_changed = PETSC_FALSE
  ! find appropriate time interval
  do
    if (time_storage%cur_time < time_storage%times(next_time_index) .or. &
        cur_time_index == next_time_index) &
      exit
    
    if (cur_time_index /= next_time_index) &
      ! toggle flag indicating a change in index
      time_storage%cur_time_index_changed = PETSC_TRUE
    cur_time_index = next_time_index
    ! ensure that time index does not go beyond end of array
    if (next_time_index < time_storage%max_time_index) then
      next_time_index = next_time_index + 1
    ! this conditional enable the code to find the correct
    ! time index for a cyclic dataset
    else if (time_storage%is_cyclic .and. time_storage%max_time_index > 1) then
      do cur_time_index = 1, time_storage%max_time_index
        time_storage%times(cur_time_index) = &
          time_storage%times(cur_time_index) + time_storage%time_shift
      enddo
      cur_time_index = 1
      next_time_index = 2
    endif
  enddo
    
  time_storage%cur_time_index = cur_time_index
  if (cur_time_index < 1) then
    return
  else if (cur_time_index < time_storage%max_time_index) then
    time_storage%cur_time_fraction_changed = PETSC_TRUE
    ! fraction = (t-t1)/(t2-t1)
    time_storage%cur_time_fraction = (time_storage%cur_time- &
                                      time_storage%times(cur_time_index)) / &
                                     (time_storage%times(next_time_index) - &
                                      time_storage%times(cur_time_index))
  else
    if (dabs(time_storage%cur_time_fraction - 1.d0) < 1.d-10) then
      ! essentially zero change
      time_storage%cur_time_fraction_changed = PETSC_FALSE
    else
      time_storage%cur_time_fraction_changed = PETSC_TRUE
      time_storage%cur_time_fraction = 1.d0
    endif
  endif

  if (time_storage%force_update) then
    time_storage%cur_time_fraction_changed = PETSC_TRUE
    time_storage%cur_time_index_changed = PETSC_TRUE
    time_storage%force_update = PETSC_FALSE
  endif

end subroutine TimeStorageUpdate

! ************************************************************************** !

subroutine TimeStorageDestroy(time_storage)
  ! 
  ! Destroys a time storage associated with a sub_condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/03/13
  ! 

  implicit none
  
  type(time_storage_type), pointer :: time_storage
  
  if (.not.associated(time_storage)) return
  
  if (associated(time_storage%times)) deallocate(time_storage%times)
  nullify(time_storage%times)
  
  deallocate(time_storage)
  nullify(time_storage)

end subroutine TimeStorageDestroy

end module Time_Storage_module
