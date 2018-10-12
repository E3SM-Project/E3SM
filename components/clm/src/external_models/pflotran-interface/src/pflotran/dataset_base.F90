module Dataset_Base_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Time_Storage_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: dataset_base_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    type(time_storage_type), pointer :: time_storage ! stores transient times
    character(len=MAXSTRINGLENGTH) :: header
    PetscInt :: rank  ! size of dims(:)
    PetscInt, pointer :: dims(:)    ! dimensions of arrays (excludes time)
    PetscInt :: data_type
    PetscInt, pointer :: iarray(:)
    PetscInt, pointer :: ibuffer(:)
    PetscReal, pointer :: rarray(:)
    PetscReal, pointer :: rbuffer(:)
    PetscInt :: buffer_slice_offset ! index of the first time slice in the buffer
    PetscInt :: buffer_nslice ! # of time slices stored in buffer
    class(dataset_base_type), pointer :: next
  end type dataset_base_type

  ! dataset types
  PetscInt, public, parameter :: DATASET_INTEGER = 1
  PetscInt, public, parameter :: DATASET_REAL = 2
  
  public :: DatasetBaseCreate, &
            DatasetBaseInit, &
            DatasetBaseCopy, &
            DatasetBaseVerify, &
            DatasetBaseInterpolateTime, &
            DatasetBaseReorder, &
            DatasetBaseGetPointer, &
            DatasetBaseAddToList, &
            DatasetBaseGetNameInfo, &
            DatasetBasePrint, &
            DatasetBaseStrip, &
            DatasetBaseDestroy
contains

! ************************************************************************** !

function DatasetBaseCreate()
  ! 
  ! Creates members of base database class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
  implicit none
  
  class(dataset_base_type), pointer :: dataset

  class(dataset_base_type), pointer :: DatasetBaseCreate
  
  allocate(dataset)
  call DatasetBaseInit(dataset)

  DatasetBaseCreate => dataset
    
end function DatasetBaseCreate

! ************************************************************************** !

subroutine DatasetBaseInit(this)
  ! 
  ! Initializes members of base database class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
  implicit none
  
  class(dataset_base_type) :: this
  
  this%name = ''
  this%filename = ''
  this%header = ''
  this%rank = 0
  this%data_type = 0
  nullify(this%time_storage)
  nullify(this%dims)
  nullify(this%iarray)
  nullify(this%ibuffer)
  nullify(this%rarray)
  nullify(this%rbuffer)
  this%buffer_slice_offset = 0
  this%buffer_nslice = 0
  nullify(this%next)
    
end subroutine DatasetBaseInit

! ************************************************************************** !

subroutine DatasetBaseCopy(this, that)
  ! 
  ! Copies members of base database class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
  implicit none
  
  class(dataset_base_type) :: this
  class(dataset_base_type) :: that
  
  that%name = this%name
  that%filename = this%filename
  that%rank = this%rank
  that%data_type = this%data_type
  if (associated(this%time_storage)) then
    that%time_storage => this%time_storage
    nullify(this%time_storage)
  endif
  if (associated(this%dims)) then
    that%dims => this%dims
    nullify(this%dims)
  endif
  if (associated(this%iarray)) then
    that%iarray => this%iarray
    nullify(this%iarray)
  endif
  if (associated(this%ibuffer)) then
    that%ibuffer => this%ibuffer
    nullify(this%ibuffer)
  endif
  if (associated(this%rarray)) then
    that%rarray => this%rarray
    nullify(this%rarray)
  endif
  if (associated(this%rbuffer)) then
    that%rbuffer => this%rbuffer
    nullify(this%rbuffer)
  endif
  that%buffer_slice_offset = this%buffer_slice_offset
  that%buffer_nslice = this%buffer_nslice
  that%next => this%next
  nullify(this%next)
    
end subroutine DatasetBaseCopy

! ************************************************************************** !

subroutine DatasetBaseVerify(this,dataset_error,option)
  ! 
  ! Verifies that data structure is properly set up.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/08/13
  ! 
  
  use Option_module
  
  implicit none
  
  class(dataset_base_type) :: this
  PetscBool :: dataset_error
  type(option_type) :: option

  if (associated(this%time_storage) .or. associated(this%iarray) .or. &
      associated(this%iarray) .or. associated(this%ibuffer) .or. &
      associated(this%rarray) .or. associated(this%rbuffer)) then
    if (len_trim(this%name) < 1) then
      this%name = 'Unnamed Dataset'
    endif
    if (this%data_type == 0) then
      option%io_buffer = '"data_type" not defined in dataset: ' // &
                         trim(this%name)
      call printMsg(option)
      dataset_error = PETSC_TRUE
    endif
    if (associated(this%ibuffer) .or. associated(this%rbuffer)) then
      if (.not.associated(this%dims)) then
        option%io_buffer = '"dims" not allocated in dataset: ' // &
                           trim(this%name)
        call printMsg(option)
        dataset_error = PETSC_TRUE
      endif
      if (this%dims(1) == 0) then
        option%io_buffer = '"dims" not defined in dataset: ' // &
                           trim(this%name)
        call printMsg(option)
        dataset_error = PETSC_TRUE
      endif
      if (size(this%dims) /= this%rank) then
        option%io_buffer = 'Size of "dims" not match "rank" in dataset: ' // &
                            trim(this%name)
        dataset_error = PETSC_TRUE
        call printMsg(option)
      endif
    endif
    if (associated(this%ibuffer) .and. .not.associated(this%iarray)) then
      allocate(this%iarray(this%dims(1)))
      this%iarray = 0
    endif
    if (associated(this%rbuffer) .and. .not.associated(this%rarray)) then
      allocate(this%rarray(this%dims(1)))
      this%rarray = 0.d0
    endif
  else if (len_trim(this%name) < 1) then
    option%io_buffer = 'ERROR: No value or dataset specified.'
    call printMsg(option)
    dataset_error = PETSC_TRUE
  endif
    
end subroutine DatasetBaseVerify

! ************************************************************************** !

subroutine DatasetBaseInterpolateTime(this)
  ! 
  ! Interpolates dataset between two buffer times
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  
  PetscInt :: array_size,i
  PetscInt :: time_interpolation_method
  PetscReal :: weight2
  PetscInt :: time1_start, time1_end, time2_start, time2_end
  
  if (.not.associated(this%rbuffer)) return
  
  time_interpolation_method = this%time_storage%time_interpolation_method
  
  if (this%time_storage%cur_time_index >= &
      this%time_storage%max_time_index) then
    ! dataset has reached the end of the time array and is not cyclic.
    ! no interpolation needed
    if (this%time_storage%cur_time_index_changed) then
      array_size = size(this%rarray)
      time1_end = array_size*(this%time_storage%max_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      this%rarray = this%rbuffer(time1_start:time1_end)
    endif
    return
  endif
  
  array_size = size(this%rarray)
  select case(time_interpolation_method)
    case(INTERPOLATION_NULL)
      ! it possible that the time_interpolation_method is null during
      ! the load process.  in that case, we need to set the array to an
      ! uninitialized value (UNINITIALIZED_DOUBLE).
      this%rarray = UNINITIALIZED_DOUBLE
    case(INTERPOLATION_STEP)
      ! if time index has not changed skip
      if (.not.this%time_storage%cur_time_index_changed) return
      time1_end = array_size*(this%time_storage%cur_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      this%rarray = this%rbuffer(time1_start:time1_end)
    case(INTERPOLATION_LINEAR)
      ! if fraction has not changed skip
      if (.not.this%time_storage%cur_time_fraction_changed) return
      time1_end = array_size*(this%time_storage%cur_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      time2_end = time1_end + array_size
      time2_start = time1_start + array_size
      this%rarray = (1.d0-this%time_storage%cur_time_fraction) * &
                        this%rbuffer(time1_start:time1_end) + &
                        this%time_storage%cur_time_fraction * &
                        this%rbuffer(time2_start:time2_end)
  end select

end subroutine DatasetBaseInterpolateTime

! ************************************************************************** !

subroutine DatasetBaseInterpolateSpace(this,xx,yy,zz,time,real_value,option)
  ! 
  ! Interpolates data from the dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  class(dataset_base_type) :: this
  PetscReal, intent(in) :: xx, yy, zz
  PetscReal :: time
  PetscReal :: real_value
  type(option_type) :: option
  
end subroutine DatasetBaseInterpolateSpace

! ************************************************************************** !

subroutine DatasetBaseReorder(this,option)
  ! 
  ! If a dataset is loaded from an HDF5 file, and it was
  ! multidimensional in the HDF5 file, the array needs to be
  ! reordered fro Fortran -> C indexing.  This subroutine
  ! takes care of the reordering.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module
  
  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  
  PetscReal, allocatable :: temp_real(:)
  PetscInt :: i, j, k, l
  PetscInt :: length_1d
  PetscInt :: dims(4), n1, n1Xn2, n1Xn2Xn3
  PetscInt :: count, index
  PetscReal, pointer :: rarray(:)
  
  if (this%data_type == DATASET_INTEGER) then
    option%io_buffer = 'Reordering of integer data sets not yet supported.'
    call printErrMsg(option)
  endif
  
  ! set each dim to 1 by default for loop below
  dims(:) = 1
  dims(1:this%rank) = this%dims(1:this%rank)
  if (associated(this%rbuffer)) then
    rarray => this%rbuffer
    dims(this%rank+1) = this%buffer_nslice
  else
    rarray => this%rarray
  endif
  
  ! Not necessary for 1D arrays
  if (maxval(dims(2:)) == 1) return
  
  length_1d = 1
  do i = 1, size(dims)
    length_1d = length_1d*dims(i)
  enddo
  allocate(temp_real(length_1d))
  
  n1 = dims(1)
  n1Xn2 = n1*dims(2)
  n1Xn2Xn3 = n1Xn2*dims(3)
  count = 0
  do i = 1, dims(1)
    do j = 0, dims(2)-1
      do k = 0, dims(3)-1
        do l = 0, dims(4)-1
          index = l*n1Xn2Xn3+k*n1Xn2+j*n1+i
          count = count+1
          temp_real(index) = rarray(count)
        enddo
      enddo
    enddo
  enddo  

  !geh: had to add 1:length_1d since in some cases, rarray is larger
  !     than the specified size based on "dims". not sure why....
  rarray(1:length_1d) = temp_real
  deallocate(temp_real)
  
end subroutine DatasetBaseReorder

! ************************************************************************** !

subroutine DatasetBaseGetTimes(this, option, max_sim_time, time_array)
  ! 
  ! Fills an array of times based on a dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  PetscReal :: max_sim_time
  PetscReal, pointer :: time_array(:)
  
  
  if (associated(this%time_storage)) then
    call TimeStorageGetTimes(this%time_storage, option, max_sim_time, &
                             time_array)
  endif
 
end subroutine DatasetBaseGetTimes

! ************************************************************************** !

subroutine DatasetBaseAddToList(dataset,list)
  ! 
  ! Adds a dataset to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  implicit none
  
  class(dataset_base_type), pointer :: dataset
  class(dataset_base_type), pointer :: list

  class(dataset_base_type), pointer :: cur_dataset
  
  if (associated(list)) then
    cur_dataset => list
    ! loop to end of list
    do
      if (.not.associated(cur_dataset%next)) exit
      cur_dataset => cur_dataset%next
    enddo
    cur_dataset%next => dataset
  else
    list => dataset
  endif
  
end subroutine DatasetBaseAddToList

! ************************************************************************** !

function DatasetBaseGetPointer(dataset_list, dataset_name, debug_string, &
                               option)
  ! 
  ! Returns the pointer to the dataset named "name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  class(dataset_base_type), pointer :: dataset_list
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: debug_string
  type(option_type) :: option

  class(dataset_base_type), pointer :: DatasetBaseGetPointer
  PetscBool :: found
  class(dataset_base_type), pointer :: cur_dataset

  found = PETSC_FALSE
  cur_dataset => dataset_list
  do 
    if (.not.associated(cur_dataset)) exit
    if (StringCompare(dataset_name, &
                      cur_dataset%name,MAXWORDLENGTH)) then
      found = PETSC_TRUE
      DatasetBaseGetPointer => cur_dataset
      return
    endif
    cur_dataset => cur_dataset%next
  enddo
  if (.not.found) then
    option%io_buffer = 'Dataset "' // trim(dataset_name) // '" in "' // &
             trim(debug_string) // '" not found among available datasets.'
    call printErrMsgByRank(option)    
  endif

end function DatasetBaseGetPointer

! ************************************************************************** !

function DatasetBaseGetNameInfo(this)
  ! 
  ! Returns naming information for dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/18
  ! 
  implicit none

  class(dataset_base_type) :: this

  character(len=MAXSTRINGLENGTH) :: DatasetBaseGetNameInfo

  character(len=MAXSTRINGLENGTH) :: string

  string = ''
  if (len_trim(this%name) > 0) then
    string = 'NAME: "' // trim(this%name) // '"'
  endif
  if (len_trim(this%filename) > 0) then
    string = trim(string) // ' FILENAME: "' // trim(this%filename) // '"'
  endif
  DatasetBaseGetNameInfo = string

end function DatasetBaseGetNameInfo

! ************************************************************************** !

subroutine DatasetBasePrint(this,option)
  ! 
  ! Prints dataset info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/13
  ! 

  use Option_module

  implicit none

  class(dataset_base_type) :: this
  type(option_type) :: option

  if (len_trim(this%filename) > 0) then
    write(option%fid_out,'(10x,''Filename: '',a)') trim(this%filename)
  endif
  if (associated(this%time_storage)) then
    write(option%fid_out,'(10x,''Is transient?: yes'')')
    write(option%fid_out,'(10x,''Number of times: '',i6)') &
      this%time_storage%max_time_index
    if (this%time_storage%is_cyclic) then
      write(option%fid_out,'(10x,''Is cyclic?: yes'')')
    else
      write(option%fid_out,'(10x,''Is cyclic?: no'')')
    endif
  else
    write(option%fid_out,'(10x,''Transient: no'')')
  endif
  if (associated(this%rbuffer)) then
    write(option%fid_out,'(10x,''Buffer:'')')
    write(option%fid_out,'(12x,''Rank: '',i2)') this%rank
    if (associated(this%dims)) then
      write(option%fid_out,'(12x,''Dims: '',10i4)') this%dims
    endif
    write(option%fid_out,'(12x,''Buffer Slice Size: '',i3)') this%buffer_nslice
  endif

end subroutine DatasetBasePrint

! ************************************************************************** !

subroutine DatasetBaseStrip(this)
  ! 
  ! Strips allocated objects within base dataset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_base_type) :: this
  
  call DeallocateArray(this%iarray)
  call DeallocateArray(this%rarray)
  call DeallocateArray(this%ibuffer)
  call DeallocateArray(this%rbuffer)
  call DeallocateArray(this%dims)
  
  call TimeStorageDestroy(this%time_storage)
  
end subroutine DatasetBaseStrip

! ************************************************************************** !

subroutine DatasetBaseDestroy(dataset)
  ! 
  ! Destroys a dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11, 05/03/13
  ! 

  implicit none
  
  class(dataset_base_type), pointer :: dataset
  
  if (.not.associated(dataset)) return
  
  call DatasetBaseStrip(dataset)
  
  deallocate(dataset)
  nullify(dataset)
  
end subroutine DatasetBaseDestroy

end module Dataset_Base_class
