module Dataset_Map_HDF5_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Dataset_Common_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(dataset_common_hdf5_type) :: dataset_map_hdf5_type
    character(len=MAXSTRINGLENGTH) :: h5_dataset_map_name
    character(len=MAXSTRINGLENGTH) :: map_filename
    PetscInt, pointer :: mapping(:,:)
    PetscInt :: map_dims_global(2)
    PetscInt :: map_dims_local(2)
    PetscInt, pointer :: datatocell_ids(:)
    PetscInt, pointer :: cell_ids_local(:)
    PetscBool :: first_time
  end type dataset_map_hdf5_type
  
  PetscInt, parameter :: MAX_NSLICE = 100
  
  public :: DatasetMapHDF5Create, &
            DatasetMapHDF5Init, &
            DatasetMapHDF5Cast, &
            DatasetMapHDF5Read, &
            DatasetMapHDF5Load, &
            DatasetMapHDF5Print, &
            DatasetMapHDF5Destroy
  
contains

! ************************************************************************** !

function DatasetMapHDF5Create()
  ! 
  ! Creates global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 
  
  implicit none
  
  class(dataset_map_hdf5_type), pointer :: dataset

  class(dataset_map_hdf5_type), pointer :: DatasetMapHDF5Create
  
  allocate(dataset)
  call DatasetMapHDF5Init(dataset)

  DatasetMapHDF5Create => dataset
    
end function DatasetMapHDF5Create

! ************************************************************************** !

subroutine DatasetMapHDF5Init(this)
  ! 
  ! Initializes members of global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 
  
  implicit none
  
  class(dataset_map_hdf5_type) :: this
  
  call DatasetCommonHDF5Init(this)
  this%h5_dataset_map_name = ''
  this%map_filename = ''
  nullify(this%mapping)
  this%map_dims_global = 0
  this%map_dims_local = 0
  nullify(this%datatocell_ids)
  nullify(this%cell_ids_local)
  this%first_time = PETSC_TRUE
    
end subroutine DatasetMapHDF5Init

! ************************************************************************** !

function DatasetMapHDF5Cast(this)
  ! 
  ! Casts a dataset_base_type to dataset_map_hdf5_type
  ! 
  ! Date: 05/03/13
  ! 

  use Dataset_Base_class
  
  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_map_hdf5_type), pointer :: DatasetMapHDF5Cast
  
  nullify(DatasetMapHDF5Cast)
  select type (this)
    class is (dataset_map_hdf5_type)
      DatasetMapHDF5Cast => this
  end select
    
end function DatasetMapHDF5Cast

! ************************************************************************** !

subroutine DatasetMapHDF5Read(this,input,option)
  ! 
  ! Reads in contents of a dataset card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11, 06/04/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(dataset_map_hdf5_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DATASET')
    call StringToUpper(keyword)   
      
    call DatasetCommonHDF5ReadSelectCase(this,input,keyword,found,option)
    if (.not.found) then
      select case(trim(keyword))
        case('MAP_HDF5_DATASET_NAME') 
          call InputReadWord(input,option,this%h5_dataset_map_name,PETSC_TRUE)
          call InputErrorMsg(input,option,'map dataset name','DATASET')
        case('MAP_FILENAME_NAME') 
          call InputReadWord(input,option,this%map_filename,PETSC_TRUE)
          call InputErrorMsg(input,option,'map filename','DATASET')
        case default
          call InputKeywordUnrecognized(keyword,'dataset',option)
      end select
    endif
  
  enddo
  
  if (len_trim(this%hdf5_dataset_name) < 1) then
    this%hdf5_dataset_name = this%name
  endif

  if (len_trim(this%map_filename) < 1) then
    this%map_filename = this%filename
  endif
  
end subroutine DatasetMapHDF5Read

! ************************************************************************** !

subroutine DatasetMapHDF5Load(this,option)
  ! 
  ! Load new data into dataset buffer
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 
  
  use Option_module
  use Time_Storage_module
  use Dataset_Base_class

  implicit none
  
  class(dataset_map_hdf5_type) :: this
  type(option_type) :: option
  
  if (DatasetCommonHDF5Load(this,option)) then
#if defined(PETSC_HAVE_HDF5)    
    if (.not.associated(this%mapping)) then
      call DatasetMapHDF5ReadMap(this,option)
    endif
    call DatasetMapHDF5ReadData(this,option)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetMapHDF5Load

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine DatasetMapHDF5ReadData(this,option)
  ! 
  ! Read an hdf5 data into a array
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/11, 05/29/13
  ! 

  use hdf5
  use Option_module
  use Units_module
  use Logging_module
  use HDF5_Aux_module
  
  implicit none
  
  class(dataset_map_hdf5_type) :: this
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attribute_id
  integer(HID_T) :: atype_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(3)
  integer(HSIZE_T) :: offset(4), length(4), stride(4)
  integer(HSIZE_T) :: num_data_values
  integer :: ndims_h5
  PetscInt :: i, temp_array(4)
  PetscInt :: time_dim, num_times
  PetscInt :: num_dims_in_h5_file, num_times_in_h5_file
  PetscMPIInt :: array_rank_mpi
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: dataset_name, word

  call PetscLogEventBegin(logging%event_dataset_map_hdf5_read, &
                          ierr);CHKERRQ(ierr)

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(this%filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(this%filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! the dataset is actually stored in a group.  the group contains
  ! a "data" dataset and optionally a "time" dataset.
  option%io_buffer = 'Opening group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gopen_f(file_id,this%hdf5_dataset_name,grp_id,hdf5_err)

  time_dim = -1
  num_times = 1
  if (associated(this%time_storage)) then
    num_times = this%time_storage%max_time_index
    time_dim = 2
  endif
  
  ! open the "data" dataset
  dataset_name = 'Data'
  if (this%realization_dependent) then
    write(word,'(i9)') option%id
    dataset_name = trim(dataset_name) // trim(adjustl(word))
  endif
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! get dataset dimensions
  if (.not.associated(this%dims)) then
    call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)
    allocate(dims_h5(ndims_h5))
    allocate(max_dims_h5(ndims_h5))
    call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
    if (associated(this%time_storage)) then
      num_times_in_h5_file = int(dims_h5(1))
    else
      num_times_in_h5_file = 0
    endif
    if ((ndims_h5 == 2 .and. .not.associated(this%time_storage)) .or. &
        (ndims_h5 == 1 .and. associated(this%time_storage))) then
      option%io_buffer = 'Inconsistent dimensions in dataset file.'
      call printErrMsg(option)
    endif
    ! rank in space will always be 1
    this%rank = 1
    allocate(this%dims(this%rank))
    ! have to invert dimensions, but do not count the time dimension
    this%dims(1) = int(dims_h5(ndims_h5))
    deallocate(dims_h5)
    deallocate(max_dims_h5) 

    ! check to ensure that dataset is properly sized
    call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
    if (num_data_values/num_times /= this%dims(1)) then
      option%io_buffer = &
        'Number of values in dataset does not match dimensions.'
      call printErrMsg(option)
    endif
    if (associated(this%time_storage) .and. &
        num_times_in_h5_file /= num_times) then
      option%io_buffer = &
        'Number of times does not match last dimension of data array.'
      call printErrMsg(option)
    endif
    if (.not.associated(this%rarray)) then
      allocate(this%rarray(this%dims(1)))
    endif
    this%rarray = 0.d0
    if (associated(this%time_storage) .and. .not.associated(this%rbuffer)) then
      ! buffered array
      allocate(this%rbuffer(size(this%rarray)*MAX_NSLICE))
      this%rbuffer = 0.d0
    endif
  endif

  call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)

  if (associated(this%time_storage)) then
    num_dims_in_h5_file = this%rank + 1
  else
    num_dims_in_h5_file = this%rank
  endif
  
  array_rank_mpi = 1
  length = 1
  ! length (or size) must be adjusted according to the size of the 
  ! remaining data in the file
  this%buffer_nslice = min(MAX_NSLICE,(num_times-this%buffer_slice_offset))
  if (time_dim > 0) then
    length(1) = size(this%rarray) * this%buffer_nslice
  else
    length(1) = size(this%rarray)
  endif
  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    

  length = 1
  stride = 1
  offset = 0
  length(1) = this%dims(1)
  ! cannot read beyond end of the buffer
  if (time_dim > 0) then
    length(time_dim) = this%buffer_nslice
    offset(time_dim) = this%buffer_slice_offset
  endif
  
  !geh: for some reason, we have to invert here.  Perhaps because the
  !     dataset was generated in C???
  temp_array(1:num_dims_in_h5_file) = int(length(1:num_dims_in_h5_file))
  do i = 1, num_dims_in_h5_file
    length(i) = temp_array(num_dims_in_h5_file-i+1)
  enddo
  temp_array(1:num_dims_in_h5_file) = int(offset(1:num_dims_in_h5_file))
  do i = 1, num_dims_in_h5_file
    offset(i) = temp_array(num_dims_in_h5_file-i+1)
  enddo
  ! stride is fine
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                             hdf5_err,stride,stride) 
  if (associated(this%rbuffer)) then
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%rbuffer,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
  else
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%rarray,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
!    this%rmax = maxval(this%rarray)
!    this%rmin = minval(this%rarray)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  

  call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)

  option%io_buffer = 'Closing group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  
  call PetscLogEventEnd(logging%event_dataset_map_hdf5_read, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine DatasetMapHDF5ReadData

! ************************************************************************** !

subroutine DatasetMapHDF5ReadMap(this,option)
  ! 
  ! Read an hdf5 array
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/11, 05/29/13
  ! 

  use hdf5
  use Option_module
  use Units_module
  use Logging_module
  use HDF5_Aux_module
  
  implicit none
  
  class(dataset_map_hdf5_type) :: this
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(2), length(2)
  integer :: ndims_hdf5
  PetscInt :: i
  PetscMPIInt :: array_rank_mpi
  PetscMPIInt :: hdf5_err
  PetscInt :: nids_local, remainder, istart, iend
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, allocatable :: tempint_array(:,:)

  call PetscLogEventBegin(logging%event_dataset_map_hdf5_read, &
                          ierr);CHKERRQ(ierr)

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(this%map_filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(this%map_filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! the dataset is actually stored in a group.  the group contains
  ! a "data" dataset and optionally a "time" dataset.
  option%io_buffer = 'Opening group: ' // trim(this%h5_dataset_map_name)
  call printMsg(option)  
  call h5gopen_f(file_id,this%h5_dataset_map_name,grp_id,hdf5_err)
  
  ! Open the "data" dataset
  dataset_name = 'Data'
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  if (ndims_hdf5 /= 2) then
    option%io_buffer='Dimension of '// trim(this%h5_dataset_map_name) // &
      '/Data dataset in ' // trim(this%filename) // ' is not equal to 2.'
    call printErrMsg(option)
  endif

  ! Get dimensions of dataset
  allocate(dims_h5(ndims_hdf5))
  allocate(max_dims_h5(ndims_hdf5))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  
  nids_local=int(dims_h5(2)/option%mycommsize)
  remainder =int(dims_h5(2))-nids_local*option%mycommsize
  if (option%myrank<remainder) nids_local=nids_local+1

  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(nids_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(nids_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
  ! Determine the length and offset of data to be read by each processor
  length(1) = dims_h5(1)
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart

  ! TOdo: Only read part of data
!  option%io_buffer='Gautam: Modify code for reading HDF5 to generate mapping '//&
!   'of dataset'
!  call printMsg(option)
!  nids_local=dims_h5(2)
!  length(:) = dims_h5(:)
!  offset(:) = 0
  
  ! Save dimension size
  this%map_dims_global(:) = int(dims_h5(:))
  this%map_dims_local(:) = int(length(:))
  
  ! Create data space for dataset
  array_rank_mpi=2
  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err)

  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif

  ! Select hyperslab
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset,length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(this%mapping(length(1), length(2)))
  allocate(tempint_array(length(1), length(2)))

  ! Read the dataset collectively
  call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, tempint_array, &
                 dims_h5, hdf5_err, memory_space_id, file_space_id,prop_id)
  this%mapping = tempint_array

  call h5pclose_f(prop_id,hdf5_err)

  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  
  
  option%io_buffer = 'Closing group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)  
  
  call PetscLogEventEnd(logging%event_dataset_map_hdf5_read, &
                        ierr);CHKERRQ(ierr)
  
end subroutine DatasetMapHDF5ReadMap
#endif

! ************************************************************************** !

subroutine DatasetMapHDF5Print(this,option)
  ! 
  ! Prints dataset info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/13
  ! 

  use Option_module

  implicit none
  
  class(dataset_map_hdf5_type), target :: this
  type(option_type) :: option
  
  class(dataset_common_hdf5_type), pointer :: dataset_hdf5

  dataset_hdf5 => this
  call DatasetCommonHDF5Print(this,option)

  if (len_trim(this%h5_dataset_map_name) > 0) then
    write(option%fid_out,'(10x,''HDF5 Dataset Map Name: '',a)') &
      trim(this%h5_dataset_map_name)
  endif
  if (len_trim(this%map_filename) > 0) then
    write(option%fid_out,'(10x,''Map Filename: '',a)') &
      trim(this%map_filename)
  endif
  write(option%fid_out,'(10x,''Global Dimensions: '',2i8)') &
    this%map_dims_global(:)
  
end subroutine DatasetMapHDF5Print

! ************************************************************************** !

subroutine DatasetMapHDF5Strip(this)
  ! 
  ! Strips allocated objects within Map dataset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_map_hdf5_type) :: this
  
  call DatasetCommonHDF5Strip(this)
  
  call DeallocateArray(this%mapping)
  call DeallocateArray(this%datatocell_ids)
  call DeallocateArray(this%cell_ids_local)
  
end subroutine DatasetMapHDF5Strip

! ************************************************************************** !

subroutine DatasetMapHDF5Destroy(this)
  ! 
  ! Destroys a dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 

  implicit none
  
  class(dataset_map_hdf5_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetMapHDF5Strip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetMapHDF5Destroy

end module Dataset_Map_HDF5_class
