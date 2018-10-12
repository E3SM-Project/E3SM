module Dataset_Global_HDF5_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Dataset_Common_HDF5_class
  use DM_Kludge_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(dataset_common_hdf5_type) :: dataset_global_hdf5_type
    PetscInt :: local_size    ! local number of entries on this process
    PetscInt :: global_size   ! global number of entries
    type(dm_ptr_type), pointer :: dm_wrapper ! pointer 
  end type dataset_global_hdf5_type

  PetscInt, parameter :: default_max_buffer_size = 10
  
  public :: DatasetGlobalHDF5Create, &
            DatasetGlobalHDF5Init, &
            DatasetGlobalHDF5Cast, &
            DatasetGlobalHDF5Load, &
            DatasetGlobalHDF5Print, &
            DatasetGlobalHDF5Destroy
  
contains

! ************************************************************************** !

function DatasetGlobalHDF5Create()
  ! 
  ! Creates global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
  implicit none
  
  class(dataset_global_hdf5_type), pointer :: dataset

  class(dataset_global_hdf5_type), pointer :: DatasetGlobalHDF5Create
  
  allocate(dataset)
  call DatasetGlobalHDF5Init(dataset)

  DatasetGlobalHDF5Create => dataset
    
end function DatasetGlobalHDF5Create

! ************************************************************************** !

subroutine DatasetGlobalHDF5Init(this)
  ! 
  ! Initializes members of global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
  implicit none
  
  class(dataset_global_hdf5_type) :: this
  
  call DatasetCommonHDF5Init(this)
  this%local_size = 0
  this%global_size = 0
  nullify(this%dm_wrapper)
    
end subroutine DatasetGlobalHDF5Init

! ************************************************************************** !

function DatasetGlobalHDF5Cast(this)
  ! 
  ! Casts a dataset_base_type to database_global_type
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 

  use Dataset_Base_class
  
  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_global_hdf5_type), pointer :: DatasetGlobalHDF5Cast
  
  nullify(DatasetGlobalHDF5Cast)
  select type (this)
    class is (dataset_global_hdf5_type)
      DatasetGlobalHDF5Cast => this
  end select
    
end function DatasetGlobalHDF5Cast

! ************************************************************************** !

subroutine DatasetGlobalHDF5Load(this,option)
  ! 
  ! Load new data into dataset buffer
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 
  
#if defined(PETSC_HAVE_HDF5)    
  use hdf5, only : H5T_NATIVE_DOUBLE
#endif
  use Option_module
  use Time_Storage_module
  use Dataset_Base_class  

  implicit none
  
  class(dataset_global_hdf5_type) :: this
  type(option_type) :: option
  
  if (.not.associated(this%dm_wrapper)) then
    option%io_buffer = 'dm_wrapper not associated in Global Dataset: ' // &
      trim(this%name)
    call printErrMsg(option)
  endif
  
  if (DatasetCommonHDF5Load(this,option)) then
    if (.not.associated(this%rarray)) then
      if (this%local_size == 0) then
        option%io_buffer = 'Local size of Global Dataset has not been set.'
        call printErrMsg(option)
      endif
      allocate(this%rarray(this%local_size))
      this%rarray = 0.d0
    endif
    if (.not.associated(this%rbuffer)) then ! not initialized
      if (this%max_buffer_size < 0) then
        this%max_buffer_size = default_max_buffer_size
      endif
      this%buffer_nslice = min(this%max_buffer_size, &
                               this%time_storage%max_time_index)
      allocate(this%rbuffer(this%local_size*this%buffer_nslice))
      this%rbuffer = 0.d0
    endif
#if defined(PETSC_HAVE_HDF5)    
    call DatasetGlobalHDF5ReadData(this,option,H5T_NATIVE_DOUBLE)
#endif  
    ! no need to reorder since it is 1D in the h5 file.
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetGlobalHDF5Load

#if defined(PETSC_HAVE_HDF5)    

! ************************************************************************** !

subroutine DatasetGlobalHDF5ReadData(this,option,data_type)
  ! 
  ! Read an hdf5 array into a Petsc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use hdf5
  use Logging_module
  use Option_module
  use HDF5_Aux_module
  
  implicit none

! Default HDF5 Mechanism 
 
  class(dataset_global_hdf5_type) :: this
  type(option_type) :: option
  integer(HID_T) :: data_type 
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3), max_dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: ndims
  PetscMPIInt :: rank_mpi
  Vec :: natural_vec
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i, istart
  PetscInt :: buffer_size, buffer_rank2_size, file_rank1_size, file_rank2_size
  integer, allocatable :: integer_buffer_i4(:)
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err

  PetscViewer :: viewer
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)

  if (this%dm_wrapper%dm /= PETSC_NULL_DM) then
    call DMCreateGlobalVector(this%dm_wrapper%dm,global_vec, &
                              ierr);CHKERRQ(ierr)
    call DMDACreateNaturalVector(this%dm_wrapper%dm,natural_vec, &
                                 ierr);CHKERRQ(ierr)
  else
    option%io_buffer = 'ugdm not yet supported in DatasetGlobalHDF5ReadData()'
    call printErrMsg(option)
  endif
  
  call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
  call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)

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

  string = trim(this%hdf5_dataset_name) // '/Data'
  if (this%realization_dependent) then
    write(word,'(i9)') option%id
    string = trim(string) // trim(adjustl(word))
  endif
  option%io_buffer = 'Opening data set: ' // trim(string)
  call printMsg(option)  
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  
  call h5sget_simple_extent_ndims_f(file_space_id,ndims,hdf5_err)
  call h5sget_simple_extent_dims_f(file_space_id,dims,max_dims,hdf5_err)
  
  buffer_size = size(this%rbuffer)
  buffer_rank2_size = buffer_size / this%local_size
  file_rank1_size = int(dims(1))
  if (ndims > 1) then
    file_rank2_size = int(dims(2))
  else
    if (option%mycommsize > 1) then
      option%io_buffer = 'Dataset "' // trim(this%hdf5_dataset_name) // &
        '" in file "' // trim (this%filename) // &
        '" must be a 2D dataset (time,cell) if PFLOTRAN is run in parallel.'
      call printErrMsg(option)
    endif
    file_rank2_size = 1
  endif

  allocate(this%dims(1))
  this%dims = size(this%rarray)
  this%rank = 1

  if (mod(buffer_size,this%local_size) /= 0) then
    write(option%io_buffer, &
          '(a," buffer dimension (",i9,") is not a multiple of local domain",&
           &" dimension (",i9,").")') trim(this%hdf5_dataset_name), &
           size(this%rbuffer,1), this%local_size
    call printErrMsg(option)   
  endif

  if (mod(file_rank1_size,this%global_size) /= 0) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") is not a multiple of domain",&
           &" dimension (",i9,").")') trim(this%hdf5_dataset_name), &
           file_rank1_size, this%global_size
    call printErrMsg(option)   
  endif

  istart = 0
  call MPI_Exscan(this%local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr)
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif

  ! must initialize here to avoid error below when closing memory space
  memory_space_id = -1
  
  ! offset is zero-based
  offset = 0
  length = 0
  stride = 1
  
  offset(1) = istart ! istart is offset in the first dimension
  length(1) = file_rank1_size
  if (ndims > 1) then
    offset(2) = this%buffer_slice_offset
    length(1) = this%local_size
    length(2) = min(buffer_rank2_size,file_rank2_size-this%buffer_slice_offset)
  else
    offset(1) = this%buffer_slice_offset*this%local_size
    length(1) = min(buffer_size, &
                    file_rank1_size-this%buffer_slice_offset*this%local_size)
  endif
  call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                             length,hdf5_err,stride,stride) 

  dims = 0
  rank_mpi = 1
  if (ndims > 1) then
    dims(1) = min(buffer_size, &
                  this%local_size*(file_rank2_size-this%buffer_slice_offset))
  else
    dims(1) = min(buffer_size, &
                  file_rank1_size-this%buffer_slice_offset*this%local_size)
  endif
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! initialize to UNINITIALIZED_INTEGER to catch errors
  this%rbuffer = UNINITIALIZED_DOUBLE
  
  if (data_type == H5T_NATIVE_DOUBLE) then
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,this%rbuffer,dims, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
  else if (data_type == H5T_NATIVE_INTEGER) then
    allocate(integer_buffer_i4(size(this%rbuffer)))
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer_i4,dims, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    do i = 1, min(buffer_size, &
                  this%local_size*(file_rank2_size-this%buffer_slice_offset))
      this%rbuffer(i) = real(integer_buffer_i4(i))
    enddo
    deallocate(integer_buffer_i4)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  string = trim(this%hdf5_dataset_name) // '/Data'
  option%io_buffer = 'Closing data set: ' // trim(string)
  call printMsg(option)  
  call h5dclose_f(data_set_id,hdf5_err)
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err) 
  
  istart = 0
  do i = 1, min(buffer_rank2_size,file_rank2_size-this%buffer_slice_offset)
    call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
    vec_ptr(:) = this%rbuffer(istart+1:istart+this%local_size)
    call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

    !geh: for debugging purposes
    !write(string,*) i
    !string = trim(adjustl(this%hdf5_dataset_name)) // '_' // &
    !         trim(adjustl(string)) // '_natural.txt'
    !call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
    !call VecView(natural_vec,viewer,ierr)
    !call PetscViewerDestroy(viewer,ierr)

    call DMDANaturalToGlobalBegin(this%dm_wrapper%dm,natural_vec, &
                                  INSERT_VALUES,global_vec,ierr);CHKERRQ(ierr)
    call DMDANaturalToGlobalEnd(this%dm_wrapper%dm,natural_vec, &
                                INSERT_VALUES,global_vec,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
    this%rbuffer(istart+1:istart+this%local_size) = vec_ptr(:)
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)

    !geh: for debugging purposes
    !write(string,*) i
    !string = trim(adjustl(this%hdf5_dataset_name)) // '_' // &
    !         trim(adjustl(string)) // '_global.txt'
    !call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
    !call VecView(global_vec,viewer,ierr)
    !call PetscViewerDestroy(viewer,ierr)

    istart = istart + this%local_size
  enddo

  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
! End of Default HDF5 Mechanism

end subroutine DatasetGlobalHDF5ReadData
#endif

! ************************************************************************** !

subroutine DatasetGlobalHDF5Print(this,option)
  ! 
  ! Prints dataset info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/13
  ! 

  use Option_module

  implicit none
  
  class(dataset_global_hdf5_type), target :: this
  type(option_type) :: option
  
  class(dataset_common_hdf5_type), pointer :: dataset_hdf5

  dataset_hdf5 => this
  call DatasetCommonHDF5Print(this,option)

  ! no need to print local_size as it varies by process
  write(option%fid_out,'(10x,''Global Size: '',i2)') this%global_size
  
end subroutine DatasetGlobalHDF5Print

! ************************************************************************** !

subroutine DatasetGlobalHDF5Strip(this)
  ! 
  ! Strips allocated objects within Global dataset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_global_hdf5_type) :: this
  
  call DatasetCommonHDF5Strip(this)
  nullify(this%dm_wrapper) ! do not deallocate, as this is solely a pointer
  
end subroutine DatasetGlobalHDF5Strip

! ************************************************************************** !

subroutine DatasetGlobalHDF5Destroy(this)
  ! 
  ! Destroys a dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  implicit none
  
  class(dataset_global_hdf5_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetGlobalHDF5Strip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetGlobalHDF5Destroy

end module Dataset_Global_HDF5_class
