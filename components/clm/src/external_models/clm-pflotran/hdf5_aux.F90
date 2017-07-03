module HDF5_Aux_module

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  use Logging_module

  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  PetscInt, parameter, public :: HDF5_READ_BUFFER_SIZE = 1000000
!#define HDF5_BROADCAST

  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: io_rank_mpi
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE  
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  public :: HDF5ReadNDimRealArray, &
#ifdef SCORPIO
            HDF5ReadDatasetInteger2D, &
            HDF5ReadDatasetReal2D, &
            HDF5GroupExists, &
            HDF5DatasetExists, &
#else
            HDF5GroupExists, &
            HDF5DatasetExists, &
#endif
! SCORPIO
            HDF5MakeStringCompatible, &
            HDF5ReadDbase, &
            HDF5OpenFileReadOnly

contains

#endif
  

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine HDF5ReadNDimRealArray(option,file_id,dataset_name,ndims,dims, &
                                 real_array)
  ! 
  ! Read in an n-dimensional array from an hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/10
  ! 

  use hdf5
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  integer(HID_T) :: file_id
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: ndims_hdf5
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  integer(HSIZE_T) :: num_reals_in_dataset
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  call PetscLogEventBegin(logging%event_read_ndim_real_array_hdf5, &
                          ierr);CHKERRQ(ierr)
                          
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  ndims = ndims_hdf5
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  allocate(dims(ndims))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  dims = int(dims_h5)
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_dataset,hdf5_err)
  temp_int = dims(1)
  do i = 2, ndims
    temp_int = temp_int * dims(i)
  enddo

  rank_mpi = 1
  offset = 0
  length = num_reals_in_dataset
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err,length)

  allocate(real_array(num_reals_in_dataset))
  real_array = 0.d0
#ifdef HDF5_BROADCAST
  if (option%myrank == option%io_rank) then                           
#endif
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_array,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
  endif
  if (option%mycommsize > 1) then
    int_mpi = num_reals_in_dataset
    call MPI_Bcast(real_array,int_mpi,MPI_DOUBLE_PRECISION, &
                   option%io_rank,option%mycomm,ierr)
  endif
#endif

  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  
  deallocate(dims_h5)
  deallocate(max_dims_h5) 

  call PetscLogEventEnd(logging%event_read_ndim_real_array_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine HDF5ReadNDimRealArray

#if defined(SCORPIO)

! ************************************************************************** !

subroutine HDF5ReadDatasetInteger2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)
  ! 
  ! Author: Gautam Bisht
  ! Date: 05/13/2010
  ! 

  use hdf5
  use Option_module
  
  implicit none
  
#if defined(SCORPIO)
  include "scorpiof.h"  
#endif

  ! in
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  integer :: read_option
  type(option_type) :: option
  
  ! out
  PetscInt,pointer :: data(:,:)
  PetscInt :: data_dims(2)
  PetscInt :: dataset_dims(2)
  
  ! local
  PetscInt :: file_id
  PetscInt :: ndims
  PetscInt :: ii, remainder

  PetscErrorCode :: ierr
  
  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  call fscorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call fscorpio_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call fscorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call fscorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call fscorpio_read_dataset( data, SCORPIO_INTEGER, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)
  
  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call fscorpio_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetInteger2D
#endif
! SCORPIO

#if defined(SCORPIO)

! ************************************************************************** !

subroutine HDF5ReadDatasetReal2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)
  use hdf5
  use Option_module
  
  implicit none
  
#if defined(SCORPIO)
  include "scorpiof.h"  
#endif

  ! in
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  integer :: read_option
  type(option_type) :: option
  
  ! out
  PetscReal,pointer :: data(:,:)
  PetscInt :: data_dims(2)
  PetscInt :: dataset_dims(2)
  
  ! local
  integer :: file_id
  integer :: ndims
  PetscInt :: ii, remainder

  PetscErrorCode :: ierr
  
  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  call fscorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call fscorpio_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call fscorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call fscorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call fscorpio_read_dataset( data, SCORPIO_DOUBLE, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)

  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call fscorpio_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetReal2D
#endif

! ************************************************************************** !

function HDF5GroupExists(filename,group_name,option)
  ! 
  ! SCORPIO
  ! Returns true if a group exists
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/26/2012
  ! 

  use hdf5
  use Option_module
  
  implicit none

#if defined(SCORPIO)
  include "scorpiof.h"  
#endif
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  type(option_type) :: option

  integer(HID_T) :: file_id  
  integer(HID_T) :: grp_id  
  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0 
  PetscBool :: group_exists
  
  PetscBool :: HDF5GroupExists

#if defined(SCORPIO)

  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  group_name = trim(group_name) // CHAR(0)
  call fscorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)
  call fscorpio_group_exists(group_name, file_id, option%ioread_group_id, ierr)
  group_exists = (ierr == 1);

  if (group_exists) then
    HDF5GroupExists = PETSC_TRUE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" found in file.'
  else
    HDF5GroupExists = PETSC_FALSE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" not found in file.  Therefore, assuming a ' // &
      'cell-indexed dataset.'
  endif
  call printMsg(option)

  call fscorpio_close_file( file_id, option%ioread_group_id, ierr)  
#else
  ! open the file
  call h5open_f(hdf5_err)
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  option%io_buffer = 'Testing group: ' // trim(group_name)
  call printMsg(option)
  ! I turn off error messaging since if the group does not exist, an error
  ! will be printed, but the user does not need to see this.
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  group_exists = .not.(hdf5_err < 0)
  call h5eset_auto_f(ON,hdf5_err)  

  if (group_exists) then
    HDF5GroupExists = PETSC_TRUE
    call h5gclose_f(grp_id,hdf5_err)  
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" found in file.'
  else
    HDF5GroupExists = PETSC_FALSE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" not found in file.  Therefore, assuming a ' // &
      'cell-indexed dataset.'
  endif
  call printMsg(option)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)  
#endif 
!SCORPIO

end function HDF5GroupExists

! ************************************************************************** !

function HDF5DatasetExists(filename,group_name,dataset_name,option)
  !
  ! SCORPIO
  ! Returns true if a dataset exists
  !
  ! Author: Gautam Bisht
  ! Date: 04/30/2015
  !

  use hdf5
  use Option_module

  implicit none

#if defined(SCORPIO)
  include "scorpiof.h"
#endif

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  character(len=MAXWORDLENGTH) :: dataset_name
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: group_name_local
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: dataset_id
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscBool :: group_exists
  PetscBool :: dataset_exists

  PetscBool :: HDF5DatasetExists

#if defined(SCORPIO)

  option%io_buffer = 'Need to extend HDF5DatasetExists() for SCORPIO.'
  call printErrMsg(option)

#else

  if (len_trim(group_name) == 0) then
    group_name_local = "/" // CHAR(0)
  else
    group_name_local = trim(group_name) // "/" // CHAR(0)
  endif

  ! open the file
  call h5open_f(hdf5_err)
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! I turn off error messaging since if the group does not exist, an error
  ! will be printed, but the user does not need to see this.
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,group_name_local,grp_id,hdf5_err)
  group_exists = .not.(hdf5_err < 0)
  call h5eset_auto_f(ON,hdf5_err)

  if (.not.group_exists) then
    HDF5DatasetExists = PETSC_FALSE
  endif

  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  dataset_exists = .not.(hdf5_err < 0)

  if (.not.dataset_exists) then
    HDF5DatasetExists = PETSC_FALSE
  else
    HDF5DatasetExists = PETSC_TRUE
    call h5dclose_f(dataset_id,hdf5_err)
  endif

  if (group_exists) call h5gclose_f(grp_id,hdf5_err)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif
!SCORPIO

end function HDF5DatasetExists

! ************************************************************************** !

subroutine HDF5MakeStringCompatible(name)
  ! 
  ! Replaces '/' in string with '_' for hdf5 names
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  character(len=*) :: name
  
  PetscInt :: len, ichar
  
  len = len_trim(name)
  do ichar = 1, len
    if (name(ichar:ichar) == '/') then
      name(ichar:ichar) = '_'
    endif
  enddo
  
  name = trim(name)

end subroutine HDF5MakeStringCompatible

! ************************************************************************** !

subroutine HDF5ReadDbase(filename,option)
  ! 
  ! Read in an ASCII database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  use String_module
  use Input_Aux_module, only : dbase
  use h5lt
  
  implicit none
  
  character(len=*) :: filename
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH), allocatable :: wbuffer(:)
  character(len=MAXWORDLENGTH) :: wbuffer_word
  PetscReal, allocatable :: rbuffer(:)
  PetscInt, allocatable :: ibuffer(:)
  PetscInt :: dummy_int
  PetscInt :: value_index
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: object_name
  character(len=MAXWORDLENGTH) :: word
#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: num_objects
  integer(HID_T) :: i_object
  integer(HID_T) :: object_type
  integer(HID_T) :: prop_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: class_id
  integer(HID_T) :: datatype_id
  integer(HID_T) :: datatype_id2
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: num_values_in_dataset
  integer(SIZE_T) size_t_int
!  integer(HSIZE_T) :: dims(1)
  integer(SIZE_T) :: type_size
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: int_mpi
  PetscMPIInt :: hdf5_err
#endif
  PetscInt :: num_ints
  PetscInt :: num_reals
  PetscInt :: num_words

  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)
  call h5gn_members_f(file_id, '.',num_objects, hdf5_err)
  num_ints = 0
  num_reals = 0
  num_words = 0
  ! index is zero-based
  do i_object = 0, num_objects-1
    ! read in to string in case the name is too large.
    call h5gget_obj_info_idx_f(file_id,'.',i_object,string, &
                               object_type,hdf5_err)
    if (len_trim(string) > MAXWORDLENGTH) then
      option%io_buffer = 'HDF5 DBASE object names must be shorter than &
        &32 characters: ' // trim(string)
      call printErrMsg(option)
    endif
    object_name = trim(string)
    if (object_type == H5G_DATASET_F) then
      call h5dopen_f(file_id,object_name,dataset_id,hdf5_err)
      call h5dget_type_f(dataset_id, datatype_id, hdf5_err)
      call h5tget_class_f(datatype_id, class_id, hdf5_err)
      ! cannot use a select case statement since the H5T definitions are not
      ! guaranteed to be constant.  the preprocessor throws an error
      if (class_id == H5T_INTEGER_F) then
        num_ints = num_ints + 1
      else if (class_id == H5T_FLOAT_F) then
        num_reals = num_reals + 1
      else if (class_id == H5T_STRING_F) then
        num_words = num_words + 1
      else
        option%io_buffer = 'Unrecognized HDF5 datatype in Dbase: ' // &
          trim(object_name)
        call printErrMsg(option)
      endif
      call h5tclose_f(datatype_id, hdf5_err)
      call h5dclose_f(dataset_id,hdf5_err)
    endif
  enddo
  allocate(dbase)
  if (num_ints > 0) then
    allocate(dbase%icard(num_ints))
    dbase%icard = ''
    allocate(dbase%ivalue(num_ints))
    dbase%ivalue = UNINITIALIZED_INTEGER
  endif
  if (num_reals > 0) then
    allocate(dbase%rcard(num_reals))
    dbase%rcard = ''
    allocate(dbase%rvalue(num_reals))
    dbase%rvalue = UNINITIALIZED_DOUBLE
  endif
  if (num_words > 0) then
    allocate(dbase%ccard(num_words))
    dbase%ccard = ''
    allocate(dbase%cvalue(num_words))
    dbase%cvalue = '-999'
  endif
  value_index = 1
  if (option%id > 0) then
    value_index = option%id
  endif
  num_ints = 0
  num_reals = 0
  num_words = 0
  do i_object = 0, num_objects-1
    call h5gget_obj_info_idx_f(file_id,'.',i_object,object_name, &
                               object_type,hdf5_err)
    if (object_type == H5G_DATASET_F) then
! use once HDF5 lite is linked in PETSc      
!      call h5ltget_dataset_info_f(file_id,object_name,dims,dummy_int, &
!                                  type_size,hdf5_err)
!      allocate(buffer(dims(1)))
!      buffer = 0.d0
!      call h5ltread_dataset_double_f(file_id,object_name,buffer, &
!                                     dims,hdf5_err)
!      dbase%card(icount) = trim(object_name)
!      if (option%id > 0) then
!        if (option%id > dims(1)) then
!          write(word,*) dims(1)
!          option%io_buffer = 'DBASE dataset "' // trim(object_name) // &
!            '" is too small (' // trim(adjustl(word)) // &
!            ') for number of realizations.'
!          call printErrMsg(option)
!        endif
!        dbase%value(icount) = buffer(option%id)
!      else
!        dbase%value(icount) = buffer(1)
!      endif
!      deallocate(buffer)

      call h5dopen_f(file_id,object_name,dataset_id,hdf5_err)
      call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
      ! should be a rank=1 data space
      call h5sget_simple_extent_npoints_f(file_space_id, &
                                          num_values_in_dataset,hdf5_err)
      if (option%id > 0) then
        if (option%id > num_values_in_dataset) then
          write(word,*) num_values_in_dataset
          option%io_buffer = 'Data in DBASE_FILENAME "' // &
            trim(object_name) // &
            '" is too small (' // trim(adjustl(word)) // &
            ') for number of realizations.'
          call printErrMsg(option)
        endif
      endif
      rank_mpi = 1
      offset = 0
      length = num_values_in_dataset
      stride = 1
      call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
      call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err, &
                              length)

      call h5dget_type_f(dataset_id, datatype_id, hdf5_err)
      call h5tget_class_f(datatype_id, class_id, hdf5_err)
      call h5tclose_f(datatype_id, hdf5_err)
#ifdef HDF5_BROADCAST
      if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      if (class_id == H5T_INTEGER_F) then
        allocate(ibuffer(num_values_in_dataset))
        ibuffer = 0
        call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,ibuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
      else if (class_id == H5T_FLOAT_F) then
        allocate(rbuffer(num_values_in_dataset))
        rbuffer = 0.d0
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,rbuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
      else if (class_id == H5T_STRING_F) then
        call h5tcopy_f(H5T_NATIVE_CHARACTER,datatype_id2,hdf5_err)
        size_t_int = MAXWORDLENGTH
        call h5tset_size_f(datatype_id2,size_t_int,hdf5_err)
        allocate(wbuffer(num_values_in_dataset))
        wbuffer = ''
        call h5dread_f(dataset_id,datatype_id2,wbuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
        wbuffer_word = wbuffer(value_index)
        deallocate(wbuffer)
        call h5tclose_f(datatype_id2,hdf5_err)
      endif
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
      endif
      if (option%mycommsize > 1) then
        int_mpi = num_values_in_dataset
        if (class_id == H5T_INTEGER_F) then
          call MPI_Bcast(ibuffer,int_mpi,MPI_INTEGER, &
                         option%io_rank,option%mycomm,ierr)
        else if (class_id == H5T_FLOAT_F) then
          call MPI_Bcast(rbuffer,int_mpi,MPI_DOUBLE_PRECISION, &
                         option%io_rank,option%mycomm,ierr)
        else if (class_id == H5T_STRING_F) then
          int_mpi = MAXWORDLENGTH
          call MPI_Bcast(wbuffer_word,int_mpi,MPI_CHARACTER, &
                         option%io_rank,option%mycomm,ierr)
        endif
      endif
#endif
      call h5pclose_f(prop_id,hdf5_err)
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      call h5sclose_f(file_space_id,hdf5_err)
      call h5dclose_f(dataset_id,hdf5_err)
      call StringToUpper(object_name)
      ! these conditionals must come after the bcasts above!!!
      if (class_id == H5T_INTEGER_F) then
        num_ints = num_ints + 1
        dbase%icard(num_ints) = trim(object_name)
        dbase%ivalue(num_ints) = ibuffer(value_index)
        deallocate(ibuffer)
      else if (class_id == H5T_FLOAT_F) then
        num_reals = num_reals + 1
        dbase%rcard(num_reals) = trim(object_name)
        dbase%rvalue(num_reals) = rbuffer(value_index)
        deallocate(rbuffer)
      else if (class_id == H5T_STRING_F) then
        num_words = num_words + 1
        dbase%ccard(num_words) = trim(object_name)
        dbase%cvalue(num_words) = wbuffer_word
      endif
    endif
  enddo
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
      
end subroutine HDF5ReadDbase

! ************************************************************************** !

subroutine HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  ! 
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  ! 
  use hdf5
  use Option_module
  
  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  type(option_type) :: option
  
  PetscMPIInt :: hdf5_err

  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  if (hdf5_err /= 0) then
    option%io_buffer = 'HDF5 file "' // trim(filename) // '" not found.'
    call printErrMsg(option)
  endif
  
end subroutine HDF5OpenFileReadOnly

#endif
! defined(PETSC_HAVE_HDF5)

end module HDF5_Aux_module
