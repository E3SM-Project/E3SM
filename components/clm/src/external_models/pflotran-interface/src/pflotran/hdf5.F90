module HDF5_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module

  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscErrorCode :: ierr

  PetscBool, public :: trick_hdf5 = PETSC_FALSE

#if defined(SCORPIO)
  include "scorpiof.h"  
#endif

#if defined(PETSC_HAVE_HDF5)
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: io_rank
      
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#endif
      
#if defined(PETSC_HAVE_HDF5)
  public :: HDF5WriteStructDataSetFromVec, &
            HDF5WriteDataSetFromVec, &
            HDF5ReadDataSetInVec, &
            HDF5WriteStructuredDataSet
#endif            

  public :: HDF5ReadRegionFromFile, &
            HDF5ReadCellIndexedIntegerArray, &
            HDF5ReadCellIndexedRealArray, &
            HDF5QueryRegionDefinition, &
            HDF5ReadRegionDefinedByVertex

contains

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine HDF5ReadIntegerArraySplit(option,file_id,dataset_name,local_size, &
                                     integer_array)
  ! 
  ! Read in local integer values from hdf5 global file
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/21/07
  ! 

  use hdf5
  
  use Grid_module
  use Option_module
  use HDF5_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

#if defined(SCORPIO)
  type(option_type) :: option
  integer(HID_T) :: file_id
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: local_size
  PetscInt, pointer :: integer_array(:)
  
  option%io_buffer = 'HDF5ReadIntegerArraySplit() not yet setup for SCORPIO'
  call printErrMsg(option)

#else
! SCORPIO is not defined

  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscInt :: istart, iend, remainder
  PetscInt :: temp_int
  integer(HSIZE_T) :: num_integers_in_file
  integer(SIZE_T) :: string_size
  
  integer, allocatable :: integer_buffer_i4(:)
  
  PetscInt :: read_block_size

! Default & Glenn's HDF5 Broadcast Mechanism (uses HDF5 Independent I/O mode)

  call PetscLogEventBegin(logging%event_read_int_array_hdf5, &
                          ierr);CHKERRQ(ierr)

  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  if (hdf5_err /= 0) then
    string_size = MAXSTRINGLENGTH
    call h5fget_name_f(file_id,string,string_size,hdf5_err)
    option%io_buffer = 'HDF5 dataset "' // trim(dataset_name) // '" not found &
      &in file "' // trim(string) // '".'
    call printErrMsg(option)
  endif
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
  ! divide across all processes
  temp_int = int(num_integers_in_file)
  read_block_size = temp_int / option%mycommsize
  remainder = temp_int - read_block_size*option%mycommsize
  if (option%myrank < temp_int - read_block_size*option%mycommsize) &
    read_block_size = read_block_size + 1
  if (local_size > 0 .and. local_size /= read_block_size) then
    write(string,*) local_size, read_block_size
    option%io_buffer = 'Array mismatch in HDF5ReadIntegerArraySplit(): ' // &
      trim(adjustl(string))
    call printErrMsgByRank(option)
  endif
  istart = 0
  iend   = 0
  call MPI_Exscan(read_block_size, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                  MPI_SUM, option%mycomm, ierr)
  
  rank_mpi = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  if (read_block_size > 0) then
    dims = 0
    memory_space_id = -1
    dims(1) = read_block_size
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
    ! offset is zero-based
    offset(1) = istart
    length(1) = dims(1)
    allocate(integer_buffer_i4(read_block_size))
    integer_buffer_i4 = UNINITIALIZED_INTEGER
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)   
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    allocate(integer_array(read_block_size))
    integer_array = integer_buffer_i4
    deallocate(integer_buffer_i4)
    if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_int_array_hdf5,ierr);CHKERRQ(ierr)

#endif
! SCORPIO

end subroutine HDF5ReadIntegerArraySplit

! ************************************************************************** !

subroutine HDF5WriteStructuredDataSet(name,array,file_id,data_type,option, &
                                      nx_global,ny_global,nz_global, &
                                      nx_local,ny_local,nz_local, &
                                      istart_local,jstart_local,kstart_local)
  ! 
  ! Writes data from an array into HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use hdf5
  use Option_module
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  character(len=*) :: name
  PetscReal :: array(:)
  type(option_type) :: option
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: istart_local, jstart_local, kstart_local

#if defined(SCORPIO_WRITE)    
  integer :: file_id
  integer :: data_type
  integer :: file_space_id
  integer :: memory_space_id
  integer :: data_set_id
  integer :: prop_id
  integer :: dims(3),mem_dims(3)
  integer :: start(3), length(3), stride(3)
  integer :: nlmax ! Sarat added nlmax 
  integer :: rank_mpi,file_space_rank_mpi
  integer :: i, j, k, count, id
  integer :: ny_local_X_nz_local
  integer :: num_to_write_mpi
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3),mem_dims(3)
  PetscMPIInt :: rank_mpi,file_space_rank_mpi  
  PetscInt :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: ny_local_X_nz_local
  PetscMPIInt :: num_to_write_mpi
#endif
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscMPIInt :: hdf5_flag

  integer, pointer :: int_array_i4(:)
  PetscReal, pointer :: double_array(:)
  character(len=128) :: fscorpio_string

  name = trim(name) // CHAR(0)

#if defined(SCORPIO_WRITE)

  !geh: kludge to get scorpio to write name properly without garbage appended.
! fscorpio_string(1:len_trim(name)) = name

  fscorpio_string = trim(name)

  fscorpio_string = trim(fscorpio_string) // CHAR(0)

!  write(option%io_buffer,'(" Writing dataset block: ", a, " type - ", i, ".")') trim(name), data_type
!  call printMsg(option)
!  write(option%io_buffer,'(" HDF_NATIVE_INTEGER is ", i, " and H5T_NATIVE_DOUBLE is ", i, " and H5T_NATIVE_INTEGER is ", i, ".")') HDF_NATIVE_INTEGER, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER
!  call printMsg(option)

  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5, &
                          ierr);CHKERRQ(ierr)
  
  ny_local_X_nz_local = ny_local*nz_local
  num_to_write_mpi = nx_local*ny_local_X_nz_local
  
  ! file space which is a 3D block
  rank_mpi = 3

  dims(1) = nx_global
  dims(2) = ny_global
  dims(3) = nz_global

  !write(option%io_buffer,'(" Writing dataset block with dimensions: ", i,i,i, ".")') dims(1), dims(2), dims(3)
  !call printMsg(option)

  ! create the hyperslab
  start(1) = istart_local
  start(2) = jstart_local
  start(3) = kstart_local
  length(1) =  nx_local
  length(2) =  ny_local
  length(3) =  nz_local

  ! Sarat added this to eliminate grid%nlmax dependency in the non-INVERT case
  ! below
  nlmax = nx_local * ny_local * nz_local

  if (num_to_write_mpi == 0) length(1) = 1
  stride = 1

    ! write the data
  if (num_to_write_mpi > 0) then
    if (data_type == H5T_NATIVE_DOUBLE) then
      allocate(double_array(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            double_array(id) = array(count)
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
       !write(option%io_buffer, &
       !   '(a," Writing double dataset1: dimensions: ",i9,i9,i9, " Data type and ndims: ",i9, i9)') & 
       !trim(name), dims(1), dims(2), dims(3), SCORPIO_DOUBLE, rank_mpi
       !call printMsg(option)   

       call fscorpio_write_dataset_block(double_array, SCORPIO_DOUBLE, rank_mpi, &
              dims, length, start, file_id, trim(fscorpio_string), &
              option%iowrite_group_id, ierr)
      !call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      !hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call DeallocateArray(double_array)
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(int_array_i4(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array_i4(id) = int(array(count))
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
       !write(option%io_buffer, &
       !   '(a," Writing integer dataset1: dimensions: ",i9,i9,i9, " Data type and ndims: ",i9, i9)') & 
       !trim(name), dims(1), dims(2), dims(3), SCORPIO_INTEGER, rank_mpi
       !call printMsg(option)   
      call fscorpio_write_dataset_block(int_array_i4, SCORPIO_INTEGER, rank_mpi, &
              dims, length, start, file_id, trim(fscorpio_string), &
              option%iowrite_group_id, ierr)
      !!call h5dwrite_f(data_set_id,data_type,int_array_i4,dims, &
      !                hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call DeallocateArray(int_array_i4)
    endif
  endif

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5, &
                        ierr);CHKERRQ(ierr)

#else
! SCORPIO_WRITE is not defined  

! Default HDF5 Write

  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5, &
                          ierr);CHKERRQ(ierr)
  
  ny_local_X_nz_local = ny_local*nz_local
  num_to_write_mpi = nx_local*ny_local_X_nz_local
  
  ! memory space which is a 1D vector  
  rank_mpi = 1
  dims = 0
  dims(1) = num_to_write_mpi
  if (num_to_write_mpi == 0) dims(1) = 1
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank_mpi = 3
  ! have to trick hdf5 for now with inverted ordering
  dims(3) = nx_global
  dims(2) = ny_global
  dims(1) = nz_global
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,name,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,name,data_type,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  ! create the hyperslab
  start(3) = istart_local
  start(2) = jstart_local
  start(1) = kstart_local
  ! recall these are inverted
  length(3) =  nx_local
  length(2) =  ny_local
  length(1) =  nz_local
  if (num_to_write_mpi == 0) length(1) = 1
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err) 
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif
  if (num_to_write_mpi > 0) then
    if (data_type == HDF_NATIVE_INTEGER) then
      allocate(int_array_i4(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array_i4(id) = int(array(count))
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call h5dwrite_f(data_set_id,data_type,int_array_i4,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      ! cannot use DeallocateArray since int_array_i4 may not be PetscInt
      deallocate(int_array_i4)
      nullify(int_array_i4)
    else
      allocate(double_array(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            double_array(id) = array(count)
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call DeallocateArray(double_array)
    endif
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
#endif
! SCORPIO_WRITE vs previous

end subroutine HDF5WriteStructuredDataSet

! ************************************************************************** !

subroutine HDF5ReadIndices(grid,option,file_id,dataset_name,dataset_size, &
                           indices)
  ! 
  ! End of Default HDF5 Write
  ! GEH - Structured Grid Dependence - End
  ! Reads cell indices from an hdf5 dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 

  use hdf5
  
  use Option_module
  use Grid_module
  use Utility_module, only : DeallocateArray
  
  implicit none

#if defined(SCORPIO)

  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
    
#if defined(SCORPIO)    
  integer :: file_id
  integer :: file_space_id
  integer :: memory_space_id
  integer :: data_set_id
  integer :: prop_id
  integer :: dims(3)
  integer :: offset(3), length(3), stride(3), globaldims(3)
  integer :: num_data_in_file
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3), globaldims(3)
  integer(HSIZE_T) :: num_data_in_file
#endif
 
  PetscMPIInt :: rank_mpi
  ! seeting to MPIInt to ensure i4
  integer, allocatable :: indices_i4(:)
  
  PetscInt :: istart, iend

  call PetscLogEventBegin(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)
                        
  istart = 0  ! this will be zero-based
  iend = 0
  
  ! first determine upper and lower bound on PETSc global array
  call MPI_Scan(grid%nlmax,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                option%mycomm,ierr)
  istart = iend - grid%nlmax
  
  ! should be a rank=1 data space
  call fscorpio_get_dataset_size( num_data_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  globaldims(1) = num_data_in_file 
!???? get size of dataset  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)
  !if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
  !  write(option%io_buffer, &
  !        '(a," data space dimension (",i9,") does not match the dimensions",&
  !         &" of the domain (",i9,").")') trim(dataset_name), &
  !         num_data_in_file,dataset_size
  !  call printErrMsg(option)   
  !else
  !  dataset_size = num_data_in_file
  !endif  
  
  dataset_size = int(num_data_in_file)

  if (istart < num_data_in_file) then
  
    allocate(indices_i4(-1:iend-istart))
    allocate(indices(-1:iend-istart))
    indices_i4(-1) = istart
    indices_i4(0) = iend
  
    rank_mpi = 1
    offset = 0
    length = 0
    stride = 1
  
    dims = 0
    dims(1) = iend-istart
    memory_space_id = -1

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
                               
    call fscorpio_read_dataset(indices_i4(1:iend-istart), SCORPIO_INTEGER, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
    !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,indices_i4(1:iend-istart), &
                   !dims,hdf5_err,memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    indices(-1:iend-istart) = indices_i4(-1:iend-istart)                
    deallocate(indices_i4)
  endif
  
  call PetscLogEventEnd(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)

#else
! SCORPIO is not defined

! Default HDF5 Mechanism 

  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: ndims_h5
  PetscMPIInt :: rank_mpi
  ! seeting to MPIInt to ensure i4
  integer, allocatable :: indices_i4(:)
  integer(HSIZE_T) :: num_data_in_file
  integer(SIZE_T) :: string_size
  
  PetscInt :: istart, iend

  call PetscLogEventBegin(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)
                        
  istart = 0  ! this will be zero-based
  iend = 0
  
  ! first determine upper and lower bound on PETSc global array
  call MPI_Scan(grid%nlmax,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                option%mycomm,ierr)
  istart = iend - grid%nlmax
  
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  if (hdf5_err /= 0) then
    string_size = MAXSTRINGLENGTH
    call h5fget_name_f(file_id,string,string_size,hdf5_err)
    option%io_buffer = 'HDF5 dataset "' // trim(dataset_name) // '" not found &
      &in file "' // trim(string) // '".'
    call printErrMsg(option)
  endif
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)
  if (ndims_h5 /= 1) then
    write(option%io_buffer, &
          '(a," data space dimension (",i2,"D) must be 1D.")') &
          trim(dataset_name), ndims_h5
    call printErrMsg(option)
  endif
  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)
  if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_data_in_file,dataset_size
    call printErrMsg(option)   
  else
    dataset_size = int(num_data_in_file)
  endif  
  
  if (istart < num_data_in_file) then
  
    allocate(indices_i4(-1:iend-istart))
    allocate(indices(-1:iend-istart))
    indices_i4(-1) = istart
    indices_i4(0) = iend
  
    rank_mpi = 1
    offset = 0
    length = 0
    stride = 1
  
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
    dims = 0
    dims(1) = iend-istart
    memory_space_id = -1
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
                               
    call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,indices_i4(1:iend-istart), &
                   dims,hdf5_err,memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    indices(-1:iend-istart) = indices_i4(-1:iend-istart)                
    deallocate(indices_i4)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_indices_hdf5,ierr);CHKERRQ(ierr)
! End of Default HDF5 Mechanism  
  
#endif
! SCORPIO
  
end subroutine HDF5ReadIndices

! ************************************************************************** !

subroutine HDF5ReadArray(discretization,grid,option,file_id,dataset_name, &
                         dataset_size, &
                         indices,global_vec,data_type)
  ! 
  ! Read an hdf5 array into a Petsc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  
  use Option_module
  use Grid_module
  use Discretization_module
  use Utility_module, only : DeallocateArray
  
  implicit none

#if defined(SCORPIO)
 
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  Vec :: global_vec
  integer :: data_type 
  
  integer :: file_space_id
  integer :: memory_space_id
  integer :: data_set_id
  integer :: prop_id
  integer :: dims(3), globaldims(3)
  integer :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  integer :: num_data_in_file
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt, allocatable :: indices0(:)
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
                          
  istart = 0
  iend = 0
  
  ! should be a rank=1 data space

  call fscorpio_get_dataset_size( num_data_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  globaldims(1) = num_data_in_file
!???? get size   call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)

  !if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
  !  write(option%io_buffer, &
  !        '(a," data space dimension (",i9,") does not match the dimensions",&
  !         &" of the domain (",i9,").")') trim(dataset_name), &
  !         num_data_in_file,grid%nmax
  !  call printErrMsg(option)   
  !endif

  rank_mpi = 1
  offset = 0
  length = 0
  stride = 1
  
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)
  call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)

  ! must initialize here to avoid error below when closing memory space
  memory_space_id = -1

  if (associated(indices)) then

    istart = indices(-1)
    iend = indices(0)

    dims = 0
    dims(1) = iend-istart

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    allocate(real_buffer(iend-istart))
    if (data_type == H5T_NATIVE_DOUBLE) then
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    
      call fscorpio_read_dataset(real_buffer, SCORPIO_DOUBLE, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
      !call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     !hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(integer_buffer_i4(iend-istart))
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)

      call fscorpio_read_dataset(integer_buffer_i4, SCORPIO_INTEGER, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
      !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                     !hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      do i=1,iend-istart
        real_buffer(i) = real(integer_buffer_i4(i))
      enddo
      deallocate(integer_buffer_i4)
    endif
    ! must convert indices to zero based for VecSetValues
    allocate(indices0(iend-istart))
    indices0 = indices(1:iend-istart)-1
    call VecSetValues(natural_vec,iend-istart,indices0, &
                      real_buffer,INSERT_VALUES,ierr);CHKERRQ(ierr)
    deallocate(indices0)
    deallocate(real_buffer)

  endif

  call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
  call DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec, &
                                     ONEDOF)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)

#else
! SCORPIO is not defined

! Default HDF5 Mechanism 
 
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  Vec :: global_vec
  integer(HID_T) :: data_type 
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: num_data_in_file
  integer(SIZE_T) :: string_size
  integer :: ndims_h5
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt, allocatable :: indices0(:)
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
                          
  istart = 0
  iend = 0
  
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  if (hdf5_err /= 0) then
    string_size = MAXSTRINGLENGTH
    call h5fget_name_f(file_id,string,string_size,hdf5_err)
    option%io_buffer = 'HDF5 dataset "' // trim(dataset_name) // '" not found &
      &in file "' // trim(string) // '".'
    call printErrMsg(option)
  endif
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)
  if (ndims_h5 /= 1) then
    write(option%io_buffer, &
          '(a," data space dimension (",i2,"D) must be 1D.")') &
          trim(dataset_name), ndims_h5
    call printErrMsg(option)
  endif
  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)

  if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_data_in_file,grid%nmax
    call printErrMsg(option)   
  endif

  rank_mpi = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif

  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)
  call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)

  ! must initialize here to avoid error below when closing memory space
  memory_space_id = -1

  if (associated(indices)) then

    istart = indices(-1)
    iend = indices(0)

    dims = 0
    dims(1) = iend-istart
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    allocate(real_buffer(iend-istart))
    if (data_type == H5T_NATIVE_DOUBLE) then
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(integer_buffer_i4(iend-istart))
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      do i=1,iend-istart
        real_buffer(i) = real(integer_buffer_i4(i))
      enddo
      deallocate(integer_buffer_i4)
    endif
    ! must convert indices to zero based for VecSetValues
    allocate(indices0(iend-istart))
    indices0 = indices(1:iend-istart)-1
    call VecSetValues(natural_vec,iend-istart,indices0, &
                      real_buffer,INSERT_VALUES,ierr);CHKERRQ(ierr)
    deallocate(indices0)
    deallocate(real_buffer)

  endif

  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
  call DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec, &
                                     ONEDOF)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
! End of Default HDF5 Mechanism

#endif
! SCORPIO

end subroutine HDF5ReadArray

#endif

! ************************************************************************** !

subroutine HDF5QueryRegionDefinition(region, filename, option, &
     cell_ids_exists, face_ids_exists, vert_ids_exists)

  !
  ! Queries HDF5 to determine with region definition includes which groups:
  !
  ! cell_ids_exits = true if "Regions/<Region Name>/Cell Ids" exists
  ! face_ids_exits = true if "Regions/<Region Name>/Face Ids" exists
  ! vert_ids_exits = true if "Regions/<Region Name>/Vertex Ids" exists
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/21/2015
  !

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  use HDF5_Aux_module

  implicit none

  type(region_type) :: region
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool, intent (out) :: cell_ids_exists
  PetscBool, intent (out) :: face_ids_exists
  PetscBool, intent (out) :: vert_ids_exists

  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: string

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
#endif

  PetscBool :: grp_exists

#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! Open the Regions group
  string = 'Regions'
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  ! Open the Regions group
  string = trim(region%name)
  call h5gopen_f(grp_id,string,grp_id2,hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'HDF5 group "' // trim(region%name) // '" not found.'
    call printErrMsg(option)
  endif

  ! Querry region definition
  string = "Cell Ids"
  call h5lexists_f(grp_id2, string, cell_ids_exists, hdf5_err)

  string = "Face Ids"
  call h5lexists_f(grp_id2, string, face_ids_exists, hdf5_err)

  string = "Vertex Ids"
  call h5lexists_f(grp_id2, string, vert_ids_exists, hdf5_err)

  call h5gclose_f(grp_id2,hdf5_err)
  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif

end subroutine HDF5QueryRegionDefinition

! ************************************************************************** !

subroutine HDF5ReadRegionFromFile(grid,region,filename,option)
  ! 
  ! PETSC_HAVE_HDF5
  ! Reads a region from an hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/3/08
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  use HDF5_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(option_type), pointer :: option
  type(region_type) :: region
  character(len=MAXSTRINGLENGTH) :: filename

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
#endif

  PetscInt :: num_integers
  PetscBool :: grp_exists

#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  call PetscLogEventBegin(logging%event_region_read_hdf5,ierr);CHKERRQ(ierr)
                          
  ! initialize fortran hdf5 interface 
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! Open the Regions group
  string = 'Regions' 
  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)  
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  ! Open the Regions group
  string = trim(region%name)
  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)  
  call h5gopen_f(grp_id,string,grp_id2,hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'HDF5 group "' // trim(region%name) // '" not found.'
    call printErrMsg(option)
  endif

  ! Read Cell Ids
  string = "Cell Ids"

  ! Check if the region dataset has "Cell Ids" group
  call h5lexists_f(grp_id2,string,grp_exists,hdf5_err)
  if (.not.grp_exists) then
    option%io_buffer = 'HDF5 group: "Regions/' // trim(region%name) // &
      '/Cell Ids" not found.'
    call printErrMsg(option)
  endif

  call HDF5ReadIntegerArraySplit(option,grp_id2,string,ZERO_INTEGER, &
                                 region%cell_ids)
  region%def_type = DEFINED_BY_CELL_IDS
                            
  ! can't use size(region%cell_ids) alone as a null array returns a size of 1
  if (associated(region%cell_ids)) then
    region%num_cells = size(region%cell_ids)
  endif
  string = "Face Ids"
  ! Check if the region dataset has "Face Ids" group
  call h5lexists_f(grp_id2,string,grp_exists,hdf5_err)
  if (grp_exists) then
    option%io_buffer = 'Reading dataset: ' // trim(string)
    call printMsg(option)
    call HDF5ReadIntegerArraySplit(option,grp_id2,string, &
                                   region%num_cells,region%faces)
    region%def_type = DEFINED_BY_CELL_AND_FACE_IDS
  endif

  option%io_buffer = 'Closing group: ' // trim(region%name)
  call printMsg(option)  
  call h5gclose_f(grp_id2,hdf5_err)
  option%io_buffer = 'Closing group: Regions'
  call printMsg(option)   
  call h5gclose_f(grp_id,hdf5_err)
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

#endif
!PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_region_read_hdf5,ierr);CHKERRQ(ierr)

end subroutine HDF5ReadRegionFromFile

! ************************************************************************** !

subroutine HDF5ReadRegionDefinedByVertex(option,region,filename)
  ! 
  ! Reads a region from an hdf5 file defined by Vertex Ids
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/21/11
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  use HDF5_Aux_module
  use Grid_Unstructured_Cell_module
  use Utility_module, only : DeallocateArray

  implicit none

  !class(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(region_type) :: region
  type(region_sideset_type),pointer :: sideset
  character(len=MAXSTRINGLENGTH) :: filename

  ! local
  !type(option_type), pointer :: option
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: rank_mpi
  PetscInt :: remainder
  PetscInt :: istart, iend, ii, jj
  PetscInt, pointer :: int_buffer_1d(:)
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_buffer_2d(:,:)
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: string2

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: data_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: length(2), offset(2)
  integer(SIZE_T) :: string_size
#endif
  integer :: ndims_h5


#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted unstructured grids.")')
  call printErrMsg(option)
#else
  ! Initialize FORTRAN predefined datatypes
  call h5open_f(hdf5_err)

  ! Setup file access property with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)

#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

  ! Open the file collectively
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! Open dataset
  string = 'Regions/' // trim(region%name) // '/Vertex Ids'
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  if (hdf5_err /= 0) then
    string_size = MAXSTRINGLENGTH
    call h5fget_name_f(file_id,string2,string_size,hdf5_err)
    option%io_buffer = 'HDF5 dataset "' // trim(string) // '" not found &
      &in file "' // trim(string2) // '".'
    call printErrMsg(option)
  endif

  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id,data_space_id,hdf5_err)

  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id,ndims_h5,hdf5_err)
  if (ndims_h5 /= 2) then
    option%io_buffer='Dimension of '//string//' dataset in ' // trim(filename) // &
     ' is /= 2.'
  call printErrMsg(option)
  endif

  ! Allocate memory
  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))

  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5,hdf5_err)

  ! Create storage for sideset
  region%sideset => RegionCreateSideset()
  sideset => region%sideset

  sideset%nfaces = int(dims_h5(2)/option%mycommsize)
  remainder = int(dims_h5(2)) - sideset%nfaces*option%mycommsize
  if (option%myrank < remainder) sideset%nfaces = sideset%nfaces + 1

  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(sideset%nfaces,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                  MPI_SUM,option%mycomm,ierr)
  call MPI_Scan(sideset%nfaces,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                MPI_SUM,option%mycomm,ierr)

  ! Determine the length and offset of data to be read by each processor
  length(1) = dims_h5(1)
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart

  !
  rank_mpi = 2
  memory_space_id = -1

  ! Create data space for dataset
  call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)

  ! Select hyperslab
  call h5sselect_hyperslab_f(data_space_id,H5S_SELECT_SET_F,offset,length,hdf5_err)

  ! Initialize data buffer
  allocate(int_buffer_2d(length(1),length(2)))

  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F,hdf5_err)
#endif

  ! Read the dataset collectively
  call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,int_buffer_2d,&
       dims_h5,hdf5_err,memory_space_id,data_space_id)

  !
  ! Input data is list of Vertices
  !
  ! allocate array to store vertices for each cell
  region%def_type = DEFINED_BY_SIDESET_UGRID
  allocate(sideset%face_vertices(MAX_VERT_PER_FACE,sideset%nfaces))
  sideset%face_vertices = UNINITIALIZED_INTEGER

  do ii = 1,sideset%nfaces
     do jj = 1,int(dims_h5(1))
        sideset%face_vertices(jj,ii) = int_buffer_2d(jj,ii)
     enddo
  enddo
  ! cannot use DeallocateArray since int_array_i4 may not be PetscInt
  deallocate(int_buffer_2d)
  nullify(int_buffer_2d)

  deallocate(dims_h5)
  deallocate(max_dims_h5)

  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(data_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

#endif
! if defined(PETSC_HAVE_HDF5)

end subroutine HDF5ReadRegionDefinedByVertex

! ************************************************************************** !

subroutine HDF5ReadCellIndexedIntegerArray(realization,global_vec,filename, &
                                           group_name, &
                                           dataset_name,append_realization_id)
  ! 
  ! Reads an array of integer values from an
  ! hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/3/08; 02/18/09
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_Subsurface_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use HDF5_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(realization_subsurface_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscInt, allocatable :: integer_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_int_read_hdf5, &
                          ierr);CHKERRQ(ierr)
  
#if defined(SCORPIO)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
     option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
     call printMsg(option) 
  end if   
  filename = trim(filename) //CHAR(0)
  call fscorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, &
          file_id, ierr)

  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)

  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // "/Cell Ids" // CHAR(0)
  else
    string = "/Cell Ids" // CHAR(0)
  endif  
    
  call HDF5ReadIndices(grid,option,file_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif

  string = trim(dataset_name) // adjustl(trim(string))
  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // '/' // trim(string) // CHAR(0)
  else
    string = trim(string) // CHAR(0)
  endif  

  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   

  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
  
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read integer array.")') &
    tend-tstart
  call printMsg(option)  

  call DeallocateArray(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)   
    call fscorpio_close_file(file_id, option%ioread_group_id, ierr)
  endif

#else
! if SCORPIO is not defined

 ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option) 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open group if necessary
  if (len_trim(group_name) > 1) then
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)   
    call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  else
    grp_id = file_id
  endif

  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)


  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,grp_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
  
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read integer array.")') &
    tend-tstart
  call printMsg(option)  

  call DeallocateArray(indices)

  if (file_id /= grp_id) then
    option%io_buffer = 'Closing group: ' // trim(group_name)
    call printMsg(option)   
    call h5gclose_f(grp_id,hdf5_err)
  endif
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif  
! if SCORPIO is not defined

#endif
! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_int_read_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine HDF5ReadCellIndexedIntegerArray

! ************************************************************************** !

subroutine HDF5ReadCellIndexedRealArray(realization,global_vec,filename, &
                                        group_name, &
                                        dataset_name,append_realization_id)
  ! 
  ! Reads an array of real values from an hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/09, 02/18/09
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_Subsurface_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use HDF5_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(realization_subsurface_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscReal, allocatable :: real_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_real_read_hdf5, &
                          ierr);CHKERRQ(ierr)

#if defined(SCORPIO)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
     option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
     call printMsg(option) 
  end if   
  filename = trim(filename) //CHAR(0)
  call fscorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, &
          file_id, ierr)


! only new approach (old approach is removed)
  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)

  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // "/Cell Ids" // CHAR(0)
  else
    string = "/Cell Ids" // CHAR(0)
  endif  
    
  call HDF5ReadIndices(grid,option,file_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif

  string = trim(dataset_name) // adjustl(trim(string))
  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // '/' // trim(string) // CHAR(0)
  else
    string = trim(string) // CHAR(0)
  endif  

  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   

  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,H5T_NATIVE_DOUBLE)
  
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read real array.")') &
    tend-tstart
  call printMsg(option)  

  call DeallocateArray(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)   
    call fscorpio_close_file(file_id, option%ioread_group_id, ierr)
  endif

#else
! if SCORPIO is not defined

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option) 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open group if necessary
  if (len_trim(group_name) > 1) then
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)   
    call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  else
    grp_id = file_id
  endif

  ! Read Cell Ids
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,H5T_NATIVE_DOUBLE)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read real array")') &
    tend-tstart
  call printMsg(option)  

  call DeallocateArray(indices)

  if (file_id /= grp_id) then
    option%io_buffer = 'Closing group: ' // trim(group_name)
    call printMsg(option)   
    call h5gclose_f(grp_id,hdf5_err)
  endif
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

#endif  
! if SCORPIO is not defined

#endif
! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_real_read_hdf5, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine HDF5ReadCellIndexedRealArray

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine HDF5WriteStructDataSetFromVec(name,realization_base,vec,file_id,data_type)
  ! 
  ! Writes data from a PetscVec to HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  character(len=*) :: name
  class(realization_base_type) :: realization_base
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  
  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
!GEH - Structured Grid Dependence - Begin
  call HDF5WriteStructuredDataSet(trim(name), &
                                 vec_ptr,file_id,data_type,option, &
                                  grid%structured_grid%nx, &
                                  grid%structured_grid%ny, &
                                  grid%structured_grid%nz, &
                                  grid%structured_grid%nlx, &
                                  grid%structured_grid%nly, &
                                  grid%structured_grid%nlz, &
                                  grid%structured_grid%lxs, &
                                  grid%structured_grid%lys, &
                                  grid%structured_grid%lzs)
!GEH - Structured Grid Dependence - End
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine HDF5WriteStructDataSetFromVec

! ************************************************************************** !

subroutine HDF5WriteDataSetFromVec(name,option,vec,file_id,data_type)
  ! 
  ! This routine writes data from a PETSc Vec to HDF5 file for unstructured
  ! grids.
  ! subroutine HDF5WriteDataSetFromVec(name,realization,vec,file_id,data_type)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/31/12
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Subsurface_class
  use Grid_module
  use Option_module
  use Patch_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  character(len=32) :: name
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscReal, pointer :: vec_ptr(:)
  
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: istart
  PetscInt :: local_size,global_size,i
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: double_array(:)

  call VecGetLocalSize(vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetSize(vec,global_size,ierr);CHKERRQ(ierr)
  
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 1
  dims = 0
  dims(1) = global_size
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,name,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,name,data_type,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = local_size
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif


  if (data_type == H5T_NATIVE_DOUBLE) then
    allocate(double_array(local_size))  
    call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    do i=1,local_size
      double_array(i) = vec_ptr(i)
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    call DeallocateArray(double_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  if (data_type == H5T_NATIVE_INTEGER) then
    allocate(int_array(local_size))
    call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    do i=1,local_size
      int_array(i) = int(vec_ptr(i))
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    ! cannot use DeallocateArray since int_array_i4 may not be PetscInt
    deallocate(int_array)
    nullify(int_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  
end subroutine HDF5WriteDataSetFromVec

! ************************************************************************** !

subroutine HDF5ReadDataSetInVec(name, option, vec, file_id, data_type)
  ! 
  ! This routine reads data from HDF5 file in a PETSc Vec
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/16/2015
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Subsurface_class
  use Grid_module
  use Option_module
  use Patch_module
  use Utility_module, only : DeallocateArray

  implicit none

  character(len=32) :: name
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscReal, pointer :: vec_ptr(:)

  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: istart
  PetscInt :: local_size,global_size,i
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: double_array(:)

  call VecGetLocalSize(vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetSize(vec,global_size,ierr);CHKERRQ(ierr)

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank_mpi = 1
  dims = 0
  dims(1) = global_size
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,name,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,name,data_type,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = local_size
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif


  if (data_type == H5T_NATIVE_DOUBLE) then
    allocate(double_array(local_size))

    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,data_type,double_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    do i=1,local_size
      vec_ptr(i) = double_array(i)
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

    call DeallocateArray(double_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  if (data_type == H5T_NATIVE_INTEGER) then
    allocate(int_array(local_size))

    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    do i=1,local_size
      vec_ptr(i) = real(int_array(i))
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

    ! cannot use DeallocateArray since int_array_i4 may not be PetscInt
    deallocate(int_array)
    nullify(int_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine HDF5ReadDataSetInVec

#endif

end module HDF5_module



