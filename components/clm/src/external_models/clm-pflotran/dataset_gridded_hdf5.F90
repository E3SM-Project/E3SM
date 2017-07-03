module Dataset_Gridded_HDF5_class
 
  use Dataset_Common_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(dataset_common_hdf5_type) :: dataset_gridded_hdf5_type
    PetscBool :: is_cell_centered
    PetscInt :: data_dim
    PetscInt :: interpolation_method
    PetscReal, pointer :: origin(:)
    PetscReal, pointer :: extent(:)
    PetscReal, pointer :: discretization(:)
  end type dataset_gridded_hdf5_type
  
  PetscInt, parameter, public :: DIM_NULL = 0
  PetscInt, parameter, public :: DIM_X = 1
  PetscInt, parameter, public :: DIM_Y = 2
  PetscInt, parameter, public :: DIM_Z = 3
  PetscInt, parameter, public :: DIM_XY = 4
  PetscInt, parameter, public :: DIM_XZ = 5
  PetscInt, parameter, public :: DIM_YZ = 6
  PetscInt, parameter, public :: DIM_XYZ = 7
  
  PetscInt, parameter :: default_max_buffer_size = 10
  
  public :: DatasetGriddedHDF5Create, &
            DatasetGriddedHDF5Init, &
            DatasetGriddedHDF5Cast, &
            DatasetGriddedHDF5Load, &
            DatasetGriddedHDF5InterpolateReal, &
            DatasetGriddedHDF5Print, &
            DatasetGriddedHDF5Strip, &
            DatasetGriddedHDF5Destroy
  
contains

! ************************************************************************** !

function DatasetGriddedHDF5Create()
  ! 
  ! Creates global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 
  
  implicit none
  
  class(dataset_gridded_hdf5_type), pointer :: dataset

  class(dataset_gridded_hdf5_type), pointer :: DatasetGriddedHDF5Create
  
  allocate(dataset)
  call DatasetGriddedHDF5Init(dataset)

  DatasetGriddedHDF5Create => dataset
    
end function DatasetGriddedHDF5Create

! ************************************************************************** !

subroutine DatasetGriddedHDF5Init(this)
  ! 
  ! Initializes members of global dataset class
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 

  use Dataset_Base_class
  
  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  
  call DatasetCommonHDF5Init(this)
  this%is_cell_centered = PETSC_FALSE
  this%interpolation_method = INTERPOLATION_LINEAR
  this%data_dim = DIM_NULL
  nullify(this%origin)
  nullify(this%extent)
  nullify(this%discretization)
    
end subroutine DatasetGriddedHDF5Init

! ************************************************************************** !

function DatasetGriddedHDF5Cast(this)
  ! 
  ! Casts a dataset_base_type to dataset_gridded_hdf5_type
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/13
  ! 
  
  use Dataset_Base_class

  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_gridded_hdf5_type), pointer :: DatasetGriddedHDF5Cast
  
  nullify(DatasetGriddedHDF5Cast)
  select type (this)
    class is (dataset_gridded_hdf5_type)
      DatasetGriddedHDF5Cast => this
    class default
      ! to catch a class that is not gridded.
  end select
    
end function DatasetGriddedHDF5Cast

! ************************************************************************** !

subroutine DatasetGriddedHDF5Load(this,option)
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
  
  class(dataset_gridded_hdf5_type) :: this
  type(option_type) :: option
  
  if (DatasetCommonHDF5Load(this,option)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetGriddedHDF5ReadData(this,option)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetGriddedHDF5Load

#if defined(PETSC_HAVE_HDF5)    

! ************************************************************************** !

subroutine DatasetGriddedHDF5ReadData(this,option)
  ! 
  ! Read an hdf5 data into arrays
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/11, 05/29/13
  ! 

  use hdf5
  use Option_module
  use Units_module
  use Logging_module
  use HDF5_Aux_module
  use String_module
  
  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attribute_id
  integer(HID_T) :: ndims_h5
  integer(HID_T) :: atype_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(3)
  integer(HSIZE_T) :: offset(4), length(4), stride(4)
  integer(HSIZE_T) :: num_data_values
  integer(SIZE_T) size_t_int
  PetscInt :: i, temp_int, temp_array(4)
  PetscInt :: num_spatial_dims, time_dim, num_times
  PetscInt :: num_dims_in_h5_file, num_times_in_h5_file
  PetscMPIInt :: array_rank_mpi, mpi_int
  PetscBool :: attribute_exists
  PetscBool :: first_time
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  character(len=MAXWORDLENGTH) :: attribute_name, dataset_name, word

  call PetscLogEventBegin(logging%event_dataset_gridded_hdf5_read, &
                          ierr);CHKERRQ(ierr)

  first_time = (this%data_dim == DIM_NULL)

!#define TIME_DATASET
#ifdef TIME_DATASET
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
#endif

#define BROADCAST_DATASET
#ifdef BROADCAST_DATASET
  if (first_time .or. option%myrank == option%io_rank) then
#endif

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(this%filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
!geh: I don't believe that we ever need this
!  if (first_time) then
!  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
!  endif
#endif
  call HDF5OpenFileReadOnly(this%filename,file_id,prop_id,option)
  call h5pclose_f(prop_id,hdf5_err)

  ! the dataset is actually stored in a group.  the group contains
  ! a "data" dataset and optionally a "time" dataset.
  option%io_buffer = 'Opening group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gopen_f(file_id,this%hdf5_dataset_name,grp_id,hdf5_err)

  if (hdf5_err < 0) then
    option%io_buffer = 'A group named "' // trim(this%hdf5_dataset_name) // &
      '" not found in HDF5 file "' // trim(this%filename) // '".'
    call printErrMsg(option)  
  endif

  ! only want to read on first time through
  if (this%data_dim == DIM_NULL) then
    ! read in attributes if they exist
    attribute_name = "Dimension"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim = 1
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdf5_err)
      size_t_int = MAXWORDLENGTH
      call h5tset_size_f(atype_id,size_t_int,hdf5_err)
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,atype_id,word,attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
      call h5tclose_f(atype_id,hdf5_err)
      ! set dimensionality of dataset
      call DatasetGriddedHDF5SetDimension(this,word)
    else
      option%io_buffer = &
        'Dimension attribute must be included in hdf5 dataset file.'
      call printErrMsg(option)
    endif
    attribute_name = "Discretization"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGriddedHDF5GetNDimensions(this)
      allocate(this%discretization(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,this%discretization, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    else
      option%io_buffer = &
        '"Discretization" attribute must be included in GRIDDED hdf5 ' // &
        'dataset file.'
      call printErrMsg(option)
    endif
    attribute_name = "Origin"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGriddedHDF5GetNDimensions(this)
      allocate(this%origin(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,this%origin, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    endif
    attribute_name = "Cell Centered"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      this%is_cell_centered = PETSC_TRUE
    endif
    attribute_name = "Interpolation Method"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdf5_err)
      size_t_int = MAXWORDLENGTH
      call h5tset_size_f(atype_id,size_t_int,hdf5_err)
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,atype_id,word,attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
      call h5tclose_f(atype_id,hdf5_err)
      call StringToUpper(word)
      select case(trim(word))
        case('STEP')
          this%interpolation_method = INTERPOLATION_STEP
        case('LINEAR')
          this%interpolation_method = INTERPOLATION_LINEAR
        case default
          option%io_buffer = '"Interpolation Method" not recognized in ' // &
            'Gridded HDF5 Dataset "' // trim(this%name) // '".'
          call printErrMsg(option)
      end select
    endif
    ! this%max_buffer_size is initially set to UNINITIALIZED_INTEGER to force initializaion
    ! either here, or in the reading of the dataset block.
    attribute_name = "Max Buffer Size"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists .and. this%max_buffer_size < 0) then
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      attribute_dim(1) = 1
      call h5aread_f(attribute_id,H5T_NATIVE_INTEGER,this%max_buffer_size, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    else if (this%max_buffer_size < 0) then
      this%max_buffer_size = default_max_buffer_size
    endif
  endif ! this%data_dim == DIM_NULL

#ifdef BROADCAST_DATASET
  endif
#endif

  num_spatial_dims = DatasetGriddedHDF5GetNDimensions(this)
  
  ! num_times and time_dim must be calcualted by all processes; does not 
  ! require communication  
  time_dim = -1
  num_times = 1
  if (associated(this%time_storage)) then
    num_times = this%time_storage%max_time_index
    time_dim = num_spatial_dims + 1
  endif
  
#ifdef BROADCAST_DATASET
  if (first_time .or. option%myrank == option%io_rank) then
#endif
  ! open the "data" dataset
  dataset_name = 'Data'
  if (this%realization_dependent) then
    write(word,'(i9)') option%id
    dataset_name = trim(dataset_name) // trim(adjustl(word))
  endif
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  if (hdf5_err < 0) then
    option%io_buffer = 'A dataset named "Data" not found in HDF5 file "' // &
      trim(this%filename) // '".'
    call printErrMsg(option)  
  endif
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! get dataset dimensions
  if (.not.associated(this%dims)) then
    call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)
    allocate(dims_h5(ndims_h5))
    allocate(max_dims_h5(ndims_h5))
    call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
    if (associated(this%time_storage)) then
      ! if dataset is time dependent, need to remove that time index from 
      ! dimensions and decrement rank (we don't want to include time)
      this%rank = ndims_h5-1
      ! the first dimension of dims_h5 is the time dimension
      num_times_in_h5_file = int(dims_h5(1))
    else
      this%rank = ndims_h5
      num_times_in_h5_file = 0
    endif
    allocate(this%dims(this%rank))
    ! have to invert dimensions, but do not count the time dimension
    do i = 1, this%rank
      this%dims(i) = int(dims_h5(ndims_h5-i+1))
    enddo
    deallocate(dims_h5)
    deallocate(max_dims_h5) 
    
    allocate(this%extent(num_spatial_dims))
    do i = 1, num_spatial_dims
      temp_int = this%dims(i)
      if (.not.this%is_cell_centered) then
        temp_int = temp_int - 1
      endif
      this%extent(i) = this%origin(i) + this%discretization(i) * temp_int
    enddo
  
    call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
  
    temp_int = this%dims(1)
    do i = 2, num_spatial_dims
      temp_int = temp_int * this%dims(i)
    enddo
    if (num_data_values/num_times /= temp_int) then
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
      allocate(this%rarray(temp_int))
    endif
    this%rarray = 0.d0
    if (associated(this%time_storage) .and. .not.associated(this%rbuffer)) then
      ! buffered array
      allocate(this%rbuffer(size(this%rarray)*this%max_buffer_size))
      this%rbuffer = 0.d0
    endif
  endif

#ifdef BROADCAST_DATASET
  endif
#endif

#ifdef TIME_DATASET
  call MPI_Barrier(option%mycomm,ierr)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up dataset ",a,".")') &
    tend-tstart, trim(this%hdf5_dataset_name) // ' (' // &
    trim(option%group_prefix) // ')'
  if (option%myrank == option%io_rank) then
    print *, trim(option%io_buffer)
  endif
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
#endif

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
  this%buffer_nslice = min(this%max_buffer_size,(num_times-this%buffer_slice_offset))
  if (time_dim > 0) then
    length(1) = size(this%rarray) * this%buffer_nslice
  else
    length(1) = size(this%rarray)
  endif

#ifdef BROADCAST_DATASET
  if (option%myrank == option%io_rank) then
#endif

  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    

  length = 1
  stride = 1
  offset = 0
  do i = 1, num_spatial_dims
    length(i) = this%dims(i)
  enddo
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

#ifdef BROADCAST_DATASET
  endif !if (option%myrank == option%io_rank) then
  if (associated(this%rbuffer)) then
    mpi_int = size(this%rbuffer)
    call MPI_Bcast(this%rbuffer,mpi_int,MPI_DOUBLE_PRECISION,option%io_rank, &
                   option%mycomm,ierr)
  else
    mpi_int = size(this%rarray)
    call MPI_Bcast(this%rarray,mpi_int,MPI_DOUBLE_PRECISION,option%io_rank, &
                   option%mycomm,ierr)
  endif
#endif
  
  call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)

#ifdef BROADCAST_DATASET
  if (first_time .or. option%myrank == option%io_rank) then
#endif  
  option%io_buffer = 'Closing group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#ifdef BROADCAST_DATASET
  endif
#endif

#ifdef TIME_DATASET
  call MPI_Barrier(option%mycomm,ierr)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer,'(f6.2," Seconds to read dataset ",a,".")') &
    tend-tstart, trim(this%hdf5_dataset_name) // ' (' // &
    trim(option%group_prefix) // ')'
  if (option%myrank == option%io_rank) then
    print *, trim(option%io_buffer)
  endif
#endif

  call PetscLogEventEnd(logging%event_dataset_gridded_hdf5_read, &
                        ierr);CHKERRQ(ierr)
                          
end subroutine DatasetGriddedHDF5ReadData

#endif

! ************************************************************************** !

subroutine DatasetGriddedHDF5SetDimension(this,word)
  ! 
  ! Sets the dimension of the dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/11, 05/29/13
  ! 

  use String_module

  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  character(len=MAXWORDLENGTH) :: word

  call StringToUpper(word)
  select case(word)
    case('X')
      this%data_dim = DIM_X
    case('Y')
      this%data_dim = DIM_Y
    case('Z')
      this%data_dim = DIM_Z
    case('XY')
      this%data_dim = DIM_XY
    case('XZ')
      this%data_dim = DIM_XZ
    case('YZ')
      this%data_dim = DIM_YZ
    case('XYZ')
      this%data_dim = DIM_XYZ
  end select
      
end subroutine DatasetGriddedHDF5SetDimension

! ************************************************************************** !

function DatasetGriddedHDF5GetDimensionString(this)
  ! 
  ! Returns a string describing dimension
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/11, 05/29/13, 10/22/13
  ! 

  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  
  character(len=MAXWORDLENGTH) :: DatasetGriddedHDF5GetDimensionString

  select case(this%data_dim)
    case(DIM_X)
      DatasetGriddedHDF5GetDimensionString = 'X'
    case(DIM_Y)
      DatasetGriddedHDF5GetDimensionString = 'Y'
    case(DIM_Z)
      DatasetGriddedHDF5GetDimensionString = 'Z'
    case(DIM_XY)
      DatasetGriddedHDF5GetDimensionString = 'XY'
    case(DIM_XZ)
      DatasetGriddedHDF5GetDimensionString = 'XZ'
    case(DIM_YZ)
      DatasetGriddedHDF5GetDimensionString = 'YZ'
    case(DIM_XYZ)
      DatasetGriddedHDF5GetDimensionString = 'XYZ'
  end select
      
end function DatasetGriddedHDF5GetDimensionString

! ************************************************************************** !

function DatasetGriddedHDF5GetNDimensions(this)
  ! 
  ! Returns the number of dimensions
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/11, 05/29/13
  ! 

  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  
  PetscInt :: DatasetGriddedHDF5GetNDimensions

  select case(this%data_dim)
    case(DIM_X,DIM_Y,DIM_Z)
      DatasetGriddedHDF5GetNDimensions = ONE_INTEGER
    case(DIM_XY,DIM_XZ,DIM_YZ)
      DatasetGriddedHDF5GetNDimensions = TWO_INTEGER
    case(DIM_XYZ)
      DatasetGriddedHDF5GetNDimensions = THREE_INTEGER
    case default
      DatasetGriddedHDF5GetNDimensions = ZERO_INTEGER
  end select
      
end function DatasetGriddedHDF5GetNDimensions

! ************************************************************************** !

subroutine DatasetGriddedHDF5InterpolateReal(this,xx,yy,zz,real_value,option)
  ! 
  ! Interpolates data from the dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/29/13
  ! 

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  PetscReal, intent(in) :: xx, yy, zz
  PetscReal :: real_value
  type(option_type) :: option
  
  PetscInt :: spatial_interpolation_method
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  PetscReal :: x1, x2, y1, y2, z1
  PetscReal :: v1, v2, v3, v4
  PetscInt :: index
  PetscInt :: ii, jj, kk
  PetscInt :: i_upper, j_upper, k_upper
  PetscReal :: dx, dy, dz
  PetscInt :: nx, ny
  PetscBool :: lerr
  character(len=MAXWORDLENGTH) :: word
  
  call DatasetGriddedHDF5GetIndices(this,xx,yy,zz,i,j,k,x,y,z)
  
  ! in the below, i,j,k,xx,yy,zz to not reflect the 
  ! coordinates of the problem domain in 3D.  They
  ! are transfored to the dimensions of the dataset
  lerr = PETSC_FALSE
  select case(this%interpolation_method)
    case(INTERPOLATION_STEP)
      select case(this%data_dim)
        case(DIM_X,DIM_Y,DIM_Z)
          if (this%is_cell_centered) then
            i_upper = i
          else
            i_upper = i+1
          endif
          if (i < 1 .or. i_upper > this%dims(1)) then 
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_X)
                write(option%io_buffer,*) 'X value (', xx, &
                  ') outside of X bounds (', this%origin(1), this%extent(1), &
                  '), i =' 
              case(DIM_Y)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
              case(DIM_Z)
                write(option%io_buffer,*) 'Z value (', zz, &
                  ') outside of Z bounds (', this%origin(3), this%extent(3), &
                  '), k =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(word) // ' for gridded dataset "' // trim(this%name) // '".'
            call printErrMsgByRank(option)
          endif
          index = i
          if (.not.this%is_cell_centered) then
            dx = this%discretization(1)
            x1 = this%origin(1) + (i-1)*dx
            if ((x-x1) / dx > 0.5) then
              index = i+1
            endif
          endif
        case(DIM_XY,DIM_XZ,DIM_YZ)
          if (this%is_cell_centered) then
            i_upper = i
            j_upper = j
          else
            i_upper = i+1
            j_upper = j+1
          endif
          if (i < 1 .or. i_upper > this%dims(1)) then
            lerr = PETSC_TRUE
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY,DIM_XZ)
                write(option%io_buffer,*) 'X value (', xx, &
                  ') outside of X bounds (', this%origin(1), this%extent(1), &
                  '), i =' 
              case(DIM_YZ)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // word
            call printMsgByRank(option)
          endif
          if (j < 1 .or. j_upper > this%dims(2)) then
            lerr = PETSC_TRUE
            write(word,*) j
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
              case(DIM_YZ,DIM_XZ)
                write(option%io_buffer,*) 'Z value (', zz, &
                  ') outside of Z bounds (', this%origin(3), this%extent(3), &
                  '), k =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // word
            call printMsgByRank(option)
          endif
          if (lerr) then
            word = this%name
            option%io_buffer = 'gridded dataset "' // trim(word) // &
                               '" out of bounds.'
            call printErrMsgByRank(option)
          endif
          ii = i
          jj = j
          nx = this%dims(1)
          if (.not.this%is_cell_centered) then
            dx = this%discretization(1)
            dy = this%discretization(2)
            x1 = this%origin(1) + (i-1)*dx
            y1 = this%origin(2) + (j-1)*dy
            if ((x-x1) / dx > 0.5d0) then
              ii = i+1
            endif
            if ((y-y1) / dy > 0.5d0) then
              jj = j+1
            endif
          endif
          index = ii + (jj-1)*nx
        case(DIM_XYZ)
          if (this%is_cell_centered) then
            i_upper = i
            j_upper = j
            k_upper = k
          else
            i_upper = i+1
            j_upper = j+1
            k_upper = k+1
          endif          
          if (i < 1 .or. i_upper > this%dims(1)) then
            lerr = PETSC_TRUE
            write(word,*) i
            word = adjustl(word)
            write(option%io_buffer,*) 'X value (', xx, &
              ') outside of X bounds (', this%origin(1), this%extent(1), &
              '), i = ' // trim(word)
            call printMsgByRank(option)
          endif
          if (j < 1 .or. j_upper > this%dims(2)) then
            lerr = PETSC_TRUE
            write(word,*) j
            word = adjustl(word)
            write(option%io_buffer,*) 'Y value (', yy, &
              ') outside of Y bounds (', this%origin(2), this%extent(2), &
              '), j = ' // trim(word)
            call printMsgByRank(option)
          endif
          if (k < 1 .or. k_upper > this%dims(3)) then
            lerr = PETSC_TRUE
            write(word,*) k
            word = adjustl(word)
            write(option%io_buffer,*) 'Z value (', zz, &
              ') outside of Z bounds (', this%origin(3), this%extent(3), &
              '), k = ' // trim(word)
            call printMsgByRank(option)
          endif
          if (lerr) then
            word = this%name
            option%io_buffer = 'gridded dataset "' // trim(word) // &
                               '" out of bounds.'
            call printErrMsgByRank(option)
          endif
          ii = i
          jj = j
          kk = k
          nx = this%dims(1)
          ny = this%dims(2)
          if (.not.this%is_cell_centered) then
            dx = this%discretization(1)
            dy = this%discretization(2)
            dz = this%discretization(3)
            x1 = this%origin(1) + (i-1)*dx
            y1 = this%origin(2) + (j-1)*dy
            z1 = this%origin(3) + (k-1)*dz
            if ((x-x1) / dx > 0.5d0) then
              ii = i+1
            endif
            if ((y-y1) / dy > 0.5d0) then
              jj = j+1
            endif
            if ((z-z1) / dz > 0.5d0) then
              kk = k+1
            endif
          endif
          index = ii + (jj-1)*nx + (kk-1)*nx*ny
      end select
      real_value = this%rarray(index)
    case(INTERPOLATION_LINEAR)
      select case(this%data_dim)
        case(DIM_X,DIM_Y,DIM_Z)
          if (i < 1 .or. i+1 > this%dims(1)) then 
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_X)
                write(option%io_buffer,*) 'X value (', xx, &
                  ') outside of X bounds (', this%origin(1), this%extent(1), &
                  '), i =' 
              case(DIM_Y)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
              case(DIM_Z)
                write(option%io_buffer,*) 'Z value (', zz, &
                  ') outside of Z bounds (', this%origin(3), this%extent(3), &
                  '), k =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(word) // ' for gridded dataset "' // trim(this%name) // '".'
            call printErrMsgByRank(option)
          endif
          dx = this%discretization(1)
          x1 = this%origin(1) + (i-1)*dx
          if (this%is_cell_centered) x1 = x1 + 0.5d0*dx
          v1 = this%rarray(i)
          v2 = this%rarray(i+1)
          real_value = v1 + (x-x1)/dx*(v2-v1)
        case(DIM_XY,DIM_XZ,DIM_YZ)
          if (i < 1 .or. i+1 > this%dims(1)) then
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY,DIM_XZ)
                write(option%io_buffer,*) 'X value (', xx, &
                  ') outside of X bounds (', this%origin(1), this%extent(1), &
                  '), i =' 
              case(DIM_YZ)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(word) // ' for gridded dataset "' // trim(this%name) // '".'
            call printErrMsgByRank(option)
          endif
          if (j < 1 .or. j+1 > this%dims(2)) then
            write(word,*) j
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY)
                write(option%io_buffer,*) 'Y value (', yy, &
                  ') outside of Y bounds (', this%origin(2), this%extent(2), &
                  '), j =' 
              case(DIM_YZ,DIM_XZ)
                write(option%io_buffer,*) 'Z value (', zz, &
                  ') outside of Z bounds (', this%origin(3), this%extent(3), &
                  '), k =' 
            end select
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(word) // ' for gridded dataset "' // trim(this%name) // '".'
            call printErrMsgByRank(option)
          endif
          dx = this%discretization(1)
          dy = this%discretization(2)
          nx = this%dims(1)

          x1 = this%origin(1) + (i-1)*dx
          if (this%is_cell_centered) x1 = x1 + 0.5d0*dx
          x2 = x1 + dx
          
          index = i + (j-1)*nx
          v1 = this%rarray(index)
          v2 = this%rarray(index+1)
          
          y1 = this%origin(2) + (j-1)*dy
          if (this%is_cell_centered) y1 = y1 + 0.5d0*dy
          y2 = y1 + dy
          
           ! really (j1-1+1)
          index = i + j*nx
          v3 = this%rarray(index)
          v4 = this%rarray(index+1)
          
          real_value = InterpolateBilinear(x,y,x1,x2,y1,y2,v1,v2,v3,v4)
        case(DIM_XYZ)
          option%io_buffer = 'Trilinear interpolation not yet supported'
          call printErrMsgByRank(option)
      end select
  end select
  
end subroutine DatasetGriddedHDF5InterpolateReal

! ************************************************************************** !

subroutine DatasetGriddedHDF5GetIndices(this,xx,yy,zz,i,j,k,x,y,z)
  ! 
  ! Returns bounding indices for point in dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11, 05/29/13
  ! 

  implicit none

  class(dataset_gridded_hdf5_type) :: this
  PetscReal, intent(in) :: xx, yy, zz
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  
  PetscReal :: discretization_offset
  PetscInt :: upper_index_offset
  PetscReal :: tol
  PetscReal, parameter :: tolerance_scale = 1.d-3

  select case(this%data_dim)
    ! since these are 1D array, always use first dimension
    case(DIM_X)
      x = xx
    case(DIM_Y)
      x = yy
    case(DIM_XY)
      x = xx
      y = yy
    case(DIM_XYZ)
      x = xx
      y = yy
      z = zz
    case(DIM_Z)
      x = zz
    case(DIM_XZ)
      x = xx
      y = zz
    case(DIM_YZ)
      x = yy
      y = zz
  end select

  upper_index_offset = 1 
  if (this%is_cell_centered) then
    select case(this%interpolation_method)
      case(INTERPOLATION_STEP)
        discretization_offset = 1.d0
        upper_index_offset = 0
      case(INTERPOLATION_LINEAR)
        discretization_offset = 0.5d0
    end select
    i = int((x - this%origin(1))/ &
            this%discretization(1) + discretization_offset)
!    i = max(1,min(i,this%dims(1)-1))
    if (this%data_dim > DIM_Z) then ! at least 2D
      j = int((y - this%origin(2))/ &
              this%discretization(2) + discretization_offset)
!      j = max(1,min(j,this%dims(2)-1))
    endif
    if (this%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - this%origin(3))/ &
              this%discretization(3) + discretization_offset)
!      k = max(1,min(k,this%dims(3)-1))
    endif
  else
    select case(this%interpolation_method)
      case(INTERPOLATION_STEP)
        discretization_offset = 1.5d0
      case(INTERPOLATION_LINEAR)
        discretization_offset = 1.d0
    end select
    i = int((x - this%origin(1))/ &
            this%discretization(1) + discretization_offset)
    if (this%data_dim > DIM_Z) then ! at least 2D
      j = int((y - this%origin(2))/ &
              this%discretization(2) + discretization_offset)
    endif
    if (this%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - this%origin(3))/ &
              this%discretization(3) + discretization_offset)
    endif
  endif
  
  ! if indices are out of bounds, check if on boundary and reset index
  !geh: the tolerance allows one to go outside the bounds 
  if (i < 1 .or. i+upper_index_offset > this%dims(1)) then
    tol = this%discretization(1) * tolerance_scale
    if (x >= this%origin(1)-tol .and. x <= this%extent(1)+tol) then
      i = min(max(i,1),this%dims(1)-upper_index_offset)
    endif
  endif
  if (this%data_dim > DIM_Z) then ! at least 2D
    if (j < 1 .or. j+upper_index_offset > this%dims(2)) then
      tol = this%discretization(2) * tolerance_scale
      if (y >= this%origin(2)-tol .and. y <= this%extent(2)+tol) then
        j = min(max(j,1),this%dims(2)-upper_index_offset)
      endif
    endif
  endif  
  if (this%data_dim > DIM_YZ) then ! at least 2D
    if (k < 1 .or. k+upper_index_offset > this%dims(3)) then
      tol = this%discretization(3) * tolerance_scale
      if (z >= this%origin(3)-tol .and. z <= this%extent(3)+tol) then
        k = min(max(k,1),this%dims(3)-upper_index_offset)
      endif
    endif
  endif  
  
end subroutine DatasetGriddedHDF5GetIndices

! ************************************************************************** !

subroutine DatasetGriddedHDF5Print(this,option)
  ! 
  ! Prints dataset info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/13
  ! 

  use Option_module

  implicit none
  
  class(dataset_gridded_hdf5_type), target :: this
  type(option_type) :: option

  class(dataset_common_hdf5_type), pointer :: dataset_hdf5

  dataset_hdf5 => this
  call DatasetCommonHDF5Print(this,option)
  
  write(option%fid_out,'(10x,''Grid Dimension: '',a)') &
    trim(DatasetGriddedHDF5GetDimensionString(this))
  if (this%is_cell_centered) then
    write(option%fid_out,'(10x,''Is cell-centered?: yes'')') 
  else
    write(option%fid_out,'(10x,''Is cell-centered?: no'')')
  endif
  write(option%fid_out,'(10x,''Origin: '',3es12.4)') this%origin(:)
  write(option%fid_out,'(10x,''Discretization: '',3es12.4)') &
    this%discretization(:)
  
end subroutine DatasetGriddedHDF5Print

! ************************************************************************** !

subroutine DatasetGriddedHDF5Strip(this)
  ! 
  ! Strips allocated objects within XYZ dataset object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_gridded_hdf5_type) :: this
  
  call DatasetCommonHDF5Strip(this)
  
  call DeallocateArray(this%origin)
  call DeallocateArray(this%extent)
  call DeallocateArray(this%discretization)
  
end subroutine DatasetGriddedHDF5Strip

! ************************************************************************** !

subroutine DatasetGriddedHDF5Destroy(this)
  ! 
  ! Destroys a dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/13
  ! 

  implicit none
  
  class(dataset_gridded_hdf5_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetGriddedHDF5Strip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetGriddedHDF5Destroy

end module Dataset_Gridded_HDF5_class
