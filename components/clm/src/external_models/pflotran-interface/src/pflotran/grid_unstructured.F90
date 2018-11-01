module Grid_Unstructured_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private 

#if defined(SCORPIO)
  include "scorpiof.h"
#endif

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: UGridRead, &
#if defined(PETSC_HAVE_HDF5)
            UGridReadHDF5, &
#endif
#if defined(SCORPIO)
            UGridReadHDF5PIOLib, &
#endif
            UGridReadSurfGrid, &
#if defined(PETSC_HAVE_HDF5)
            UGridReadHDF5SurfGrid, &
#endif
            UGridDecompose, &
            UGridComputeInternConnect, &
            UGridPopulateConnection, &
            UGridComputeCoord, &
            UGridComputeVolumes, &
            UGridComputeAreas, &
            UGridComputeQuality, &
            UGridGetCellFromPoint, &
            UGridGetCellsInRectangle, &
            UGridEnsureRightHandRule, &
            UGridMapSideSet, &
            UGridMapSideSet2, &
            UGridGrowStencilSupport, &
            UGridMapBoundFacesInPolVol

contains

! ************************************************************************** !

subroutine UGridRead(unstructured_grid,filename,option)
  ! 
  ! Reads an unstructured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/09
  ! 

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: num_cells_local_save
  PetscInt :: num_cells_local
  PetscInt :: num_vertices_local_save
  PetscInt :: num_vertices_local
  PetscInt :: num_to_read
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscReal, allocatable :: vertex_coordinates(:,:)

  PetscInt :: icell, ivertex, idir, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename,option)

  ! initial guess is 8 vertices per cell
  unstructured_grid%max_nvert_per_cell = 8

! Format of unstructured grid file
! type: H=hexahedron, T=tetrahedron, W=wedge, P=pyramid
! vertn(H) = 8
! vertn(T) = 4
! vertn(W) = 6
! vertn(P) = 5
! -----------------------------------------------------------------
! num_cells num_vertices  (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 1 (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 2
! ...
! ...
! type vert1 vert2 vert3 ... vertn  ! for cell num_cells
! xcoord ycoord zcoord ! coordinates of vertex 1 (real)
! xcoord ycoord zcoord ! coordinates of vertex 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
! -----------------------------------------------------------------

  hint = 'Unstructured Grid'

  call InputReadPflotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,hint)  

  ! read num_cells
  call InputReadInt(input,option,unstructured_grid%nmax)
  call InputErrorMsg(input,option,'number of cells',hint)
  ! read num_vertices
  call InputReadInt(input,option,unstructured_grid%num_vertices_global)
  call InputErrorMsg(input,option,'number of vertices',hint)

  ! divide cells across ranks
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices(unstructured_grid%max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = UNINITIALIZED_INTEGER

  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(unstructured_grid%max_nvert_per_cell, &
                            num_cells_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = UNINITIALIZED_INTEGER
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do icell = 1, num_to_read
        ! read in the vertices defining the grid cell
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'element type',hint)
        call StringToUpper(word)
        select case(word)
          case('H')
            num_vertices = 8
          case('W')
            num_vertices = 6
          case('P')
            num_vertices = 5
          case('T')
            num_vertices = 4
          case('Q')
            num_vertices = 4
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,icell))
          call InputErrorMsg(input,option,'vertex id',hint)
        enddo
      enddo
      
      ! if the cells reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_cells_local
        string = trim(adjustl(string)) // ' cells stored on p0'
        print *, trim(string)
#endif
        unstructured_grid%cell_vertices(:,1:num_cells_local) = &
          temp_int_array(:,1:num_cells_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' cells sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*unstructured_grid%max_nvert_per_cell
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_cells_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' cells received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_cells_local*unstructured_grid%max_nvert_per_cell
    call MPI_Recv(unstructured_grid%cell_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif


  ! divide vertices across ranks
  num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                 num_vertices_local + 1

  allocate(vertex_coordinates(3,num_vertices_local))
  vertex_coordinates = 0.d0

  ! just like above, but this time for vertex coordinates
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(3,num_vertices_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_vertices_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do ivertex = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        do idir = 1, 3
          call InputReadDouble(input,option,temp_real_array(idir,ivertex))
          call InputErrorMsg(input,option,'vertex coordinate',hint)
        enddo
      enddo
      
      if (irank == option%io_rank) then
        vertex_coordinates(:,1:num_vertices_local) = &
          temp_real_array(:,1:num_vertices_local)
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    int_mpi = num_vertices_local*3
    call MPI_Recv(vertex_coordinates, &
                  int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif
  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ivertex = 1, num_vertices_local
    unstructured_grid%vertices(ivertex)%id = 0
    unstructured_grid%vertices(ivertex)%x = vertex_coordinates(1,ivertex)
    unstructured_grid%vertices(ivertex)%y = vertex_coordinates(2,ivertex)
    unstructured_grid%vertices(ivertex)%z = vertex_coordinates(3,ivertex)
  enddo
  deallocate(vertex_coordinates)

  unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local

  call InputDestroy(input)

end subroutine UGridRead

! ************************************************************************** !

subroutine UGridReadSurfGrid(unstructured_grid,filename,surf_filename,option)
  ! 
  ! UGridRead: Reads an unstructured grid
  ! 
  ! Author: Gautam Bisht
  ! Date: 01/09/2012
  ! 

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: surf_filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: num_cells_local_save
  PetscInt :: num_cells_local
  PetscInt :: num_vertices_local_save
  PetscInt :: num_vertices_local
  PetscInt :: num_to_read
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscReal, allocatable :: vertex_coordinates(:,:)

  PetscInt :: icell, ivertex, idir, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename,option)

  ! initial guess is 8 vertices per cell
  unstructured_grid%max_nvert_per_cell = 8

! Format of unstructured grid file
! type: H=hexahedron, T=tetrahedron, W=wedge, P=pyramid
! vertn(H) = 8
! vertn(T) = 4
! vertn(W) = 6
! vertn(P) = 5
! -----------------------------------------------------------------
! num_cells num_vertices  (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 1 (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 2
! ...
! ...
! type vert1 vert2 vert3 ... vertn  ! for cell num_cells
! xcoord ycoord zcoord ! coordinates of vertex 1 (real)
! xcoord ycoord zcoord ! coordinates of vertex 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
! -----------------------------------------------------------------

  hint = 'Unstructured Grid'

  call InputReadPflotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,hint)  

  ! read num_cells
  call InputReadInt(input,option,unstructured_grid%nmax)
  call InputErrorMsg(input,option,'number of cells',hint)
  ! read num_vertices
  call InputReadInt(input,option,unstructured_grid%num_vertices_global)
  call InputErrorMsg(input,option,'number of vertices',hint)

  ! divide cells across ranks
  !num_cells_local = unstructured_grid%nmax/option%mycommsize 
  !num_cells_local_save = num_cells_local
  !remainder = unstructured_grid%nmax - &
  !            num_cells_local*option%mycommsize
  !if (option%myrank < remainder) num_cells_local = &
  !                               num_cells_local + 1

  ! allocate array to store vertices for each cell
  !allocate(unstructured_grid%cell_vertices(unstructured_grid%max_nvert_per_cell, &
  !                                           num_cells_local))
  !unstructured_grid%cell_vertices = UNINITIALIZED_INTEGER

  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(unstructured_grid%max_nvert_per_cell, &
                            unstructured_grid%nmax))
    ! read for other processors
    temp_int_array = UNINITIALIZED_INTEGER
    num_to_read = unstructured_grid%nmax
    do icell = 1, num_to_read
      ! read in the vertices defining the grid cell
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,hint)  
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'element type',hint)
      call StringToUpper(word)
      select case(word)
        case('H')
          num_vertices = 8
        case('W')
          num_vertices = 6
        case('P')
          num_vertices = 5
        case('T')
          num_vertices = 4
        case('Q')
          num_vertices = 4
      end select
      do ivertex = 1, num_vertices
        call InputReadInt(input,option,temp_int_array(ivertex,icell))
        call InputErrorMsg(input,option,'vertex id',hint)
      enddo
    enddo
  endif


  ! divide vertices across ranks
  num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                 num_vertices_local + 1

  allocate(vertex_coordinates(3,num_vertices_local))
  vertex_coordinates = 0.d0

  ! just like above, but this time for vertex coordinates
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(3,num_vertices_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_vertices_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do ivertex = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        do idir = 1, 3
          call InputReadDouble(input,option,temp_real_array(idir,ivertex))
          call InputErrorMsg(input,option,'vertex coordinate',hint)
        enddo
      enddo
      
      if (irank == option%io_rank) then
        vertex_coordinates(:,1:num_vertices_local) = &
          temp_real_array(:,1:num_vertices_local)
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    int_mpi = num_vertices_local*3
    call MPI_Recv(vertex_coordinates, &
                  int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif
  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ivertex = 1, num_vertices_local
    unstructured_grid%vertices(ivertex)%id = 0
    unstructured_grid%vertices(ivertex)%x = vertex_coordinates(1,ivertex)
    unstructured_grid%vertices(ivertex)%y = vertex_coordinates(2,ivertex)
    unstructured_grid%vertices(ivertex)%z = vertex_coordinates(3,ivertex)
  enddo
  deallocate(vertex_coordinates)

  !unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local
  
  call InputDestroy(input)
  if (option%myrank == option%io_rank) deallocate(temp_int_array)


  input => InputCreate(fileid,surf_filename,option)
  call InputReadPflotranString(input,option)
  string = 'unstructured sideset'
  call InputReadStringErrorMsg(input,option,hint)  

  ! read num_cells
  call InputReadInt(input,option,unstructured_grid%nmax)
  call InputErrorMsg(input,option,'number of cells',hint)

  ! divide cells across ranks
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  ! allocate array to store vertices for each faces
  allocate(unstructured_grid%cell_vertices(unstructured_grid%max_nvert_per_cell, &
                                 num_cells_local))
  unstructured_grid%cell_vertices = UNINITIALIZED_INTEGER

  ! for now, read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(unstructured_grid%max_nvert_per_cell, &
                            num_cells_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = UNINITIALIZED_INTEGER
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1

      do icell = 1, num_to_read
        ! read in the vertices defining the cell face
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'element type',hint)
        call StringToUpper(word)
        select case(word)
          case('Q')
            num_vertices = 4
          case('T')
            num_vertices = 3
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,icell))
          call InputErrorMsg(input,option,'vertex id',hint)
        enddo
      enddo

      ! if the faces reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_cells_local
        string = trim(adjustl(string)) // ' cells stored on p0'
        print *, trim(string)
#endif
        unstructured_grid%cell_vertices(:,1:num_cells_local) = &
          temp_int_array(:,1:num_cells_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' cells sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*unstructured_grid%max_nvert_per_cell
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_cells_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' cells received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_cells_local*unstructured_grid%max_nvert_per_cell
    call MPI_Recv(unstructured_grid%cell_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif

  unstructured_grid%nlmax = num_cells_local

  call InputDestroy(input)

end subroutine UGridReadSurfGrid

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine UGridReadHDF5SurfGrid(unstructured_grid,filename,option)
  ! 
  ! This routine reads unstructured grid from HDF5 for surface mesh.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 06/01/12
  ! 

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use Input_Aux_module
  use Option_module
  use HDF5_Aux_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscMPIInt :: hdf5_err
  PetscMPIInt :: rank_mpi
  PetscInt :: istart, iend, ii, jj
  PetscInt :: num_cells_local
  PetscInt :: num_cells_local_save
  PetscInt :: num_vertices_local
  PetscInt :: num_vertices_local_save
  PetscInt :: remainder
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_buffer(:,:)
  PetscReal, pointer :: double_buffer(:,:)
  PetscInt, parameter :: max_nvert_per_cell = 8  
  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: num_data_in_file
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(2), length(2), stride(2), block(2), dims(2)
#endif
  integer :: ndims_h5

  ! Initialize FORTRAN predefined datatypes
  call h5open_f(hdf5_err)

  ! Setup file access property with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, hdf5_err)

#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm, MPI_INFO_NULL, hdf5_err)
#endif

  ! Open the file collectively
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id, hdf5_err)
  
  !
  ! Regions/top
  !
  
  ! Open group
  group_name = "/Regions/top/Vertex Ids"
  option%io_buffer = 'Opening group: ' // trim(group_name)
  call printMsg(option)

  ! Open dataset
  call h5dopen_f(file_id, "/Regions/top/Vertex Ids", data_set_id, hdf5_err)

  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims_h5, hdf5_err)
  if (ndims_h5 /= 2) then
    option%io_buffer='Dimension of Domain/Cells dataset in ' // trim(filename) // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%nmax = INT(dims_h5(2))
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                  num_cells_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_cells_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_cells_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
  ! Determine the length and offset of data to be read by each processor
  length(1) = INT(dims_h5(1))
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart
  
  !
  rank_mpi = 2
  memory_space_id = -1
  
  ! Create data space for dataset
  call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)
  
  ! Select hyperslab
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(int_buffer(length(1), length(2)))
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices(max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -1
  
  do ii = 1, num_cells_local
    do jj = 1, INT(dims_h5(1))
      if (int_buffer(jj, ii) > 0) then
        unstructured_grid%cell_vertices(jj, ii) = int_buffer(jj, ii)
      end if
    enddo
  enddo
  
  call h5dclose_f(data_set_id, hdf5_err)
  
  deallocate(int_buffer)
  nullify(int_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)

  !
  ! Domain/Vertices
  !
  
  ! Open dataset
  call h5dopen_f(file_id, "Domain/Vertices", data_set_id, hdf5_err)
  
  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims_h5, hdf5_err)
  if (ndims_h5 /= 2) then
    option%io_buffer='Dimension of Domain/Vertices dataset in ' // trim(filename) // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%num_vertices_global = INT(dims_h5(2))
  num_vertices_local  = &
                                       unstructured_grid%num_vertices_global/ &
                                       option%mycommsize 
  num_cells_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                  num_vertices_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_vertices_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_vertices_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
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
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(double_buffer(length(1), length(2)))
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, double_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  call h5dclose_f(data_set_id, hdf5_err)
  !call h5gclose_f(grp_id, hdf5_err)
  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)

  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ii = 1, num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo
  
  
  deallocate(double_buffer)
  nullify(double_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)
  
  unstructured_grid%max_nvert_per_cell = max_nvert_per_cell
  unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local
  
end subroutine UGridReadHDF5SurfGrid

#endif
! End PETSC_HAVE_HDF5

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine UGridReadHDF5(unstructured_grid,filename,option)
  ! 
  ! Reads an unstructured grid from HDF5
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/25/11
  ! 

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use Input_Aux_module
  use Option_module
  use HDF5_Aux_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: rank_mpi
  PetscInt :: istart, iend, ii, jj
  PetscInt :: num_cells_local
  PetscInt :: num_cells_local_save
  PetscInt :: num_vertices_local
  PetscInt :: num_vertices_local_save
  PetscInt :: remainder
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_buffer(:,:)
  PetscReal, pointer :: double_buffer(:,:)
  PetscInt, parameter :: max_nvert_per_cell = 8  
  PetscInt :: error_count
  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: num_data_in_file
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(2), length(2), stride(2), block(2), dims(2)
#endif
  integer :: ndims_h5

  ! Initialize FORTRAN predefined datatypes
  call h5open_f(hdf5_err)

  ! Setup file access property with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, hdf5_err)

#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm, MPI_INFO_NULL, hdf5_err)
#endif

  ! Open the file collectively
  call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
  call h5pclose_f(prop_id, hdf5_err)
  
  !
  ! Domain/Cells
  !
  
  ! Open group
  group_name = "Domain"
  option%io_buffer = 'Opening group: ' // trim(group_name)
  call printMsg(option)

  ! Open dataset
  call h5dopen_f(file_id, "Domain/Cells", data_set_id, hdf5_err)

  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims_h5, hdf5_err)
  if (ndims_h5 /= 2) then
    option%io_buffer='Dimension of Domain/Cells dataset in ' // trim(filename) // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%nmax = INT(dims_h5(2))
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                  num_cells_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_cells_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_cells_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
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
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(int_buffer(length(1), length(2)))
  int_buffer = UNINITIALIZED_INTEGER
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices(max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -1

  ! check for incorrectly assigned cell types
  error_count = 0
  do ii = 1, num_cells_local
    select case(int_buffer(1,ii))
      case(4,5,6,8)
      case default
        write(string,*) int_buffer(1,ii)
        option%io_buffer = 'Unknown cell type : ' // trim(adjustl(string))
        error_count = error_count + 1
        if (error_count < 10) then
          call printMsgByRank(option)
        endif
    end select
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,error_count,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (error_count > 0) then
    option%io_buffer = 'Unknown cell types in ' // trim(filename) // '.'
    call printErrMsg(option)
  endif
  
  do ii = 1, num_cells_local
    do jj = 2, int_buffer(1,ii) + 1
      unstructured_grid%cell_vertices(jj-1, ii) = int_buffer(jj, ii)
    enddo
  enddo
  
  call h5dclose_f(data_set_id, hdf5_err)
  
  deallocate(int_buffer)
  nullify(int_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)

  !
  ! Domain/Vertices
  !
  
  ! Open dataset
  call h5dopen_f(file_id, "Domain/Vertices", data_set_id, hdf5_err)
  
  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims_h5, hdf5_err)
  if (ndims_h5 /= 2) then
    option%io_buffer='Dimension of Domain/Vertices dataset in ' // trim(filename) // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%num_vertices_global = INT(dims_h5(2))
  num_vertices_local  = &
                                       unstructured_grid%num_vertices_global/ &
                                       option%mycommsize 
  num_cells_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                  num_vertices_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_vertices_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_vertices_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
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
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(double_buffer(length(1), length(2)))
  double_buffer = UNINITIALIZED_DOUBLE
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, double_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  call h5dclose_f(data_set_id, hdf5_err)
  !call h5gclose_f(grp_id, hdf5_err)
  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)

  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ii = 1, num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo
  
  deallocate(double_buffer)
  nullify(double_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)
  
  unstructured_grid%max_nvert_per_cell = max_nvert_per_cell
  unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local
  
end subroutine UGridReadHDF5

#endif
! End PETSC_HAVE_HDF5

#if defined(SCORPIO)

! ************************************************************************** !

subroutine UGridReadHDF5PIOLib(unstructured_grid, filename, &
                                          option)
!
! UGridReadHDF5PIOLib: Reads an unstructured grid from HDF5
! Author: Gautam Bisht
! Date: 05/13/11
!
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

!#include "petsc/finclude/petscsys.h"

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use Input_Aux_module
  use Option_module
  use HDF5_Aux_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscInt,pointer :: int_buffer(:,:)
  PetscReal,pointer :: double_buffer(:,:)
  PetscInt :: ii, jj
  PetscInt :: dims(2), dataset_dims(2)
  PetscInt, parameter :: max_nvert_per_cell = 8
  PetscInt :: num_cells_local

  character(len=MAXSTRINGLENGTH) :: cell_dataset_name = &
                                                       '/Domain/Cells'//CHAR(0)
  character(len=MAXSTRINGLENGTH) :: vert_dataset_name = &
                                                    '/Domain/Vertices'//CHAR(0)

  ! Read Domain/Cells
  call HDF5ReadDatasetInteger2D(filename, &
                                cell_dataset_name, &
                                SCORPIO_NONUNIFORM_CONTIGUOUS_READ, &
                                option, &
                                int_buffer, &
                                dims, &
                                dataset_dims)

  ! Allocate array to store vertices for each cell
  num_cells_local  = dims(2)
  unstructured_grid%nmax = dataset_dims(2)
  allocate(unstructured_grid%cell_vertices(max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -1

  ! Fill the cell data structure
  do ii = 1, num_cells_local
    do jj = 2, int_buffer(1, ii) + 1
      unstructured_grid%cell_vertices(jj-1, ii) = int_buffer(jj, ii)
    enddo
  enddo
  deallocate(int_buffer)
  nullify(int_buffer)

  ! Read Vertices
  call HDF5ReadDatasetReal2D(filename, &
                             vert_dataset_name, &
                             SCORPIO_NONUNIFORM_CONTIGUOUS_READ, &
                             option, &
                             double_buffer, &
                             dims, &
                             dataset_dims)

  unstructured_grid%num_vertices_local = dims(2)
  unstructured_grid%num_vertices_global= dataset_dims(2)
  allocate(unstructured_grid%vertices(unstructured_grid%num_vertices_local))
  ! fill the vertices data structure
  do ii = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo
  deallocate(double_buffer)
  nullify(double_buffer)

  unstructured_grid%max_nvert_per_cell = max_nvert_per_cell
  unstructured_grid%nlmax = num_cells_local

end subroutine UGridReadHDF5PIOLib

#endif

! ************************************************************************** !

subroutine UGridDecompose(unstructured_grid,option)
  ! 
  ! Decomposes an unstructured grid across ranks
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/09
  ! 
  
#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  use Utility_module, only: reallocateIntArray, SearchOrderedArray
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  
  PetscInt :: local_id, local_id2
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  PetscInt :: count, vertex_count
  PetscInt :: vertex_offset, global_vertex_offset
  PetscInt :: stride
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscInt :: num_rows, num_cols, istart, iend, icol
  PetscBool :: success
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  PetscViewer :: viewer
  Mat :: Adj_mat
  Mat :: Dual_mat
  MatPartitioning :: Part
  Vec :: elements_natural
  Vec :: elements_local
  Vec :: elements_old
  Vec :: vertices_old
  Vec :: vertices_new
  IS :: is_new
  IS :: is_scatter
  IS :: is_gather

  VecScatter :: vec_scatter
  
  PetscInt :: vertex_ids_offset
  PetscInt :: dual_offset
  PetscInt :: natural_id_offset

  PetscInt :: max_int_count
  PetscInt :: temp_int
  PetscInt :: min_value
  PetscInt :: num_cells_local_new
  PetscInt :: num_cells_local_old  
  PetscInt :: global_offset_old
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, pointer :: int_array_pointer(:)
  
  PetscInt :: idual, dual_id
  PetscInt :: iflag
  PetscBool :: found

!  cell distribution across processors (size = num_cores + 1)
!  core i owns cells cell_distribution(i):cell_distribution(i+1), note
!  the zero-based indexing
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%nlmax,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)

  num_cells_local_old = unstructured_grid%nlmax

  ! recalculate maximum number of vertices for any given cell
  temp_int = 0
  min_value = 2 ! min value should be either 0 or 1 after global reduction
  do local_id = 1, num_cells_local_old
    vertex_count = 0
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point, cell vertex can be 0
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) exit
      if (unstructured_grid%cell_vertices(ivertex,local_id) < min_value) then
        min_value = unstructured_grid%cell_vertices(ivertex,local_id)
      endif
      vertex_count = vertex_count+1
    enddo
    if (vertex_count > temp_int) temp_int = vertex_count
  enddo
  call MPI_Allreduce(temp_int,unstructured_grid%max_nvert_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(min_value,index_format_flag, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN,option%mycomm,ierr)

  ! let's make it Fortran indexing
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point we may be zero-based
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) then
        ! change no_value (UNINITIALIZED_INTEGER) to '0'
        unstructured_grid%cell_vertices(ivertex,local_id) = 0
      else
        if (index_format_flag == 0) then
          ! let's make it Fortran indexing
          unstructured_grid%cell_vertices(ivertex,local_id) = &
            unstructured_grid%cell_vertices(ivertex,local_id) + 1
        endif
      endif
    enddo
  enddo

#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_nvert_per_cell
  option%io_buffer = 'Maximum number of vertices per cell: ' // adjustl(string)
  call printMsg(option)
  write(string,*) index_format_flag
  option%io_buffer = 'Vertex indexing starts at: ' // adjustl(string)
  call printMsg(option)
  if (index_format_flag == 0) then
    option%io_buffer = 'Changing vertex indexing to 1-based.'
    call printMsg(option)
  endif
#endif

  num_cells_local_old = unstructured_grid%nlmax 
  allocate(local_vertices(unstructured_grid%max_nvert_per_cell* &
                          num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  local_vertices = 0
  local_vertex_offset = 0
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      if (unstructured_grid%cell_vertices(ivertex,local_id) == 0) exit
      count = count + 1
      ! local vertices must be zero-based for MatCreateMPIAdj; thus subtract 1
      local_vertices(count) = &
        unstructured_grid%cell_vertices(ivertex,local_id) - 1
    enddo
    local_vertex_offset(local_id+1) = count 
  enddo
    
  select case (unstructured_grid%grid_type)
    case(TWO_DIM_GRID)
      num_common_vertices = 2 ! cells must share at least this number of vertices
    case(THREE_DIM_GRID)
      num_common_vertices = 3 ! cells must share at least this number of vertices
    case default
        option%io_buffer = 'Grid type not recognized '
        call printErrMsg(option)
    end select

  ! determine the global offset from 0 for cells on this rank
  global_offset_old = 0
  call MPI_Exscan(num_cells_local_old,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! create an adjacency matrix for calculating the duals (connnections)
#if UGRID_DEBUG
  call printMsg(option,'Adjacency matrix')
#endif

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat, &
                       ierr);CHKERRQ(ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'Adj_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'Adj_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call MatView(Adj_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

#if UGRID_DEBUG
  call printMsg(option,'Dual matrix')
#endif

#if defined(PETSC_HAVE_PARMETIS)
  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat, &
                          ierr);CHKERRQ(ierr)
#endif
  call MatDestroy(Adj_mat,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'Dual_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'Dual_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call MatView(Dual_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call UGridPartition(unstructured_grid,option,Dual_mat,is_new, &
                      num_cells_local_new)
  
  if (allocated(local_vertices)) deallocate(local_vertices)
  if (allocated(local_vertex_offset)) deallocate(local_vertex_offset)
  
  ! second argument of ZERO_INTEGER means to use 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  if (.not.success .or. num_rows /= num_cells_local_old) then
    print *, option%myrank, num_rows, success, num_cells_local_old
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  ! calculate maximum number of connections for any given cell
  unstructured_grid%max_ndual_per_cell = 0
  do local_id = 1, num_cells_local_old
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) &
      unstructured_grid%max_ndual_per_cell = num_cols
  enddo
  temp_int = unstructured_grid%max_ndual_per_cell
  call MPI_Allreduce(temp_int,unstructured_grid%max_ndual_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  
#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_ndual_per_cell
  option%io_buffer = 'Maximum number of duals per cell: ' // adjustl(string)
  call printMsg(option)
#endif
  
  if (unstructured_grid%max_ndual_per_cell > 0) then
    call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                            num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  
  ! in order to redistributed vertex/cell data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  vertex_ids_offset = 1 + 1 ! +1 for -777
  dual_offset = vertex_ids_offset + unstructured_grid%max_nvert_per_cell + 1 ! +1 for -888
  stride = dual_offset+ unstructured_grid%max_ndual_per_cell + 1 ! +1 for -999999
  natural_id_offset = 1

  ! Information for each cell is packed in a strided petsc vec
  ! The information is ordered within each stride as follows:
  ! -cell_N   ! global cell id (negative indicates 1-based)
  ! -777      ! separator between cell id and vertex ids for cell_N
  ! vertex1   ! in cell_N
  ! vertex2
  ! ...
  ! vertexN   
  ! -888      ! separator between vertex and dual ids
  ! dual1     ! dual ids between cell_N and others
  ! dual2
  ! ...
  ! dualN     
  ! -999999   ! separator indicating end of information for cell_N
  
  ! the purpose of -777, -888, and -999999 is to allow one to use cells of 
  ! various geometry.  Currently, the max # vertices = 8 and max # duals = 6.
  ! But this will be generalized in the future.
  
  call UGridCreateOldVec(unstructured_grid,option,elements_old, &
                                num_cells_local_old, &
                                is_new,is_scatter,stride)

  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(elements_old,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  vertex_count = 0
  do local_id = 1, num_cells_local_old
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(global_offset_old+local_id)
    count = count + 1
    ! add the separator
    vec_ptr(count) = -777  ! help differentiate
    ! add the vertex ids
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      count = count + 1
      vertex_count = vertex_count + 1
      ! increment for 1-based ordering
      vec_ptr(count) = unstructured_grid%cell_vertices(ivertex,local_id)
    enddo


    count = count + 1 
    ! another vertex/dual separator
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Dual matrix is larger then max_ndual_per_cell.'
      call printErrMsgByRank(option)
    endif
    do icol = 1, unstructured_grid%max_ndual_per_cell
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo
    count = count + 1 
    ! final separator
    vec_ptr(count) = -999999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr);CHKERRQ(ierr)
  
  if (unstructured_grid%max_ndual_per_cell > 0) then
    call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                            num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  call MatDestroy(Dual_mat,ierr);CHKERRQ(ierr)
 
  
  call UGridNaturalToPetsc(unstructured_grid,option, &
                           elements_old,elements_local, &
                           num_cells_local_new,stride,dual_offset, &
                           natural_id_offset,is_scatter)
  
  ! make a list of local vertices
  max_int_count = 2*unstructured_grid%ngmax
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
  ! note that the vertices are still in natural numbering
  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id=1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride))
      if (vertex_id < 1) exit
      vertex_count = vertex_count + 1
      if (vertex_count > max_int_count) then
        call reallocateIntArray(int_array_pointer,max_int_count)
      endif
      vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride) = vertex_count
      int_array_pointer(vertex_count) = vertex_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)

  ! sort the vertex ids
  allocate(int_array(vertex_count))
  int_array(1:vertex_count) = int_array_pointer(1:vertex_count)
  allocate(int_array2(vertex_count))
  do ivertex = 1, vertex_count
    int_array2(ivertex) = ivertex 
  enddo
  deallocate(int_array_pointer)
  nullify(int_array_pointer)
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(vertex_count,int_array,int_array2, &
                                   ierr);CHKERRQ(ierr)
  int_array2 = int_array2+1

  ! remove duplicates
  allocate(int_array3(vertex_count))
  allocate(int_array4(vertex_count))
  int_array3 = 0
  int_array4 = 0
  int_array3(1) = int_array(int_array2(1))
  count = 1
  int_array4(int_array2(1)) = count
  do ivertex = 2, vertex_count
    vertex_id = int_array(int_array2(ivertex))
    if (vertex_id > int_array3(count)) then
      count = count + 1
      int_array3(count) = vertex_id
    endif
    int_array4(int_array2(ivertex)) = count
  enddo
  vertex_count = count
  deallocate(int_array)

  allocate(unstructured_grid%vertex_ids_natural(vertex_count))
  unstructured_grid%vertex_ids_natural = int_array3(1:vertex_count)

  ! now load all the vertices needed to define all the local cells
  ! on the processor
  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  ! allocate the array that will store the vertex ids for each cell.
  ! remember that max_nvert_per_cell is the max # of vertices in a cell
  ! currently hardwired to 8.
  deallocate(unstructured_grid%cell_vertices)
  allocate(unstructured_grid%cell_vertices( &
             0:unstructured_grid%max_nvert_per_cell,unstructured_grid%ngmax))
  unstructured_grid%cell_vertices = 0
  
  ! permute the local ids calculated earlier in the int_array4
  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! extract the original vertex id
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride))
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices(0,ghosted_id)+1
      unstructured_grid%cell_vertices(count,ghosted_id) = &
        int_array4(vertex_id)
      unstructured_grid%cell_vertices(0,ghosted_id) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride) = &
        int_array4(vertex_id)
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(int_array2)
  deallocate(int_array3)
  deallocate(int_array4)

#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'elements_vert_local' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'elements_vert_local' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(elements_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  call VecDestroy(elements_local,ierr);CHKERRQ(ierr)

  ! now we need to work on aligning the original vertex coordinates with 
  ! the current ordering or permuted/rearranged ordering.

  ! IS for gather operation - need local numbering
  allocate(int_array(vertex_count))
  ! vertex_count = # of local vertices (I believe ghosted+non-ghosted)
  do ivertex = 1, vertex_count
    int_array(ivertex) = ivertex-1
  enddo

  ! include cell ids (use block ids, not indices)
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  ! create a parallel petsc vector with a stride of 3.
  !call VecCreateMPI(option%mycomm,unstructured_grid%num_vertices_local*3, &
  !                  PETSC_DETERMINE,vertices_old,ierr)
  call VecCreate(option%mycomm,vertices_old,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_old,unstructured_grid%num_vertices_local*3, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_old,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_old,ierr);CHKERRQ(ierr)

  ! create serial petsc vector with a stride of 3
  !call VecCreateSeq(PETSC_COMM_SELF,vertex_count*3,vertices_new,ierr)
  call VecCreate(PETSC_COMM_SELF,vertices_new,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_new,vertex_count*3,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_new,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_new,ierr);CHKERRQ(ierr)

!  call VecCreate(option%mycomm,vertices_new,ierr)
!  call VecSetSizes(vertices_new, &
!                   vertex_count*3,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_new,ierr)
  
!  call VecCreate(option%mycomm,vertices_old,ierr)
!  call VecSetSizes(vertices_old, &
!                   3*unstructured_grid%num_vertices_local,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_old,ierr)
! load up the coordinates
  call VecGetArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    vec_ptr((ivertex-1)*3+1) = unstructured_grid%vertices(ivertex)%x
    vec_ptr((ivertex-1)*3+2) = unstructured_grid%vertices(ivertex)%y
    vec_ptr((ivertex-1)*3+3) = unstructured_grid%vertices(ivertex)%z
  enddo
  call VecRestoreArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)

  ! IS for scatter - provide petsc global numbering
  allocate(int_array(vertex_count))
  do ivertex = 1, vertex_count
    int_array(ivertex) = (needed_vertices_petsc(ivertex)-1)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_scatter, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

  ! resize vertex array to new size
  unstructured_grid%num_vertices_natural = unstructured_grid%num_vertices_local
  unstructured_grid%num_vertices_local = vertex_count
  allocate(unstructured_grid%vertices(vertex_count))
  do ivertex = 1, vertex_count
    unstructured_grid%vertices(ivertex)%x = 0.d0
    unstructured_grid%vertices(ivertex)%y = 0.d0
    unstructured_grid%vertices(ivertex)%z = 0.d0
  enddo

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_vert_old_to_new_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_vert_old_to_new_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_gather,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecScatterCreate(vertices_old,is_scatter,vertices_new,is_gather, &
                        vec_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_gather,ierr);CHKERRQ(ierr)
  call VecScatterBegin(vec_scatter,vertices_old,vertices_new, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,vertices_old,vertices_new, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call VecView(vertices_old,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecDestroy(vertices_old,ierr);CHKERRQ(ierr)


  call VecGetArrayF90(vertices_new,vec_ptr,ierr);CHKERRQ(ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ivertex)%id = needed_vertices_petsc(ivertex)
    unstructured_grid%vertices(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    unstructured_grid%vertices(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    unstructured_grid%vertices(ivertex)%z = vec_ptr((ivertex-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_new,vec_ptr,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'vertex_coord_new' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'vertex_coord_new' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(vertices_new,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecDestroy(vertices_new,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call printMsg(option,'Setting cell types')
#endif

  allocate(unstructured_grid%cell_type(unstructured_grid%ngmax))

  select case(unstructured_grid%grid_type)
    case(THREE_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        ! Determine number of faces and cell-type of the current cell
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(8)
            unstructured_grid%cell_type(ghosted_id) = HEX_TYPE
          case(6)
            unstructured_grid%cell_type(ghosted_id) = WEDGE_TYPE
          case(5)
            unstructured_grid%cell_type(ghosted_id) = PYR_TYPE
          case(4)
            unstructured_grid%cell_type(ghosted_id) = TET_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call printErrMsg(option)
        end select      
      enddo
    case(TWO_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(4)
            unstructured_grid%cell_type = QUAD_TYPE
          case(3)
            unstructured_grid%cell_type = TRI_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call printErrMsg(option)
        end select
      end do
    case default
      option%io_buffer = 'Grid type not recognized: '
      call printErrMsg(option)
  end select
  
end subroutine UGridDecompose

! ************************************************************************** !

function UGridComputeInternConnect(unstructured_grid,grid_x,grid_y,grid_z, &
                                   option)
  ! 
  ! computes internal connectivity of an
  ! unstructured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/09
  ! 

  use Connection_module
  use Option_module
  use Utility_module, only : DotProduct, CrossProduct
  use Geometry_module  

  implicit none

  type(connection_set_type), pointer :: UGridComputeInternConnect
  type(option_type) :: option
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(grid_unstructured_type) :: unstructured_grid

  type(connection_set_type), pointer :: connections
  PetscInt :: nconn, iconn
  PetscInt :: idual, dual_id

  PetscInt, allocatable :: face_to_vertex(:,:)
  PetscInt, allocatable :: cell_to_face(:,:)
  PetscInt, allocatable :: face_to_cell(:,:)
  PetscInt, allocatable :: vertex_to_cell(:,:)
  PetscInt, allocatable :: temp_int(:)
  PetscInt, allocatable :: temp_int_2d(:,:)
  PetscBool, allocatable :: local_boundary_face(:)
  PetscInt :: num_match
  PetscInt :: found_count
  PetscBool :: found
  PetscBool :: match_found
  PetscInt :: face_count
  PetscInt :: count, i
  PetscInt :: iface, iface2, iside
  PetscInt :: face_id, face_id2
  PetscInt :: ghosted_id, ghosted_id2
  PetscInt :: local_id, local_id2
  PetscInt :: cell_id, cell_id2
  PetscInt :: dual_local_id
  PetscInt :: ivertex, ivertex2
  PetscInt :: vertex_id, vertex_id2
  PetscInt :: vertex_ids4(4)
  PetscInt :: nfaces, nfaces2, nvertices, nvertices2, cell_type, cell_type2
  PetscInt :: face_type, face_type2
  PetscBool :: face_found, vertex_found
  
  PetscReal :: v1(3), v2(3), v3(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: vcross(3), magnitude
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  PetscInt :: ivert
  
  type(plane_type) :: plane1, plane2
  type(point3d_type) :: point1, point2, point3, point4
  type(point3d_type) :: point_up, point_dn
  type(point3d_type) :: intercept1, intercept2, intercept

  character(len=MAXSTRINGLENGTH) :: string  
  
  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL* &
           unstructured_grid%ngmax))
  face_to_vertex = 0
  allocate(cell_to_face(MAX_FACE_PER_CELL, &
                        unstructured_grid%ngmax))
  cell_to_face = 0
  allocate(face_to_cell(2,MAX_FACE_PER_CELL* &
                        unstructured_grid%ngmax))
  face_to_cell = 0
  allocate(vertex_to_cell(0:unstructured_grid%max_cells_sharing_a_vertex, &
                          unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  allocate(unstructured_grid%face_to_vertex_natural(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL*unstructured_grid%ngmax))
  unstructured_grid%face_to_vertex_natural = 0

  face_count = 0
  do ghosted_id = 1, unstructured_grid%ngmax
    cell_type = unstructured_grid%cell_type(ghosted_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do iface = 1, nfaces
      face_count = face_count + 1
      cell_to_face(iface,ghosted_id) = face_count
      face_to_cell(1,face_count) = ghosted_id
      call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                      vertex_ids4)
      do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices(vertex_ids4(ivertex),ghosted_id)
          if (face_to_vertex(ivertex,face_count) > 0) then
            unstructured_grid%face_to_vertex_natural(ivertex,face_count) = &
              unstructured_grid%vertex_ids_natural(face_to_vertex(ivertex,face_count))
          endif
      enddo
    enddo
  enddo

  !
  ! Remove duplicate faces:
  !
  ! A cell (cell_id) and Neighboring-Cell (cell_id2) will only share ONE face.
  ! Find the face that cell_id ane cell_id2 share and remove it.
  !
  ! Method:
  !        - Pick i-th face (iface) of cell_id and check if ALL the vertices of
  !          the iface are present in cell_id2. If all the vertices of iface are
  !          not present in cell_id2, move to the next face.
  !        - After finding the iface, now find iface2 in cell_id2 that
  !          corresponds to iface.
  !        - Check to ensure that atleast on face of cell_id is shared
  !          with cell_id2.
  !
  !
  !
  ! NOTE: For a cell_type = WEDGE_TYPE, faces 1-3 have 4 vertices; while
  !       faces 4-5 have 3 vertices
  !
  do local_id = 1, unstructured_grid%nlmax
    ! Selet a cell and find number of vertices
    cell_id = local_id
    ! cell_type is ghosted, but local cells are in the first nlmax entries
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      ! Select a neighboring cell
      ! ghosted neighbors have a negative id
      cell_id2 = &
        abs(unstructured_grid%cell_neighbors_local_ghosted(idual,local_id))
      cell_type2 = unstructured_grid%cell_type(cell_id2)
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      ! Find the number of vertices for neighboring cell
      nfaces2 = UCellGetNFaces(cell_type2,option)
      ! Initialize
      face_found = PETSC_FALSE
      do iface = 1, nfaces
        ! Select a face and find number of vertices forming the face
        face_id = cell_to_face(iface,cell_id)
        nvertices = UCellGetNFaceVertices(cell_type,iface,option)
        do ivertex = 1, nvertices
          ! Select a vertex and initialize vertex_found
          vertex_id = face_to_vertex(ivertex,face_id) ! face_to_vertex is 1-based indexing
          vertex_found = PETSC_FALSE
          do ivertex2 = 1, unstructured_grid%cell_vertices(0,cell_id2)
            vertex_id2 = unstructured_grid%cell_vertices(ivertex2,cell_id2)
            if (vertex_id == vertex_id2) then
              vertex_found = PETSC_TRUE
              exit
            endif
          enddo
          !
          ! If ivertex of iface of the Cell is not present as vertices of the
          ! Neighboring-Cell, then iface is not the shared face. Skip iterating
          ! over the remaing vertices of iface
          if (.not.vertex_found) exit
        enddo
        
        if (vertex_found) then
          ! All the vertices of iface are present in the Neighboring cells.
          ! Thus, iface is the shared face.
          face_found = PETSC_TRUE
          
          ! Now, we have to find iface2 that corresponds to iface
          do iface2 = 1, nfaces2
            face_id2 = cell_to_face(iface2,cell_id2)
            !geh nvertices2 = 4
            !gehcomment: I believe that cell_type and iface on next line shoudl be the "2" versions
            !geh if ((cell_type == WEDGE_TYPE).and.(iface.gt.3)) nvertices2 = 3
            nvertices2 = UCellGetNFaceVertices(cell_type2,iface2,option)
            ! Both iface and iface2 need to have same number of vertices
            if (nvertices == nvertices2) then
              ! Count the number of vertices of iface which match vertices
              ! of iface2
              num_match = 0
              do ivertex = 1,nvertices
                vertex_id = face_to_vertex(ivertex,face_id)
                vertex_found = PETSC_FALSE ! gehbug - used to be PETSC_TRUE
                
                do ivertex2 = 1, nvertices2 ! gehbug - used to be nvertices
                  vertex_id2 = face_to_vertex(ivertex2,face_id2)
                  if (vertex_id == vertex_id2) then
                    vertex_found = PETSC_TRUE 
                    num_match = num_match + 1
                    vertex_ids4(num_match) = vertex_id
                    exit
                  endif
                enddo
                !
                ! If vertex_id of face_id not found as one of vertices of face_id2,
                ! face_id2 is not shared between cells
                if (.not.vertex_found) exit
              enddo
              if (num_match == nvertices) then
                ! remove duplicate face
                !geh: I believe that face_id2 will always be removed
                if (face_id2 > face_id) then
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id2, ' -> ', face_id
                  option%io_buffer = 'Duplicated face removed:' // trim(string)
                  call printMsg(option)
#endif
                  cell_to_face(iface2,cell_id2) = face_id
                  ! flag face_id2 as removed
                  face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
                  ! add cell_id2 to face_ids list
                  face_to_cell(2,face_id) = cell_id2
                else
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id, ' -> ', face_id2
                  option%io_buffer = 'Duplicated face removed:' // trim(string)
                  call printMsg(option)
#endif
                  cell_to_face(iface,cell_id) = face_id2
                  ! flag face_id as removed  
                  face_to_cell(1,face_id) = -face_to_cell(1,face_id)
                  ! add cell_id to face_ids2 list
                  face_to_cell(2,face_id2) = cell_id
                endif
                exit
              endif
            endif
          enddo
          exit
        endif
      enddo ! iface-loop
      
      ! Check that one shared face was found between the Cell and Neighboring-Cell
      if (.not.face_found) then
        write(string,*) option%myrank
        string = '(' // trim(adjustl(string)) // ')'
        write(*,'(a,'' local_id = '',i6,'' natural_id = '',i6,&
                  &''  vertices: '',8i6)') &
                   trim(string), &
                   cell_id,unstructured_grid%cell_ids_natural(cell_id), &
                   (unstructured_grid%vertex_ids_natural( &
                     unstructured_grid%cell_vertices(ivertex,cell_id)), &
                     ivertex=1,unstructured_grid%cell_vertices(0,cell_id))
        write(*,'(a,'' local_id2 = '',i6,'' natural_id2 = '',i6,&
                  &''  vertices2: '',8i6)') &
                   trim(string), &
                   cell_id2,unstructured_grid%cell_ids_natural(cell_id2), &
                   (unstructured_grid%vertex_ids_natural( &
                     unstructured_grid%cell_vertices(ivertex2,cell_id2)), &
                     ivertex2=1,unstructured_grid%cell_vertices(0,cell_id2))
        option%io_buffer='No shared face found.'
        call printErrMsgByRank(option)
      endif
    enddo ! idual-loop
  enddo  ! local_id-loop

  ! count up the # of faces
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) &
      face_count = face_count + 1
  enddo
  allocate(unstructured_grid%face_to_vertex(MAX_VERT_PER_FACE,face_count))
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      unstructured_grid%face_to_vertex(:,face_count) = face_to_vertex(:,iface)
    endif
  enddo
  deallocate(face_to_vertex)
  ! reallocate face_to_cell to proper size
  allocate(temp_int_2d(2,face_count))
  allocate(temp_int(size(face_to_cell,2)))
  temp_int = 0
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      temp_int_2d(:,face_count) = face_to_cell(:,iface)
      temp_int(iface) = face_count
    endif
  enddo
  deallocate(face_to_cell)
  allocate(face_to_cell(2,face_count))
  face_to_cell = temp_int_2d
  deallocate(temp_int_2d)

  
  ! remap faces in cells using temp_int from above
  do iface = 1, size(face_to_cell,2)
    face_id = iface
    do i = 1,2
      cell_id = face_to_cell(i,face_id)
      ! check for exterior face
      if (cell_id < 1) cycle
      found = PETSC_FALSE
      cell_type = unstructured_grid%cell_type(cell_id)
      nfaces = UCellGetNFaces(cell_type,option)
      do iface2 = 1, nfaces
        face_id2 = cell_to_face(iface2,cell_id)
        if (face_id < 0) cycle
        if (face_id == temp_int(face_id2)) then
          found = PETSC_TRUE
          cell_to_face(iface2,cell_id) = face_id
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Remapping of cell face id unsuccessful'
        call printErrMsg(option)
      endif
    enddo
  enddo
  deallocate(temp_int)
  
  do ghosted_id = 1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      if ( vertex_id <= 0) cycle 
      count = vertex_to_cell(0,vertex_id) + 1
      if (count > unstructured_grid%max_cells_sharing_a_vertex) then
        write(string,*) 'Vertex can be shared by at most by ', &
              unstructured_grid%max_cells_sharing_a_vertex, &
              ' cells. Rank = ', option%myrank, ' vertex_id = ', vertex_id, ' exceeds it.'
        option%io_buffer = string
        call printErrMsg(option)
      endif
      vertex_to_cell(count,vertex_id) = ghosted_id
      vertex_to_cell(0,vertex_id) = count
    enddo
  enddo
  
  nconn = 0
  do local_id = 1, unstructured_grid%nlmax
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_id = unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      ! count all ghosted connections (dual_id < 0)
      ! only count connection with cells of larger ids to avoid double counts
!geh: we need to cound all local connection, but just once (local_id < dual_id) and all
!      ghosted connections (dual_id < 0)
      if (dual_id < 0 .or. local_id < dual_id) then
!geh: Nope      if (dual_id > 0 .and. local_id < dual_id) then !sp 
        nconn = nconn + 1
      endif
    enddo
  enddo


  connections => ConnectionCreate(nconn,INTERNAL_CONNECTION_TYPE)
  
  allocate(unstructured_grid%face_area(face_count))
  allocate(unstructured_grid%connection_to_face(nconn))
  unstructured_grid%connection_to_face = 0

  ! loop over connection again
  iconn = 0
  do local_id = 1, unstructured_grid%nlmax
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_local_id = &
        unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      ! abs(dual_local_id) to accommodate connections to ghost cells where 
      ! the dual id is < 0.
      if (local_id < abs(dual_local_id)) then 
        iconn = iconn + 1
        ! find face
        found = PETSC_FALSE
        do iface = 1, unstructured_grid%cell_vertices(0,local_id)
          face_id = cell_to_face(iface,local_id)
          do iside = 1,2
            cell_id2 = face_to_cell(iside,face_id)
            if (cell_id2 == abs(dual_local_id)) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (found) exit
        enddo
        if (found) then
          unstructured_grid%connection_to_face(iconn) = face_id
        else
          write(string,*) option%myrank,local_id,dual_local_id 
          option%io_buffer = 'face not found in connection loop' // trim(string)
          call printErrMsg(option)
        endif
        face_type = &
          UCellGetFaceType(unstructured_grid%cell_type(local_id),iface,option)
        found = PETSC_FALSE
        do iface2 = 1, unstructured_grid%cell_vertices(0,cell_id2)
          if (cell_to_face(iface,local_id) == &
              cell_to_face(iface2,cell_id2)) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) then
          face_type2 = &
            UCellGetFaceType(unstructured_grid%cell_type(cell_id2), &
                                                                 iface2,option)
          if (face_type /= face_type2) then
            write(string,*) option%myrank, local_id, cell_id2 
            option%io_buffer = 'face types do not match' // trim(string)
            call printErrMsg(option)
          endif
        else
          write(string,*) option%myrank, iface, cell_id2
          option%io_buffer = 'global face not found' // trim(string)
          call printErrMsg(option)
        endif
        connections%id_up(iconn) = local_id
        connections%id_dn(iconn) = abs(dual_local_id)
        connections%face_id(iconn) = cell_to_face(iface,local_id)
        if (face_type == LINE_FACE_TYPE) then

          point_up%x = grid_x(local_id)
          point_up%y = grid_y(local_id)
          point_up%z = grid_z(local_id)
          point_dn%x = grid_x(abs(dual_local_id))
          point_dn%y = grid_y(abs(dual_local_id))
          point_dn%z = grid_z(abs(dual_local_id))
          point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
          point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))

          call UcellGetLineIntercept(point1,point2,point_up,intercept1)
          call UcellGetLineIntercept(point1,point2,point_dn,intercept2)
          intercept%x = 0.5d0*(intercept1%x + intercept2%x)
          intercept%y = 0.5d0*(intercept1%y + intercept2%y)
          intercept%z = 0.5d0*(intercept1%z + intercept2%z)

          !v1(1) = point_dn%x-point_up%x
          !v1(2) = point_dn%y-point_up%y
          !v1(3) = point_dn%z-point_up%z
          v1(1) = point1%x-point2%x
          v1(2) = point1%y-point2%y
          v1(3) = point1%z-point2%z

          area1 = sqrt(DotProduct(v1,v1))
          area2 = 0.d0
        else
        
          ! need to add the surface areas, distance, etc.
          point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
          point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
          point3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
          if (face_type == QUAD_FACE_TYPE) then
            point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(4,face_id))
          endif
          
          call GeometryComputePlaneWithPoints(plane1,point1,point2,point3)
         
          point_up%x = grid_x(local_id)
          point_up%y = grid_y(local_id)
          point_up%z = grid_z(local_id)
          point_dn%x = grid_x(abs(dual_local_id))
          point_dn%y = grid_y(abs(dual_local_id))
          point_dn%z = grid_z(abs(dual_local_id))
          v1(1) = point_dn%x-point_up%x
          v1(2) = point_dn%y-point_up%y
          v1(3) = point_dn%z-point_up%z
          n_up_dn = v1 / sqrt(DotProduct(v1,v1))
          call GeometryGetPlaneIntercept(plane1,point_up,point_dn,intercept1)
          
          v1(1) = point3%x-point2%x
          v1(2) = point3%y-point2%y
          v1(3) = point3%z-point2%z
          v2(1) = point1%x-point2%x
          v2(2) = point1%y-point2%y
          v2(3) = point1%z-point2%z
          !geh: area = 0.5 * |v1 x v2|
          vcross = CrossProduct(v1,v2)
          !geh: but then we have to project the area onto the vector between
          !     the cell centers (n_up_dn)
          magnitude = sqrt(DotProduct(vcross,vcross))
          n1 = vcross/magnitude
          area1 = 0.5d0*magnitude
          area1 = dabs(area1*DotProduct(n1,n_up_dn))
          !geh: The below does not project onto the vector between cell centers.
          !gehbug area1 = 0.5d0*sqrt(DotProduct(n1,n1))
          
          if (face_type == QUAD_FACE_TYPE) then
            call GeometryComputePlaneWithPoints(plane2,point3,point4,point1)
            call GeometryGetPlaneIntercept(plane2,point_up,point_dn,intercept2)
            v1(1) = point1%x-point4%x
            v1(2) = point1%y-point4%y
            v1(3) = point1%z-point4%z
            v2(1) = point3%x-point4%x
            v2(2) = point3%y-point4%y
            v2(3) = point3%z-point4%z
            magnitude = sqrt(DotProduct(vcross,vcross))
            n2 = vcross/magnitude
            area2 = 0.5d0*magnitude
            area2 = dabs(area2*DotProduct(n2,n_up_dn))
          else 
            area2 = 0.d0
          endif

          if (face_type == QUAD_FACE_TYPE) then
            intercept%x = 0.5d0*(intercept1%x + intercept2%x)
            intercept%y = 0.5d0*(intercept1%y + intercept2%y)
            intercept%z = 0.5d0*(intercept1%z + intercept2%z)
          else
            intercept%x = intercept1%x
            intercept%y = intercept1%y
            intercept%z = intercept1%z
          endif
        endif
        
        !geh: this is very crude, but for now use average location of intercept
        v1(1) = intercept%x-point_up%x
        v1(2) = intercept%y-point_up%y
        v1(3) = intercept%z-point_up%z
        v2(1) = point_dn%x-intercept%x
        v2(2) = point_dn%y-intercept%y
        v2(3) = point_dn%z-intercept%z
        dist_up = sqrt(DotProduct(v1,v1))
        dist_dn = sqrt(DotProduct(v2,v2))
        
        connections%dist(-1:3,iconn) = 0.d0
        connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connections%dist(0,iconn) = dist_up + dist_dn
        v3 = v1 + v2
        connections%dist(1:3,iconn) = v3/sqrt(DotProduct(v3,v3))
        connections%area(iconn) = area1 + area2
        connections%intercp(1,iconn) = intercept%x
        connections%intercp(2,iconn) = intercept%y
        connections%intercp(3,iconn) = intercept%z
       
      endif
    enddo
  enddo
  
  ! Save area and centroid of faces
  allocate(unstructured_grid%face_centroid(face_count))
  do iface = 1,face_count
    unstructured_grid%face_centroid(iface)%id = 0
  enddo
  
  do local_id = 1, unstructured_grid%nlmax
    do iface = 1,MAX_FACE_PER_CELL
      face_id = cell_to_face(iface, local_id)
      if (face_id == 0) cycle
      if ( unstructured_grid%face_centroid(face_id)%id == 0) then
        count = 0
        unstructured_grid%face_centroid(face_id)%x = 0.d0
        unstructured_grid%face_centroid(face_id)%y = 0.d0
        unstructured_grid%face_centroid(face_id)%z = 0.d0

        if (unstructured_grid%face_to_vertex(3,face_id) == 0) then
          face_type = LINE_FACE_TYPE
        else
          if (unstructured_grid%face_to_vertex(4,face_id) == 0) then
            face_type = TRI_FACE_TYPE
          else
            face_type = QUAD_FACE_TYPE
          endif
        endif

        if (face_type == LINE_FACE_TYPE) then

          point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
          point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))

          v1(1) = point1%x-point2%x
          v1(2) = point1%y-point2%y
          v1(3) = point1%z-point2%z

          area1 = sqrt(DotProduct(v1,v1))
          area2 = 0.d0
        else

          point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
          point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
          point3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
          if (face_type == QUAD_FACE_TYPE) then
            point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(4,face_id))
          else
            point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
          endif
          v1(1) = point3%x-point2%x
          v1(2) = point3%y-point2%y
          v1(3) = point3%z-point2%z
          v2(1) = point1%x-point2%x
          v2(2) = point1%y-point2%y
          v2(3) = point1%z-point2%z
          n1 = CrossProduct(v1,v2)
          area1 = 0.5d0*sqrt(DotProduct(n1,n1))
        
          v1(1) = point1%x-point4%x
          v1(2) = point1%y-point4%y
          v1(3) = point1%z-point4%z
          v2(1) = point3%x-point4%x
          v2(2) = point3%y-point4%y
          v2(3) = point3%z-point4%z
          n2 = CrossProduct(v1,v2)
          area2 = 0.5d0*sqrt(DotProduct(n2,n2))
        endif
        
        unstructured_grid%face_area(face_id) = area1 + area2
        
        do ivert = 1,MAX_VERT_PER_FACE
          vertex_id = unstructured_grid%face_to_vertex(ivert,face_id)
          if (vertex_id.ne.0) then
            unstructured_grid%face_centroid(face_id)%x = &
              unstructured_grid%face_centroid(face_id)%x + &
              unstructured_grid%vertices(vertex_id)%x
            unstructured_grid%face_centroid(face_id)%y = &
              unstructured_grid%face_centroid(face_id)%y + &
              unstructured_grid%vertices(vertex_id)%y
            unstructured_grid%face_centroid(face_id)%z = &
              unstructured_grid%face_centroid(face_id)%z + &
              unstructured_grid%vertices(vertex_id)%z
            count = count +1
          endif
        enddo
        unstructured_grid%face_centroid(face_id)%id = face_id
        unstructured_grid%face_centroid(face_id)%x  = &
          unstructured_grid%face_centroid(face_id)%x/count
        unstructured_grid%face_centroid(face_id)%y  = &
          unstructured_grid%face_centroid(face_id)%y/count
        unstructured_grid%face_centroid(face_id)%z  = &
          unstructured_grid%face_centroid(face_id)%z/count
      endif
    enddo
  enddo

  allocate(unstructured_grid%face_to_cell_ghosted(size(face_to_cell,1), &
                                                  size(face_to_cell,2)))
  unstructured_grid%face_to_cell_ghosted = face_to_cell
  allocate(unstructured_grid%cell_to_face_ghosted(size(cell_to_face,1), &
                                                  size(cell_to_face,2)))
  unstructured_grid%cell_to_face_ghosted(:,:) = cell_to_face(:,:)

  deallocate(cell_to_face)
  deallocate(face_to_cell)
  deallocate(vertex_to_cell)

  UGridComputeInternConnect => connections

end function UGridComputeInternConnect

! ************************************************************************** !

subroutine UGridPopulateConnection(unstructured_grid, connection, iface_cell, &
                                   iconn, ghosted_id, option)
  ! 
  ! Computes details of connection (area, dist, etc)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/30/09
  ! 

  use Connection_module
  use Utility_module, only : DotProduct
  use Option_module
  use Grid_Unstructured_Cell_module
  use Geometry_module
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface_cell
  PetscInt :: iconn
  PetscInt :: ghosted_id
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  PetscInt :: face_id
  PetscInt :: ivert,vert_id
  PetscInt :: face_type
  PetscReal :: v1(3),v2(3),n_dist(3), dist
  type(point3d_type) :: vertex_8(8)
  type(plane_type) :: plane
  type(point3d_type) :: point, vertex1, vertex2, vertex3, intercept
  character(len=MAXWORDLENGTH) :: word
  
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      if (iface_cell == 0) then
        write(word,*) ghosted_id
        option%io_buffer = 'Face id undefined for cell ' // &
          trim(adjustl(word)) // &
          ' in boundary condition.  Should this be a source/sink?'
        call printErrMsgByRank(option)
      endif
      ! Compute cell centeroid
      v2 = 0.d0
      do ivert = 1, unstructured_grid%cell_vertices(0, ghosted_id)
        vert_id = unstructured_grid%cell_vertices(ivert, ghosted_id)
        vertex_8(ivert)%x = unstructured_grid%vertices(vert_id)%x
        vertex_8(ivert)%y = unstructured_grid%vertices(vert_id)%y
        vertex_8(ivert)%z = unstructured_grid%vertices(vert_id)%z
      enddo
      v2 = UCellComputeCentroid(unstructured_grid%cell_type(ghosted_id), &
                                vertex_8,option)
! Instead of connecting centroid with face center, calculate the shortest
! distance between the centroid and face and use that distance - geh
#if 0
      ! Get face-centroid vector
      face_id = unstructured_grid%cell_to_face_ghosted(iface_cell, ghosted_id)
      v1(1) = unstructured_grid%face_centroid(face_id)%x
      v1(2) = unstructured_grid%face_centroid(face_id)%y
      v1(3) = unstructured_grid%face_centroid(face_id)%z
      
#endif
      !TODO(geh): add support for a quad face
      !TODO(geh): replace %face_to_vertex array with function that returns vertices
      !           based on cell type and iface
      point%x = v2(1)
      point%y = v2(2)
      point%z = v2(3)
      face_id = unstructured_grid%cell_to_face_ghosted(iface_cell, ghosted_id)
      face_type = UCellGetFaceType(unstructured_grid%cell_type(ghosted_id),face_id,option)
      vertex1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
      vertex2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
      if (face_type == LINE_FACE_TYPE) then
        call UCellGetLineIntercept(vertex1,vertex2,point,intercept)
      else
        vertex3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
        call GeometryComputePlaneWithPoints(plane,vertex1,vertex2,vertex3)
        call GeometryProjectPointOntoPlane(plane,point,intercept)
      endif
      
      ! Compute distance vector: cell_center - face_centroid
      v1(1) = v2(1) - intercept%x
      v1(2) = v2(2) - intercept%y
      v1(3) = v2(3) - intercept%z
      
      dist = sqrt(DotProduct(v1, v1))
      n_dist = v1/dist
      connection%dist(0, iconn) = dist
      connection%dist(1, iconn) = n_dist(1)
      connection%dist(2, iconn) = n_dist(2)
      connection%dist(3, iconn) = n_dist(3)
      connection%area(iconn)    = unstructured_grid%face_area(face_id)
      connection%intercp(1,iconn)= intercept%x
      connection%intercp(2,iconn)= intercept%y
      connection%intercp(3,iconn)= intercept%z
      connection%face_id(iconn)  = face_id
      
  end select
  
end subroutine UGridPopulateConnection

! ************************************************************************** !

subroutine UGridComputeCoord(unstructured_grid,option, &
                             grid_x,grid_y,grid_z, &
                             x_min,x_max,y_min,y_max,z_min,z_max)
  ! 
  ! Computes coordinates in x,y,z of unstructured grid cells
  ! 11/2/10 Major rewrite to extend coordinates to ghost cells SP and GEH
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_8(8)
  PetscReal :: centroid(3)
  PetscErrorCode :: ierr 

  do ghosted_id = 1, unstructured_grid%ngmax 
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    centroid = UCellComputeCentroid(unstructured_grid%cell_type(ghosted_id), &
                                    vertex_8,option)
    grid_x(ghosted_id) = centroid(1)
    grid_y(ghosted_id) = centroid(2)
    grid_z(ghosted_id) = centroid(3)
  enddo

  do ivertex = 1, unstructured_grid%num_vertices_local
    if (x_max < unstructured_grid%vertices(ivertex)%x) &
      x_max = unstructured_grid%vertices(ivertex)%x
    if (x_min > unstructured_grid%vertices(ivertex)%x) &
      x_min = unstructured_grid%vertices(ivertex)%x
    if (y_max < unstructured_grid%vertices(ivertex)%y) &
      y_max = unstructured_grid%vertices(ivertex)%y
    if (y_min > unstructured_grid%vertices(ivertex)%y) &
      y_min = unstructured_grid%vertices(ivertex)%y
    if (z_max < unstructured_grid%vertices(ivertex)%z) &
      z_max = unstructured_grid%vertices(ivertex)%z
    if (z_min > unstructured_grid%vertices(ivertex)%z) &
      z_min = unstructured_grid%vertices(ivertex)%z
  enddo
      
end subroutine UGridComputeCoord

! ************************************************************************** !

subroutine UGridComputeVolumes(unstructured_grid,option,volume)
  ! 
  ! Computes volume of unstructured grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/06/09
  ! 

  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  Vec :: volume
  

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_8(8)
  PetscReal, pointer :: volume_p(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(volume,volume_p,ierr);CHKERRQ(ierr)

  do local_id = 1, unstructured_grid%nlmax
    ! ghosted_id = local_id on unstructured grids
    ghosted_id = local_id
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    volume_p(local_id) = UCellComputeVolume(unstructured_grid%cell_type( &
                           ghosted_id),vertex_8,option)
  enddo
      
  call VecRestoreArrayF90(volume,volume_p,ierr);CHKERRQ(ierr)

end subroutine UGridComputeVolumes

! ************************************************************************** !

subroutine UGridComputeAreas(unstructured_grid,option,area)
  ! 
  ! Computes area of unstructured grid cells
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/07/2012
  ! 

  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  Vec :: area
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_4(4)
  PetscReal, pointer :: area_p(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(area,area_p,ierr);CHKERRQ(ierr)

  do local_id = 1, unstructured_grid%nlmax
    ! ghosted_id = local_id on unstructured grids
    ghosted_id = local_id
    if (unstructured_grid%cell_vertices(0,ghosted_id) > 4 ) then
      option%io_buffer = 'ERROR: In UGridComputeAreas the no. of vertices > 4'
      call printErrMsg(option)
    endif
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_4(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_4(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_4(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    area_p(local_id) = UCellComputeArea(unstructured_grid%cell_type( &
                           ghosted_id),vertex_4,option)
  enddo
      
  call VecRestoreArrayF90(area,area_p,ierr);CHKERRQ(ierr)

end subroutine UGridComputeAreas

! ************************************************************************** !

subroutine UGridComputeQuality(unstructured_grid,option)
  ! 
  ! Computes quality of unstructured grid cells
  ! geh: Yes, this is very primitive as mesh quality can be based on any
  ! number of metrics (e.g., see http://cubit.sandia.gov/help-version8/
  ! Chapter_5/Mesh_Quality_Assessment.html).  However, the current edge
  ! length-based formula gives a ballpark estimate.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/17/12
  ! 

  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_8(8)
  PetscReal :: quality, mean_quality, max_quality, min_quality
  PetscErrorCode :: ierr

  mean_quality = 0.d0
  max_quality = -1.d20
  min_quality = 1.d20
  
  do local_id = 1, unstructured_grid%nlmax
    ! ghosted_id = local_id on unstructured grids
    ghosted_id = local_id
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    quality = UCellQuality(unstructured_grid%cell_type( &
                           ghosted_id),vertex_8,option)
    if (quality < min_quality) min_quality = quality
    if (quality > max_quality) max_quality = quality
    mean_quality = mean_quality + quality
  enddo

  call MPI_Allreduce(MPI_IN_PLACE,mean_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  mean_quality = mean_quality / unstructured_grid%nmax

  call MPI_Allreduce(MPI_IN_PLACE,max_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  call MPI_Allreduce(MPI_IN_PLACE,min_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)

  if (OptionPrintToScreen(option)) then
    write(*,'(/," ---------- Mesh Quality ----------", &
            & /,"   Mean Quality: ",es10.2, &
            & /,"   Max Quality : ",es10.2, &
            & /,"   Min Quality : ",es10.2, &
            & /," ----------------------------------",/)') &
              mean_quality, max_quality, min_quality
  endif

end subroutine UGridComputeQuality

! ************************************************************************** !

subroutine UGridEnsureRightHandRule(unstructured_grid,x,y,z,nG2A,nl2G,option)
  ! 
  ! Rearranges order of vertices within each cell
  ! so that when the right hand rule is applied to a
  ! face, the thumb points away from the centroid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/11
  ! 

  use Option_module
  use Utility_module, only : DotProduct, CrossProduct
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscReal :: x(:), y(:), z(:)
  PetscInt :: nG2A(:)
  PetscInt :: nL2G(:)
  type(option_type) :: option

  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(point3d_type) :: point, point1, point2, point3
  type(plane_type) :: plane1
  PetscReal :: v1(3),v2(3),vcross(3),magnitude
  PetscReal :: distance
  PetscInt :: cell_vertex_ids_before(8), cell_vertex_ids_after(8)
  PetscInt :: face_vertex_ids(4)
  type(point3d_type) :: vertex_8(8)
  PetscInt :: ivertex, vertex_id
  PetscInt :: num_vertices, iface, cell_type, num_faces, face_type, i
  PetscInt :: num_face_vertices
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: error_found

  error_found = PETSC_FALSE
  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id
    cell_type = unstructured_grid%cell_type(local_id)
    num_vertices = UCellGetNVertices(cell_type,option)
    cell_vertex_ids_before(1:num_vertices) = &
      unstructured_grid%cell_vertices(1:num_vertices,ghosted_id)
    cell_vertex_ids_after = cell_vertex_ids_before
    ! point is the centroid of cell
    point%x = x(ghosted_id)
    point%y = y(ghosted_id)
    point%z = z(ghosted_id)
    num_faces = UCellGetNFaces(cell_type,option)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface,option)
      num_face_vertices = UCellGetNFaceVertices(cell_type,iface,option)
      call UCellGetFaceVertices(option,cell_type,iface,face_vertex_ids)
      ! Need to find distance of a point (centroid) from a line (formed by
      ! joining vertices of a line)
      point1 = &
        unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(1)))
      point2 = &
        unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(2)))
      if (face_type == LINE_FACE_TYPE) then
        point3%x = point2%x
        point3%y = point2%y
        point3%z = point2%z + 1.d0
      else
        point3 = &
          unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(3)))
      endif

      call GeometryComputePlaneWithPoints(plane1,point1,point2,point3)
      distance = GeomComputeDistanceFromPlane(plane1,point)

      if (distance > 0.d0) then
        ! need to swap so that distance is negative (point lies below plane)
        if (cell_type == TRI_TYPE .or. cell_type == QUAD_TYPE) then
          ! Error message for 2D cell type
          option%io_buffer = 'Cell:'
          write(string,'(i13)') nG2A(nL2G(local_id))
          option%io_buffer = trim(option%io_buffer) // ' ' // &
            trim(adjustl(string)) // ' of type "' // &
            trim(UCellTypeToWord(cell_type,option)) // '" with vertices:'
          do i = 1, num_vertices
            write(string,'(i13)') &
              unstructured_grid%vertex_ids_natural(cell_vertex_ids_before(i))
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(adjustl(string))
          enddo
          option%io_buffer = trim(option%io_buffer) // &
            ' violates right hand rule at face "' // &
            trim(UCellFaceTypeToWord(face_type,option)) // &
            '" based on face vertices:'
          do i = 1, num_face_vertices
            write(string,'(i13)') face_vertex_ids(i)
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(adjustl(string))
          enddo
          do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
            vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
            vertex_8(ivertex)%x = &
              unstructured_grid%vertices(vertex_id)%x
            vertex_8(ivertex)%y = &
              unstructured_grid%vertices(vertex_id)%y
            vertex_8(ivertex)%z = &
              unstructured_grid%vertices(vertex_id)%z
          enddo
          write(string,'(es12.4)') &
            UCellComputeArea(cell_type,vertex_8,option)
          option%io_buffer = trim(option%io_buffer) // ' and area: ' // &
            trim(adjustl(string)) // '.'
          call printMsgAnyRank(option)
          error_found = PETSC_TRUE
        else
          ! Error message for 3D cell type
          option%io_buffer = 'Cell:'
          write(string,'(i13)') nG2A(nL2G(local_id))
          option%io_buffer = trim(option%io_buffer) // ' ' // &
            trim(adjustl(string)) // ' of type "' // &
            trim(UCellTypeToWord(cell_type,option)) // '" with vertices:'
          do i = 1, num_vertices
            write(string,'(i13)') &
              unstructured_grid%vertex_ids_natural(cell_vertex_ids_before(i))
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(adjustl(string))
          enddo
          option%io_buffer = trim(option%io_buffer) // &
            ' violates right hand rule at face "' // &
            trim(UCellFaceTypeToWord(face_type,option)) // &
            '" based on face vertices:'
          do i = 1, num_face_vertices
            write(string,'(i13)') face_vertex_ids(i)
            option%io_buffer = trim(option%io_buffer) // ' ' // &
              trim(adjustl(string))
          enddo
          do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
            vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
            vertex_8(ivertex)%x = &
              unstructured_grid%vertices(vertex_id)%x
            vertex_8(ivertex)%y = &
              unstructured_grid%vertices(vertex_id)%y
            vertex_8(ivertex)%z = &
              unstructured_grid%vertices(vertex_id)%z
          enddo
          write(string,'(es12.4)') &
            UCellComputeVolume(cell_type,vertex_8,option)
          option%io_buffer = trim(option%io_buffer) // ' and volume: ' // &
            trim(adjustl(string)) // '.'
          call printMsgAnyRank(option)
          error_found = PETSC_TRUE
        endif
      endif
    enddo
  enddo
  
  if (error_found) then
    option%io_buffer = 'Cells founds that violate right hand rule.'
    call printErrMsgByRank(option)
  endif

end subroutine UGridEnsureRightHandRule

! ************************************************************************** !

subroutine UGridGetCellFromPoint(x,y,z,unstructured_grid,option,icell)
  ! 
  ! Returns the cell that encompasses a point in space
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/09
  ! 

  use Option_module
  use Geometry_module  

  implicit none
  
  PetscReal :: x, y, z
  PetscInt :: icell
  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  
  PetscInt :: cell_type, num_faces, iface, face_type
  PetscInt :: vertex_ids(4)
  type(plane_type) :: plane1, plane2
  type(point3d_type) :: point, point1, point2, point3, point4
  PetscInt :: local_id, ghosted_id
  PetscReal :: distance
  PetscBool :: inside
  
  icell = 0
  
  point%x = x
  point%y = y
  point%z = z
  
  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type(ghosted_id)
    num_faces = UCellGetNFaces(cell_type,option)
 
    ! vertices should be ordered counter-clockwise so that a cross product
    ! of the two vectors v1-v2 and v1-v3 points outward.
    ! if the distance from the point to the planes making up the faces is always
    ! negative using counter-clockwise ordering, the point is within the volume
    ! encompassed by the faces.
    inside = PETSC_TRUE
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface,option)
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)
      point1 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(1),ghosted_id))
      point2 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(2),ghosted_id))
      point3 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(3),ghosted_id))
      call GeometryComputePlaneWithPoints(plane1,point1,point2,point3)
      distance = GeomComputeDistanceFromPlane(plane1,point)
      if (distance > 0.d0) then
        inside = PETSC_FALSE
        exit
      endif
      if (face_type == QUAD_FACE_TYPE) then
        point4 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(4),ghosted_id))
        call GeometryComputePlaneWithPoints(plane2,point3,point4,point1)
        distance = GeomComputeDistanceFromPlane(plane2,point)
        if (distance > 0.d0) then
          inside = PETSC_FALSE
          exit
        endif
      endif
    enddo
    
    if (inside) then
      icell = local_id
      exit
    endif

  enddo
  
end subroutine UGridGetCellFromPoint

! ************************************************************************** !

subroutine UGridGetCellsInRectangle(x_min,x_max,y_min,y_max,z_min,z_max, &
                                    unstructured_grid,option,num_cells, &
                                    cell_ids,cell_face_ids)
  ! 
  ! Returns the cell that encompasses a point in space
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/09
  ! 
  use Option_module
  use Utility_module, only : reallocateIntArray
  use Geometry_module  
  
  implicit none
                  
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  PetscInt :: num_cells
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: cell_face_ids(:)
  
  PetscInt :: cell_type, num_faces, iface, face_type
  PetscInt :: vertex_ids(4)
  PetscInt :: num_vertices, ivertex
  PetscInt :: local_id, ghosted_id
  type(point3d_type) :: point
  
  PetscReal :: x_min_adj, x_max_adj, y_min_adj, y_max_adj, z_min_adj, z_max_adj
  PetscReal :: pert
  PetscBool :: in_rectangle
  
  PetscInt, pointer :: temp_cell_array(:), temp_face_array(:)
  PetscInt :: temp_array_size
  
  temp_array_size = 100
  allocate(temp_cell_array(temp_array_size))
  allocate(temp_face_array(temp_array_size))
  temp_cell_array = 0
  temp_face_array = 0
  
  ! enlarge box slightly
  pert = max(1.d-8*(x_max-x_min),1.d-8)
  x_min_adj = x_min - pert 
  x_max_adj = x_max + pert 
  pert = max(1.d-8*(y_max-y_min),1.d-8)
  y_min_adj = y_min - pert 
  y_max_adj = y_max + pert 
  pert = max(1.d-8*(z_max-z_min),1.d-8)
  z_min_adj = z_min - pert 
  z_max_adj = z_max + pert 
  
  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type(ghosted_id)
    num_faces = UCellGetNFaces(cell_type,option)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface,option)
      num_vertices = UCellGetNFaceVertices(cell_type,iface,option)
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)
      in_rectangle = PETSC_TRUE
      do ivertex = 1, num_vertices
        point = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(ivertex),ghosted_id))
        if (point%x < x_min_adj .or. &
            point%x > x_max_adj .or. &
            point%y < y_min_adj .or. &
            point%y > y_max_adj .or. &
            point%z < z_min_adj .or. &
            point%z > z_max_adj) then
          in_rectangle = PETSC_FALSE
          exit
        endif
      enddo
     
      if (in_rectangle) then
        num_cells = num_cells + 1
        if (num_cells > temp_array_size) then
          call reallocateIntArray(temp_cell_array,temp_array_size)
          temp_array_size = temp_array_size / 2 ! convert back for next call
          call reallocateIntArray(temp_face_array,temp_array_size)
        endif
        temp_cell_array(num_cells) = local_id
        temp_face_array(num_cells) = iface
      endif

    enddo
  enddo
  
  allocate(cell_ids(num_cells))
  allocate(cell_face_ids(num_cells))
  cell_ids = temp_cell_array(1:num_cells)
  cell_face_ids = temp_face_array(1:num_cells)
  deallocate(temp_cell_array)
  nullify(temp_cell_array)
  deallocate(temp_face_array)
  nullify(temp_face_array)
  
end subroutine UGridGetCellsInRectangle

! ************************************************************************** !

subroutine UGridMapSideSet(unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)
  ! 
  ! Maps a global boundary side set to the faces of local
  ! ghosted cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/16/11
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscInt :: face_vertices(:,:)
  PetscInt :: n_ss_faces
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)
  
  Mat :: Mat_vert_to_face
  Vec :: Vertex_vec, Face_vec
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscReal :: real_array4(4)
  PetscInt, allocatable :: boundary_faces(:)
  PetscInt, allocatable :: temp_int(:,:)
  PetscInt :: boundary_face_count
  PetscInt :: mapped_face_count
  PetscInt :: nfaces, nvertices
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: local_id
  PetscInt :: cell_type
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ivertex, cell_id, vertex_id_local
  PetscErrorCode :: ierr
  PetscReal :: min_verts_req
  PetscInt :: largest_vert_id, v_id_n
  Vec :: sideset_vert_vec
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  IS :: is_tmp1, is_tmp2
  VecScatter :: scatter_gton


  ! fill matrix with boundary faces of local cells
  ! count up the number of boundary faces
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    nfaces = UCellGetNFaces(unstructured_grid%cell_type(local_id),option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
      endif
    enddo
  enddo

  call MatCreateAIJ(option%mycomm, &
                       boundary_face_count, &
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       unstructured_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face, &
                       ierr);CHKERRQ(ierr)
  call MatZeroEntries(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  real_array4 = 1.d0

  offset=0
  call MPI_Exscan(boundary_face_count,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(boundary_faces(boundary_face_count))
  boundary_faces = 0
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
        boundary_faces(boundary_face_count) = face_id
        call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                        int_array4)

        ! For this matrix:
        !   irow = local face id
        !   icol = natural (global) vertex id
        do ivertex = 1, nvertices
          vertex_id_local = &
            unstructured_grid%cell_vertices(int_array4(ivertex),local_id)
          int_array4_0(ivertex) = &
            unstructured_grid%vertex_ids_natural(vertex_id_local)-1
        enddo
        call MatSetValues(Mat_vert_to_face,1,boundary_face_count-1+offset, &
                          nvertices,int_array4_0,real_array4, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
  enddo

  call MatAssemblyBegin(Mat_vert_to_face,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecCreateMPI(option%mycomm,PETSC_DETERMINE, &
                    unstructured_grid%num_vertices_global, &
                    Vertex_vec,ierr);CHKERRQ(ierr)
  call VecZeroEntries(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyBegin(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(Vertex_vec,ierr);CHKERRQ(ierr)

  ! For this vector:
  !   irow = natural (global) vertex id
  nvertices = 0
  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        nvertices = nvertices + 1
      endif
    enddo
  enddo

  offset=0
  call MPI_Exscan(nvertices,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array(nvertices))
  do local_id = 1, nvertices 
    int_array(local_id) = (local_id-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nvertices, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'is_tmp1_' // trim(region_name) // '_subsurf.out'
  else
    string = 'is_tmp1_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call ISView(is_tmp1,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  nvertices = 0
  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        nvertices = nvertices + 1
        int_array(nvertices) = face_vertices(ivertex,iface)-1
      endif
    enddo
  enddo

  call ISCreateGeneral(option%mycomm,nvertices, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'is_tmp2_' // trim(region_name) // '_subsurf.out'
  else
    string = 'is_tmp2_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call ISView(is_tmp2,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call VecCreateMPI(option%mycomm,nvertices, PETSC_DETERMINE, &
                    sideset_vert_vec,ierr);CHKERRQ(ierr)
  call VecSet(sideset_vert_vec,1.d0,ierr);CHKERRQ(ierr)

  call VecScatterCreate(sideset_vert_vec,is_tmp1, &
                        Vertex_vec,is_tmp2,scatter_gton,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp1,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp2,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'scatter_gton_' // trim(region_name) // '_subsurf.out'
  else
    string = 'scatter_gton_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(scatter_gton,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call VecScatterBegin(scatter_gton,sideset_vert_vec,Vertex_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_gton,sideset_vert_vec,Vertex_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_gton,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Vertex_vec_' // trim(region_name) // '_global' // &
              '_subsurf.out'
  else
    string = 'Vertex_vec_' // trim(region_name) // '_global' // &
              '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(Vertex_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  call VecCreateMPI(option%mycomm,boundary_face_count,PETSC_DETERMINE,Face_vec, &
                    ierr);CHKERRQ(ierr)
  call MatMult(Mat_vert_to_face,Vertex_vec,Face_vec,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Face_vec_' // trim(region_name) // '_global_subsurf.out'
  else
    string = 'Face_vec_' // trim(region_name) // '_global_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(Face_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  allocate(temp_int(MAX_FACE_PER_CELL,boundary_face_count))
  temp_int = 0
  
  mapped_face_count = 0
  if ( unstructured_grid%grid_type == THREE_DIM_GRID) then
    min_verts_req = 3.d0
  else
    min_verts_req = 2.d0
  endif
  
  call VecGetArrayF90(Face_vec,vec_ptr,ierr);CHKERRQ(ierr)
  ! resulting vec contains the number of natural vertices in the sideset that
  ! intersect a local face
  do iface = 1, boundary_face_count
    face_id = boundary_faces(iface)
    if (vec_ptr(iface) >= min_verts_req) then ! 3 or more vertices in sideset
      ! need to ensure that the right number of vertices are included
      cell_id = unstructured_grid%face_to_cell_ghosted(1,face_id)
      cell_type = unstructured_grid%cell_type(cell_id)
      nfaces = UCellGetNFaces(cell_type,option)
      nvertices = 0
      do iface2 = 1, nfaces
        face_id2 = unstructured_grid%cell_to_face_ghosted(iface2,cell_id)
        if (face_id == face_id2) then
          nvertices = UCellGetNFaceVertices(cell_type,iface2,option)
          exit
        endif
      enddo
      if (nvertices == 0) then ! the case if not found 
        option%io_buffer = 'Face not found in UGridMapSideSet'
        call printErrMsgByRank(option)
      endif
      if (abs(nvertices - vec_ptr(iface)) < 0.5d0) then
        mapped_face_count = mapped_face_count + 1
        temp_int(1,mapped_face_count) = cell_id
        temp_int(2,mapped_face_count) = iface2
      endif
    endif
  enddo
  call VecRestoreArrayF90(Face_vec,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(boundary_faces)
  
  allocate(cell_ids(mapped_face_count))
  allocate(face_ids(mapped_face_count))
  
  cell_ids(:) = temp_int(1,1:mapped_face_count)
  face_ids(:) = temp_int(2,1:mapped_face_count)

  call MatDestroy(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  call VecDestroy(Face_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(sideset_vert_vec,ierr);CHKERRQ(ierr)
  
end subroutine UGridMapSideSet

! ************************************************************************** !

subroutine UGridMapSideSet2(unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)
  ! 
  ! Maps a global boundary side set to the faces of local
  ! ghosted cells
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/21/17
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscInt :: face_vertices(:,:)
  PetscInt :: n_ss_faces
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)

  Mat :: Mat_vert_to_face
  Mat :: Mat_region_vert_to_face
  Mat :: Mat_region_face_to_vert
  Mat :: Mat_face
  Mat :: Mat_face_loc
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscInt, allocatable :: boundary_faces(:)
  PetscInt, allocatable :: temp_int(:,:)
  PetscInt :: boundary_face_count
  PetscInt :: mapped_face_count
  PetscInt :: nfaces, nvertices
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: local_id
  PetscInt :: cell_type
  PetscInt :: ivertex, cell_id, vertex_id_local
  PetscInt :: largest_vert_id, v_id_n
  PetscInt :: offset
  PetscInt, pointer :: cell_ids_for_all_boundary_faces(:)
  PetscInt, pointer :: face_ids_for_all_boundary_faces(:)
  PetscInt, pointer :: ia_p(:),ja_p(:)
  PetscInt :: nrow
  PetscInt :: row, col
  PetscInt :: min_nverts
  PetscInt :: ii, jj

  PetscReal :: max_value
  PetscReal :: real_array4(4)
  PetscReal :: min_verts_req

  PetscErrorCode :: ierr

  PetscBool :: done
  PetscScalar, pointer :: aa_v(:)

  ! fill matrix with boundary faces of local cells
  ! count up the number of boundary faces
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    nfaces = UCellGetNFaces(unstructured_grid%cell_type(local_id),option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
      endif
    enddo
  enddo

  call MatCreateAIJ(option%mycomm, &
                       boundary_face_count, &
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       unstructured_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face, &
                       ierr);CHKERRQ(ierr)
  call MatZeroEntries(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  real_array4 = 1.d0

  offset=0
  call MPI_Exscan(boundary_face_count,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(face_ids_for_all_boundary_faces(boundary_face_count))
  allocate(cell_ids_for_all_boundary_faces(boundary_face_count))
  face_ids_for_all_boundary_faces = 0
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
        face_ids_for_all_boundary_faces(boundary_face_count) = iface
        cell_ids_for_all_boundary_faces(boundary_face_count) = local_id
        call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                        int_array4)

        ! For this matrix:
        !   irow = local face id
        !   icol = natural (global) vertex id
        do ivertex = 1, nvertices
          vertex_id_local = &
            unstructured_grid%cell_vertices(int_array4(ivertex),local_id)
          int_array4_0(ivertex) = &
            unstructured_grid%vertex_ids_natural(vertex_id_local)-1
        enddo
        call MatSetValues(Mat_vert_to_face,1,boundary_face_count-1+offset, &
                          nvertices,int_array4_0,real_array4, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
  enddo

  call MatAssemblyBegin(Mat_vert_to_face,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call MatCreateAIJ(option%mycomm, &
                       n_ss_faces, &
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       unstructured_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_region_vert_to_face, &
                       ierr);CHKERRQ(ierr)
  call MatZeroEntries(Mat_region_vert_to_face,ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(n_ss_faces,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        call MatSetValue(Mat_region_vert_to_face, &
                         iface-1+offset, &
                         face_vertices(ivertex,iface)-1, &
                         1.d0, &
                         INSERT_VALUES,ierr);CHKERRQ(ierr)

      endif
    enddo
  enddo

  call MatAssemblyBegin(Mat_region_vert_to_face,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_region_vert_to_face,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_region_vert_to_face_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_region_vert_to_face_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_region_vert_to_face,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call MatTranspose(Mat_region_vert_to_face, MAT_INITIAL_MATRIX, &
                    Mat_region_face_to_vert, ierr); CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_region_face_to_vert_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_region_face_to_vert_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_region_face_to_vert,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call MatMatMult(Mat_vert_to_face, Mat_region_face_to_vert, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                           Mat_face, ierr); CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_face_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_face_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_face,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  if (option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(Mat_face,MAT_INITIAL_MATRIX,Mat_face_loc, &
                              ierr);CHKERRQ(ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(Mat_face_loc,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(Mat_face_loc,aa_v,ierr);CHKERRQ(ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(Mat_face,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(Mat_face,aa_v,ierr);CHKERRQ(ierr)
  endif

  min_nverts = 3
  if (unstructured_grid%grid_type == THREE_DIM_GRID) min_nverts = 3
  if (unstructured_grid%grid_type == TWO_DIM_GRID  ) min_nverts = 2

  ! Determine the total number of faces mapped
  mapped_face_count = 0
  row = 1
  col = 0
  do ii = 1,nrow
    max_value = 0.d0
    do jj = ia_p(ii),ia_p(ii + 1) - 1
      if (aa_v(jj) > max_value) then
        max_value = aa_v(jj)
      endif
    enddo
    if (max_value >= min_nverts) then
      mapped_face_count = mapped_face_count + 1
    endif
  enddo

  allocate(cell_ids(mapped_face_count))
  allocate(face_ids(mapped_face_count))

  ! Determine the total number of faces mapped
  mapped_face_count = 0
  row = 1
  col = 0
  do ii = 1,nrow
    max_value = 0.d0
    do jj = ia_p(ii),ia_p(ii + 1) - 1
      if (aa_v(jj) > max_value) then
        max_value = aa_v(jj)
      endif
    enddo
    if (max_value >= min_nverts) then
      mapped_face_count = mapped_face_count + 1
      cell_ids(mapped_face_count) = cell_ids_for_all_boundary_faces(ii)
      face_ids(mapped_face_count) = face_ids_for_all_boundary_faces(ii)
    endif
  enddo

  if (option%mycommsize>1) then
    call MatSeqAIJRestoreArrayF90(Mat_face_loc,aa_v,ierr);CHKERRQ(ierr)
    call MatDestroy(Mat_face_loc,ierr);CHKERRQ(ierr)
  else
    call MatSeqAIJRestoreArrayF90(Mat_face,aa_v,ierr);CHKERRQ(ierr)
  endif

  call MatDestroy(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_region_vert_to_face,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_region_face_to_vert,ierr);CHKERRQ(ierr)
  deallocate(face_ids_for_all_boundary_faces)
  deallocate(cell_ids_for_all_boundary_faces)

end subroutine UGridMapSideSet2

! ************************************************************************** !

subroutine UGridMapBoundFacesInPolVol(unstructured_grid,polygonal_volume, &
                                      region_name,option, &
                                      cell_ids,face_ids)
  ! 
  ! Maps all global boundary cell faces within a
  ! polygonal volume to a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/16/11
  ! 
  use Option_module
  use Geometry_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(polygonal_volume_type) :: polygonal_volume
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)

  PetscInt :: ivertex
  PetscInt :: iface, face_id
  PetscInt :: iface2, face_id2
  PetscInt :: nfaces
  PetscInt :: cell_id, cell_type
  PetscInt :: vertex_id
  type(point3d_type) :: vertex
  PetscInt :: mapped_face_count
  PetscBool :: found
  PetscInt :: boundary_face_count
  PetscInt, pointer :: boundary_faces(:)
  
  nullify(boundary_faces)
  
  call UGridGetBoundaryFaces(unstructured_grid,option,boundary_faces)
  
  if (associated(boundary_faces)) then
  
    boundary_face_count = size(boundary_faces)
    
    mapped_face_count = 0
    do iface = 1, boundary_face_count
      face_id = boundary_faces(iface)
      found = GeometryPointInPolygonalVolume( &
                unstructured_grid%face_centroid(face_id)%x, &
                unstructured_grid%face_centroid(face_id)%y, &
                unstructured_grid%face_centroid(face_id)%z, &
                polygonal_volume,option)
      if (found) then
        mapped_face_count = mapped_face_count + 1
        ! if inside, shift the face earlier in the array to same array space
        boundary_faces(mapped_face_count) = boundary_faces(iface)
      endif
    enddo

    if (mapped_face_count > 0) then
      allocate(cell_ids(mapped_face_count))
      cell_ids = 0
      allocate(face_ids(mapped_face_count))
      face_ids = 0
      do iface = 1, mapped_face_count
        face_id = boundary_faces(iface)
        cell_id = &
          unstructured_grid%face_to_cell_ghosted(1,face_id)
        cell_type = unstructured_grid%cell_type(cell_id)
        nfaces = UCellGetNFaces(cell_type,option)
        found = PETSC_FALSE
        do iface2 = 1, nfaces
          face_id2 = unstructured_grid%cell_to_face_ghosted(iface2,cell_id)
          if (face_id == face_id2) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (.not.found) then
          option%io_buffer = &
            'Boundary face mismatch in UGridMapBoundFacesInPolVol()'
          call printErrMsg(option)
        else
          cell_ids(iface) = cell_id
          face_ids(iface) = iface2
        endif
      enddo
    endif
  
    deallocate(boundary_faces)
    nullify(boundary_faces)
  
  endif  
  
end subroutine UGridMapBoundFacesInPolVol

! ************************************************************************** !

subroutine UGridGetBoundaryFaces(unstructured_grid,option,boundary_faces)
  ! 
  ! Returns an array of ids for cell faces on boundary
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscInt, pointer :: boundary_faces(:)
  type(option_type) :: option
  
  PetscInt :: boundary_face_count
  PetscInt :: nfaces
  PetscInt :: iface
  PetscInt :: face_id
  PetscInt :: local_id
  PetscInt :: cell_type
  PetscErrorCode :: ierr
    
  ! fill matrix with boundary faces of local cells
  ! count up the number of boundary faces
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    nfaces = UCellGetNFaces(unstructured_grid%cell_type(local_id),option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
      endif
    enddo
  enddo

  if (boundary_face_count > 0) then
    allocate(boundary_faces(boundary_face_count))
    boundary_faces = 0
    boundary_face_count = 0
    do local_id = 1, unstructured_grid%nlmax
      cell_type = unstructured_grid%cell_type(local_id)
      nfaces = UCellGetNFaces(cell_type,option)
      do iface = 1, nfaces
        face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
        if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
          ! boundary face, since not connected to 2 cells
          boundary_face_count = boundary_face_count + 1
          boundary_faces(boundary_face_count) = face_id
        endif
      enddo
    enddo
  endif
  
end subroutine UGridGetBoundaryFaces

! ************************************************************************** !

subroutine UGridGrowStencilSupport(unstructured_grid,stencil_width, &
                                   ghosted_level,option)
  ! 
  ! This routine will update the mesh to accomodate larger stencil width.
  ! -1) Stencil support will be increased by one cell at a time.
  ! -2) Find updated list of local+ghost cells (Note: Only the list of
  ! ghost cells get updated).
  ! -3) Find the 'new' ghost cells from the updated list found in (2)
  ! -4) Lastly update the mesh
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/17/12
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  PetscInt :: stencil_width
  PetscInt, pointer :: ghosted_level(:)
  
  Mat :: Mat_vert_to_cell  !
  Mat :: Mat_vert_to_proc  !
  Mat :: Mat_proc_to_vert  !
  
  PetscInt :: offset
  PetscInt :: local_id,ghosted_id
  PetscInt :: ivertex
  PetscInt :: cell_type
  PetscInt :: nvertices
  PetscInt :: vertex_id_local
  PetscInt :: vertex_id_nat
  PetscInt :: ngmax_new
  PetscInt :: swidth

  PetscInt, pointer :: ia_p(:), ja_p(:)
  PetscInt :: nrow,rstart,rend,icol(1)
  PetscOffset :: iia,jja,aaa,iicol,jj
  PetscBool :: done
  PetscScalar :: aa(1)

  PetscReal, allocatable :: real_arrayV(:)
  PetscInt, allocatable :: int_arrayV(:)
  PetscInt, allocatable :: cell_ids_natural(:)
  PetscInt, allocatable :: cell_ids_petsc(:)
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: cids_new(:)
  PetscInt, pointer :: ghosted_level_new(:)
  
  Vec :: Vec_cids_local
  PetscReal, pointer :: vec_ptr(:)

  PetscErrorCode :: ierr
  PetscViewer :: viewer

  IS :: is_from
  IS :: is_to
  
  VecScatter :: vec_scatter

  PetscInt :: nghost_new
  PetscInt,allocatable :: ghost_cids_new(:)
  PetscInt,allocatable :: ghost_cids_new_petsc(:)

  ! There are no ghost cells when running with a single processor, so get out
  ! of here
  if (option%mycommsize==1) return
  
  allocate(real_arrayV(unstructured_grid%max_nvert_per_cell))
  allocate(int_arrayV(unstructured_grid%max_nvert_per_cell))
  real_arrayV=1.d0

  ! Allocate memory for a matrix to saves mesh connectivity
  ! size(Mat_vert_to_cell) = global_num_cell x global_num_vertices
  call MatCreateAIJ(option%mycomm, &
                    unstructured_grid%nlmax, &
                    PETSC_DETERMINE, &
                    PETSC_DETERMINE, &
                    unstructured_grid%num_vertices_global, &
                    unstructured_grid%max_nvert_per_cell, &
                    PETSC_NULL_INTEGER, &
                    unstructured_grid%max_nvert_per_cell, &
                    PETSC_NULL_INTEGER, &
                    Mat_vert_to_cell, &
                    ierr);CHKERRQ(ierr)

  call MatZeroEntries(Mat_vert_to_cell,ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(unstructured_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! Create the mesh connectivity matrix
  do local_id=1,unstructured_grid%nlmax
    cell_type = unstructured_grid%cell_type(local_id)
    nvertices=UCellGetNVertices(cell_type,option)
    do ivertex=1,nvertices
      vertex_id_local=unstructured_grid%cell_vertices(ivertex,local_id)
      vertex_id_nat=unstructured_grid%vertex_ids_natural(vertex_id_local)
      int_arrayV(ivertex)=vertex_id_nat-1
    enddo
    call MatSetValues(Mat_vert_to_cell,1,local_id-1+offset, &
                      nvertices,int_arrayV,real_arrayV, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(Mat_vert_to_cell,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_cell,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  
  ! Create a vector which has natural cell ids in PETSc order
  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax, &
                    PETSC_DETERMINE, &
                    Vec_cids_local,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(Vec_cids_local,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1,unstructured_grid%nlmax
    vec_ptr(local_id) = unstructured_grid%cell_ids_natural(local_id)
  enddo
  call VecRestoreArrayF90(Vec_cids_local,vec_ptr,ierr);CHKERRQ(ierr)

  ! Now begin expanding stencil support
  do swidth = 1,stencil_width

    ! Create a matrix that saves natural id of vertices present on each
    ! processor
    call MatCreateAIJ(option%mycomm, &
                      1, &
                      PETSC_DETERMINE, &
                      PETSC_DETERMINE, &
                      unstructured_grid%num_vertices_global, &
                      unstructured_grid%num_vertices_global, &
                      PETSC_NULL_INTEGER, &
                      unstructured_grid%num_vertices_global, &
                      PETSC_NULL_INTEGER, &
                      Mat_vert_to_proc, &
                      ierr);CHKERRQ(ierr)

    call MatZeroEntries(Mat_vert_to_proc,ierr);CHKERRQ(ierr)

    if (swidth==1) then
      ! When the stencil width counter = 1, loop over only local cells present
      do local_id=1,unstructured_grid%nlmax
        cell_type = unstructured_grid%cell_type(local_id)
        nvertices=UCellGetNVertices(cell_type,option)
        do ivertex=1,nvertices
          vertex_id_local=unstructured_grid%cell_vertices(ivertex,local_id)
          vertex_id_nat=unstructured_grid%vertex_ids_natural(vertex_id_local)
          call MatSetValues(Mat_vert_to_proc,1,option%myrank, &
                            1,vertex_id_nat-1,1.d0,INSERT_VALUES, &
                            ierr);CHKERRQ(ierr)
        enddo
      enddo
    else
      ! When the stencil width counter is > 1, loop over ghosted cells
      do ghosted_id=1,unstructured_grid%ngmax
        cell_type = unstructured_grid%cell_type(ghosted_id)
        nvertices=UCellGetNVertices(cell_type,option)
        do ivertex=1,nvertices
          vertex_id_local=unstructured_grid%cell_vertices(ivertex,ghosted_id)
          vertex_id_nat=unstructured_grid%vertex_ids_natural(vertex_id_local)
          call MatSetValues(Mat_vert_to_proc,1,option%myrank, &
                            1,vertex_id_nat-1,1.d0,INSERT_VALUES, &
                            ierr);CHKERRQ(ierr)
        enddo
      enddo
    endif

    ! Assemble the matrix
    call MatAssemblyBegin(Mat_vert_to_proc,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(Mat_vert_to_proc,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)

    ! Transpose the matrix
    call MatTranspose(Mat_vert_to_proc,MAT_INITIAL_MATRIX, &
                      Mat_proc_to_vert,ierr);CHKERRQ(ierr)
    call MatDestroy(Mat_vert_to_proc,ierr);CHKERRQ(ierr)

    ! Find the number and natural ids of cells (local+ghost) when stencil width
    ! is increased by one.
    call UGridFindCellIDsAfterGrowingStencilWidthByOne( &
                                    Mat_vert_to_cell, &
                                    Mat_proc_to_vert, &
                                    Vec_cids_local, &
                                    cids_new, &
                                    ngmax_new, &
                                    option)

    ! Find additional ghost cells
    call UGridFindNewGhostCellIDsAfterGrowingStencilWidth(unstructured_grid,&
                      cids_new, &
                      ngmax_new, &
                      ghost_cids_new, &
                      ghost_cids_new_petsc, &
                      nghost_new, &
                      option)
                                          
    ! Update the mesh by adding the new ghost cells
    call UGridUpdateMeshAfterGrowingStencilWidth(unstructured_grid,&
          ghost_cids_new,ghost_cids_new_petsc,nghost_new,option)

    ! Update ghosted_level array
    if (swidth==1) then
      ! In this case, ghosted_level will have only two values: 
      !   0 - local cells
      !   1 - ghost cells
      allocate(ghosted_level(unstructured_grid%ngmax))
      do local_id=1,unstructured_grid%nlmax
        ghosted_level(local_id)=0
      enddo
      
      do ghosted_id=unstructured_grid%nlmax+1,unstructured_grid%ngmax
        ghosted_level(ghosted_id)=1
      enddo
    else
    
      ! ghosted_level of all new ghost cells will be 'swidth' 
      allocate(ghosted_level_new(unstructured_grid%ngmax))
      do ghosted_id=1,unstructured_grid%ngmax-nghost_new
        ghosted_level_new(ghosted_id)=ghosted_level(ghosted_id)
      enddo
      
      do ghosted_id=unstructured_grid%ngmax-nghost_new+1,unstructured_grid%ngmax
        ghosted_level_new(ghosted_id)=swidth
      enddo
      
      deallocate(ghosted_level)
      allocate(ghosted_level(unstructured_grid%ngmax))
      ghosted_level=ghosted_level_new
      deallocate(ghosted_level_new)
    endif
    
    ! Free up the memory
    call MatDestroy(Mat_vert_to_proc,ierr);CHKERRQ(ierr)
    call MatDestroy(Mat_proc_to_vert,ierr);CHKERRQ(ierr)
    deallocate(ghost_cids_new)
    deallocate(ghost_cids_new_petsc)
    deallocate(cids_new)

  enddo

  call MatDestroy(Mat_vert_to_cell,ierr);CHKERRQ(ierr)
  call VecDestroy(Vec_cids_local,ierr);CHKERRQ(ierr)

end subroutine UGridGrowStencilSupport

! ************************************************************************** !

subroutine UGridFindCellIDsAfterGrowingStencilWidthByOne(Mat_vert_to_cell, &
                                      Mat_proc_to_vert, &
                                      Vec_cids_local, &
                                      cids_new, &
                                      ngmax_new, &
                                      option)
  ! 
  ! This routine finds the cells that are required on a given processor, if
  ! stencil width is increased by one.
  ! - It used the same algorithm used in UGridMapSidesets, but instead of a
  ! matrix-vector product, matrix-matrix product is used in this subroutine.
  ! - Returns a list of natural ids of all cells (local+ghost)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/17/12
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(option_type) :: option
  Mat :: Mat_vert_to_cell
  Vec :: Vec_cids_local
  !PetscInt, intent(out) :: ngmax_new
  PetscInt :: ngmax_new
  
  Mat :: Mat_proc_to_vert  !
  Mat :: Mat_cell_to_proc  !
  Mat :: Mat_proc_to_cell  !
  Mat :: Mat_cell_to_proc_loc
  
  PetscInt :: offset
  PetscInt :: local_id
  PetscInt :: ivertex
  PetscInt :: cell_type
  PetscInt :: nvertices
  PetscInt :: vertex_id_local
  PetscInt :: vertex_id_nat

  PetscInt, pointer :: ia_p(:), ja_p(:)
  PetscInt :: nrow,rstart,rend,icol(1)
  PetscOffset :: iia,jja,aaa,iicol,jj
  PetscBool :: done
  PetscScalar :: aa(1)

  PetscReal, allocatable :: real_arrayV(:)
  PetscInt, allocatable :: int_arrayV(:)
  PetscInt, allocatable :: cell_ids_natural(:)
  PetscInt, allocatable :: cell_ids_petsc(:)
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: cids_new(:)
  
  Vec :: Vec_cids_ghosted
  PetscReal, pointer :: vec_ptr(:)

  PetscErrorCode :: ierr

  IS :: is_from
  IS :: is_to

  VecScatter :: vec_scatter
  
  ! Perform a matrix-matrix multiplication
  call MatMatMult(Mat_vert_to_cell,Mat_proc_to_vert, &
                    MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,Mat_proc_to_cell, &
                  ierr);CHKERRQ(ierr)

  ! Transpose of the result gives: cell ids that are needed after growing stencil
  ! width by one
  call MatTranspose(Mat_proc_to_cell,MAT_INITIAL_MATRIX, &
                    Mat_cell_to_proc,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_proc_to_cell,ierr);CHKERRQ(ierr)

  if (option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(Mat_cell_to_proc,MAT_INITIAL_MATRIX,Mat_cell_to_proc_loc, &
                              ierr);CHKERRQ(ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(Mat_cell_to_proc_loc, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(Mat_cell_to_proc_loc, aa, aaa, ierr);CHKERRQ(ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(Mat_cell_to_proc, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(Mat_cell_to_proc, aa, aaa, ierr);CHKERRQ(ierr)
  endif

  ! Obtain the PETSc index of all cells required.
  ! Note: We get PETSc index because the rows of Mat_vert_to_cell correspond to
  !       PETSc index of cells.
  ngmax_new = ia_p(2)-ia_p(1)
  allocate(cell_ids_petsc(ngmax_new))
  do jj=ia_p(1),ia_p(2)-1
    cell_ids_petsc(jj)=ja_p(jj)
  enddo

  ! Now, find natural ids of all cells required from PETSc index. This is done
  ! by scattering the 'Vec_cids_local'
  
  ! Create MPI vector to save natural ids of cells required
  call VecCreateMPI(option%mycomm,ngmax_new,PETSC_DETERMINE,Vec_cids_ghosted, &
                    ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(ngmax_new,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array(ngmax_new))
  do jj=1,ngmax_new
    int_array(jj)=INT(jj+offset)
  enddo
  int_array=int_array-1
  
  ! Create a index set to scatter to
  call ISCreateGeneral(option%mycomm,ngmax_new, &
                       int_array,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  ! Create a index set to scatter from
  cell_ids_petsc = cell_ids_petsc - 1
  call ISCreateGeneral(option%mycomm,ngmax_new, &
                       cell_ids_petsc,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)
  cell_ids_petsc = cell_ids_petsc + 1

  ! Create a vec-scatter contex
  call VecScatterCreate(Vec_cids_local,is_from,Vec_cids_ghosted,is_to, &
                        vec_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  ! Scatter the data
  call VecScatterBegin(vec_scatter,Vec_cids_local,Vec_cids_ghosted, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,Vec_cids_local,Vec_cids_ghosted, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Save the natural ids of all cells required after growning the stencil
  ! width by one.
  allocate(cids_new(ngmax_new))
  call VecGetArrayF90(Vec_cids_ghosted,vec_ptr,ierr);CHKERRQ(ierr)
  do jj=1,ngmax_new
    cids_new(jj) = INT(vec_ptr(jj))
  enddo
  call VecRestoreArrayF90(Vec_cids_ghosted,vec_ptr,ierr);CHKERRQ(ierr)
  
  call MatDestroy(Mat_cell_to_proc,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_cell_to_proc_loc,ierr);CHKERRQ(ierr)
  call VecDestroy(Vec_cids_ghosted,ierr);CHKERRQ(ierr)

end subroutine UGridFindCellIDsAfterGrowingStencilWidthByOne

! ************************************************************************** !

subroutine UGridFindNewGhostCellIDsAfterGrowingStencilWidth(unstructured_grid, &
                      cids_new, &
                      ngmax_new, &
                      ghost_cids_new, &
                      ghost_cids_new_petsc, &
                      nghost_new, &
                      option)
  ! 
  ! This routine finds new ghosts cells needed to be saved on a local processor
  ! after stencil width is increased.
  ! - Returns the natural index of new ghosts cells.
  ! - Also, returns the PETSc index of new ghosts cells. (Required for creating
  ! gather/scater contexts in UGridCreateUGDM)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/17/12
  ! 
#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscInt :: cids_new(:)
  type(option_type) :: option
  PetscInt :: ngmax_new

  ! local
  PetscInt :: count
  PetscInt :: ii
  PetscInt :: ghosted_id,local_id,nat_id
  PetscInt :: nghost_new
  PetscInt :: offset
  PetscInt,allocatable :: ghost_cids_new(:)
  PetscInt,allocatable :: ghost_cids_new_petsc(:)
  
  PetscInt,allocatable :: int_array1(:)
  PetscInt,allocatable :: int_array2(:)
  PetscScalar,pointer :: tmp_scl_array(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  Vec :: cids_petsc
  Vec :: ghosts_petsc
  Vec :: cells_on_proc
  Vec :: cids_on_proc
  IS :: is_from
  IS :: is_to
  VecScatter :: vec_scatter

  ! Step-1: Find additional ghost cells

  ! 1.0) Find which cells in 'cids_new' are local or ghost
  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax,PETSC_DETERMINE,cells_on_proc, &
                    ierr);CHKERRQ(ierr)
  
  allocate(int_array1(unstructured_grid%nlmax))
  allocate(tmp_scl_array(unstructured_grid%nlmax))
  do ii=1,unstructured_grid%nlmax
    int_array1(ii)=unstructured_grid%cell_ids_natural(ii)
    tmp_scl_array(ii)=option%myrank
  enddo
  int_array1=int_array1-1
  
  call VecSetValues(cells_on_proc,unstructured_grid%nlmax,int_array1,tmp_scl_array,INSERT_VALUES, &
                    ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  deallocate(tmp_scl_array)
  call VecAssemblyBegin(cells_on_proc,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(cells_on_proc,ierr);CHKERRQ(ierr)
  
  allocate(int_array1(ngmax_new))
  int_array1=cids_new-1
  call ISCreateGeneral(option%mycomm,ngmax_new, &
                       int_array1,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  
  call VecCreateMPI(option%mycomm,ngmax_new,PETSC_DETERMINE,cids_on_proc, &
                    ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(ngmax_new,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array1(ngmax_new))
  do ii=1,ngmax_new
    int_array1(ii)=ii-1+offset
  enddo
  call ISCreateGeneral(option%mycomm,ngmax_new, &
                       int_array1,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)

  call VecScatterCreate(cells_on_proc,is_from,cids_on_proc,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,cells_on_proc,cids_on_proc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,cells_on_proc,cids_on_proc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)
  call VecDestroy(cells_on_proc,ierr);CHKERRQ(ierr)
  
  deallocate(int_array1)

  
  ! 1.1) Create an array containing cell-ids of 'exisiting' ghost cells + 
  !      ghost cells from 'cids_new'
  !
  count = ngmax_new-unstructured_grid%nlmax + &
          unstructured_grid%ngmax - unstructured_grid%nlmax
  allocate(int_array1(count))
  allocate(int_array2(count))

  count=0
  do ii=1,unstructured_grid%ngmax-unstructured_grid%nlmax
    count=count+1
    ghosted_id=ii+unstructured_grid%nlmax
    int_array1(count)=unstructured_grid%cell_ids_natural(ghosted_id)
    int_array2(count)=count
  enddo
  
  call VecGetArrayF90(cids_on_proc,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,ngmax_new
    if (vec_ptr(ii)/=option%myrank) then
      count=count+1
      int_array1(count)=cids_new(ii)
      int_array2(count)=count
    endif
  enddo
  call VecRestoreArrayF90(cids_on_proc,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(cids_on_proc,ierr);CHKERRQ(ierr)

  ! 1.2) Sort the array
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(count,int_array1, &
                                   int_array2,ierr);CHKERRQ(ierr)
  int_array2 = int_array2+1

  ! 1.3) Count the entries in the sorted array which appear only once.
  nghost_new=0
  ii=1
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) nghost_new=nghost_new+1
  
  do ii=2,count-1
    if ((int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))).and. &
       (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) ) nghost_new=nghost_new+1
  enddo

  ii=count
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))) nghost_new=nghost_new+1
  
  ! 1.4) Save the entries in the sorted array which appear only once.
  allocate(ghost_cids_new(nghost_new))
  nghost_new=0
  ii=1
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) then
    nghost_new=nghost_new+1
    ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
  endif
  
  do ii=2,count-1
    if ((int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))).and. &
       (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) ) then
      nghost_new=nghost_new+1
      ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
    endif
  enddo

  ii=count
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))) then
    nghost_new=nghost_new+1
    ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
  endif
  
  deallocate(int_array1)
  deallocate(int_array2)

  ! Step-2: Find PETSc index of additional ghost cells
  call VecCreateMPI(option%mycomm, &
                    unstructured_grid%nlmax, &
                    PETSC_DETERMINE, &
                    cids_petsc,ierr);CHKERRQ(ierr)
  
  offset=0
  call MPI_Exscan(unstructured_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array1(unstructured_grid%nlmax))
  allocate(tmp_scl_array(unstructured_grid%nlmax))
  
  do local_id=1,unstructured_grid%nlmax
    nat_id=unstructured_grid%cell_ids_natural(local_id)
    int_array1(local_id)=nat_id - 1
    tmp_scl_array(local_id)=local_id+offset+0.d0
  enddo
  
  call VecSetValues(cids_petsc,unstructured_grid%nlmax,int_array1,tmp_scl_array,INSERT_VALUES, &
                    ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  deallocate(tmp_scl_array)
  
  call VecAssemblyBegin(cids_petsc,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(cids_petsc,ierr);CHKERRQ(ierr)

  call VecCreateMPI(option%mycomm,nghost_new,PETSC_DETERMINE,ghosts_petsc, &
                    ierr);CHKERRQ(ierr)
  allocate(int_array1(nghost_new))

  int_array1=ghost_cids_new-1
  call ISCreateGeneral(option%mycomm,nghost_new, &
                       int_array1,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(nghost_new,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii=1,nghost_new
    int_array1(ii)=ii-1+offset
  enddo
  call ISCreateGeneral(option%mycomm,nghost_new, &
                       int_array1,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  
  call VecScatterCreate(cids_petsc,is_from,ghosts_petsc,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,cids_petsc,ghosts_petsc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,cids_petsc,ghosts_petsc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)
  
  allocate(ghost_cids_new_petsc(nghost_new))
  call VecGetArrayF90(ghosts_petsc,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,nghost_new
    ghost_cids_new_petsc(ii)=INT(vec_ptr(ii))
  enddo
  call VecRestoreArrayF90(ghosts_petsc,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(cids_petsc,ierr);CHKERRQ(ierr)
  call VecDestroy(ghosts_petsc,ierr);CHKERRQ(ierr)

end subroutine UGridFindNewGhostCellIDsAfterGrowingStencilWidth

! ************************************************************************** !

subroutine UGridUpdateMeshAfterGrowingStencilWidth(unstructured_grid, &
              ghost_cids_new,ghost_cids_new_petsc,nghost_new,option)
  ! 
  ! This routine updates the mesh after additional ghost cells have be found
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/17/12
  ! 


#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  PetscInt :: ngmax_new

  ! local
  PetscInt :: count,count2
  PetscInt :: ii,jj
  PetscInt :: ivertex
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: nghost_new
  PetscInt :: vertex_id_nat
  PetscInt :: vertex_id_loc
  PetscInt :: offset
  PetscInt :: nverts
  PetscInt :: nverts_new
  PetscInt :: cell_type
  PetscInt :: nvertices
    
  PetscReal, pointer :: vec_ptr(:)
  PetscViewer :: viewer

  PetscErrorCode :: ierr
  
  Vec :: elements_petsc
  Vec :: elements_ghost_cells
!  Vec :: cids_nat
!  Vec :: cids_nat2petsc
!  Vec :: needed_ghosts_cids_petsc
  Vec :: vertices_nat
  Vec :: vertices_loc
  
  IS :: is_from
  IS :: is_to
  VecScatter :: vec_scatter
  PetscInt,allocatable :: ghost_cids_new(:)
  PetscInt,allocatable :: ghost_cids_new_petsc(:)
  PetscInt,allocatable :: int_array1(:)
  PetscInt,allocatable :: int_array2(:)
  PetscInt,allocatable :: int_array3(:)
  PetscInt,allocatable :: int_array4(:)
  
  PetscInt,allocatable :: cell_vertices(:,:)
  PetscInt,allocatable :: cell_ids_natural(:)
  PetscInt,allocatable :: ghost_cell_ids_petsc(:)

  PetscScalar,pointer :: tmp_scl_array(:)
  
  ! Step-1: Find the natural ids for vertices forming new ghost cells
  
  ! Create a vector listing the vertices forming each cell in PETSc-order
  call VecCreateMPI(option%mycomm, &
                    unstructured_grid%nlmax*unstructured_grid%max_nvert_per_cell, &
                    PETSC_DETERMINE, &
                    elements_petsc,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(elements_petsc,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr = -9999
  do local_id=1,unstructured_grid%nlmax
    do ivertex = 1, unstructured_grid%cell_vertices(0,local_id)
      vertex_id_loc=unstructured_grid%cell_vertices(ivertex,local_id)
      vertex_id_nat=unstructured_grid%vertex_ids_natural(vertex_id_loc)
      vec_ptr((local_id-1)*unstructured_grid%max_nvert_per_cell+ivertex)=vertex_id_nat
    enddo
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(nghost_new,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array1(nghost_new*unstructured_grid%max_nvert_per_cell))
  
  do ii=1,nghost_new
    do jj=1,unstructured_grid%max_nvert_per_cell
      int_array1((ii-1)*unstructured_grid%max_nvert_per_cell + jj) = &
        (ghost_cids_new_petsc(ii)-1)*unstructured_grid%max_nvert_per_cell + jj-1
    enddo
  enddo

  call ISCreateGeneral(option%mycomm,nghost_new*unstructured_grid%max_nvert_per_cell, &
                       int_array1,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)
  deallocate(int_array1)

  call VecCreateMPI(option%mycomm, &
                    nghost_new*unstructured_grid%max_nvert_per_cell, &
                    PETSC_DETERMINE, &
                    elements_ghost_cells,ierr);CHKERRQ(ierr)
  
  allocate(int_array1(nghost_new*unstructured_grid%max_nvert_per_cell))
  do ii=1,nghost_new
    do jj=1,unstructured_grid%max_nvert_per_cell
      int_array1((ii-1)*unstructured_grid%max_nvert_per_cell + jj) = &
        (ii-1)*unstructured_grid%max_nvert_per_cell + jj-1 + &
        offset*unstructured_grid%max_nvert_per_cell
    enddo
  enddo
  call ISCreateGeneral(option%mycomm,nghost_new*unstructured_grid%max_nvert_per_cell, &
                       int_array1,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)

  deallocate(int_array1)
  
  call VecScatterCreate(elements_petsc,is_from,elements_ghost_cells,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,elements_petsc,elements_ghost_cells, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,elements_petsc,elements_ghost_cells, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Step-2: Given already existing vertices + vertices of new ghost
  !         cells, find additional vertices that need to be now saved.
  !         Once, the new vertices are found, all vertices will be reorder
  !
  ! Note: Algorithm is similar to the one used in UGridDecompose
  
  allocate(int_array1((unstructured_grid%ngmax+nghost_new)*unstructured_grid%max_nvert_per_cell))
  allocate(int_array2((unstructured_grid%ngmax+nghost_new)*unstructured_grid%max_nvert_per_cell))

  ! save vertices of local+ghost cells
  count=0
  do ghosted_id=1,unstructured_grid%ngmax
    cell_type = unstructured_grid%cell_type(ghosted_id)
    nvertices=UCellGetNVertices(cell_type,option)
    do ivertex=1,nvertices
      vertex_id_loc=unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_id_nat=unstructured_grid%vertex_ids_natural(vertex_id_loc)
      
      count=count+1
      int_array1(count) = vertex_id_nat
      int_array2(count) = count
    enddo
  enddo
  
  ! save vertices of new ghost cells
  call VecGetArrayF90(elements_ghost_cells,vec_ptr,ierr);CHKERRQ(ierr)
  do ii =1,nghost_new
    do jj=1,unstructured_grid%max_nvert_per_cell
      if (vec_ptr((ii-1)*unstructured_grid%max_nvert_per_cell+jj)/=-9999) then
        count=count+1
        int_array1(count)=INT(vec_ptr((ii-1)*unstructured_grid%max_nvert_per_cell+jj))
        int_array2(count)=count
      endif
    enddo
  enddo
  call VecRestoreArrayF90(elements_ghost_cells,vec_ptr,ierr);CHKERRQ(ierr)

  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(count,int_array1,int_array2, &
                                   ierr);CHKERRQ(ierr)
  int_array2 = int_array2+1

  ! remove duplicates
  allocate(int_array3(count))
  allocate(int_array4(count))

  int_array3 = 0
  int_array4 = 0
  int_array3(1) = int_array1(int_array2(1))
  count2 = 1
  int_array4(int_array2(1)) = count2
  do ii = 2, count
    jj = int_array1(int_array2(ii))
    if (jj > int_array3(count2)) then
      count2 = count2 + 1
      int_array3(count2) = jj
    endif
    int_array4(int_array2(ii)) = count2
  enddo

  deallocate(int_array1)
  
  allocate(cell_vertices(0:unstructured_grid%max_nvert_per_cell,unstructured_grid%ngmax+nghost_new))
  cell_vertices=0
  
  ! Update vertices for local+ghost cells
  count=0
  do ghosted_id=1,unstructured_grid%ngmax
    cell_type = unstructured_grid%cell_type(ghosted_id)
    nvertices = UCellGetNVertices(cell_type,option)
    cell_vertices(0,ghosted_id)=nvertices
    do ivertex=1,nvertices
      count=count+1
      cell_vertices(ivertex,ghosted_id)=int_array4(count)
    enddo
  enddo

  call VecGetArrayF90(elements_ghost_cells,vec_ptr,ierr);CHKERRQ(ierr)
  do ii =1,nghost_new
    do jj=1,unstructured_grid%max_nvert_per_cell
      if (vec_ptr((ii-1)*unstructured_grid%max_nvert_per_cell+jj)/=-9999) then
        count=count+1
        cell_vertices(jj,ii+unstructured_grid%ngmax)=int_array4(count)
        cell_vertices(0 ,ii+unstructured_grid%ngmax)=cell_vertices(0 ,ii+unstructured_grid%ngmax)+1
      endif
    enddo
  enddo
  call VecRestoreArrayF90(elements_ghost_cells,vec_ptr,ierr);CHKERRQ(ierr)

  ! Make local copies of array which need to be updated.
  
  ! cells
  allocate(cell_ids_natural(unstructured_grid%ngmax+nghost_new))
  allocate(ghost_cell_ids_petsc(unstructured_grid%num_ghost_cells+nghost_new))

  cell_ids_natural(1:unstructured_grid%ngmax)=unstructured_grid%cell_ids_natural(:)
  ghost_cell_ids_petsc(1:unstructured_grid%num_ghost_cells)=unstructured_grid%ghost_cell_ids_petsc(:)
  
  do ii=1,nghost_new
    cell_ids_natural(ii+unstructured_grid%ngmax)=ghost_cids_new(ii)
    ghost_cell_ids_petsc(ii+unstructured_grid%num_ghost_cells)=ghost_cids_new_petsc(ii)
  enddo
  
  ! Save location of vertices needed on a given processor
  call VecCreateMPI(option%mycomm, &
                    PETSC_DETERMINE, &
                    unstructured_grid%num_vertices_global*3, &
                    vertices_nat,ierr);CHKERRQ(ierr)

  allocate(tmp_scl_array(unstructured_grid%num_vertices_local*3))
  allocate(int_array1(unstructured_grid%num_vertices_local*3))
  count=0
  do ivertex=1,unstructured_grid%num_vertices_local
    count=count+1
    tmp_scl_array(count)=unstructured_grid%vertices(ivertex)%x
    int_array1(count)=(unstructured_grid%vertex_ids_natural(ivertex)-1)*3+0

    count=count+1
    tmp_scl_array(count)=unstructured_grid%vertices(ivertex)%y
    int_array1(count)=(unstructured_grid%vertex_ids_natural(ivertex)-1)*3+1

    count=count+1
    tmp_scl_array(count)=unstructured_grid%vertices(ivertex)%z
    int_array1(count)=(unstructured_grid%vertex_ids_natural(ivertex)-1)*3+2
  enddo
  
  call VecSetValues(vertices_nat,count,int_array1, &
                      tmp_scl_array,INSERT_VALUES,ierr);CHKERRQ(ierr)
  call VecAssemblyBegin(vertices_nat,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(vertices_nat,ierr);CHKERRQ(ierr)
  
  deallocate(int_array1)
  deallocate(tmp_scl_array)

  allocate(int_array1(count2*3))
  do ii=1,count2
    int_array1((ii-1)*3+1)=(int_array3(ii)-1)*3
    int_array1((ii-1)*3+2)=(int_array3(ii)-1)*3+1
    int_array1((ii-1)*3+3)=(int_array3(ii)-1)*3+2
  enddo

  call ISCreateGeneral(option%mycomm,count2*3, &
                       int_array1,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  
  offset=0
  call MPI_Exscan(count2*3,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array1(count2*3))
  do ii=1,count2*3
    int_array1(ii)=ii-1+offset
  enddo
  call ISCreateGeneral(option%mycomm,count2*3, &
                       int_array1,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)

  call VecCreateMPI(option%mycomm, &
                    count2*3, &
                    PETSC_DETERMINE, &
                    vertices_loc,ierr);CHKERRQ(ierr)

  call VecScatterCreate(vertices_nat,is_from,vertices_loc,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,vertices_nat,vertices_loc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,vertices_nat,vertices_loc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)
  call VecDestroy(vertices_nat,ierr);CHKERRQ(ierr)

  ! Update the mesh

  ! cell update
  unstructured_grid%ngmax = unstructured_grid%ngmax+nghost_new

  deallocate(unstructured_grid%cell_vertices)
  allocate(unstructured_grid%cell_vertices(0:unstructured_grid%max_nvert_per_cell,unstructured_grid%ngmax))
  unstructured_grid%cell_vertices = cell_vertices
  deallocate(cell_vertices)

  deallocate(unstructured_grid%cell_ids_natural)
  allocate(unstructured_grid%cell_ids_natural(unstructured_grid%ngmax))
  unstructured_grid%cell_ids_natural=cell_ids_natural
  deallocate(cell_ids_natural)
  
  ! vertex update
  deallocate(unstructured_grid%vertex_ids_natural)
  allocate(unstructured_grid%vertex_ids_natural(count2))
  unstructured_grid%vertex_ids_natural(:)=int_array3(1:count2)
  unstructured_grid%num_vertices_local=count2
  
  deallocate(unstructured_grid%vertices)
  allocate(unstructured_grid%vertices(count2))
  
  call VecGetArrayF90(vertices_loc,vec_ptr,ierr);CHKERRQ(ierr)
  unstructured_grid%num_vertices_local=count2
  do ii=1,count2
    unstructured_grid%vertices(ii)%id = int_array3(ii)
    unstructured_grid%vertices(ii)%x = vec_ptr((ii-1)*3+1)
    unstructured_grid%vertices(ii)%y = vec_ptr((ii-1)*3+2)
    unstructured_grid%vertices(ii)%z = vec_ptr((ii-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_loc,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(vertices_loc,ierr);CHKERRQ(ierr)
  
  ! ghost cell update
  unstructured_grid%num_ghost_cells=unstructured_grid%ngmax-unstructured_grid%nlmax
  deallocate(unstructured_grid%ghost_cell_ids_petsc)
  allocate(unstructured_grid%ghost_cell_ids_petsc(unstructured_grid%num_ghost_cells))
  unstructured_grid%ghost_cell_ids_petsc=ghost_cell_ids_petsc
  deallocate(ghost_cell_ids_petsc)
  
  ! cell type update
  deallocate(unstructured_grid%cell_type)
  allocate(unstructured_grid%cell_type(unstructured_grid%ngmax))

  select case(unstructured_grid%grid_type)
    case(THREE_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        ! Determine number of faces and cell-type of the current cell
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(8)
            unstructured_grid%cell_type(ghosted_id) = HEX_TYPE
          case(6)
            unstructured_grid%cell_type(ghosted_id) = WEDGE_TYPE
          case(5)
            unstructured_grid%cell_type(ghosted_id) = PYR_TYPE
          case(4)
            unstructured_grid%cell_type(ghosted_id) = TET_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call printErrMsg(option)
        end select      
      enddo
    case(TWO_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(4)
            unstructured_grid%cell_type = QUAD_TYPE
          case(3)
            unstructured_grid%cell_type = TRI_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call printErrMsg(option)
        end select
      end do
    case default
      option%io_buffer = 'Grid type not recognized: '
      call printErrMsg(option)
  end select

  call VecDestroy(elements_petsc,ierr);CHKERRQ(ierr)
  call VecDestroy(elements_ghost_cells,ierr);CHKERRQ(ierr)

end subroutine UGridUpdateMeshAfterGrowingStencilWidth


end module Grid_Unstructured_module
