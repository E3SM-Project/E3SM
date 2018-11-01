module Grid_Unstructured_Polyhedra_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Geometry_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use PFLOTRAN_Constants_module

  implicit none

  private

#if defined(SCORPIO)
  include "scorpiof.h"
#endif

  public :: UGridPolyhedraRead, &
            UGridPolyhedraDecompose, &
            UGridPolyhedraSetCellCentroids, &
            UGridPolyhedraComputeInternConnect, &
            UGridPolyhedraComputeVolumes, &
            UGridPolyhedraPopulateConnection, &
            UGridPolyhedraGetCellsInRectangle, &
            UGridPolyhedraComputeOutputInfo

contains

! ************************************************************************** !

subroutine UGridPolyhedraRead(ugrid, filename, option)
  ! 
  ! This routine reads unstructured polyhedra grid in ASCII.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 09/29/13
  ! 

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none

  type(grid_unstructured_type) :: ugrid 
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(unstructured_polyhedra_type), pointer :: pgrid
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word, card
  PetscInt :: fileid, icell, iface, irank, ivert
  PetscInt :: remainder, temp_int, num_to_read

  PetscInt :: num_cells
  PetscInt :: num_cells_local
  PetscInt :: num_cells_local_save
  PetscInt :: num_faces
  PetscInt :: num_faces_local
  PetscInt :: num_faces_local_save
  PetscInt :: num_vertices
  PetscInt :: num_vertices_local
  PetscInt :: num_vertices_local_save
  PetscInt :: max_nvert_per_cell
  PetscInt :: max_nface_per_cell
  PetscInt :: max_nvert_per_face
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscErrorCode :: ierr
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscInt, allocatable :: nfaces_per_proc(:)

  pgrid => ugrid%polyhedra_grid
  allocate(nfaces_per_proc(option%mycommsize))
  nfaces_per_proc = 0
  num_faces_local_save = 0

  max_nvert_per_cell = -1
  if (option%myrank == option%io_rank) then
    fileid = 86
    input => InputCreate(fileid,filename,option)

    call InputReadPflotranString(input,option)
    ! read CELL card, though we already know the
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'CELLS'
    call InputErrorMsg(input,option,word,card)
    if (.not.StringCompare(word,card)) then
      option%io_buffer = 'Unrecognized keyword "' // trim(card) // &
        '" in explicit grid file.'
      call printErrMsgByRank(option)
    endif
  
    hint = 'Polyhedra Unstructured Grid CELLS'
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of cells',hint)
  endif

  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
  num_cells = temp_int
  pgrid%num_cells_global = num_cells
  ugrid%nmax = num_cells

  ! divide cells across ranks
  num_cells_local = num_cells/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = num_cells - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  ! allocate memory
  allocate(pgrid%cell_ids(num_cells_local))
  allocate(pgrid%cell_nfaces(num_cells_local))
  allocate(pgrid%cell_nverts(num_cells_local))
  allocate(pgrid%cell_volumes(num_cells_local))
  allocate(pgrid%cell_centroids(num_cells_local))
  pgrid%cell_ids = 0
  pgrid%cell_nfaces = 0
  pgrid%cell_nverts = 0
  pgrid%cell_volumes = 0.d0
  do icell = 1, num_cells_local
    pgrid%cell_centroids(icell)%x = 0.d0
    pgrid%cell_centroids(icell)%y = 0.d0
    pgrid%cell_centroids(icell)%z = 0.d0
  enddo

  ! Read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  max_nface_per_cell = 0
  if (option%myrank == option%io_rank) then

    allocate(temp_real_array(7,num_cells_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do icell = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'cell id',hint)
        temp_real_array(1,icell) = dble(temp_int)

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'num faces',hint)
        temp_real_array(2,icell) = dble(temp_int)
        nfaces_per_proc(irank+1) = nfaces_per_proc(irank+1) + &
          temp_int
        if (temp_int > max_nface_per_cell) &
          max_nface_per_cell = temp_int

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'num vertices',hint)
        temp_real_array(3,icell) = dble(temp_int)
        if (temp_int>max_nvert_per_cell) max_nvert_per_cell = temp_int

        call InputReadDouble(input,option,temp_real_array(4,icell))
        call InputErrorMsg(input,option,'cell x coordinate',hint)

        call InputReadDouble(input,option,temp_real_array(5,icell))
        call InputErrorMsg(input,option,'cell y coordinate',hint)

        call InputReadDouble(input,option,temp_real_array(6,icell))
        call InputErrorMsg(input,option,'cell z coordinate',hint)

        call InputReadDouble(input,option,temp_real_array(7,icell))
        call InputErrorMsg(input,option,'cell volume',hint)
      enddo
      if (nfaces_per_proc(irank+1)>num_faces_local_save) &
        num_faces_local_save = nfaces_per_proc(irank+1)

      if (irank == option%io_rank) then
        ! cells reside on io_rank
        num_faces_local = 0
        pgrid%num_cells_local = num_cells_local
        do icell = 1, num_cells_local
          pgrid%cell_ids(icell) = int(temp_real_array(1,icell))
          pgrid%cell_nfaces(icell) = int(temp_real_array(2,icell))
          pgrid%cell_nverts(icell) = int(temp_real_array(3,icell))
          pgrid%cell_centroids(icell)%x = temp_real_array(4,icell)
          pgrid%cell_centroids(icell)%y = temp_real_array(5,icell)
          pgrid%cell_centroids(icell)%z = temp_real_array(6,icell)
          pgrid%cell_volumes(icell) = temp_real_array(7,icell)

          num_faces_local = num_faces_local + &
            pgrid%cell_nfaces(icell)
        enddo
        pgrid%num_faces_local = num_faces_local
      else
        ! otherwise communicate to other ranks
        int_mpi = num_to_read*7
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo

  else

      ! other ranks pos the recv
      allocate(temp_real_array(7,num_cells_local))
      int_mpi = num_cells_local*7
      call MPI_Recv(temp_real_array,int_mpi, &
                    MPI_DOUBLE_PRECISION,option%io_rank, &
                    MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
      num_faces_local = 0
      pgrid%num_cells_local = num_cells_local
      do icell = 1, num_cells_local
        pgrid%cell_ids(icell) = int(temp_real_array(1,icell))
        pgrid%cell_nfaces(icell) = int(temp_real_array(2,icell))
        pgrid%cell_nverts(icell) = int(temp_real_array(3,icell))
        pgrid%cell_centroids(icell)%x = temp_real_array(4,icell)
        pgrid%cell_centroids(icell)%y = temp_real_array(5,icell)
        pgrid%cell_centroids(icell)%z = temp_real_array(6,icell)
        pgrid%cell_volumes(icell) = temp_real_array(7,icell)

        num_faces_local = num_faces_local + &
          pgrid%cell_nfaces(icell)
      enddo
      num_faces_local_save = num_faces_local
      pgrid%num_faces_local = num_faces_local
  endif
  deallocate(temp_real_array)

  call MPI_Bcast(max_nvert_per_cell,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
  ugrid%max_nvert_per_cell = max_nvert_per_cell
  pgrid%max_nvert_per_cell = max_nvert_per_cell
  allocate(pgrid%cell_vertids(max_nvert_per_cell,num_cells_local))

  if (option%myrank == option%io_rank) then

    call InputReadPflotranString(input,option)
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'FACES'
    call InputErrorMsg(input,option,word,card)
    if (.not.StringCompare(word,card)) then
      option%io_buffer = 'Unrecongnized keyword "' // trim(card) // &
        '" in polyhedra grid file.'
      call printErrMsgByRank(option)
    endif

    hint = 'Polyhedra Unstructured Grid FACES'
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of faces',hint)
  endif

  int_mpi = 1
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                  option%mycomm,ierr)
  num_faces = temp_int
  pgrid%num_faces_global = num_faces

  ! divide faces across ranks
  allocate(pgrid%face_ids(num_faces_local))
  allocate(pgrid%face_cellids(num_faces_local))
  allocate(pgrid%face_nverts(num_faces_local))
  allocate(pgrid%face_vertids(max_nvert_per_cell,num_faces_local))
  allocate(pgrid%face_areas(num_faces_local))
  allocate(pgrid%face_centroids(num_faces_local))
  do iface = 1, num_faces_local
    pgrid%face_centroids(iface)%x = 0.0d0
    pgrid%face_centroids(iface)%y = 0.0d0
    pgrid%face_centroids(iface)%z = 0.0d0
  enddo

  ! read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  max_nvert_per_face = 0
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(7+max_nvert_per_cell,num_faces_local_save))
    do irank = 0, option%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_read = nfaces_per_proc(irank+1)
      do iface = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'face id',hint)
        temp_real_array(1,iface) = dble(temp_int)

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'cell id',hint)
        temp_real_array(2,iface) = dble(temp_int)

        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'number of vertices',hint)
        temp_real_array(3,iface) = dble(temp_int)
        if (temp_int > max_nvert_per_face) &
          max_nvert_per_face = temp_int

        do ivert = 1, int(temp_real_array(3,iface))
          call InputReadInt(input,option,temp_int)
          call InputErrorMsg(input,option,'face - vertex id',hint)
          temp_real_array(3+ivert,iface) = dble(temp_int)
        enddo

        call InputReadDouble(input,option,temp_real_array(4+max_nvert_per_cell,iface))
        call InputErrorMsg(input,option,'face x coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(5+max_nvert_per_cell,iface))
        call InputErrorMsg(input,option,'face y coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(6+max_nvert_per_cell,iface))
        call InputErrorMsg(input,option,'face z coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(7+max_nvert_per_cell,iface))
        call InputErrorMsg(input,option,'face area',hint)

      enddo

      if (irank == option%io_rank) then
        do iface = 1, num_faces_local
          pgrid%face_ids(iface) = int(temp_real_array(1,iface))
          pgrid%face_cellids(iface) = int(temp_real_array(2,iface))
          pgrid%face_nverts(iface) = int(temp_real_array(3,iface))
          pgrid%face_vertids(:,iface) = UNINITIALIZED_INTEGER
          do ivert = 1, pgrid%face_nverts(iface)
            pgrid%face_vertids(ivert,iface) = &
              int(temp_real_array(3+ivert,iface))
          enddo
          pgrid%face_centroids(iface)%x = &
            temp_real_array(4+max_nvert_per_cell,iface)
          pgrid%face_centroids(iface)%y = &
            temp_real_array(5+max_nvert_per_cell,iface)
          pgrid%face_centroids(iface)%z = &
            temp_real_array(6+max_nvert_per_cell,iface)
          pgrid%face_areas(iface) = temp_real_array(7+max_nvert_per_cell,iface)
        enddo
      else

        int_mpi = num_to_read*(7+max_nvert_per_cell)
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif

    enddo
  else
    ! other ranks post the recv
    allocate(temp_real_array(7+max_nvert_per_cell,num_faces_local+1))
    int_mpi = num_faces_local*(7+max_nvert_per_cell)
    call MPI_Recv(temp_real_array,int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)

    do iface = 1, num_faces_local
      pgrid%face_ids(iface) = int(temp_real_array(1,iface))
      pgrid%face_cellids(iface) = int(temp_real_array(2,iface))
      pgrid%face_nverts(iface) = int(temp_real_array(3,iface))
      pgrid%face_vertids(:,iface) = UNINITIALIZED_INTEGER
      do ivert = 1, pgrid%face_nverts(iface)
        pgrid%face_vertids(ivert,iface) = &
          int(temp_real_array(3+ivert,iface))
      enddo
      pgrid%face_centroids(iface)%x = &
        temp_real_array(4+max_nvert_per_cell,iface)
      pgrid%face_centroids(iface)%y = &
        temp_real_array(5+max_nvert_per_cell,iface)
      pgrid%face_centroids(iface)%z = &
        temp_real_array(6+max_nvert_per_cell,iface)
      pgrid%face_areas(iface) = temp_real_array(7+max_nvert_per_cell,iface)

    enddo
  endif
  deallocate(temp_real_array)

  if (option%myrank == option%io_rank) then
    call InputReadPflotranString(input,option)
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'VERTICES'
    call InputErrorMsg(input,option,word,card)
    if (.not.StringCompare(word,card)) then
      option%io_buffer = 'Unrecongnized keyword "' // trim(card) // &
        '" in polyhedra grid file.'
      call printErrMsgByRank(option)
    endif
    hint = 'Polyhedra Unstructured Grid VERTICES'
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of vertices',hint)
  endif

  int_mpi = 1
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                  option%mycomm,ierr)
  num_vertices = temp_int
  pgrid%num_vertices_global = num_vertices
  ugrid%num_vertices_global = num_vertices

   ! divide cells across ranks
  num_vertices_local = num_vertices/option%mycommsize 
  num_vertices_local_save = num_vertices_local
  remainder = num_vertices - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                 num_vertices_local + 1

  allocate(pgrid%vertex_coordinates(num_vertices_local))
  do ivert = 1, num_vertices_local
    pgrid%vertex_coordinates(ivert)%x = 0.d0
    pgrid%vertex_coordinates(ivert)%y = 0.d0
    pgrid%vertex_coordinates(ivert)%z = 0.d0
  enddo

  ! read all vertices from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(3,num_vertices_local_save+1))
    ! read for all processors
    do irank = 0, option%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_read = num_vertices_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do ivert = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadDouble(input,option,temp_real_array(1,ivert))
        call InputErrorMsg(input,option,'vertex x coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(2,ivert))
        call InputErrorMsg(input,option,'vertex y coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(3,ivert))
        call InputErrorMsg(input,option,'vertex z coordinate',hint)
      enddo

      if (irank == option%io_rank) then
        pgrid%num_vertices_local = num_vertices_local
        do ivert = 1, num_vertices_local
          pgrid%vertex_coordinates(ivert)%x = temp_real_array(1,ivert)
          pgrid%vertex_coordinates(ivert)%y = temp_real_array(2,ivert)
          pgrid%vertex_coordinates(ivert)%z = temp_real_array(3,ivert)
        enddo
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif

    enddo
  else

    pgrid%num_vertices_local = num_vertices_local
    allocate(temp_real_array(3,num_vertices_local))
    int_mpi = num_vertices_local*3
    call MPI_Recv(temp_real_array,int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    do ivert = 1, num_vertices_local
      pgrid%vertex_coordinates(ivert)%x = temp_real_array(1,ivert)
      pgrid%vertex_coordinates(ivert)%y = temp_real_array(2,ivert)
      pgrid%vertex_coordinates(ivert)%z = temp_real_array(3,ivert)
    enddo
  endif

  deallocate(temp_real_array)
  deallocate(nfaces_per_proc)

  temp_int = max_nface_per_cell
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                  option%mycomm,ierr)
  pgrid%max_nface_per_cell = temp_int

  temp_int = max_nvert_per_face
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                  option%mycomm,ierr)
  pgrid%max_nvert_per_face = temp_int

  if (option%myrank == option%io_rank) then
    call InputDestroy(input)
  endif

end subroutine UGridPolyhedraRead

! ************************************************************************** !

subroutine UGridPolyhedraDecompose(ugrid, option)
  ! 
  ! This routine decomposes a polyhedra grid across ranks.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 09/29/13
  ! 

#include "petsc/finclude/petscdm.h"
  use petscdm
  use Input_Aux_module
  use Option_module
  use String_module
  use Grid_Unstructured_Cell_module
  use Utility_module, only: reallocateIntArray, SearchOrderedArray

  implicit none

  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option

  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt :: icell
  PetscInt :: iface
  PetscInt :: ivertex
  PetscInt :: ivertex2
  PetscInt :: iface_beg
  PetscInt :: iface_end
  PetscInt :: local_id
  PetscInt :: min_vertex_id
  PetscInt :: index_format_flag
  PetscInt :: num_cells_local_old
  PetscInt :: global_offset_old
  PetscInt :: num_common_vertices
  PetscInt :: num_cells_local_new
  PetscInt :: count
  PetscInt :: max_nvert_per_cell
  PetscInt :: max_ndual_per_cell
  PetscInt :: max_nvert_per_face
  PetscInt :: max_nface_per_cell
  PetscInt :: vertex_ids_offset
  PetscInt :: temp_int
  PetscInt :: dual_offset
  PetscInt :: face_offset
  PetscInt :: natural_id_offset
  PetscInt :: cell_stride
  PetscInt :: face_count
  PetscInt :: max_int_count
  PetscInt :: vertex_count
  PetscInt :: vertex_id
  PetscInt :: ghosted_id
  PetscInt :: nface
  PetscInt :: iface_vert
  PetscInt :: num_faces_local
  PetscInt :: idx
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array1(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, allocatable :: vert_n2g(:,:)
  PetscInt, pointer :: int_array_pointer(:)
  PetscErrorCode :: ierr
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscBool :: success
  PetscInt :: num_rows, num_cols, istart, iend, icol
  character(len=MAXSTRINGLENGTH) :: string

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

  pgrid => ugrid%polyhedra_grid
  max_nvert_per_cell = ugrid%max_nvert_per_cell

  pgrid%cell_vertids = UNINITIALIZED_INTEGER

  min_vertex_id = 2 ! min value should be either 0 or 1 after global
                    ! reduction to idenitfy if vertex ids is 0-based
                    ! or 1-based.
  iface_beg = 1
  iface_end = 0
  ! find vertex ids forming a cell
  do icell = 1, pgrid%num_cells_local
    iface_end = iface_beg + pgrid%cell_nfaces(icell) - 1
    do iface = iface_beg, iface_end
      if (pgrid%face_cellids(iface) /= &
          pgrid%cell_ids(icell)) then
        option%io_buffer = 'Face ID does not correspond to cell'
        call printErrMsgByRank(option)
      endif
      do ivertex = 1, pgrid%face_nverts(iface)
        do ivertex2 = 1, max_nvert_per_cell
          if (pgrid%cell_vertids(ivertex2,icell) == UNINITIALIZED_INTEGER) then
              pgrid%cell_vertids(ivertex2,icell) = &
              pgrid%face_vertids(ivertex,iface)

              if (pgrid%cell_vertids(ivertex2,icell) < min_vertex_id) &
                min_vertex_id = pgrid%cell_vertids(ivertex2,icell)
            exit
          endif
          if (pgrid%cell_vertids(ivertex2,icell) == &
              pgrid%face_vertids(ivertex,iface)) exit
        enddo
      enddo

    enddo
    iface_beg = iface_end + 1
  enddo

  allocate(int_array1(max_nvert_per_cell))
  allocate(int_array2(max_nvert_per_cell))

  ! for a given cell, sort vertices in ascending order
  do icell = 1, pgrid%num_cells_local
    do ivertex = 1, pgrid%cell_nverts(icell)
      int_array1(ivertex) = pgrid%cell_vertids(ivertex,icell)
      int_array2(ivertex) = ivertex - 1
    enddo
    call PetscSortIntWithPermutation(pgrid%cell_nverts(icell), &
            int_array1,int_array2,ierr);CHKERRQ(ierr)
    int_array2 = int_array2 + 1
    do ivertex = 1, pgrid%cell_nverts(icell)
      pgrid%cell_vertids(ivertex,icell) = &
        int_array1(int_array2(ivertex))
    enddo
  enddo
  deallocate(int_array1)
  deallocate(int_array2)

  call MPI_Allreduce(min_vertex_id,index_format_flag, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN,option%mycomm,ierr)

  if (index_format_flag /= 0 .and. index_format_flag /= 1) then
    call printErrMsg(option,'Min. vertex id is neither 0 nor 1. Check input mesh.')
  endif

  num_cells_local_old = pgrid%num_cells_local
  ! let's make it Fortran indexing (i.e. 1-based)
  do local_id = 1, num_cells_local_old
    do ivertex = 1, max_nvert_per_cell
      ! at this point we may be zero-based
      if (pgrid%cell_vertids(ivertex,local_id) < 0) then
        ! change no_value (UNINITIALIZED_INTEGER) to '0'
        pgrid%cell_vertids(ivertex,local_id) = 0
      else
        if (index_format_flag == 0) then
          ! let's make it Fortran indexing
          pgrid%cell_vertids(ivertex,local_id) = &
            pgrid%cell_vertids(ivertex,local_id) + 1
        endif
      endif
    enddo
  enddo

#if UGRID_DEBUG
  write(string,*) ugrid%max_nvert_per_cell
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

  num_cells_local_old = pgrid%num_cells_local
  allocate(local_vertices(max_nvert_per_cell*num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  local_vertices = 0
  local_vertex_offset = 0
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, ugrid%max_nvert_per_cell
      if (pgrid%cell_vertids(ivertex,local_id) == 0) exit
      count = count + 1
      ! local vertices must be zero-based for MatCreateMPIAdj; thus subtract 1
      local_vertices(count) = &
        pgrid%cell_vertids(ivertex,local_id) - 1
    enddo
    local_vertex_offset(local_id+1) = count
  enddo

  select case (ugrid%grid_type)
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
                       pgrid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat, &
                       ierr);CHKERRQ(ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
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
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'Dual_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'Dual_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call MatView(Dual_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call UGridPartition(ugrid,option,Dual_mat,is_new, &
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


  if (.not.success .or. num_rows /= num_cells_local_old) then
    print *, option%myrank, num_rows, success, num_cells_local_old
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  ! calculate maximum number of connections for any given cell
  max_ndual_per_cell = 0
  do local_id = 1, num_cells_local_old
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > max_ndual_per_cell) &
      max_ndual_per_cell = num_cols
  enddo
  temp_int = max_ndual_per_cell
  call MPI_Allreduce(temp_int,max_ndual_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  ugrid%max_ndual_per_cell = max_ndual_per_cell

#if UGRID_DEBUG
  write(string,*) max_ndual_per_cell
  option%io_buffer = 'Maximum number of duals per cell: ' // adjustl(string)
  call printMsg(option)
#endif
  
  call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                          num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  
  ! in order to redistributed vertex/cell data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  max_nface_per_cell = pgrid%max_nface_per_cell
  max_nvert_per_face = pgrid%max_nvert_per_face

  face_offset = 8 + 1 ! +1 for -666
  vertex_ids_offset = face_offset + &
    ( &
      1                  + & ! natural face id
      1                  + & ! natural cell id
      1                  + & ! number of vertices
      max_nvert_per_face + & ! maximum vertices forming a face
      3                  + & ! face centroid (xc,yc,zc)
      1                    & ! face area
    ) * max_nface_per_cell + &
   1 ! for -777
  dual_offset = vertex_ids_offset + max_nvert_per_cell + 1 ! +1 for -888
  cell_stride = dual_offset+ max_ndual_per_cell + 1 ! +1 for -999999
  natural_id_offset = 2

  ! Information for each cell is packed in a strided petsc vec
  ! The information is ordered within each stride as follows:
  ! -cell_N   ! global cell id (negative indicates 1-based)
  ! natural cell id
  ! number of faces
  ! number of vertices
  ! cell x coordinate
  ! cell y coordinate
  ! cell z coordinate
  ! cell volume
  ! -666      ! separator between cell info and face info
  ! face1_id
  ! face1_cellid
  ! face1_nverts
  ! face1_vert1
  ! face1_vert2
  ! ...
  ! face1_vertN
  ! face1_xc
  ! face1_yc
  ! face1_zc
  ! face1_area
  ! face2
  ! face2_...
  ! ...
  ! face2_area
  ! ...
  ! ...
  ! faceM_area
  ! -777      ! separator between face info and vertices
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
  call UGridCreateOldVec(ugrid,option, elements_old, &
                         num_cells_local_old, is_new, is_scatter, cell_stride)
  
  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(elements_old,vec_ptr,ierr);CHKERRQ(ierr)
  count =  0
  face_count = 0
  do local_id = 1, num_cells_local_old
    count = count + 1
    vec_ptr(count) = -(global_offset_old+local_id)
    count = count + 1
    vec_ptr(count) = pgrid%cell_ids(local_id)
    count = count + 1
    vec_ptr(count) = pgrid%cell_nfaces(local_id)
    count = count + 1
    vec_ptr(count) = pgrid%cell_nverts(local_id)
    count = count + 1
    vec_ptr(count) = pgrid%cell_centroids(local_id)%x
    count = count + 1
    vec_ptr(count) = pgrid%cell_centroids(local_id)%y
    count = count + 1
    vec_ptr(count) = pgrid%cell_centroids(local_id)%z
    count = count + 1
    vec_ptr(count) = pgrid%cell_volumes(local_id)
    count = count + 1
    vec_ptr(count) = -666

    do iface = 1, pgrid%cell_nfaces(local_id)
      face_count = face_count + 1
      count = count + 1
      vec_ptr(count) = pgrid%face_ids(face_count)
      count = count + 1
      vec_ptr(count) = pgrid%face_cellids(face_count)
      count = count + 1
      vec_ptr(count) = pgrid%face_nverts(face_count)

      do ivertex = 1, pgrid%face_nverts(face_count)
        count = count + 1
        vec_ptr(count) = pgrid%face_vertids(ivertex,face_count)
      enddo

      do ivertex = pgrid%face_nverts(face_count)+1, max_nvert_per_face
        count = count + 1
        vec_ptr(count) = 0
      enddo

      count = count + 1
      vec_ptr(count) = pgrid%face_centroids(face_count)%x
      count = count + 1
      vec_ptr(count) = pgrid%face_centroids(face_count)%y
      count = count + 1
      vec_ptr(count) = pgrid%face_centroids(face_count)%z
      count = count + 1
      vec_ptr(count) = pgrid%face_areas(face_count)

    enddo

    do iface = pgrid%cell_nfaces(local_id)+1, max_nface_per_cell

      count = count + 1
      vec_ptr(count) = 0 ! face_id
      count = count + 1
      vec_ptr(count) = 0 ! face_cellid
      count = count + 1
      vec_ptr(count) = 0 ! face_nverts

      do ivertex = 1, max_nvert_per_face
        count = count + 1
        vec_ptr(count) = 0 ! face_vertids
      enddo

      count = count + 1
      vec_ptr(count) = 0 ! face_centoid_x
      count = count + 1
      vec_ptr(count) = 0 ! face_centoid_y
      count = count + 1
      vec_ptr(count) = 0 ! face_centoid_z
      count = count + 1
      vec_ptr(count) = 0 ! face_area

    enddo

    count = count + 1
    vec_ptr(count) = -777

    do ivertex = 1, max_nvert_per_cell
      count = count + 1
      if (ivertex <= pgrid%cell_nverts(local_id)) then
        vec_ptr(count) = pgrid%cell_vertids(ivertex,local_id)
      else
        vec_ptr(count) = 0
      endif
    enddo

    count = count + 1
    vec_ptr(count) = -888

    ! add the dual ids
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Dual matrix is larger then max_ndual_per_cell.'
      call printErrMsgByRank(option)
    endif
    do icol = 1, max_ndual_per_cell
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

  call MatRestoreRowIJF90(Dual_mat, ZERO_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                          num_rows, ia_ptr, ja_ptr, success,  &
                          ierr);CHKERRQ(ierr)
  call MatDestroy(Dual_mat,ierr);CHKERRQ(ierr)

  call UGridNaturalToPetsc(ugrid, option, &
                           elements_old, elements_local, &
                           num_cells_local_new, cell_stride, dual_offset, &
                           natural_id_offset, is_scatter)

  ! make a list of local vertices
  max_int_count = 2*ugrid%ngmax
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
  ! note that the vertices are still in natural numbering
  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id=1, ugrid%ngmax
    do ivertex = 1, ugrid%max_nvert_per_cell
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*cell_stride))
      if (vertex_id < 1) exit
      vertex_count = vertex_count + 1
      if (vertex_count > max_int_count) then
        call reallocateIntArray(int_array_pointer,max_int_count)
      endif
      vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*cell_stride) = vertex_count
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

  allocate(ugrid%vertex_ids_natural(vertex_count))
  ugrid%vertex_ids_natural = int_array3(1:vertex_count)

  ! now load all the vertices needed to define all the local cells
  ! on the processor
  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  ! allocate the array that will store the vertex ids for each cell.
  ! remember that max_nvert_per_cell is the max # of vertices in a cell
  ! currently hardwired to 8.
  !deallocate(ugrid%cell_vertices)
  allocate(ugrid%cell_vertices( &
             0:ugrid%max_nvert_per_cell,ugrid%ngmax))
  ugrid%cell_vertices = 0

  ! permute the local ids calculated earlier in the int_array4
  allocate(vert_n2g(ugrid%max_nvert_per_cell,2))

  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, ugrid%ngmax
     vert_n2g = 0
    do ivertex = 1, ugrid%max_nvert_per_cell
      ! extract the original vertex id
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*cell_stride))
      if (vertex_id < 1) exit
      count = ugrid%cell_vertices(0,ghosted_id)+1
      ugrid%cell_vertices(count,ghosted_id) = &
        int_array4(vertex_id)
      ugrid%cell_vertices(0,ghosted_id) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*cell_stride) = &
        int_array4(vertex_id)

      vert_n2g(ivertex,1) = int_array(vertex_id)
      vert_n2g(ivertex,2) = int_array4(vertex_id)
    enddo

    nface = int(vec_ptr((ghosted_id-1)*cell_stride+3))
    do iface = 1,nface
       do iface_vert = 1,max_nvert_per_face
          vertex_id = (ghosted_id-1)*cell_stride + &
                      face_offset + &
                      (iface-1)*(7 + max_nvert_per_face) + &
                      3 + iface_vert
          if (vec_ptr(vertex_id) < 1) exit
          do ivertex = 1, max_nvert_per_cell
            if (vec_ptr(vertex_id) == vert_n2g(ivertex,1)) then
              vec_ptr(vertex_id) = vert_n2g(ivertex,2)
              exit
            endif
          enddo
       enddo
    enddo

  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(int_array)
  deallocate(int_array2)
  deallocate(int_array3)
  deallocate(int_array4)

#if UGRID_DEBUG
  write(string,*) option%myrank
  if (ugrid%grid_type == THREE_DIM_GRID) then
    string = 'elements_vert_local' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'elements_vert_local' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(elements_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  deallocate(pgrid%cell_ids)
  deallocate(pgrid%cell_nfaces)
  deallocate(pgrid%cell_vertids)
  deallocate(pgrid%cell_nverts)
  deallocate(pgrid%cell_volumes)
  deallocate(pgrid%cell_centroids)

  deallocate(pgrid%face_ids)
  deallocate(pgrid%face_cellids)
  deallocate(pgrid%face_nverts)
  deallocate(pgrid%face_vertids)
  deallocate(pgrid%face_areas)
  deallocate(pgrid%face_centroids)

  allocate(pgrid%cell_ids(ugrid%ngmax))
  allocate(pgrid%cell_nfaces(ugrid%ngmax))
  allocate(pgrid%cell_vertids(max_nvert_per_cell,ugrid%ngmax))
  allocate(pgrid%cell_faceids(max_nface_per_cell,ugrid%ngmax))
  allocate(pgrid%cell_nverts(ugrid%ngmax))
  allocate(pgrid%cell_volumes(ugrid%ngmax))
  allocate(pgrid%cell_centroids(ugrid%ngmax))

  pgrid%cell_faceids = 0
  pgrid%cell_vertids = 0

  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)

  num_faces_local = 0
  do ghosted_id = 1, ugrid%ngmax
    idx = (ghosted_id-1)*cell_stride
    pgrid%cell_ids(ghosted_id) = int(vec_ptr(idx + 2))
    pgrid%cell_nfaces(ghosted_id) = int(vec_ptr(idx + 3))
    pgrid%cell_nverts(ghosted_id) = int(vec_ptr(idx + 4))
    pgrid%cell_centroids(ghosted_id)%x = vec_ptr(idx + 5)
    pgrid%cell_centroids(ghosted_id)%y = vec_ptr(idx + 6)
    pgrid%cell_centroids(ghosted_id)%z = vec_ptr(idx + 7)
    pgrid%cell_volumes(ghosted_id) = vec_ptr(idx + 8)

    pgrid%cell_vertids(1:max_nvert_per_cell,ghosted_id) = &
      ugrid%cell_vertices(1:max_nvert_per_cell,ghosted_id)

    do iface = 1,pgrid%cell_nfaces(ghosted_id)
      num_faces_local = num_faces_local + 1
      pgrid%cell_faceids(iface,ghosted_id) = num_faces_local
    enddo
  enddo

  allocate(pgrid%face_ids(num_faces_local))
  allocate(pgrid%face_cellids(num_faces_local))
  allocate(pgrid%face_nverts(num_faces_local))
  allocate(pgrid%face_vertids(max_nvert_per_face,num_faces_local))
  allocate(pgrid%face_areas(num_faces_local))
  allocate(pgrid%face_centroids(num_faces_local))
  pgrid%num_faces_local = num_faces_local
  pgrid%face_vertids = 0

  count = 0
  do ghosted_id = 1,ugrid%ngmax
    idx = (ghosted_id-1)*cell_stride + face_offset
    do iface = 1,pgrid%cell_nfaces(ghosted_id)
      count = count + 1
      pgrid%face_ids(count) = int(vec_ptr(idx + 1))
      !pgrid%face_cellids(count) = int(vec_ptr(idx + 2))
      pgrid%face_cellids(count) = ghosted_id
      pgrid%face_nverts(count) = int(vec_ptr(idx + 3))

      do ivertex = 1,pgrid%face_nverts(count)
        pgrid%face_vertids(ivertex,count) = int(vec_ptr(idx + 3 + ivertex))
      enddo

      idx = idx + 3 + max_nvert_per_face

      pgrid%face_centroids(count)%x = vec_ptr(idx + 1)
      pgrid%face_centroids(count)%y = vec_ptr(idx + 2)
      pgrid%face_centroids(count)%z = vec_ptr(idx + 3)
      pgrid%face_areas(count) = vec_ptr(idx + 4)

      idx = idx + 4
    enddo
  enddo

  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(elements_local,ierr);CHKERRQ(ierr)

  ! now we need to work on aligning the original vertex coordinates with 
  ! the current ordering or permuted/rearranged ordering.

  ugrid%num_vertices_local = pgrid%num_vertices_local

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
  call VecCreate(option%mycomm,vertices_old,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_old,ugrid%num_vertices_local*3, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_old,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_old,ierr);CHKERRQ(ierr)

  ! create serial petsc vector with a stride of 3
  call VecCreate(PETSC_COMM_SELF,vertices_new,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_new,vertex_count*3,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_new,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_new,ierr);CHKERRQ(ierr)

  call VecCreate(option%mycomm,vertices_old,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_old,ugrid%num_vertices_local*3, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_old,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_old,ierr);CHKERRQ(ierr)

  ! create serial petsc vector with a stride of 3
  call VecCreate(PETSC_COMM_SELF,vertices_new,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_new,vertex_count*3,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_new,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_new,ierr);CHKERRQ(ierr)

  ! load up the coordinates
  call VecGetArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  do ivertex = 1, pgrid%num_vertices_local
    vec_ptr((ivertex-1)*3+1) = pgrid%vertex_coordinates(ivertex)%x
    vec_ptr((ivertex-1)*3+2) = pgrid%vertex_coordinates(ivertex)%y
    vec_ptr((ivertex-1)*3+3) = pgrid%vertex_coordinates(ivertex)%z
  enddo
  call VecRestoreArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(pgrid%vertex_coordinates)
  nullify(pgrid%vertex_coordinates)

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
  ugrid%num_vertices_natural = ugrid%num_vertices_local
  ugrid%num_vertices_local = vertex_count
  allocate(ugrid%vertices(vertex_count))
  do ivertex = 1, vertex_count
    ugrid%vertices(ivertex)%x = 0.d0
    ugrid%vertices(ivertex)%y = 0.d0
    ugrid%vertices(ivertex)%z = 0.d0
  enddo
  pgrid%num_vertices_local = vertex_count
  allocate(pgrid%vertex_coordinates(vertex_count))
  do ivertex = 1, vertex_count
    pgrid%vertex_coordinates(ivertex)%x = 0.d0
    pgrid%vertex_coordinates(ivertex)%y = 0.d0
    pgrid%vertex_coordinates(ivertex)%z = 0.d0
  enddo

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  if (ugrid%grid_type == THREE_DIM_GRID) then
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
  if (ugrid%grid_type == THREE_DIM_GRID) then
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
  do ivertex = 1, ugrid%num_vertices_local
    ugrid%vertices(ivertex)%id = needed_vertices_petsc(ivertex)
    ugrid%vertices(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    ugrid%vertices(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    ugrid%vertices(ivertex)%z = vec_ptr((ivertex-1)*3+3)
    pgrid%vertex_coordinates(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    pgrid%vertex_coordinates(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    pgrid%vertex_coordinates(ivertex)%z = vec_ptr((ivertex-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_new,vec_ptr,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (ugrid%grid_type == THREE_DIM_GRID) then
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

  allocate(ugrid%cell_type(ugrid%ngmax))

  select case(ugrid%grid_type)
    case(THREE_DIM_GRID)
      do ghosted_id = 1, ugrid%ngmax
        ugrid%cell_type(ghosted_id) = POLY_TYPE
      enddo
    case default
      option%io_buffer = 'Grid type not recognized in UGridPolyhedraDecompose.'
      call printErrMsg(option)
  end select

end subroutine UGridPolyhedraDecompose

! ************************************************************************** !

subroutine UGridPolyhedraSetCellCentroids(pgrid,x,y,z, &
                                         x_min,x_max,y_min,y_max,z_min,z_max,option)
  ! 
  ! This routine set cell centroids for local+ghosted control volumes.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/28/13
  ! 

  use Option_module
  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscReal :: x(:), y(:), z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  type(option_type) :: option

  PetscInt :: icell
  PetscInt :: ivertex

  do icell = 1,size(pgrid%cell_centroids)
    x(icell) = pgrid%cell_centroids(icell)%x
    y(icell) = pgrid%cell_centroids(icell)%y
    z(icell) = pgrid%cell_centroids(icell)%z
  enddo

  do ivertex = 1, pgrid%num_vertices_local
    if (x_max < pgrid%vertex_coordinates(ivertex)%x) &
      x_max = pgrid%vertex_coordinates(ivertex)%x
    if (x_min > pgrid%vertex_coordinates(ivertex)%x) &
      x_min = pgrid%vertex_coordinates(ivertex)%x
    if (y_max < pgrid%vertex_coordinates(ivertex)%y) &
      y_max = pgrid%vertex_coordinates(ivertex)%y
    if (y_min > pgrid%vertex_coordinates(ivertex)%y) &
      y_min = pgrid%vertex_coordinates(ivertex)%y
    if (z_max < pgrid%vertex_coordinates(ivertex)%z) &
      z_max = pgrid%vertex_coordinates(ivertex)%z
    if (z_min > pgrid%vertex_coordinates(ivertex)%z) &
      z_min = pgrid%vertex_coordinates(ivertex)%z
  enddo

end subroutine UGridPolyhedraSetCellCentroids

! ************************************************************************** !

function UGridPolyhedraComputeInternConnect(ugrid, grid_x, &
                                             grid_y, grid_z, option)
  ! 
  ! This routine compute internal connectivity of an unstrucutred polyhedra
  ! grid.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/28/13
  ! 

  use Connection_module
  use Option_module
  use Grid_Unstructured_module
  use Utility_module, only : DotProduct, CrossProduct
  use Geometry_module

  implicit none

  type(connection_set_type), pointer :: UGridPolyhedraComputeInternConnect
  type(grid_unstructured_type) :: ugrid 
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(option_type) :: option

  type(connection_set_type), pointer :: connections
  type(unstructured_polyhedra_type), pointer :: pgrid

  PetscInt :: nconn, iconn
  PetscInt :: idual, dual_id

  PetscInt, allocatable :: face_to_vertex(:,:)
  PetscInt, allocatable :: cell_to_face(:,:)
  PetscInt, allocatable :: face_to_cell(:,:)
  PetscInt, allocatable :: vertex_to_cell(:,:)
  PetscInt, allocatable :: temp_int(:)
  PetscInt, allocatable :: temp_int_2d(:,:)
  PetscInt, allocatable :: face_nverts(:)
  PetscInt, allocatable :: vertex_ids(:)
  PetscInt, allocatable :: dup_face_id(:)
  PetscBool, allocatable :: local_boundary_face(:)

  PetscInt :: max_face_per_cell
  PetscInt :: max_vert_per_cell
  PetscInt :: max_vert_per_face

  PetscInt :: num_match
  PetscInt :: found_count
  PetscInt :: face_count
  PetscInt :: count
  PetscInt :: iside
  PetscInt :: icell
  PetscInt :: dual_local_id
  PetscInt :: cell_nvert

  PetscInt :: iface,      iface2
  PetscInt :: ivertex,    ivertex2
  PetscInt :: face_id,    face_id2
  PetscInt :: ghosted_id, ghosted_id2
  PetscInt :: local_id,   local_id2
  PetscInt :: cell_id,    cell_id2
  PetscInt :: vertex_id,  vertex_id2
  PetscInt :: cell_type,  cell_type2
  PetscInt :: face_type,  face_type2
  PetscInt :: nfaces,     nfaces2
  PetscInt :: nvertices,  nvertices2
  
  PetscInt :: nintercp
  PetscInt :: iintercp
  PetscInt :: idx
  PetscBool :: face_found
  PetscBool :: vertex_found
  PetscBool :: cell_found
    
  PetscReal :: v1(3), v2(3), v3(3)
  PetscReal :: n1(3), n2(3), n_up_dn(3)
  PetscReal :: dist_up, dist_dn

  type(plane_type) :: plane1
  type(point3d_type) :: point1, point2, point3
  type(point3d_type) :: point_up, point_dn
  type(point3d_type) :: intercept1, intercept2, intercept

  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string  

  pgrid => ugrid%polyhedra_grid

  max_face_per_cell = pgrid%max_nface_per_cell
  max_vert_per_face = pgrid%max_nvert_per_face
  max_vert_per_cell = ugrid%max_nvert_per_cell

  allocate(vertex_ids(max_vert_per_face))

  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(max_vert_per_face, &
           max_face_per_cell*ugrid%ngmax))
  face_to_vertex = 0

  allocate(cell_to_face(max_face_per_cell, &
                        ugrid%ngmax))
  cell_to_face = 0
  allocate(face_to_cell(2,max_face_per_cell* &
                        ugrid%ngmax))
  face_to_cell = 0
  allocate(vertex_to_cell(0:ugrid%max_cells_sharing_a_vertex, &
                          ugrid%num_vertices_local))
  vertex_to_cell = 0

  allocate(ugrid%face_to_vertex_natural(max_vert_per_face, &
           max_face_per_cell*ugrid%ngmax))
  ugrid%face_to_vertex_natural = 0

  allocate(dup_face_id(max_face_per_cell*ugrid%ngmax))
  dup_face_id = 0

  face_count = 0
  do ghosted_id = 1, ugrid%ngmax
     nfaces = pgrid%cell_nfaces(ghosted_id)
     do iface = 1, nfaces
       face_count = face_count + 1
       cell_to_face(iface,ghosted_id) = face_count
       face_to_cell(1,face_count) = ghosted_id
       nvertices = pgrid%face_nverts(face_count)
       vertex_ids = pgrid%face_vertids(:,face_count)
       do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = vertex_ids(ivertex)
          if (face_to_vertex(ivertex,face_count) > 0) then
            ugrid%face_to_vertex_natural(ivertex,face_count) = &
              ugrid%vertex_ids_natural(face_to_vertex(ivertex,face_count))
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
  iconn = 0
  do local_id = 1, ugrid%nlmax
    ! Select a cell and find number of vertices
    cell_id = local_id
    nfaces = pgrid%cell_nfaces(local_id)
    do idual = 1, ugrid%cell_neighbors_local_ghosted(0,local_id)
      cell_id2 = &
        abs(ugrid%cell_neighbors_local_ghosted(idual,local_id))
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      nfaces2 = pgrid%cell_nfaces(cell_id2)
      face_found = PETSC_FALSE

      do iface = 1, nfaces
        face_id = cell_to_face(iface,cell_id)
        nvertices = pgrid%face_nverts(face_id)

        do ivertex = 1, nvertices
          vertex_id = face_to_vertex(ivertex,face_id)
          vertex_found = PETSC_FALSE
          do ivertex2 = 1, ugrid%cell_vertices(0,cell_id2)
             vertex_id2 = ugrid%cell_vertices(ivertex2,cell_id2)
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
        enddo ! do-loop 'ivertex'

        if (vertex_found) then
          ! All the vertices of iface are present in the Neighboring cells.
          ! Thus, iface is the shared face.
          face_found = PETSC_TRUE

          ! Now, we have to find iface2 that corresponds to iface
          do iface2 = 1, nfaces2
            face_id2 = cell_to_face(iface2,cell_id2)
            nvertices2 = pgrid%face_nverts(face_id2)
            ! iface and iface2 need to have same number of vertices
            if (nvertices == nvertices2) then
              ! Count the number of vertices of iface which match vertices of iface2
              num_match = 0
              do ivertex = 1, nvertices
                vertex_id = face_to_vertex(ivertex,face_id)
                vertex_found = PETSC_FALSE
                do ivertex2 = 1, nvertices2
                  vertex_id2 = face_to_vertex(ivertex2,face_id2)
                  if (vertex_id == vertex_id2) then
                    vertex_found = PETSC_TRUE
                    num_match = num_match + 1
                    exit
                  endif
                enddo
                !
                ! If vertex_id of face_id not found as one of the vertices of face_id2,
                ! face_id2 is not shared between cells
                if (.not.vertex_found) exit
              enddo

              if (num_match == nvertices) then
                ! remove duplicate face
                if (face_id2 > face_id) then
#ifdef UGRID_DEBUG
                  write(string,*) option%myrank, face_id2, ' -> ', face_id
                  option%io_buffer = 'Duplicated face removed:' // trim(string)
                  call printMsg(option)
#endif
                  cell_to_face(iface2,cell_id2) = face_id
                  ! flag face_id2 as removed
                  face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
                  face_to_cell(2,face_id2) = cell_id
                  ! add cell_id2 to face_ids list
                  face_to_cell(2,face_id) = cell_id2
                  dup_face_id(face_id) = face_id2
                else
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id, ' -> ', face_id2
                  option%io_buffer = 'Duplicated face removed:' // trim(string)
                  call printMsg(option)
#endif
                  cell_to_face(iface,cell_id) = face_id2
                  ! flag face_id as removed  
                  face_to_cell(1,face_id) = -face_to_cell(1,face_id)
                  face_to_cell(2,face_id) = cell_id2
                  ! add cell_id to face_ids2 list
                  face_to_cell(2,face_id2) = cell_id
                  dup_face_id(face_id2)= face_id
                endif ! if (face_id2 > face_id)
                exit
              endif ! if (num_match == nvertices)
            endif ! if (nvertices == nvertices2)

            ! Check that one shared face was found between the Cell and Neighboring-Cell
            if (.not.face_found) then
               write(string,*) option%myrank
               string = '(' // trim(adjustl(string)) // ')'
               write(*,'(a,'' local_id = '',i3,'' natural_id = '',i3,''  vertices: '',8i3)') &
                     trim(string), &
                     cell_id,ugrid%cell_ids_natural(cell_id), &
                     (ugrid%vertex_ids_natural( &
                     ugrid%cell_vertices(ivertex,cell_id)), &
                     ivertex=1,ugrid%cell_vertices(0,cell_id))
               write(*,'(a,'' local_id2 = '',i3,'' natural_id2 = '',i3,''  vertices2: '',8i3)') &
                     trim(string), &
                     cell_id2,ugrid%cell_ids_natural(cell_id2), &
                     (ugrid%vertex_ids_natural( &
                     ugrid%cell_vertices(ivertex2,cell_id2)), &
                     ivertex2=1,ugrid%cell_vertices(0,cell_id2))
               option%io_buffer='No shared face found.'
               call printErrMsgByRank(option)
            endif

          enddo ! do-loop iface2

          exit
        endif ! if (vertex_found)
      enddo ! do-loop 'iface'

      ! Check that one shared face was found between the Cell and Neighboring-Cell
      if (.not.face_found) then
        write(string,*) option%myrank
        string = '(' // trim(adjustl(string)) // ')'
        write(*,'(a,'' local_id = '',i3,'' natural_id = '',i3,''  vertices: '',8i3)') &
                   trim(string), &
                   cell_id,ugrid%cell_ids_natural(cell_id), &
                   (ugrid%vertex_ids_natural( &
                     ugrid%cell_vertices(ivertex,cell_id)), &
                     ivertex=1,ugrid%cell_vertices(0,cell_id))
        write(*,'(a,'' local_id2 = '',i3,'' natural_id2 = '',i3,''  vertices2: '',8i3)') &
                   trim(string), &
                   cell_id2,ugrid%cell_ids_natural(cell_id2), &
                   (ugrid%vertex_ids_natural( &
                     ugrid%cell_vertices(ivertex2,cell_id2)), &
                     ivertex2=1,ugrid%cell_vertices(0,cell_id2))
        option%io_buffer='No shared face found.'
        call printErrMsgByRank(option)
      endif

    enddo ! do-loop 'idual'
  enddo ! do-loop 'local_id'

  ! GB: Add a check for dup_face_id

  ! count up the # of faces
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) &
      face_count = face_count + 1
  enddo
  allocate(ugrid%face_to_vertex(max_vert_per_face,face_count))
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      ugrid%face_to_vertex(:,face_count) = face_to_vertex(:,iface)
    endif
  enddo
  deallocate(face_to_vertex)
  ! reallocate face_to_cell to proper size
  allocate(temp_int_2d(2,face_count))
  allocate(face_nverts(face_count))
  allocate(temp_int(size(face_to_cell,2)))
  temp_int = 0
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      temp_int_2d(:,face_count) = face_to_cell(:,iface)
      face_nverts(face_count) = pgrid%face_nverts(iface)
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
    do icell = 1,2
      cell_id = face_to_cell(icell,face_id)
      ! check for exterior face
      if (cell_id < 1) cycle
      cell_found = PETSC_FALSE
      nfaces = pgrid%cell_nfaces(cell_id)
      do iface2 = 1, nfaces
        face_id2 = cell_to_face(iface2,cell_id)
        if (face_id < 0) cycle
        if (face_id == temp_int(face_id2)) then
          cell_found = PETSC_TRUE
          cell_to_face(iface2,cell_id) = face_id
          exit
        endif
      enddo

      if (.not.cell_found) then
        option%io_buffer = 'POLYHEDRA_UGRID: Remapping of cell face id unsuccessful'
        call printErrMsg(option)
      endif
    enddo
  enddo
  deallocate(temp_int)

  do ghosted_id = 1, ugrid%ngmax
    do ivertex = 1, ugrid%cell_vertices(0,ghosted_id)
      vertex_id = ugrid%cell_vertices(ivertex,ghosted_id)
      if ( vertex_id <= 0) cycle 
      count = vertex_to_cell(0,vertex_id) + 1
      if (count > ugrid%max_cells_sharing_a_vertex) then
        write(string,*) 'Vertex can be shared by at most by ', &
              ugrid%max_cells_sharing_a_vertex, &
              ' cells. Rank = ', option%myrank, ' vertex_id = ', vertex_id, ' exceeds it.'
        option%io_buffer = string
        call printErrMsg(option)
      endif
      vertex_to_cell(count,vertex_id) = ghosted_id
      vertex_to_cell(0,vertex_id) = count
    enddo
  enddo

  nconn = 0
  do local_id = 1, ugrid%nlmax
    do idual = 1, ugrid%cell_neighbors_local_ghosted(0,local_id)
      dual_id = ugrid%cell_neighbors_local_ghosted(idual,local_id)
      ! count all ghosted connections (dual_id < 0)
      ! only count connection with cells of larger ids to avoid double counts
      if (dual_id < 0 .or. local_id < dual_id) then
        nconn = nconn + 1
      endif
    enddo
  enddo

  connections => ConnectionCreate(nconn,INTERNAL_CONNECTION_TYPE)

  allocate(ugrid%connection_to_face(nconn))
  ugrid%connection_to_face = 0

  ! loop over connection again
  iconn = 0
  do local_id = 1, ugrid%nlmax
     do idual = 1, ugrid%cell_neighbors_local_ghosted(0,local_id)
       dual_local_id = &
         ugrid%cell_neighbors_local_ghosted(idual,local_id)
       if (local_id < abs(dual_local_id)) then
        iconn = iconn + 1
        ! find face
        face_found = PETSC_FALSE
        do iface = 1, ugrid%cell_vertices(0,local_id)
          face_id = cell_to_face(iface,local_id)
          do iside = 1,2
            cell_id2 = face_to_cell(iside,face_id)
            if (cell_id2 == abs(dual_local_id)) then
              face_found = PETSC_TRUE
              exit
            endif
          enddo
          if (face_found) exit
        enddo
        if (face_found) then
          ugrid%connection_to_face(iconn) = face_id
        else
          write(string,*) option%myrank,local_id,dual_local_id 
          option%io_buffer = 'face not found in connection loop' // trim(string)
          call printErrMsg(option)
        endif

        face_found = PETSC_FALSE
        do iface2 = 1, ugrid%cell_vertices(0,cell_id2)
          if (cell_to_face(iface,local_id) == &
              cell_to_face(iface2,cell_id2)) then
            face_found = PETSC_TRUE
            exit
          endif
        enddo

        if (.not.face_found) then
          write(string,*) option%myrank, iface, cell_id2
          option%io_buffer = 'global face not found' // trim(string)
          call printErrMsg(option)
        endif

        connections%id_up(iconn) = local_id
        connections%id_dn(iconn) = abs(dual_local_id)
        connections%face_id(iconn) = cell_to_face(iface,local_id)

        point_up%x = grid_x(local_id)
        point_up%y = grid_y(local_id)
        point_up%z = grid_z(local_id)
        point_dn%x = grid_x(abs(dual_local_id))
        point_dn%y = grid_y(abs(dual_local_id))
        point_dn%z = grid_z(abs(dual_local_id))

        ! Find intercept
        nintercp = face_nverts(cell_to_face(iface,local_id)) - 2

        intercept%x = 0.d0
        intercept%y = 0.d0
        intercept%z = 0.d0

        do iintercp = 0, nintercp - 1 
          idx = ugrid%face_to_vertex(1 + iintercp,face_id)
          point1 = ugrid%vertices(idx)
          idx = ugrid%face_to_vertex(2 + iintercp,face_id)
          point2 = ugrid%vertices(idx)
          idx = ugrid%face_to_vertex(3 + iintercp,face_id)
          point3 = ugrid%vertices(idx)

          call GeometryComputePlaneWithPoints(plane1,point1,point2,point3)
          call GeometryGetPlaneIntercept(plane1,point_up,point_dn,intercept1)

          intercept%x = intercept%x + intercept1%x
          intercept%y = intercept%y + intercept1%y
          intercept%z = intercept%z + intercept1%z

        enddo

        intercept%x = intercept%x/nintercp
        intercept%y = intercept%y/nintercp
        intercept%z = intercept%z/nintercp

        ! This is very crude, but for now use average location of intercept
        v1(1) = intercept%x-point_up%x
        v1(2) = intercept%y-point_up%y
        v1(3) = intercept%z-point_up%z
        v2(1) = point_dn%x-intercept%x
        v2(2) = point_dn%y-intercept%y
        v2(3) = point_dn%z-intercept%z
        dist_up = sqrt(DotProduct(v1,v1))
        dist_dn = sqrt(DotProduct(v2,v2))

        connections%dist(-1:3,iconn) = 0.d0
        connections%dist(-1,iconn) = dist_up/(dist_up + dist_dn)
        connections%dist(0,iconn) = dist_up + dist_dn
        v3 = v1 + v2
        connections%dist(1:3,iconn) = v3/sqrt(DotProduct(v3,v3))
        connections%area(iconn) = pgrid%face_areas(connections%face_id(iconn))
        connections%intercp(1,iconn) = intercept%x
        connections%intercp(2,iconn) = intercept%y
        connections%intercp(3,iconn) = intercept%z

       endif ! (local_id < abs(dual_local_id))

     enddo
  enddo ! 'local_id'

  allocate(ugrid%face_area(face_count))
  allocate(ugrid%face_centroid(face_count))

  do local_id = 1, ugrid%nlmax
    do iface = 1, pgrid%cell_nfaces(local_id)
      face_id = cell_to_face(iface,local_id)

      ugrid%face_centroid(face_id)%x = &
        pgrid%face_centroids(face_id)%x
      ugrid%face_centroid(face_id)%y = &
        pgrid%face_centroids(face_id)%y
      ugrid%face_centroid(face_id)%z = &
        pgrid%face_centroids(face_id)%z

      ugrid%face_area(face_id) = &
        pgrid%face_areas(face_id)

    enddo
  enddo

  allocate(ugrid%face_to_cell_ghosted(size(face_to_cell,1), &
                                                  size(face_to_cell,2)))
  ugrid%face_to_cell_ghosted = face_to_cell
  allocate(ugrid%cell_to_face_ghosted(size(cell_to_face,1), &
                                                  size(cell_to_face,2)))
  ugrid%cell_to_face_ghosted(:,:) = cell_to_face(:,:)


#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'face_to_cell' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do iface = 1, size(face_to_cell,2)
    write(86,'(i5)',advance='no') iface
    write(86,'(i5)',advance='no') face_to_cell(1,iface)
    write(86,'(i5)',advance='no') face_to_cell(2,iface)
    write(86,'(a)') ""
  enddo
  close(86)

  write(string,*) option%myrank
  string = 'cell_to_face' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, ugrid%ngmax
    write(86,'(i5)',advance='no') ghosted_id
    do iface = 1,max_face_per_cell
      write(86,'(i5)',advance='no') cell_to_face(iface,ghosted_id)
    enddo
    write(86,'(a)') ""
  enddo
  close(86)

  write(string,*) option%myrank
  string = 'poly_cell_vertids' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, ugrid%ngmax
    write(86,'(i5)',advance='no') ghosted_id
    do ivertex = 1,ugrid%max_nvert_per_cell
      write(86,'(i5)',advance='no') pgrid%cell_vertids(ivertex,ghosted_id)
    enddo
    write(86,'(a)') ""
  enddo
  close(86)

  write(string,*) option%myrank
  string = 'poly_cell_faceids' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, ugrid%ngmax
    write(86,'(i5)',advance='no') ghosted_id
    do ivertex = 1,pgrid%max_nface_per_cell
      write(86,'(i5)',advance='no') pgrid%cell_faceids(ivertex,ghosted_id)
    enddo
    write(86,'(a)') ""
  enddo
  close(86)

  write(string,*) option%myrank
  string = 'poly_face_vertids' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do iface = 1, pgrid%num_faces_local
    write(86,'(i5)',advance='no') iface
    do ivertex = 1,pgrid%max_nvert_per_face
      write(86,'(i5)',advance='no') pgrid%face_vertids(ivertex,iface)
    enddo
    write(86,'(a)') ""
  enddo
  close(86)

#endif

  deallocate(cell_to_face)
  deallocate(face_to_cell)
  deallocate(vertex_to_cell)
  deallocate(vertex_ids)
  deallocate(dup_face_id)

  UGridPolyhedraComputeInternConnect => connections

end function UGridPolyhedraComputeInternConnect

! ************************************************************************** !

subroutine UGridPolyhedraComputeVolumes(ugrid, option, volume)
  ! 
  ! This routine sets volumes of local control volumes.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/28/13
  ! 

  use Option_module

  implicit none
  
  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  Vec :: volume

  type(unstructured_polyhedra_type), pointer :: pgrid
  
  PetscInt :: icell
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid

  call VecGetArrayF90(volume,vec_ptr,ierr);CHKERRQ(ierr)
  do icell = 1, ugrid%nlmax
    vec_ptr(icell) = pgrid%cell_volumes(icell)
  enddo
  call VecRestoreArrayF90(volume,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine UGridPolyhedraComputeVolumes

! ************************************************************************** !

subroutine UGridPolyhedraPopulateConnection(ugrid, connection, iface_cell, &
                                             iconn, ghosted_id, option)
  ! 
  ! This routine computes details about boundary connections (area, dist, etc)
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/28/13
  ! 

  use Connection_module
  use Utility_module, only : DotProduct
  use Option_module
  use Grid_Unstructured_Cell_module
  use Geometry_module  
  
  implicit none
  
  type(grid_unstructured_type) :: ugrid
  type(connection_set_type) :: connection
  PetscInt :: iface_cell
  PetscInt :: iconn
  PetscInt :: ghosted_id
  type(option_type) :: option
  
  PetscInt :: face_id
  PetscInt :: ivert,vert_id
  PetscInt :: face_type
  PetscReal :: v1(3),v2(3),n_dist(3), dist
  type(point3d_type) :: vertex_8(8)
  type(plane_type) :: plane
  type(point3d_type) :: point, vertex1, vertex2, vertex3, intercept
  type(unstructured_polyhedra_type), pointer :: pgrid
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid
  
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
      v2(1) = pgrid%cell_centroids(ghosted_id)%x
      v2(2) = pgrid%cell_centroids(ghosted_id)%y
      v2(3) = pgrid%cell_centroids(ghosted_id)%z

      ! Instead of connecting centroid with face center, calculate the shortest
      ! distance between the centroid and face and use that distance - geh

      point%x = v2(1)
      point%y = v2(2)
      point%z = v2(3)

      !face_id = ugrid%cell_to_face_ghosted(iface_cell, ghosted_id)
      face_id = pgrid%cell_faceids(iface_cell, ghosted_id)
      intercept%x = pgrid%face_centroids(face_id)%x
      intercept%y = pgrid%face_centroids(face_id)%y
      intercept%z = pgrid%face_centroids(face_id)%z

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
      connection%area(iconn)    = pgrid%face_areas(face_id)
      connection%intercp(1,iconn)= intercept%x
      connection%intercp(2,iconn)= intercept%y
      connection%intercp(3,iconn)= intercept%z
      connection%face_id(iconn)  = face_id
      
  end select

end subroutine UGridPolyhedraPopulateConnection

! ************************************************************************** !

subroutine UGridPolyhedraGetCellsInRectangle(x_min, x_max, y_min, y_max, z_min, z_max, &
                                              ugrid, option, num_cells, &
                                              cell_ids, cell_face_ids)
  ! 
  ! This routine returns cells that are within a cube.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/28/13
  ! 
  use Option_module
  use Utility_module, only : reallocateIntArray
  use Geometry_module
  
  implicit none
                  
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  PetscInt :: num_cells
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: cell_face_ids(:)
  
  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt :: cell_type, num_faces, iface, face_type
  PetscInt :: vertex_id
  PetscInt :: num_vertices, ivertex
  PetscInt :: local_id, ghosted_id
  type(point3d_type) :: point
  
  PetscReal :: x_min_adj, x_max_adj, y_min_adj, y_max_adj, z_min_adj, z_max_adj
  PetscReal :: pert
  PetscBool :: in_rectangle
  
  PetscInt, pointer :: temp_cell_array(:), temp_face_array(:)
  PetscInt :: temp_array_size
  PetscInt :: face_id
  
  pgrid => ugrid%polyhedra_grid

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
  
  do local_id = 1, ugrid%nlmax
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    num_faces = pgrid%cell_nfaces(ghosted_id)
    do iface = 1, num_faces
      face_id = ugrid%cell_to_face_ghosted(iface, ghosted_id)
      num_vertices = pgrid%face_nverts(ghosted_id)
      in_rectangle = PETSC_TRUE
      do ivertex = 1, num_vertices
        !vertex_id = pgrid%face_vertids(ivertex,face_id)
        vertex_id = ugrid%face_to_vertex(ivertex,face_id)
        point = ugrid%vertices(vertex_id)
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

end subroutine UGridPolyhedraGetCellsInRectangle

! ************************************************************************** !

subroutine UGridPolyhedraComputeOutputInfo(ugrid, nL2G, nG2L, nG2A, option)
  ! 
  ! This routine computes informations later required to write tecplot output
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 12/29/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module

  implicit none

  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nG2A(:)

  type(unstructured_polyhedra_type), pointer :: pgrid

  Vec :: nat_cv_proc_rank
  Vec :: ghosted_cv_proc_rank
  VecScatter :: vec_scat
  IS :: is_scatter
  IS :: is_gather

  Vec :: ghosted_vec
  Vec :: natural_vec

  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: iface
  PetscInt :: ivertex
  PetscInt :: count

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: dual_id
  PetscInt :: face_id
  PetscInt :: vertex_id

  PetscInt :: pgridface_id

  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: real_array(:)
  PetscScalar,pointer :: v_loc_p(:)

  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid

  pgrid%ugrid_num_faces_local = size(ugrid%face_to_cell_ghosted,2)
  allocate(pgrid%ugridf2pgridf(pgrid%ugrid_num_faces_local))
  pgrid%ugridf2pgridf = 0

  ! Initialize mapping of faces between unstructured_grid and polyhedra_grid
  do ghosted_id = 1, ugrid%ngmax
    do iface = 1, pgrid%cell_nfaces(ghosted_id)
      face_id = ugrid%cell_to_face_ghosted(iface, ghosted_id)

      if (pgrid%ugridf2pgridf(face_id) == 0) then ! mapping not initialized
        pgrid%ugridf2pgridf(face_id) = pgrid%cell_faceids(iface, ghosted_id)
      else
        ! iface-th of ghosted_id-cell is an internal face shared by:
        ! - two local control volumes, or
        ! - a local and ghosted control volume

        ! Determine if ghosted_id is the upwind control volume.
        if (ugrid%face_to_cell_ghosted(1, face_id) == ghosted_id) then
          ! Map ugrid face to corresponding face of upwind control volume
          ! in pgrid. Upwind control volume is choosen because the order of
          ! vertices forming a face is required for output. Also, unit normal
          ! vector point from upwind-to-downwind control volume.
          pgrid%ugridf2pgridf(face_id) = pgrid%cell_faceids(iface, ghosted_id)
        endif
      endif
    enddo
  enddo

  ! Find number of global unique faces. This is required for output
  call VecCreateMPI(option%mycomm, ugrid%nlmax, PETSC_DETERMINE, nat_cv_proc_rank,  &
                    ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm, ugrid%ngmax, PETSC_DETERMINE, ghosted_cv_proc_rank,  &
                    ierr);CHKERRQ(ierr)

  ! Populate a vector that contains rank of procoessor on which a given
  ! control volume is active.
  allocate(int_array(ugrid%nlmax))
  allocate(real_array(ugrid%nlmax))
  do local_id = 1,ugrid%nlmax
    int_array(local_id) = nG2A(nL2G(local_id))
    real_array(local_id) = option%myrank
  enddo
  int_array = int_array - 1

  call VecSetValues(nat_cv_proc_rank, ugrid%nlmax, int_array, real_array, INSERT_VALUES, &
                    ierr);CHKERRQ(ierr)
  call VecAssemblyBegin(nat_cv_proc_rank, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(nat_cv_proc_rank, ierr);CHKERRQ(ierr)
  deallocate(int_array)
  deallocate(real_array)

  ! Find processor rank for ghost control volumes by scattering data stored in
  ! vector nat_cv_proc_rank.
  allocate(int_array(ugrid%ngmax))

  do ghosted_id = 1, ugrid%ngmax
    int_array(ghosted_id) = nG2A(ghosted_id)
  enddo
  int_array = int_array - 1
  call ISCreateBlock(option%mycomm, 1, ugrid%ngmax, int_array, PETSC_COPY_VALUES, &
                      is_scatter, ierr);CHKERRQ(ierr)

  call VecGetOwnershipRange(ghosted_cv_proc_rank, istart, iend,  &
                            ierr);CHKERRQ(ierr)
  do ghosted_id = 1, ugrid%ngmax
    int_array(ghosted_id) = ghosted_id + istart
  enddo
  int_array = int_array - 1
  call ISCreateBlock(option%mycomm, 1, ugrid%ngmax, int_array, PETSC_COPY_VALUES, &
                      is_gather, ierr);CHKERRQ(ierr)

  call VecScatterCreate(nat_cv_proc_rank, is_scatter, ghosted_cv_proc_rank, is_gather, &
                        vec_scat, ierr);CHKERRQ(ierr)
  call ISDestroy(is_scatter, ierr);CHKERRQ(ierr)
  call ISDestroy(is_gather, ierr);CHKERRQ(ierr)
  deallocate(int_array)

  call VecScatterBegin(vec_scat, nat_cv_proc_rank, ghosted_cv_proc_rank, &
                        INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scat, nat_cv_proc_rank, ghosted_cv_proc_rank, &
                      INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scat, ierr);CHKERRQ(ierr)

  ! Find the number of unique faces
  allocate(pgrid%uface_localids(pgrid%ugrid_num_faces_local))
  allocate(pgrid%uface_nverts(pgrid%ugrid_num_faces_local))
  allocate(pgrid%uface_left_natcellids(pgrid%ugrid_num_faces_local))
  allocate(pgrid%uface_right_natcellids(pgrid%ugrid_num_faces_local))

  pgrid%uface_localids = -1
  pgrid%uface_nverts = -1
  pgrid%uface_left_natcellids = -1
  pgrid%uface_right_natcellids = -1

  call VecGetArrayF90(ghosted_cv_proc_rank, v_loc_p, ierr);CHKERRQ(ierr)
  pgrid%num_ufaces_local = 0
  pgrid%uface_nverts = 0
  do iface = 1, pgrid%ugrid_num_faces_local
    ghosted_id = ugrid%face_to_cell_ghosted(1, iface)
    dual_id = ugrid%face_to_cell_ghosted(2, iface)
    local_id = nG2L(ghosted_id)

    if (ghosted_id > ugrid%nlmax) cycle

    if (dual_id == 0) then
      pgrid%num_ufaces_local = pgrid%num_ufaces_local + 1
      pgridface_id = pgrid%ugridf2pgridf(iface)
      pgrid%uface_localids(pgrid%num_ufaces_local) = pgridface_id
      pgrid%uface_nverts(pgrid%num_ufaces_local) = pgrid%face_nverts(pgridface_id)
      pgrid%uface_left_natcellids(pgrid%num_ufaces_local) = nG2A(ghosted_id)
      pgrid%uface_right_natcellids(pgrid%num_ufaces_local) = 0
      pgrid%num_verts_of_ufaces_local = pgrid%num_verts_of_ufaces_local + pgrid%face_nverts(pgridface_id)
    else
      if (nG2L(dual_id) == 0) then
        if (v_loc_p(dual_id) >= option%myrank) then
          pgrid%num_ufaces_local = pgrid%num_ufaces_local + 1
          pgridface_id = pgrid%ugridf2pgridf(iface)
          pgrid%uface_localids(pgrid%num_ufaces_local) = pgridface_id
          pgrid%uface_nverts(pgrid%num_ufaces_local) = pgrid%face_nverts(pgridface_id)
          pgrid%uface_left_natcellids(pgrid%num_ufaces_local) = nG2A(ghosted_id)
          pgrid%uface_right_natcellids(pgrid%num_ufaces_local) = nG2A(dual_id)
          pgrid%num_verts_of_ufaces_local = pgrid%num_verts_of_ufaces_local + pgrid%face_nverts(pgridface_id)
        endif
      else
        pgrid%num_ufaces_local = pgrid%num_ufaces_local + 1
        pgridface_id = pgrid%ugridf2pgridf(iface)
        pgrid%uface_localids(pgrid%num_ufaces_local) = pgridface_id
        pgrid%uface_nverts(pgrid%num_ufaces_local) = pgrid%face_nverts(pgridface_id)
        pgrid%num_verts_of_ufaces_local = pgrid%num_verts_of_ufaces_local + pgrid%face_nverts(pgridface_id)
        pgrid%uface_left_natcellids(pgrid%num_ufaces_local) = nG2A(ghosted_id)
        pgrid%uface_right_natcellids(pgrid%num_ufaces_local) = nG2A(dual_id)
      endif
    endif
  enddo
  call VecRestoreArrayF90(ghosted_cv_proc_rank, v_loc_p, ierr);CHKERRQ(ierr)

  call VecDestroy(ghosted_cv_proc_rank, ierr);CHKERRQ(ierr)


  call MPI_Allreduce(pgrid%num_ufaces_local, pgrid%num_ufaces_global, &
                      ONE_INTEGER_MPI, MPI_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Allreduce(pgrid%num_verts_of_ufaces_local, pgrid%num_verts_of_ufaces_global, &
                      ONE_INTEGER_MPI, MPI_INTEGER, MPI_SUM, option%mycomm, ierr)

  allocate(pgrid%uface_natvertids(pgrid%num_verts_of_ufaces_local))
  count = 0
  do iface = 1, pgrid%num_ufaces_local

    face_id = pgrid%uface_localids(iface)

    do ivertex = 1, pgrid%uface_nverts(iface)
      vertex_id = pgrid%face_vertids(ivertex, face_id)
      count = count + 1
      pgrid%uface_natvertids(count) = ugrid%vertex_ids_natural(vertex_id)
    enddo
  enddo

end subroutine UGridPolyhedraComputeOutputInfo

end module Grid_Unstructured_Polyhedra_module
