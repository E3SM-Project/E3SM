module Grid_Unstructured_Explicit_module
  
  use Geometry_module
  use Grid_Unstructured_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private 
  
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#if defined(SCORPIO)
  include "scorpiof.h"
#endif

  public :: UGridExplicitRead, &
            UGridExplicitDecompose, &
            UGridExplicitSetInternConnect, &
            UGridExplicitSetCellCentroids, &
            UGridExplicitComputeVolumes, &
            UGridExplicitSetBoundaryConnect, &
            UGridExplicitSetConnections

contains

! ************************************************************************** !

subroutine UGridExplicitRead(unstructured_grid,filename,option)
  ! 
  ! Reads an explicit unstructured grid in parallel
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/03/12
  ! 

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
 
  type(grid_unstructured_type) :: unstructured_grid 
  type(unstructured_explicit_type), pointer :: explicit_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word, card
  PetscInt :: fileid, icell, iconn, irank, remainder, temp_int, num_to_read
  
  PetscInt :: num_cells, num_connections, num_elems
  PetscInt :: num_cells_local, num_cells_local_save
  PetscInt :: num_connections_local, num_connections_local_save
  PetscInt :: num_elems_local, num_elems_local_save
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscErrorCode :: ierr
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscInt :: ivertex, num_vertices, num_grid_vertices 

  explicit_grid => unstructured_grid%explicit_grid 
! Format of explicit unstructured grid file
! id_, id_up_, id_dn_ = integer
! x_, y_, z_, area_, volume_ = real
! definitions
! id_ = id of grid cell
! id_up_ = id of upwind grid cell in connection
! id_dn_ = id of downwind grid cell in connection
! x_ = x coordinate of cell center
! y_ = y coordinate of cell center
! z_ = z coordinate of cell center
! volume_ = volume of grid cell
! -----------------------------------------------------------------
! CELLS <integer>    integer = # cells (N)
! id_1 x_1 y_1 z_1 volume_1
! id_2 x_2 y_2 z_2 volume_2
! ...
! ...
! id_N x_N y_N z_N volume_N
! CONNECTIONS <integer>   integer = # connections (M)
! id_up_1 id_dn_1 x_1 y_1 z_1 area_1
! id_up_2 id_dn_2 x_2 y_2 z_2 area_2
! ...
! ...
! id_up_M id_dn_M x_M y_M z_M area_M
! -----------------------------------------------------------------

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
  
    hint = 'Explicit Unstructured Grid CELLS'
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of cells',hint)
  endif
  
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
  num_cells = temp_int
  explicit_grid%num_cells_global = num_cells

   ! divide cells across ranks
  num_cells_local = num_cells/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = num_cells - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  allocate(explicit_grid%cell_ids(num_cells_local))
  explicit_grid%cell_ids = 0
  allocate(explicit_grid%cell_volumes(num_cells_local))
  explicit_grid%cell_volumes = 0
  allocate(explicit_grid%cell_centroids(num_cells_local))
  do icell = 1, num_cells_local
    explicit_grid%cell_centroids(icell)%x = 0.d0
    explicit_grid%cell_centroids(icell)%y = 0.d0
    explicit_grid%cell_centroids(icell)%z = 0.d0
  enddo
  
  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(5,num_cells_local_save+1))
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
        call InputReadDouble(input,option,temp_real_array(2,icell))
        call InputErrorMsg(input,option,'cell x coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(3,icell))
        call InputErrorMsg(input,option,'cell y coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(4,icell))
        call InputErrorMsg(input,option,'cell z coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(5,icell))
        call InputErrorMsg(input,option,'cell volume',hint)
      enddo

      ! if the cells reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_cells_local
        string = trim(adjustl(string)) // ' cells stored on p0'
        print *, trim(string)
#endif
        do icell = 1, num_cells_local
          explicit_grid%cell_ids(icell) = int(temp_real_array(1,icell))
          explicit_grid%cell_centroids(icell)%x = temp_real_array(2,icell)
          explicit_grid%cell_centroids(icell)%y = temp_real_array(3,icell)
          explicit_grid%cell_centroids(icell)%z = temp_real_array(4,icell)
          explicit_grid%cell_volumes(icell) = temp_real_array(5,icell)
        enddo
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' cells sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*5
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
  else
    ! other ranks post the recv
#if UGRID_DEBUG
    write(string,*) num_cells_local
    write(word,*) option%myrank
    string = trim(adjustl(string)) // ' cells received from p0 at p' // &
              trim(adjustl(word))
    print *, trim(string)
#endif
    allocate(temp_real_array(5,num_cells_local))
    int_mpi = num_cells_local*5
    call MPI_Recv(temp_real_array,int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    do icell = 1, num_cells_local
      explicit_grid%cell_ids(icell) = int(temp_real_array(1,icell))
      explicit_grid%cell_centroids(icell)%x = temp_real_array(2,icell)
      explicit_grid%cell_centroids(icell)%y = temp_real_array(3,icell)
      explicit_grid%cell_centroids(icell)%z = temp_real_array(4,icell)
      explicit_grid%cell_volumes(icell) = temp_real_array(5,icell)
    enddo
    
  endif
  deallocate(temp_real_array)
  
  if (option%myrank == option%io_rank) then
  
 
    call InputReadPflotranString(input,option)
    ! read CONNECTIONS card, though we already know the
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'CONNECTIONS'
    call InputErrorMsg(input,option,word,card)
    if (.not.StringCompare(word,card)) then
      option%io_buffer = 'Unrecognized keyword "' // trim(card) // &
        '" in explicit grid file.'
      call printErrMsgByRank(option)
    endif
  
    hint = 'Explicit Unstructured Grid CONNECTIONS'
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of connections',hint)
  endif
  
  int_mpi = 1
  call MPI_Bcast(temp_int,ONE_INTEGER_MPI,MPI_INTEGER,option%io_rank, &
                  option%mycomm,ierr)
  num_connections = temp_int
        
   ! divide cells across ranks
  num_connections_local = num_connections/option%mycommsize 
  num_connections_local_save = num_connections_local
  remainder = num_connections - &
              num_connections_local*option%mycommsize
  if (option%myrank < remainder) num_connections_local = &
                                 num_connections_local + 1

  allocate(explicit_grid%connections(2,num_connections_local))
  explicit_grid%connections = 0
  allocate(explicit_grid%face_areas(num_connections_local))
  explicit_grid%face_areas = 0    
  allocate(explicit_grid%face_centroids(num_connections_local))
  do iconn = 1, num_connections_local
    explicit_grid%face_centroids(iconn)%x = 0.d0
    explicit_grid%face_centroids(iconn)%y = 0.d0
    explicit_grid%face_centroids(iconn)%z = 0.d0
  enddo
        
  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(6,num_connections_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_read = num_connections_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do iconn = 1, num_to_read
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'cell id upwind',hint)
        temp_real_array(1,iconn) = dble(temp_int)
        call InputReadInt(input,option,temp_int)
        call InputErrorMsg(input,option,'cell id downwind',hint)
        temp_real_array(2,iconn) = dble(temp_int)
        call InputReadDouble(input,option,temp_real_array(3,iconn))
        call InputErrorMsg(input,option,'face x coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(4,iconn))
        call InputErrorMsg(input,option,'face y coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(5,iconn))
        call InputErrorMsg(input,option,'face z coordinate',hint)
        call InputReadDouble(input,option,temp_real_array(6,iconn))
        call InputErrorMsg(input,option,'face area',hint)
      enddo

      ! if the cells reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_connections_local
        string = trim(adjustl(string)) // ' connections stored on p0'
        print *, trim(string)
#endif
        do iconn = 1, num_connections_local
          explicit_grid%connections(1,iconn) = int(temp_real_array(1,iconn))
          explicit_grid%connections(2,iconn) = int(temp_real_array(2,iconn))
          explicit_grid%face_centroids(iconn)%x = temp_real_array(3,iconn)
          explicit_grid%face_centroids(iconn)%y = temp_real_array(4,iconn)
          explicit_grid%face_centroids(iconn)%z = temp_real_array(5,iconn)
          explicit_grid%face_areas(iconn) = temp_real_array(6,iconn)
        enddo
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' connections sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*6
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
  else
    ! other ranks post the recv
#if UGRID_DEBUG
    write(string,*) num_connections_local
    write(word,*) option%myrank
    string = trim(adjustl(string)) // ' connections received from p0 at p' // &
              trim(adjustl(word))
    print *, trim(string)
#endif
    allocate(temp_real_array(6,num_connections_local))
    int_mpi = num_connections_local*6
    call MPI_Recv(temp_real_array,int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    do iconn = 1, num_connections_local
      explicit_grid%connections(1,iconn) = int(temp_real_array(1,iconn))
      explicit_grid%connections(2,iconn) = int(temp_real_array(2,iconn))
      explicit_grid%face_centroids(iconn)%x = temp_real_array(3,iconn)
      explicit_grid%face_centroids(iconn)%y = temp_real_array(4,iconn)
      explicit_grid%face_centroids(iconn)%z = temp_real_array(5,iconn)
      explicit_grid%face_areas(iconn) = temp_real_array(6,iconn)
    enddo
    
  endif
  deallocate(temp_real_array)  
  
  if (option%myrank == option%io_rank) then
    call InputReadPflotranString(input,option)
    ! read ELEMENTS card, we only use this for tecplot output
    ! not used while solving the PDEs
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'ELEMENTS'
    if (.not.StringCompare(word,card)) return
    card = 'Explicit Unstruct. Grid ELEMENTS'
    call InputReadInt(input,option,num_elems)
    call InputErrorMsg(input,option,'number of elements',card)
        explicit_grid%num_elems = num_elems
    unstructured_grid%max_nvert_per_cell = 8 ! Initial guess
    allocate(explicit_grid%cell_vertices(0:unstructured_grid% &
                                  max_nvert_per_cell,num_elems)) 
    do iconn = 1, num_elems
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,card)  
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'element_type',card)
      call StringtoUpper(word)
      select case (word)
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
        case('TRI')
          num_vertices = 3
      end select
      explicit_grid%cell_vertices(0,iconn) = num_vertices
      do ivertex = 1, num_vertices
        call InputReadInt(input,option,explicit_grid%cell_vertices(ivertex,iconn))
        call InputErrorMsg(input,option,'vertex id',hint)
      enddo
    enddo
    call InputReadPflotranString(input,option)
    ! read VERTICES card, not used for calcuations, only tecplot output
    call InputReadWord(input,option,card,PETSC_TRUE)
    word = 'VERTICES'
    call InputErrorMsg(input,option,word,card)
    if (.not.StringCompare(word,card)) then
      option%io_buffer = 'Unrecognized keyword "' // trim(card) // &
        '" in explicit grid file.'
      call printErrMsgByRank(option)
    endif

    !at this point, as we read the grid, the output_mesh_type is not known yet 
    call InputReadInt(input,option,num_grid_vertices)

    if (InputError(input)) then
      input%ierr = 0
      !if num_grid_vertices not entered assumes vertex_centered based - default
      explicit_grid%num_vertices = explicit_grid%num_cells_global
    else   
      explicit_grid%num_vertices = num_grid_vertices
    end if

    allocate(explicit_grid%vertex_coordinates(explicit_grid%num_vertices))
    do icell = 1, explicit_grid%num_vertices
      explicit_grid%vertex_coordinates(icell)%x = 0.d0
      explicit_grid%vertex_coordinates(icell)%y = 0.d0
      explicit_grid%vertex_coordinates(icell)%z = 0.d0
    enddo
    do icell = 1, explicit_grid%num_vertices
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,card)  
      call InputReadDouble(input,option, &
                           explicit_grid%vertex_coordinates(icell)%x)
      call InputErrorMsg(input,option,'vertex 1',card)
      call InputReadDouble(input,option, &
                           explicit_grid%vertex_coordinates(icell)%y)
      call InputErrorMsg(input,option,'vertex 2',card)
      call InputReadDouble(input,option, &
                           explicit_grid%vertex_coordinates(icell)%z)
      call InputErrorMsg(input,option,'vertex 3',card)
    enddo
  endif
  
  if (option%myrank == option%io_rank) then
    call InputDestroy(input)
  endif
    
end subroutine UGridExplicitRead

! ************************************************************************** !

subroutine UGridExplicitDecompose(ugrid,option)
  ! 
  ! Decomposes an explicit unstructured grid across
  ! ranks
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/17/12
  ! 

  use Option_module
  use Utility_module, only: reallocateIntArray, SearchOrderedArray
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscdm.h" 
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscviewer.h"
  
  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option

  type(unstructured_explicit_type), pointer :: explicit_grid
  PetscViewer :: viewer
  
  Mat :: M_mat,M_mat_loc
  Vec :: M_vec
  Mat :: Adj_mat
  Mat :: Dual_mat
  MatPartitioning :: Part
  IS :: is_new
  IS :: is_scatter  
  IS :: is_gather  
  PetscInt :: num_cells_local_new, num_cells_local_old
  Vec :: cells_old, cells_local
  Vec :: connections_old, connections_local
  VecScatter :: vec_scatter
  
  PetscInt :: global_offset_old
  PetscInt :: global_offset_new
  PetscInt :: ghosted_id
  PetscInt, allocatable :: local_connections(:), local_connection_offsets(:)
  PetscInt, allocatable :: local_connections2(:), local_connection_offsets2(:)
  PetscInt, allocatable :: int_array(:), int_array2(:), int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: int_array2d(:,:)
  PetscInt :: num_connections_local_old, num_connections_local
  PetscInt :: num_connections_total 
  PetscInt :: num_connections_global, global_connection_offset
  PetscInt :: id_up, id_dn, iconn, icell, count, offset
  PetscInt :: conn_id, dual_id
  PetscBool :: found
  PetscInt :: i, temp_int, idual
  PetscReal :: temp_real
  
  PetscInt :: iflag
  PetscBool :: success
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscInt, pointer :: ia_ptr2(:), ja_ptr2(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscInt :: num_rows, num_cols, istart, iend, icol
  PetscInt :: cell_stride, dual_offset, connection_offset, connection_stride
  PetscInt :: natural_id_offset
  PetscErrorCode :: ierr
  PetscInt :: icell_up,icell_dn
  
  character(len=MAXSTRINGLENGTH) :: string

  explicit_grid => ugrid%explicit_grid
  
#if UGRID_DEBUG
  call printMsg(option,'Adjacency matrix')
#endif

  num_cells_local_old = size(explicit_grid%cell_ids)
  
  call MPI_Allreduce(num_cells_local_old,ugrid%nmax, &
                     ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr)  

  ! determine the global offset from 0 for cells on this rank
  global_offset_old = 0
  call MPI_Exscan(num_cells_local_old,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)  
  
  num_connections_local_old = size(explicit_grid%connections,2)
  
  call MPI_Allreduce(num_connections_local_old,num_connections_global, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  global_connection_offset = 0
  call MPI_Exscan(num_connections_local_old,global_connection_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  call VecCreateMPI(option%mycomm,num_cells_local_old,ugrid%nmax, &   
                    M_vec,ierr);CHKERRQ(ierr)
  call VecZeroEntries(M_vec,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local_old
    do i = 1, 2
      icell = explicit_grid%connections(i,iconn)-1
      call VecSetValue(M_vec,icell,1.d0,ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecAssemblyBegin(M_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(M_vec,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'M_vec.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(M_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecMax(M_vec,PETSC_NULL_INTEGER,temp_real,ierr);CHKERRQ(ierr)
  call VecDestroy(M_vec,ierr);CHKERRQ(ierr)
  ugrid%max_ndual_per_cell = int(temp_real+0.1d0)
  call MatCreateAIJ(option%mycomm,num_cells_local_old,PETSC_DECIDE, &
                    ugrid%nmax,num_connections_global, &
                    ugrid%max_ndual_per_cell,PETSC_NULL_INTEGER, &
                    ugrid%max_ndual_per_cell,PETSC_NULL_INTEGER, &
                    M_mat,ierr);CHKERRQ(ierr)
  call MatZeroEntries(M_mat,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local_old
    temp_int = iconn+global_connection_offset-1
    do i = 1, 2
      icell = explicit_grid%connections(i,iconn)-1
      call MatSetValue(M_mat,icell,temp_int,1.d0,INSERT_VALUES, &
                       ierr);CHKERRQ(ierr)
    enddo
  enddo
  call MatAssemblyBegin(M_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(M_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'M_mat.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(M_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! GB: When MatConvert() is used, the diagonal entries are lost in Adj_mat
  !call MatConvert(M_mat,MATMPIADJ,MAT_INITIAL_MATRIX,Adj_mat,ierr)
  !call MatDestroy(M_mat,ierr)

  ! Alternate method of creating Adj_mat
  if (option%mycommsize>1) then
    call MatMPIAIJGetLocalMat(M_mat,MAT_INITIAL_MATRIX,M_mat_loc, &
                              ierr);CHKERRQ(ierr)
    call MatGetRowIJF90(M_mat_loc,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  else
    call MatGetRowIJF90(M_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif

  count=0
  do icell = 1,num_rows
    istart = ia_ptr(icell)
    iend = ia_ptr(icell+1)-1
    num_cols = iend-istart+1
    count = count+num_cols
  enddo
  allocate(local_connections(count))
  allocate(local_connection_offsets(num_rows+1))
  local_connection_offsets(1:num_rows+1) = ia_ptr(1:num_rows+1)
  local_connections(1:count)             = ja_ptr(1:count)

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       num_connections_global, &
                       local_connection_offsets, &
                       local_connections,PETSC_NULL_INTEGER,Adj_mat, &
                       ierr);CHKERRQ(ierr)

  if (option%mycommsize>1) then
    call MatRestoreRowIJF90(M_mat_loc,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  else
    call MatRestoreRowIJF90(M_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  call MatDestroy(M_mat,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Adj.out',viewer,ierr);CHKERRQ(ierr)
  call MatView(Adj_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
!  call printErrMsg(option,'debugg')

  ! Create the Dual matrix.
  call MatCreateAIJ(option%mycomm,num_cells_local_old,PETSC_DECIDE, &
                    ugrid%nmax,ugrid%nmax, &
                    ugrid%max_ndual_per_cell,PETSC_NULL_INTEGER, &
                    ugrid%max_ndual_per_cell,PETSC_NULL_INTEGER, &
                    M_mat,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local_old
    icell_up = explicit_grid%connections(1,iconn)-1
    icell_dn = explicit_grid%connections(2,iconn)-1
    call MatSetValue(M_mat,icell_up,icell_dn,1.d0,INSERT_VALUES, &
                     ierr);CHKERRQ(ierr)
    call MatSetValue(M_mat,icell_dn,icell_up,1.d0,INSERT_VALUES, &
                     ierr);CHKERRQ(ierr)
  enddo
  
  call MatAssemblyBegin(M_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(M_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  
  !call MatConvert(M_mat,MATMPIADJ,MAT_INITIAL_MATRIX,Dual_mat,ierr)
  !call MatDestroy(M_mat,ierr)

  ! Alternate method of creating Dual_mat
  if (option%mycommsize>1) then
    call MatMPIAIJGetLocalMat(M_mat,MAT_INITIAL_MATRIX,M_mat_loc, &
                              ierr);CHKERRQ(ierr)
    call MatGetRowIJF90(M_mat_loc,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  else
    call MatGetRowIJF90(M_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif

  count=0
  do icell = 1,num_rows
    istart = ia_ptr(icell)
    iend = ia_ptr(icell+1)-1
    num_cols = iend-istart+1
    count = count+num_cols
  enddo
  allocate(local_connections2(count))
  allocate(local_connection_offsets2(num_rows+1))
  local_connection_offsets2(1:num_rows+1) = ia_ptr(1:num_rows+1)
  local_connections2(1:count)             = ja_ptr(1:count)

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       ugrid%nmax, &
                       local_connection_offsets2, &
                       local_connections2,PETSC_NULL_INTEGER,Dual_mat, &
                       ierr);CHKERRQ(ierr)

  if (option%mycommsize>1) then
    call MatRestoreRowIJF90(M_mat_loc,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  else
    call MatRestoreRowIJF90(M_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                        ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  call MatDestroy(M_mat,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Dual_mat.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(Dual_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call UGridPartition(ugrid,option,Dual_mat,is_new, &
                      num_cells_local_new)

  ! second argument of ZERO_INTEGER means to use 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  if (.not.success .or. num_rows /= num_cells_local_old) then
    print *, option%myrank, num_rows, success, num_cells_local_old
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                          num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  
  ! in order to redistributed cell/connection data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  connection_offset = 6 + 1 ! +1 for -777
  dual_offset = connection_offset + ugrid%max_ndual_per_cell + 1 ! +1 for -888
  cell_stride = dual_offset + ugrid%max_ndual_per_cell + 1 ! +1 for -999999
  natural_id_offset = 2

  ! Information for each cell is packed in a strided petsc vec
  ! The information is ordered within each stride as follows:
  ! -cell_N   ! global cell id (negative indicates 1-based)
  ! natural cell id
  ! cell x coordinate
  ! cell y coordinate
  ! cell z coordinate
  ! cell volume
  ! -777      ! separator between cell info and connection info
  ! conn1     ! connection ids between cell_N and others
  ! conn1
  ! ...
  ! connN     
  ! -888      ! separator between connection info and dual ids
  ! dual1     ! dual ids between cell_N and others
  ! dual2
  ! ...
  ! dualN     
  ! -999999   ! separator indicating end of information for cell_N
  
  ! the purpose of -888, and -999999 is to allow one to use cells of 
  ! various geometry.  
  
  call UGridCreateOldVec(ugrid,option,cells_old, &
                         num_cells_local_old, &
                         is_new,is_scatter,cell_stride)

  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  ! pointers to Dual mat
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  ! pointers to Adj mat
  call MatGetRowIJF90(Adj_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,temp_int, &
                      ia_ptr2,ja_ptr2,success,ierr);CHKERRQ(ierr)
  
  if (num_rows /= temp_int) then
    write(string,*) num_rows, temp_int
    option%io_buffer = 'Number of rows in Adj and Dual matrices inconsistent:'
    option%io_buffer = trim(option%io_buffer) // trim(adjustl(string))
    call printErrMsgByRank(option)
  endif

  call VecGetArrayF90(cells_old,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  do icell = 1, num_cells_local_old
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(global_offset_old+icell)
    count = count + 1
    vec_ptr(count) = explicit_grid%cell_ids(icell)
    count = count + 1
    vec_ptr(count) = explicit_grid%cell_centroids(icell)%x
    count = count + 1
    vec_ptr(count) = explicit_grid%cell_centroids(icell)%y
    count = count + 1
    vec_ptr(count) = explicit_grid%cell_centroids(icell)%z
    count = count + 1
    vec_ptr(count) = explicit_grid%cell_volumes(icell)
    ! add the separator
    count = count + 1
    vec_ptr(count) = -777  ! help differentiate
 
    ! add the connection ids
    istart = ia_ptr2(icell)
    iend = ia_ptr2(icell+1)-1
    num_cols = iend-istart+1
    if (num_cols > ugrid%max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Adj matrix is larger then max_ndual_per_cell.'
      call printErrMsgByRank(option)
    endif
    do icol = 1, ugrid%max_ndual_per_cell
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr2(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo
    
    ! add the separator
    count = count + 1
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
    istart = ia_ptr(icell)
    iend = ia_ptr(icell+1)-1
    num_cols = iend-istart+1
    if (num_cols > ugrid%max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Dual matrix is larger then max_ndual_per_cell.'
      call printErrMsgByRank(option)
    endif
    do icol = 1, ugrid%max_ndual_per_cell
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo

    ! final separator
    count = count + 1 
    vec_ptr(count) = -999999  ! help differentiate
  enddo
  call VecRestoreArrayF90(cells_old,vec_ptr,ierr);CHKERRQ(ierr)
  
  ! pointers to Dual mat
  call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                          num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  ! pointers to Adj mat
  call MatRestoreRowIJF90(Adj_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                          temp_int,ia_ptr2,ja_ptr2,success,ierr);CHKERRQ(ierr)
  call MatDestroy(Dual_mat,ierr);CHKERRQ(ierr)
  call MatDestroy(Adj_mat,ierr);CHKERRQ(ierr)
  deallocate(local_connections)
  deallocate(local_connection_offsets)
  deallocate(local_connections2)
  deallocate(local_connection_offsets2)


  ! is_scatter is destroyed within UGridNaturalToPetsc
  call UGridNaturalToPetsc(ugrid,option, &
                           cells_old,cells_local, &
                           num_cells_local_new,cell_stride,dual_offset, &
                           natural_id_offset,is_scatter)
  
  call VecDestroy(cells_old,ierr);CHKERRQ(ierr)

  ! set up connections
  connection_stride = 8
  ! create strided vector with the old connection distribution
  call VecCreate(option%mycomm,connections_old,ierr);CHKERRQ(ierr)
  call VecSetSizes(connections_old, &
                   connection_stride*num_connections_local_old, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(connections_old,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(connections_old,vec_ptr,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local_old
    offset = (iconn-1)*connection_stride
    vec_ptr(offset+1) = explicit_grid%connections(1,iconn)
    vec_ptr(offset+2) = explicit_grid%connections(2,iconn)
    vec_ptr(offset+3) = explicit_grid%face_centroids(iconn)%x
    vec_ptr(offset+4) = explicit_grid%face_centroids(iconn)%y
    vec_ptr(offset+5) = explicit_grid%face_centroids(iconn)%z
    vec_ptr(offset+6) = explicit_grid%face_areas(iconn)
    vec_ptr(offset+7) = 1.d0 ! flag for local connections
    vec_ptr(offset+8) = -888.d0
  enddo
  call VecRestoreArrayF90(connections_old,vec_ptr,ierr);CHKERRQ(ierr)

    
#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'connections_old.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(connections_old,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif    
  
  ! count the number of cells and their duals  
  call VecGetArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  do ghosted_id=1, ugrid%ngmax
    do iconn = 1, ugrid%max_ndual_per_cell
      conn_id = int(vec_ptr(iconn + connection_offset + (ghosted_id-1)*cell_stride))
      if (conn_id < 1) exit ! here we hit the 0 at the end of last connection
      ! yes, we will be counting them twice
      count = count + 1
    enddo
  enddo   
  num_connections_total = count ! many of these are redundant and will be removed
  ! allocate and fill an array with the natural cell and dual ids
  allocate(int_array(num_connections_total))
  count = 0
  do ghosted_id=1, ugrid%ngmax
    do iconn = 1, ugrid%max_ndual_per_cell
      conn_id = int(vec_ptr(iconn + connection_offset + (ghosted_id-1)*cell_stride))
      if (conn_id < 1) exit ! again we hit the 0 
      count = count + 1
      int_array(count) = conn_id
    enddo
  enddo
  call VecRestoreArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  
  allocate(int_array2(num_connections_total))
  do iconn = 1, num_connections_total
    int_array2(iconn) = iconn
  enddo
  
  ! sort connections - int_array2 will return the reordering while int_array 
  !                    remains the same.
  int_array2 = int_array2 - 1
  call PetscSortIntWithPermutation(num_connections_total,int_array, &
                                   int_array2,ierr);CHKERRQ(ierr)
  int_array2 = int_array2 + 1
  
  ! permute, remove duplicate connections, and renumber to local ordering
  allocate(int_array3(num_connections_total))
  allocate(int_array4(num_connections_total))
  int_array3 = UNINITIALIZED_INTEGER
  int_array4 = UNINITIALIZED_INTEGER
  int_array3(1) = int_array(int_array2(1))
  count = 1
  do iconn = 1, num_connections_total
    if (int_array(int_array2(iconn)) > int_array3(count)) then
      count = count + 1
      int_array3(count) = int_array(int_array2(iconn))
    endif
    int_array4(int_array2(iconn)) = count
  enddo
  deallocate(int_array)
  deallocate(int_array2)
  
  num_connections_local = count ! new compressed count
  allocate(int_array(num_connections_local))
  int_array = int_array3(1:num_connections_local)
  deallocate(int_array3)
  
  ! replace original connections ids (naturally numbered) with locally 
  ! numbered connection ids (int_array4)
  call VecGetArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  do ghosted_id=1, ugrid%ngmax
    do iconn = 1, ugrid%max_ndual_per_cell
      conn_id = int(vec_ptr(iconn + connection_offset + (ghosted_id-1)*cell_stride))
      if (conn_id < 1) exit ! again we hit the 0 
      count = count + 1
      vec_ptr(iconn + connection_offset + (ghosted_id-1)*cell_stride) = &
        int_array4(count)
    enddo
  enddo
  deallocate(int_array4)
  call VecRestoreArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  
  ! check to ensure that the number before/after are consistent
  if (count /= num_connections_total) then
    write(string,'(2i6)') count, num_connections_total
    option%io_buffer = 'Inconsistent values for num_connections_total: ' // &
      trim(adjustl(string))
    call printErrMsgByRank(option)
  endif
  num_connections_total = UNINITIALIZED_INTEGER ! set to uninitialized value to catch bugs
  
  call VecCreate(PETSC_COMM_SELF,connections_local,ierr);CHKERRQ(ierr)
  call VecSetSizes(connections_local,num_connections_local*connection_stride, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(connections_local,connection_stride,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(connections_local,ierr);CHKERRQ(ierr)

  int_array = int_array-1
  call ISCreateBlock(option%mycomm,connection_stride,num_connections_local, &
                     int_array,PETSC_COPY_VALUES,is_scatter, &
                     ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local
    int_array(iconn) = iconn-1
  enddo
  call ISCreateBlock(option%mycomm,connection_stride,num_connections_local, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter_conn_old_to_local.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(is_scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'is_gather_conn_old_to_local.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(is_gather,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  call printMsg(option,'Scatter/gathering local connection info')
#endif  
  
  ! scatter all the connection data from the old to local
  call VecScatterCreate(connections_old,is_scatter,connections_local, &
                        is_gather,vec_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_gather,ierr);CHKERRQ(ierr)
  call ISDestroy(is_scatter,ierr);CHKERRQ(ierr)
  call VecScatterBegin(vec_scatter,connections_old,connections_local, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,connections_old,connections_local, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

    
#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'connections_local_nat' // trim(adjustl(string)) // '.out'
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(connections_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif   
  
  ! loop over cells and change the natural ids in the duals to local ids
  allocate(int_array2d(2,num_connections_local))
  int_array2d = UNINITIALIZED_INTEGER
  call VecGetArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id=1, ugrid%ngmax
    do iconn = 1, ugrid%max_ndual_per_cell
      ! this connection id is now local
      conn_id = int(vec_ptr(iconn + connection_offset + (ghosted_id-1)*cell_stride))
      if (conn_id < 1) exit ! again we hit the 0
      do i = 1, 2
        if (int_array2d(i,conn_id) <= UNINITIALIZED_INTEGER) then
          int_array2d(i,conn_id) = ghosted_id
          exit
        endif
        if (i > 2) then
          write(string,'(2i5)') ghosted_id, conn_id
          option%io_buffer = 'Too many local cells match connection: ' // &
            trim(adjustl(string))
          call printErrMsgByRank(option)
        endif
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)

  ! map natural ids in connections to local ids
  ! negate connection ids as a flag
  int_array2d = -1*int_array2d
  call VecGetArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local
    offset = connection_stride*(iconn-1)
    ! all values should be negative at this point, unless uninitialized
    if (maxval(int_array2d(:,iconn)) >= 999) then
      ! connection is between two ghosted cells
      vec_ptr(offset+7) = 0.d0
      cycle
    endif
    id_up = int(vec_ptr(offset+1)) ! this is the natural id
    id_dn = int(vec_ptr(offset+2))
    count = 0
    found = PETSC_FALSE
    do i = 1, 2
      if (ugrid%cell_ids_natural(abs(int_array2d(i,iconn))) == id_up) then
        int_array2d(i,iconn) = abs(int_array2d(i,iconn))
        found = PETSC_TRUE
        count = count + 1
      endif
      if (ugrid%cell_ids_natural(abs(int_array2d(i,iconn))) == id_dn) then
        int_array2d(i,iconn) = abs(int_array2d(i,iconn))
        found = PETSC_TRUE
        count = count - 1
      endif
    enddo
    ! count should be zero, meaning it found the upwind and downwind cell
    ! ids
    if (count /= 0 .or. .not.found) then
      write(string,*) iconn, id_up, id_dn
      if (.not.found) then
        option%io_buffer = 'upwind/downwind cells not found: '
      else if (count < 0) then
        option%io_buffer = 'upwind cell not found: '
      else
        option%io_buffer = 'downwind cell not found: '
      endif
      option%io_buffer = trim(option%io_buffer) // trim(adjustl(string))
      call printErrMsgByRank(option)
    endif
    id_up = int_array2d(1,iconn)
    id_dn = int_array2d(2,iconn)
    if (id_up < id_dn) then
      vec_ptr(offset+1) = id_up  ! now local ids
      vec_ptr(offset+2) = id_dn 
    else
      vec_ptr(offset+1) = id_dn
      vec_ptr(offset+2) = id_up 
    endif
    if (id_up > ugrid%nlmax .and. id_dn > ugrid%nlmax) then
      ! connection is between two ghosted cells
      vec_ptr(offset+7) = 0.d0
    endif
  enddo
  call VecRestoreArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(int_array2d)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'connections_local_loc' // trim(adjustl(string)) // '.out'
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(connections_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif   

  ! deallocate/allocate grid cell info locally
  deallocate(explicit_grid%cell_ids)
  deallocate(explicit_grid%cell_volumes)
  deallocate(explicit_grid%cell_centroids)

  allocate(explicit_grid%cell_ids(ugrid%ngmax))
  explicit_grid%cell_ids = UNINITIALIZED_INTEGER
  allocate(explicit_grid%cell_volumes(ugrid%ngmax))
  explicit_grid%cell_volumes = UNINITIALIZED_DOUBLE
  allocate(explicit_grid%cell_centroids(ugrid%ngmax))
  do icell = 1, ugrid%ngmax
    explicit_grid%cell_centroids(icell)%x = UNINITIALIZED_DOUBLE
    explicit_grid%cell_centroids(icell)%y = UNINITIALIZED_DOUBLE
    explicit_grid%cell_centroids(icell)%z = UNINITIALIZED_DOUBLE
  enddo

  call VecGetArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id=1, ugrid%ngmax
    offset = cell_stride*(ghosted_id-1)
    explicit_grid%cell_ids(ghosted_id) = int(vec_ptr(offset + 2))
    explicit_grid%cell_centroids(ghosted_id)%x = vec_ptr(offset + 3)
    explicit_grid%cell_centroids(ghosted_id)%y = vec_ptr(offset + 4)
    explicit_grid%cell_centroids(ghosted_id)%z = vec_ptr(offset + 5)
    explicit_grid%cell_volumes(ghosted_id) = vec_ptr(offset + 6)
  enddo
  call VecRestoreArrayF90(cells_local,vec_ptr,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'cells_local_raw' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, ugrid%ngmax
    write(86,'(i5,4f10.3)') explicit_grid%cell_ids(ghosted_id), &
                explicit_grid%cell_centroids(ghosted_id)%x, &
                explicit_grid%cell_centroids(ghosted_id)%y, &
                explicit_grid%cell_centroids(ghosted_id)%z, &
                explicit_grid%cell_volumes(ghosted_id)
  enddo
  close(86)
#endif     
  
  ! deallocate/allocate connection info locally
  deallocate(explicit_grid%connections)
  deallocate(explicit_grid%face_areas)
  deallocate(explicit_grid%face_centroids)

  count = 0
  call VecGetArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)
  do iconn = 1, num_connections_local
    offset = connection_stride*(iconn-1)
    if (vec_ptr(offset+7) > 0.1d0) count = count + 1
  enddo
  call VecRestoreArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)

  allocate(explicit_grid%connections(2,count))
  explicit_grid%connections = 0
  allocate(explicit_grid%face_areas(count))
  explicit_grid%face_areas = 0    
  allocate(explicit_grid%face_centroids(count))
  do iconn = 1, count
    explicit_grid%face_centroids(iconn)%x = 0.d0
    explicit_grid%face_centroids(iconn)%y = 0.d0
    explicit_grid%face_centroids(iconn)%z = 0.d0
  enddo  
  call VecGetArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  do iconn = 1, num_connections_local
    offset = connection_stride*(iconn-1)
    if (vec_ptr(offset+7) > 0.1d0) then
      count = count + 1
      explicit_grid%connections(1,count) = int(vec_ptr(offset+1))
      explicit_grid%connections(2,count) = int(vec_ptr(offset+2))
      explicit_grid%face_centroids(count)%x = vec_ptr(offset+3)
      explicit_grid%face_centroids(count)%y = vec_ptr(offset+4)
      explicit_grid%face_centroids(count)%z = vec_ptr(offset+5)
      explicit_grid%face_areas(count) = vec_ptr(offset+6)
    endif
  enddo
  call VecRestoreArrayF90(connections_local,vec_ptr,ierr);CHKERRQ(ierr)
  num_connections_local = count

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'connections_local_raw' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do iconn = 1, num_connections_local
    write(86,'(2i5,4f7.3)') explicit_grid%connections(1,iconn), &
                explicit_grid%connections(2,iconn), &
                explicit_grid%face_centroids(iconn)%x, &
                explicit_grid%face_centroids(iconn)%y, &
                explicit_grid%face_centroids(iconn)%z, &
                explicit_grid%face_areas(iconn)
  enddo
  close(86)
#endif     
  
  call VecDestroy(connections_old,ierr);CHKERRQ(ierr)
  call VecDestroy(connections_local,ierr);CHKERRQ(ierr)
  call VecDestroy(cells_local,ierr);CHKERRQ(ierr)
  
end subroutine UGridExplicitDecompose

! ************************************************************************** !

subroutine UGridExplicitSetCellCentroids(explicit_grid,x,y,z, &
                                         x_min,x_max,y_min,y_max,z_min,z_max)
  ! 
  ! Sets the centroid of each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/17/12
  ! 

  use Option_module

  implicit none
  
  type(unstructured_explicit_type) :: explicit_grid
  PetscReal :: x(:), y(:), z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: icell
  
  do icell = 1, size(explicit_grid%cell_centroids)
    x(icell) = explicit_grid%cell_centroids(icell)%x
    y(icell) = explicit_grid%cell_centroids(icell)%y
    z(icell) = explicit_grid%cell_centroids(icell)%z
  enddo
  
  x_min = minval(x)
  x_max = maxval(x)
  y_min = minval(y)
  y_max = maxval(y)
  z_min = minval(z)
  z_max = maxval(z)
      
end subroutine UGridExplicitSetCellCentroids

! ************************************************************************** !

function UGridExplicitSetInternConnect(explicit_grid,option)
  ! 
  ! Sets up the internal connectivity within
  ! the connectivity object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/17/12
  ! 

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: UGridExplicitSetInternConnect
  
  type(unstructured_explicit_type) :: explicit_grid
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id_up, id_dn
  PetscReal :: v(3), v_up(3), v_dn(3)
  PetscReal :: distance
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: error 
  
  num_connections = size(explicit_grid%connections,2)
  connections => ConnectionCreate(num_connections,INTERNAL_CONNECTION_TYPE)
  
  error = PETSC_FALSE
  do iconn = 1, num_connections
    id_up = explicit_grid%connections(1,iconn)
    id_dn = explicit_grid%connections(2,iconn)
    connections%id_up(iconn) = id_up
    connections%id_dn(iconn) = id_dn
    
    v_up(1) = explicit_grid%face_centroids(iconn)%x - &
              explicit_grid%cell_centroids(id_up)%x
    v_up(2) = explicit_grid%face_centroids(iconn)%y - &
              explicit_grid%cell_centroids(id_up)%y
    v_up(3) = explicit_grid%face_centroids(iconn)%z - &
              explicit_grid%cell_centroids(id_up)%z

    v_dn(1) = explicit_grid%cell_centroids(id_dn)%x - &
              explicit_grid%face_centroids(iconn)%x
    v_dn(2) = explicit_grid%cell_centroids(id_dn)%y - &
              explicit_grid%face_centroids(iconn)%y
    v_dn(3) = explicit_grid%cell_centroids(id_dn)%z - &
              explicit_grid%face_centroids(iconn)%z

    v = v_up + v_dn
    distance = sqrt(DotProduct(v,v))
    if (dabs(distance) < 1.d-40) then
      write(string,'(2(es16.9,","),es16.9)') &
        explicit_grid%face_centroids(iconn)%x, &
        explicit_grid%face_centroids(iconn)%y, &
        explicit_grid%face_centroids(iconn)%z
      error = PETSC_TRUE
      option%io_buffer = 'Coincident cell and face centroids found at (' // &
        trim(adjustl(string)) // ') '
      call printMsgByRank(option)
    endif
    connections%dist(-1,iconn) = sqrt(DotProduct(v_up,v_up))/distance
    connections%dist(0,iconn) = distance
    connections%dist(1:3,iconn) = v/distance
    connections%area(iconn) = explicit_grid%face_areas(iconn)
  enddo
  if (error) then
    option%io_buffer = 'Coincident cell and face centroids found in ' // &
      'UGridExplicitSetInternConnect().  See details above.'
    call printErrMsgByRank(option)
  endif
  
  UGridExplicitSetInternConnect => connections

end function UGridExplicitSetInternConnect

! ************************************************************************** !

subroutine UGridExplicitComputeVolumes(ugrid,option,volume)
  ! 
  ! Sets the volume of each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/17/12
  ! 

  use Option_module

  implicit none
  
  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  Vec :: volume

  type(unstructured_explicit_type), pointer :: explicit_grid
  
  PetscInt :: icell
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  explicit_grid => ugrid%explicit_grid

  call VecGetArrayF90(volume,vec_ptr,ierr);CHKERRQ(ierr)
  do icell = 1, ugrid%nlmax
    vec_ptr(icell) = explicit_grid%cell_volumes(icell)
  enddo
  call VecRestoreArrayF90(volume,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine UGridExplicitComputeVolumes

! ************************************************************************** !

function UGridExplicitSetBoundaryConnect(explicit_grid,cell_ids, &
                                         face_centroids,face_areas, &
                                         region_name,option)
  ! 
  ! Sets up the boundary connectivity within
  ! the connectivity object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  ! 

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: UGridExplicitSetBoundaryConnect

  type(unstructured_explicit_type) :: explicit_grid
  PetscInt :: cell_ids(:)
  type(point3d_type) :: face_centroids(:)
  PetscReal :: face_areas(:)
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id
  PetscReal :: v(3)
  PetscReal :: distance
  character(len=MAXSTRINGLENGTH) :: string 
  PetscBool :: error
 
  num_connections = size(cell_ids)
  connections => ConnectionCreate(num_connections,BOUNDARY_CONNECTION_TYPE)
  
  error = PETSC_FALSE
  do iconn = 1, num_connections
    id = cell_ids(iconn)
    connections%id_dn(iconn) = id
    
    v(1) = explicit_grid%cell_centroids(id)%x - &
           face_centroids(iconn)%x
    v(2) = explicit_grid%cell_centroids(id)%y - &
           face_centroids(iconn)%y
    v(3) = explicit_grid%cell_centroids(id)%z - &
           face_centroids(iconn)%z

    distance = sqrt(DotProduct(v,v))
    if (dabs(distance) < 1.d-40) then
      write(string,'(2(es16.9,","),es16.9)') &
        face_centroids(iconn)%x, face_centroids(iconn)%y, &
        face_centroids(iconn)%z
      error = PETSC_TRUE
      option%io_buffer = 'Coincident cell and face centroids found at (' // &
        trim(adjustl(string)) // ') '
      call printMsgByRank(option)
    endif
    connections%dist(-1,iconn) = 0.d0
    connections%dist(0,iconn) = distance
    connections%dist(1:3,iconn) = v/distance
    connections%area(iconn) = face_areas(iconn)
  enddo
  if (error) then
    option%io_buffer = 'Coincident cell and face centroids found in ' // &
      'UGridExplicitSetBoundaryConnect() for region "' // trim(region_name) // &
      '".  See details above.'
    call printErrMsgByRank(option)
  endif
  
  UGridExplicitSetBoundaryConnect => connections

end function UGridExplicitSetBoundaryConnect

! ************************************************************************** !

function UGridExplicitSetConnections(explicit_grid,cell_ids,connection_type, &
                                     option)
  ! 
  ! Sets up the connectivity for a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  ! 

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: UGridExplicitSetConnections

  type(unstructured_explicit_type) :: explicit_grid
  PetscInt, pointer :: cell_ids(:)
  PetscInt :: connection_type
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id
  
  num_connections = 0
  if (associated(cell_ids)) then
    num_connections = size(cell_ids)
  endif
  connections => ConnectionCreate(num_connections,connection_type)
    
  do iconn = 1, num_connections
    id = cell_ids(iconn)
    connections%id_dn(iconn) = id
  enddo
  
  UGridExplicitSetConnections => connections

end function UGridExplicitSetConnections

end module Grid_Unstructured_Explicit_module
