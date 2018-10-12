module Grid_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Polyhedra_module
  use Connection_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: grid_type 
  
    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured_grid, implicit_unstructured_grid, etc.)
    
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt :: global_offset ! Offset of first cell on process in petsc ordering
    PetscInt :: nlmax_faces  ! Total number of non-ghosted faces in local domain.
    PetscInt :: ngmax_faces  ! Number of ghosted & non-ghosted faces in local domain.
    PetscInt :: nmax_faces  ! Number of ghosted & non-ghosted faces in local domain.
    PetscInt :: global_cell_offset, global_faces_offset  ! offsets for LP formulation
   
    ! Below, we define several arrays used for mapping between different 
    ! types of array indices.  Our terminology is as follows:
    !
    ! 'Local' indices are used to access arrays containing values that are 
    ! entirely local to the MPI process -- these arrays contain no "ghost" 
    ! entries used to hold copies of values that are owned by neighboring 
    ! processes.
    !
    ! 'Ghosted local' (or simply 'ghost') indices are used to access arrays 
    ! that contain additional entries that hold copies of values that are 
    ! owned by neighboring processes.  (These entries are filled in by 
    ! DMGlobalToLocalBegin/End() in the structured grid case.)
    !
    ! Entries of a vector created with DMCreateGlobalVector() should be 
    ! indexed using 'local' indices.  The array returned from a call to 
    ! VecGetArrayF90() on such a vector consists of local entries only and 
    ! NO ghost points.
    !
    ! Entries of a vector created with DMCreateLocalVector() should be 
    ! indexed using 'ghosted local' indices.  The array returned from a call 
    ! to VecGetArrayF90() on such a vector contains the truly3 local entries 
    ! as well as ghost points.
    !
    ! The index mapping arrays are the following:
    ! nL2G :  not collective, local processor: local  =>  ghosted local  
    ! nG2L :  not collective, local processor:  ghosted local => local  
    ! nG2A :  not collective, ghosted local => natural

    PetscInt, pointer :: nL2G(:), nG2L(:)
    PetscInt, pointer :: nG2A(:)

    PetscReal, pointer :: x(:), y(:), z(:) ! coordinates of ghosted grid cells

    PetscReal :: x_min_global, y_min_global, z_min_global
    PetscReal :: x_max_global, y_max_global, z_max_global
    PetscReal :: x_min_local, y_min_local, z_min_local
    PetscReal :: x_max_local, y_max_local, z_max_local

    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash_bins

    type(grid_structured_type), pointer :: structured_grid
    type(grid_unstructured_type), pointer :: unstructured_grid
    
    type(connection_set_list_type), pointer :: internal_connection_set_list
    type(connection_set_list_type), pointer :: boundary_connection_set_list

    ! list of connections defined over specific regions
    type(connection_set_list_type), pointer :: reg_internal_connection_set_list
    type(connection_set_list_type), pointer :: reg_boundary_connection_set_list
    
  end type grid_type
  
  type, public :: face_type
    type(connection_set_type), pointer :: conn_set_ptr
    PetscInt :: id
  end type face_type

  public :: GridCreate, &
            GridDestroy, &
            GridComputeInternalConnect, &
            GridMapIndices, &
            GridComputeSpacing, &
            GridComputeCoordinates, &
            GridComputeVolumes, &
            GridComputeAreas, &
            GridLocalizeRegions, &
            GridPopulateConnection, &
            GridCopyIntegerArrayToVec, &
            GridCopyRealArrayToVec, &
            GridCopyVecToIntegerArray, &
            GridCopyVecToRealArray, &
            GridCreateNaturalToGhostedHash, &
            GridDestroyHashTable, &
            GridGetLocalGhostedIdFromHash, &
            GridIndexToCellID, &
            GridGetGhostedNeighbors, &
            GridGetGhostedNeighborsWithCorners, &
            GridMapCellsInPolVol, &
            GridGetLocalIDFromCoordinate, &
            GridRestrictRegionalConnect, &
            GridPrintExtents
  
contains

! ************************************************************************** !

function GridCreate()
  ! 
  ! Creates a structured or unstructured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(grid_type), pointer :: GridCreate
  
  type(grid_type), pointer :: grid
  
  allocate(grid)
  grid%ctype = ''
  grid%itype = 0

  nullify(grid%structured_grid)
  nullify(grid%unstructured_grid)

  nullify(grid%internal_connection_set_list)
  nullify(grid%reg_internal_connection_set_list)

  nullify(grid%nL2G)
  nullify(grid%nG2L)
  nullify(grid%nG2A)

  nullify(grid%x)
  nullify(grid%y)
  nullify(grid%z)

  grid%x_min_global = 1.d20
  grid%x_max_global = -1.d20
  grid%y_min_global = 1.d20
  grid%y_max_global = -1.d20
  grid%z_min_global = 1.d20
  grid%z_max_global = -1.d20

  grid%x_min_local = 1.d20
  grid%x_max_local = -1.d20
  grid%y_min_local = 1.d20
  grid%y_max_local = -1.d20
  grid%z_min_local = 1.d20
  grid%z_max_local = -1.d20

  grid%nmax = 0
  grid%nlmax = 0 
  grid%ngmax = 0
  grid%global_offset = 0

  nullify(grid%hash)
  grid%num_hash_bins = 1000

  GridCreate => grid

end function GridCreate

! ************************************************************************** !

subroutine GridComputeInternalConnect(grid,option,ugdm)
  ! 
  ! computes internal connectivity of a grid
  ! sp modified December 2010
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/17/07
  ! 

  use Connection_module
  use Option_module
  use Grid_Unstructured_Explicit_module
  use Grid_Unstructured_Polyhedra_module
    
  implicit none
  
  PetscInt ierr

  type(grid_type) :: grid
  type(option_type) :: option
  type(ugdm_type), optional :: ugdm
  
  type(connection_set_type), pointer :: connection_set, connection_bound_set
  type(connection_set_type), pointer :: connection_set_2
  nullify(connection_set); nullify(connection_bound_set)
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      connection_set => &
        StructGridComputeInternConnect( grid%structured_grid, grid%x, grid%y, &
                                    grid%z, option)
    case(IMPLICIT_UNSTRUCTURED_GRID) 
      connection_set => &
        UGridComputeInternConnect(grid%unstructured_grid,grid%x,grid%y, &
                                  grid%z,option)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      connection_set => &
        UGridExplicitSetInternConnect(grid%unstructured_grid%explicit_grid, &
                                      grid%unstructured_grid% &
                                        upwind_fraction_method, &
                                      option)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      connection_set => &
        UGridPolyhedraComputeInternConnect(grid%unstructured_grid, &
                                           grid%x, grid%y, grid%z, &
                                           option)
      call UGridPolyhedraComputeOutputInfo(grid%unstructured_grid, grid%nL2G, &
                                           grid%nG2L, grid%nG2A, option)
  end select
  
  allocate(grid%internal_connection_set_list)
  call ConnectionInitList(grid%internal_connection_set_list)
  call ConnectionAddToList(connection_set,grid%internal_connection_set_list)


  select case(grid%itype)
    case(IMPLICIT_UNSTRUCTURED_GRID) 
!      connection_bound_set => &
!        UGridComputeBoundConnect(grid%unstructured_grid,option)
  end select

end subroutine GridComputeInternalConnect

! ************************************************************************** !

function ConnectionSetIntersectRegion(connection_set,region) result(reg_connection_set)
  ! 
  ! Returns a pointer to a new connection set created from the input
  ! set, where cell ids belong to the input region. Important: the
  ! cell ids of the regional set are local to the region.
  !
  ! Author: Nathan Collier
  ! Date: 09/2015
  ! 

  use Region_module
  
  implicit none
  type(connection_set_type), pointer :: connection_set,reg_connection_set  
  type(region_type),         pointer :: region

  PetscInt, allocatable :: ids(:,:)
  PetscInt              :: i,j,up,dn,nconn

  ! first pass to find connection ids and number of connections
  nconn = 0
  allocate(ids(connection_set%num_connections,3))
  do i = 1,connection_set%num_connections
     up = -1
     dn = -1
     do j = 1,region%num_cells
        if (connection_set%id_up(i) == region%cell_ids(j)) up = j
        if (connection_set%id_dn(i) == region%cell_ids(j)) dn = j
        if (up > 0 .and. dn > 0) then
           nconn = nconn + 1
           ids(nconn,1) = i
           ids(nconn,2) = up
           ids(nconn,3) = dn
           exit
        endif
     enddo
  enddo

  ! second pass to load the information
  nullify(reg_connection_set)
  if (nconn > 0) then
     reg_connection_set => ConnectionCreate(nconn,connection_set%itype)
     do i = 1,nconn
        j = ids(i,1)
        reg_connection_set%id_up  (  i) = ids(i,2)
        reg_connection_set%id_dn  (  i) = ids(i,3)
        reg_connection_set%dist   (:,i) = connection_set%dist   (:,j)
        reg_connection_set%intercp(:,i) = connection_set%intercp(:,j)
        reg_connection_set%area   (  i) = connection_set%area   (  j)
        reg_connection_set%face_id(  i) = connection_set%face_id(  j)
     enddo
  endif

  ! cleanup and return
  deallocate(ids)
  
end function ConnectionSetIntersectRegion

! ************************************************************************** !

subroutine GridRestrictRegionalConnect(grid,region)
  ! 
  ! Populates the internal regional connection list of a grid
  ! 
  ! Author: Nathan Collier
  ! Date: 09/2015
  !
  use Region_module
  
  implicit none
  
  type(grid_type)           :: grid
  type(region_type),pointer :: region

  type(connection_set_type), pointer :: cur_connection_set,reg_connection_set

  ! initialize the regional connection list
  if (.not.associated(grid%reg_internal_connection_set_list)) then
     allocate(grid%reg_internal_connection_set_list)
     call ConnectionInitList(grid%reg_internal_connection_set_list)
  endif

  ! populate the list
  cur_connection_set => grid%internal_connection_set_list%first
  do 
     if (.not.associated(cur_connection_set)) exit
     reg_connection_set => ConnectionSetIntersectRegion(cur_connection_set,region)
     if (associated(reg_connection_set)) then
        call ConnectionAddToList(reg_connection_set,grid%reg_internal_connection_set_list)
     endif
     cur_connection_set => cur_connection_set%next
  enddo
  
end subroutine GridRestrictRegionalConnect

! ************************************************************************** !

subroutine GridPopulateConnection(grid,connection,iface,iconn,cell_id_local, &
                                  option)
  ! 
  ! computes connectivity coupler to a grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  use Connection_module
  use Grid_Structured_module
  use Option_module
  
  implicit none
 
  type(grid_type) :: grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_local
  type(option_type) :: option
  
  PetscInt :: cell_id_ghosted
  
  cell_id_ghosted = grid%nL2G(cell_id_local)
  ! Use ghosted index to access dx, dy, dz because we have
  ! already done a global-to-local scatter for computing the
  ! interior node connections.
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridPopulateConnection(grid%x,grid%structured_grid,connection, &
                                        iface,iconn,cell_id_ghosted,option)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridPopulateConnection(grid%unstructured_grid,connection,iface,&
                                   iconn,cell_id_ghosted,option)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call UGridPolyhedraPopulateConnection(grid%unstructured_grid,connection,iface, &
                                            iconn,cell_id_ghosted,option)
  end select

end subroutine GridPopulateConnection

! ************************************************************************** !

subroutine GridMapIndices(grid, dm_ptr, sgrid_stencil_type,option)
  ! 
  ! maps global, local and natural indices of cells
  ! to each other
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

#include "petsc/finclude/petscdm.h"
  use petscdm

  use Option_module
  use DM_Kludge_module

  implicit none
  
  type(grid_type) :: grid
  type(dm_ptr_type) :: dm_ptr
  PetscEnum :: sgrid_stencil_type
  type(option_type) :: option

  PetscInt, allocatable :: int_tmp(:)
! PetscInt, pointer :: int_tmp(:)
  PetscInt :: n
  PetscOffset :: i_da
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridMapIndices(grid%structured_grid,sgrid_stencil_type, &
                                grid%nG2L,grid%nL2G,grid%nG2A, &
                                option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID, &
         POLYHEDRA_UNSTRUCTURED_GRID)
      call UGridMapIndices(grid%unstructured_grid, &
                           dm_ptr%ugdm, &
                           grid%nG2L,grid%nL2G,grid%nG2A,option)
  end select
 
 
end subroutine GridMapIndices

! ************************************************************************** !

subroutine GridComputeSpacing(grid,origin_global,option)
  ! 
  ! Computes grid spacing (only for structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 
  use Option_module

  implicit none
  
  type(grid_type) :: grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridComputeSpacing(grid%structured_grid,origin_global,option)
    case(IMPLICIT_UNSTRUCTURED_GRID)
  end select
  
end subroutine GridComputeSpacing

! ************************************************************************** !

subroutine GridComputeCoordinates(grid,origin_global,option,ugdm)
  ! 
  ! Computes x,y,z coordinates of grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  use Option_module
  use Grid_Unstructured_Explicit_module
  use Grid_Unstructured_Polyhedra_module
  
  implicit none

  type(grid_type) :: grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  type(ugdm_type), optional :: ugdm ! sp 
  PetscInt :: icell
  
  PetscErrorCode :: ierr  

  allocate(grid%x(grid%ngmax))
  grid%x = 0.d0
  allocate(grid%y(grid%ngmax))
  grid%y = 0.d0
  allocate(grid%z(grid%ngmax))
  grid%z = 0.d0
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridComputeCoord(grid%structured_grid,option, &
                                      origin_global, &
                                      grid%x,grid%y,grid%z, &
                                      grid%x_min_local,grid%x_max_local, &
                                      grid%y_min_local,grid%y_max_local, &
                                      grid%z_min_local,grid%z_max_local)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeCoord(grid%unstructured_grid,option, &
                             grid%x,grid%y,grid%z, &
                             grid%x_min_local,grid%x_max_local, &
                             grid%y_min_local,grid%y_max_local, &
                             grid%z_min_local,grid%z_max_local)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      call UGridExplicitSetCellCentroids(grid%unstructured_grid% &
                                         explicit_grid, &
                                         grid%x,grid%y,grid%z, &
                             grid%x_min_local,grid%x_max_local, &
                             grid%y_min_local,grid%y_max_local, &
                             grid%z_min_local,grid%z_max_local)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call UGridPolyhedraSetCellCentroids(grid%unstructured_grid%polyhedra_grid, &
                                          grid%x,grid%y,grid%z, &
                                          grid%x_min_local,grid%x_max_local, &
                                          grid%y_min_local,grid%y_max_local, &
                                          grid%z_min_local,grid%z_max_local,option)

  end select

  if (associated(grid%structured_grid)) then
    ! compute global max/min from the local max/in
    call MPI_Allreduce(grid%x_min_local,grid%x_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%y_min_local,grid%y_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%z_min_local,grid%z_min_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
    call MPI_Allreduce(grid%x_max_local,grid%x_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    call MPI_Allreduce(grid%y_max_local,grid%y_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    call MPI_Allreduce(grid%z_max_local,grid%z_max_global,ONE_INTEGER_MPI, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
 endif

  if (associated(grid%unstructured_grid)) then
     ! compute global max/min from the local max/in
     call MPI_Allreduce(grid%x_min_local,grid%x_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%y_min_local,grid%y_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%z_min_local,grid%z_min_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
     call MPI_Allreduce(grid%x_max_local,grid%x_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
     call MPI_Allreduce(grid%y_max_local,grid%y_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
     call MPI_Allreduce(grid%z_max_local,grid%z_max_global,ONE_INTEGER_MPI, &
                        MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
   !endif
 endif

end subroutine GridComputeCoordinates

! ************************************************************************** !

subroutine GridComputeVolumes(grid,volume,option)
  ! 
  ! Computes the volumes of cells in structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Grid_Unstructured_Explicit_module
  use Grid_Unstructured_Polyhedra_module
  
  implicit none

  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: volume
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridComputeVolumes(grid%x,grid%structured_grid,option, &
                                        grid%nL2G,volume)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeVolumes(grid%unstructured_grid,option,volume)
      call UGridComputeQuality(grid%unstructured_grid,option)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      call UGridExplicitComputeVolumes(grid%unstructured_grid, &
                                       option,volume)
    case(POLYHEDRA_UNSTRUCTURED_GRID)
      call UGridPolyhedraComputeVolumes(grid%unstructured_grid,option,volume)
  end select

end subroutine GridComputeVolumes

! ************************************************************************** !

subroutine GridComputeAreas(grid,area,option)
  ! 
  ! Computes the areas for 2D-mesh
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/07/2012
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: area
  
  select case(grid%itype)
    case(IMPLICIT_UNSTRUCTURED_GRID)
      call UGridComputeAreas(grid%unstructured_grid,option,area)
      call UGridComputeQuality(grid%unstructured_grid,option)
    case default
      option%io_buffer = 'ERROR: GridComputeAreas only implemented for Unstructured grid'
      call printErrMsg(option)
  end select

end subroutine GridComputeAreas

! ************************************************************************** !

subroutine GridLocalizeRegions(grid,region_list,option)
  ! 
  ! Resticts regions to cells local to processor
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

  use Option_module
  use Region_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  type(region_type), pointer :: region
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag, global_cell_count
  PetscBool :: same_point
  PetscErrorCode :: ierr
  
  iflag = 0
  region => region_list%first
  do
    if (.not.associated(region)) exit

    select case(region%def_type)
      case (DEFINED_BY_BLOCK)
        call GridLocalizeRegionFromBlock(grid,region,option)
      case (DEFINED_BY_CARTESIAN_BOUNDARY)
        call GridLocalizeRegionFromCartBound(grid,region,option)
      case (DEFINED_BY_COORD)
        call GridLocalizeRegionFromCoordinates(grid,region,option)
      case (DEFINED_BY_CELL_IDS)
        select case(grid%itype)
          case(STRUCTURED_GRID)
            call GridLocalizeRegionsFromCellIDs(grid,region,option)
          case(IMPLICIT_UNSTRUCTURED_GRID)
            call GridLocalizeRegionsFromCellIDs(grid,region,option)
          case(EXPLICIT_UNSTRUCTURED_GRID)
            call GridLocalizeRegionsFromCellIDs(grid,region,option)
        end select
      case (DEFINED_BY_CELL_AND_FACE_IDS)
        select case(grid%itype)
          case (STRUCTURED_GRID)
            call GridLocalizeRegionsFromCellIDs(grid,region,option)
          case default
            option%io_buffer = 'GridLocalizeRegions() must tbe extended ' // &
            'for unstructured region DEFINED_BY_CELL_AND_FACE_IDS'
            call printErrMsg(option)
        end select
      case (DEFINED_BY_VERTEX_IDS)
        option%io_buffer = 'GridLocalizeRegions() must tbe extended ' // &
          'for unstructured region DEFINED_BY_VERTEX_IDS'
        call printErrMsg(option)
      case (DEFINED_BY_SIDESET_UGRID)
        call UGridMapSideSet2(grid%unstructured_grid, &
                             region%sideset%face_vertices, &
                             region%sideset%nfaces,region%name, &
                             option,region%cell_ids,region%faces)
        region%num_cells = size(region%cell_ids)
      case (DEFINED_BY_FACE_UGRID_EXP)
          call GridLocalizeExplicitFaceset(grid%unstructured_grid,region, &
                                           option)
      case (DEFINED_BY_POLY_BOUNDARY_FACE)
        call UGridMapBoundFacesInPolVol(grid%unstructured_grid, &
                                        region%polygonal_volume, &
                                        region%name,option, &
                                        region%cell_ids,region%faces)
        region%num_cells = size(region%cell_ids)
      case (DEFINED_BY_POLY_CELL_CENTER)
        call GridMapCellsInPolVol(grid, &
                                  region%polygonal_volume, &
                                  region%name,option, &
                                  region%cell_ids)
        region%num_cells = size(region%cell_ids)
      case default
        option%io_buffer = 'GridLocalizeRegions: Region definition not recognized'
        call printErrMsg(option)
    end select

    if (region%num_cells == 0 .and. associated(region%cell_ids)) then
      deallocate(region%cell_ids)
      nullify(region%cell_ids)
    endif

    if (region%num_cells == 0 .and. associated(region%faces)) then
      deallocate(region%faces)
      nullify(region%faces)
    endif

    ! check to ensure that there is at least one grid cell in each region
    call MPI_Allreduce(region%num_cells,global_cell_count,ONE_INTEGER_MPI, &
                       MPI_INTEGER, MPI_SUM, option%mycomm,ierr)
    if (global_cell_count == 0) then
      option%io_buffer = 'No cells assigned to REGION "' // &
        trim(region%name) // '".'
      call printErrMsg(option)
    endif
    region => region%next

  enddo

end subroutine GridLocalizeRegions

! ************************************************************************** !

subroutine GridLocalizeRegionsFromCellIDs(grid, region, option)
  ! 
  ! Redistributed cells ids in a grid based on natural numbering to their
  ! respective global index (PETSc ordering). Sets face information too.
  ! 
  ! Author: Gautam Bisht, Glenn Hammond
  ! Date: 5/30/2011, 09/14/16

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Region_module
  use Utility_module

  implicit none

  type(grid_type) :: grid
  type(region_type) :: region
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  Vec :: vec_cell_ids
  Vec :: vec_cell_ids_loc
  PetscInt, allocatable :: tmp_int_array(:)
  PetscReal, allocatable :: tmp_scl_array(:)
  PetscInt :: count
  PetscInt :: ii
  PetscInt :: local_id, natural_id
  PetscInt :: tempint
  PetscInt :: iface
  PetscInt :: tempfacearray(20)
  PetscInt :: facecount
  PetscInt :: cell_id_max_local
  PetscInt :: cell_id_max_global
  PetscReal :: tempreal, prevreal, facereal
  PetscReal, parameter :: offset = 0.1d0
  PetscReal, pointer :: v_loc_p(:)
  PetscBool :: setup_faces
  IS :: is_from, is_to
  VecScatter :: vec_scat
  PetscInt :: istart, iend
  PetscErrorCode :: ierr

  call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DECIDE, &
                    vec_cell_ids, ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DECIDE, &
                    vec_cell_ids_loc, ierr);CHKERRQ(ierr)
    
  call VecZeroEntries(vec_cell_ids, ierr);CHKERRQ(ierr)
    
  allocate(tmp_int_array(region%num_cells))
  allocate(tmp_scl_array(region%num_cells))

  cell_id_max_local = -1
  do ii = 1, region%num_cells
    tmp_int_array(ii) = region%cell_ids(ii) - 1
    if (associated(region%faces)) then
      iface = region%faces(ii)
      if (iface > 14) then
        write(iface,*) iface
        option%io_buffer = 'Face ID (' // trim(adjustl(word)) // ') greater &
          &than 14 in GridLocalizeRegionsFromCellIDs()'
        call printErrMsg(option)
      endif
      tmp_scl_array(ii) = 10.d0**dble(region%faces(ii))
    else
      ! use 0.1 as it will take 100 connections to conflict with faces
      tmp_scl_array(ii) = 0.1d0
    endif
    cell_id_max_local = max(cell_id_max_local, region%cell_ids(ii))
  enddo

  call MPI_Allreduce(cell_id_max_local, cell_id_max_global, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_MAX, option%mycomm,ierr)
  if (cell_id_max_global > grid%nmax) then
    option%io_buffer = 'The following region includes a cell-id that is &
      &greater than number of control volumes present in the grid: ' // &
    trim(region%name)
    call printErrMsg(option)
  endif

  call VecSetValues(vec_cell_ids, region%num_cells, tmp_int_array, &
                    tmp_scl_array, ADD_VALUES, ierr);CHKERRQ(ierr)
  
  deallocate(tmp_int_array)
  deallocate(tmp_scl_array)

  call VecAssemblyBegin(vec_cell_ids, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(vec_cell_ids, ierr);CHKERRQ(ierr)

  ! create list of natural ids for local, on-process cells
  allocate(tmp_int_array(grid%nlmax))
  do local_id = 1, grid%nlmax
    natural_id = grid%nG2A(grid%nL2G(local_id))
    tmp_int_array(local_id) = natural_id
  enddo
  
  tmp_int_array = tmp_int_array - 1
  call ISCreateBlock(option%mycomm, 1, grid%nlmax, &
                     tmp_int_array, PETSC_COPY_VALUES, is_from,  &
                     ierr);CHKERRQ(ierr)
  
  call VecGetOwnershipRange(vec_cell_ids_loc,istart,iend,ierr);CHKERRQ(ierr)
  do ii=1,grid%nlmax
    tmp_int_array(ii) = ii + istart
  enddo

  tmp_int_array = tmp_int_array - 1
  call ISCreateBlock(option%mycomm, 1, grid%nlmax, &
                      tmp_int_array, PETSC_COPY_VALUES, is_to,  &
                     ierr);CHKERRQ(ierr)
  deallocate(tmp_int_array)
  
  call VecScatterCreate(vec_cell_ids,is_from,vec_cell_ids_loc,is_to, &
                        vec_scat, ierr);CHKERRQ(ierr)
  call ISDestroy(is_from, ierr);CHKERRQ(ierr)
  call ISDestroy(is_to, ierr);CHKERRQ(ierr)
  
  call VecScatterBegin(vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                       INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scat, vec_cell_ids, vec_cell_ids_loc, &
                     INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scat, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(vec_cell_ids_loc, v_loc_p, ierr);CHKERRQ(ierr)
  count = 0
  setup_faces = PETSC_FALSE
  do local_id=1, grid%nlmax
    if (v_loc_p(local_id) > 9.99d0) then
      setup_faces = PETSC_TRUE
      facereal = v_loc_p(local_id)
      tempint = int(log10(facereal))
      prevreal = mod(facereal,10.d0**dble(tempint+1)) + offset
      do ii = tempint+1, 1, -1
        tempreal = mod(facereal,10.d0**dble(ii)) + offset
        if (.not.Equal(tempreal,prevreal)) then
          count = count + 1
        endif
        prevreal = tempreal
      enddo
    elseif (v_loc_p(local_id) > 0.0999d0) then
      count = count + 1
    endif
  enddo
    
  region%num_cells = count
  call DeallocateArray(region%cell_ids)
  call DeallocateArray(region%faces)
  if (region%num_cells > 0) then
    allocate(region%cell_ids(region%num_cells))
    region%cell_ids = UNINITIALIZED_INTEGER
    if (setup_faces) then
      allocate(region%faces(region%num_cells))
      region%faces = UNINITIALIZED_INTEGER
    endif
    count = 0
    do local_id=1, grid%nlmax
      if (v_loc_p(local_id) > 9.99d0) then
        facereal = v_loc_p(local_id)
        tempint = int(log10(facereal))
        prevreal = mod(facereal,10.d0**dble(tempint+1)) + offset
        ! use a temporary array tempfacearray so that faces can be sorted
        facecount = 0
        tempfacearray = UNINITIALIZED_INTEGER
        do ii = tempint+1, 1, -1
          tempreal = mod(facereal,10.d0**dble(ii)) + offset
          if (.not.Equal(tempreal,prevreal)) then
            facecount = facecount + 1
            tempfacearray(facecount) = ii
          endif
          prevreal = tempreal
        enddo
        ! sort the faces
        call PetscSortInt(facecount,tempfacearray,ierr);CHKERRQ(ierr)
        region%cell_ids(count+1:count+facecount) = local_id
        region%faces(count+1:count+facecount) = tempfacearray(1:facecount)
        count = count + facecount
      elseif (v_loc_p(local_id) > 0.0999d0) then
        count = count + 1
        region%cell_ids(count) = local_id
      endif
    enddo
  endif
  
  call VecRestoreArrayF90(vec_cell_ids_loc,v_loc_p,ierr);CHKERRQ(ierr)
  
  call VecDestroy(vec_cell_ids,ierr);CHKERRQ(ierr)
  call VecDestroy(vec_cell_ids_loc,ierr);CHKERRQ(ierr)

end subroutine GridLocalizeRegionsFromCellIDs

! ************************************************************************** !

subroutine GridLocalizeExplicitFaceset(ugrid,region,option)
  ! 
  ! GridLocalizeExplicitFaceset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/10/12
  ! 

  use Region_module
  use Option_module

  implicit none
  
  type(grid_unstructured_type) :: ugrid
  type(region_type) :: region
  type(option_type) :: option
  Vec :: volume

  type(unstructured_explicit_type), pointer :: explicit_grid
  type(region_explicit_face_type), pointer :: faceset
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icell, count
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: real_array_2d(:,:)
  PetscErrorCode :: ierr

  explicit_grid => ugrid%explicit_grid
  faceset => region%explicit_faceset
  
  ! convert ids to petsc
  region%cell_ids = region%cell_ids - 1
  call AOApplicationToPetsc(ugrid%ao_natural_to_petsc,size(region%cell_ids), &
                            region%cell_ids,ierr);CHKERRQ(ierr)
  region%cell_ids = region%cell_ids + 1
  
  ! if petsc ids are below global_offset or above global_offset + nlmax, 
  ! they are off processor; otherwise, local

  allocate(int_array(size(region%cell_ids)))
  ! negate off processor ids
  int_array = UNINITIALIZED_INTEGER
  ! only if faceset exists
  if (associated(faceset)) then
    allocate(real_array_2d(4,size(region%cell_ids)))
    real_array_2d = UNINITIALIZED_DOUBLE
  endif
  count = 0
  do icell = 1, size(region%cell_ids)
    if (region%cell_ids(icell) > ugrid%global_offset .and. &
        region%cell_ids(icell) <= ugrid%global_offset + ugrid%nlmax) then
      count = count + 1
      ! local cell id
      int_array(count) = region%cell_ids(icell) - ugrid%global_offset
      if (associated(faceset)) then
        real_array_2d(1,count) = faceset%face_centroids(icell)%x
        real_array_2d(2,count) = faceset%face_centroids(icell)%y
        real_array_2d(3,count) = faceset%face_centroids(icell)%z
        real_array_2d(4,count) = faceset%face_areas(icell)
      endif
    endif
  enddo
  
  deallocate(region%cell_ids)
  allocate(region%cell_ids(count))
  region%cell_ids = int_array(1:count)
  deallocate(int_array)

  if (associated(faceset)) then
    deallocate(faceset%face_centroids)
    deallocate(faceset%face_areas)
    allocate(faceset%face_centroids(count))
    allocate(faceset%face_areas(count))
  
    do icell = 1, count
      faceset%face_centroids(icell)%x = real_array_2d(1,icell)
      faceset%face_centroids(icell)%y = real_array_2d(2,icell)
      faceset%face_centroids(icell)%z = real_array_2d(3,icell)
      faceset%face_areas(icell) = real_array_2d(4,icell)
    enddo
    deallocate(real_array_2d)
  endif
  

  region%num_cells = count
  
  if (region%num_cells == 0) then
    deallocate(region%cell_ids)
    nullify(region%cell_ids)
    if (associated(faceset)) then
      deallocate(faceset%face_centroids)
      nullify(faceset%face_centroids)
      deallocate(faceset%face_areas)
      nullify(faceset%face_areas)
      ! note that have to use full reference
      deallocate(region%explicit_faceset)
      nullify(region%explicit_faceset)
    endif
  endif

#if UGRID_DEBUG
  if (region%num_cells > 0) then
    write(string,*) option%myrank
    string = 'region_faceset_' // trim(region%name) // trim(adjustl(string)) // '.out'
    open(unit=86,file=trim(string))
    if (associated(faceset)) then
      do icell = 1, region%num_cells
        write(86,'(i5,4f7.3)') region%cell_ids(icell), &
                    faceset%face_centroids(icell)%x, &
                    faceset%face_centroids(icell)%y, &
                    faceset%face_centroids(icell)%z, &
                    faceset%face_areas(icell)
      enddo
    else
      do icell = 1, region%num_cells
        write(86,'(i5)') region%cell_ids(icell)
      enddo
    endif
    close(86)
  endif
#endif  

end subroutine GridLocalizeExplicitFaceset

! ************************************************************************** !

subroutine GridCopyIntegerArrayToVec(grid, array,vector,num_values)
  ! 
  ! Copies values from an integer array into a
  ! PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  type(grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90( vector,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90( vector,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine GridCopyIntegerArrayToVec

! ************************************************************************** !

subroutine GridCopyRealArrayToVec(grid,array,vector,num_values)
  ! 
  ! Copies values from an integer array into a
  ! PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none
    
  type(grid_type) :: grid
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine GridCopyRealArrayToVec

! ************************************************************************** !

subroutine GridCopyVecToIntegerArray(grid,array,vector,num_values)
  ! 
  ! Copies values from a PETSc Vec to an
  ! integer array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none
  
  type(grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  do i=1,num_values
    if (vec_ptr(i) > 0.d0) then
      array(i) = int(vec_ptr(i)+1.d-4)
    else
      array(i) = int(vec_ptr(i)-1.d-4)
    endif
  enddo
  call VecRestoreArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine GridCopyVecToIntegerArray

! ************************************************************************** !

subroutine GridCopyVecToRealArray(grid,array,vector,num_values)
  ! 
  ! Copies values from a PETSc Vec to an integer
  ! array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none
    
  type(grid_type) :: grid
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  array(1:num_values) = vec_ptr(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine GridCopyVecToRealArray

! ************************************************************************** !

subroutine GridCreateNaturalToGhostedHash(grid,option)
  ! 
  ! Creates a hash table for looking up the
  ! local ghosted id of a natural id, if it
  ! exists
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  use Option_module
  use Logging_module  
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_ghosted_id, natural_id
  PetscInt :: num_in_hash, num_ids_per_hash, hash_id, id, ierr, hash_id_2
  PetscInt :: max_num_ids_per_hash
  PetscInt, pointer :: hash(:,:,:), temp_hash(:,:,:)

  if (associated(grid%hash)) return

  call PetscLogEventBegin(logging%event_hash_create,ierr);CHKERRQ(ierr)
                          
  max_num_ids_per_hash = 0
  ! initial guess of 10% of ids per hash
  ! must be at least 5 so that reallocation (*1.2) works below
  num_ids_per_hash = max(grid%nlmax/(grid%num_hash_bins/10),5)

  allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
  hash(:,:,:) = 0

  
  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nG2A(local_ghosted_id) !nG2A is 1-based
    hash_id = mod(natural_id,grid%num_hash_bins)+1 
    num_in_hash = hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    if (num_in_hash > max_num_ids_per_hash) max_num_ids_per_hash = num_in_hash
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,1:grid%num_hash_bins) = &
                             hash(1:2,0:num_ids_per_hash,1:grid%num_hash_bins)
      deallocate(hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old to new
      do hash_id_2 = 1, grid%num_hash_bins
        do id = 1, temp_hash(1,0,hash_id_2)
          hash(1:2,id,hash_id_2) = temp_hash(1:2,id,hash_id_2)
        enddo
        hash(1,0,hash_id_2) = temp_hash(1,0,hash_id_2)
      enddo
      deallocate(temp_hash)
    endif
    hash(1,0,hash_id) = num_in_hash
    hash(1,num_in_hash,hash_id) = natural_id
    hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

  grid%hash => hash
  
!  call GridPrintHashTable(grid)
  call MPI_Allreduce(max_num_ids_per_hash,num_in_hash,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  write(option%io_buffer,'("max_num_ids_per_hash: ",i5)') num_in_hash
  call printMsg(option)

  call PetscLogEventEnd(logging%event_hash_create,ierr);CHKERRQ(ierr)

end subroutine GridCreateNaturalToGhostedHash

! ************************************************************************** !

PetscInt function GridGetLocalIdFromNaturalId(grid,natural_id)
  ! 
  ! GetLocalIdFromNaturalId: Returns the local id corresponding to a natural
  ! id or 0, if the natural id is off-processor
  ! WARNING: Extremely inefficient for large jobs
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  implicit none

  type(grid_type) :: grid

  PetscInt :: natural_id, local_id
  
  do local_id = 1, grid%nlmax
    if (natural_id == grid%nG2A(grid%nL2G(local_id))) then
      GridGetLocalIdFromNaturalId = local_id
      return
    endif
  enddo
  GridGetLocalIdFromNaturalId = 0

end function GridGetLocalIdFromNaturalId

! ************************************************************************** !

PetscInt function GridGetLocalGhostedIdFromNatId(grid,natural_id)
  ! 
  ! Returns the local ghosted id corresponding
  ! to a natural id or 0, if the natural id
  ! is off-processor
  ! WARNING: Extremely inefficient for large jobs
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    !geh: nG2A is 1-based
    if (natural_id == grid%nG2A(local_ghosted_id)) then
      GridGetLocalGhostedIdFromNatId = local_ghosted_id
      return 
    endif
  enddo
  GridGetLocalGhostedIdFromNatId = 0

end function GridGetLocalGhostedIdFromNatId

! ************************************************************************** !

PetscInt function GridGetLocalGhostedIdFromHash(grid,natural_id)
  ! 
  ! Returns the local ghosted id of a natural
  ! id, if it exists.  Otherwise 0 is returned
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: hash_id, id

  GridGetLocalGhostedIdFromHash = 0
  hash_id = mod(natural_id,grid%num_hash_bins)+1 
  do id = 1, grid%hash(1,0,hash_id)
    if (grid%hash(1,id,hash_id) == natural_id) then
      GridGetLocalGhostedIdFromHash = grid%hash(2,id,hash_id)
      return
    endif
  enddo

end function GridGetLocalGhostedIdFromHash

! ************************************************************************** !

subroutine GridDestroyHashTable(grid)
  ! 
  ! Deallocates the hash table
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(grid_type), pointer :: grid
  
  call DeallocateArray(grid%hash)

  nullify(grid%hash)
  grid%num_hash_bins = 100

end subroutine GridDestroyHashTable

! ************************************************************************** !

subroutine GridPrintHashTable(grid)
  ! 
  ! UnstructGridPrintHashTable: Prints the hashtable for viewing
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/07
  ! 

  implicit none

  type(grid_type) :: grid
  
  PetscInt :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,grid%num_hash_bins
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         grid%hash(1,0,ihash), &
                         '):', &
                         (grid%hash(1,id,ihash),id=1,grid%hash(1,0,ihash))
  enddo
  close(fid)

end subroutine GridPrintHashTable

! ************************************************************************** !

subroutine GridGetGhostedNeighbors(grid,ghosted_id,stencil_type, &
                                   stencil_width_i,stencil_width_j, &
                                   stencil_width_k,x_count,y_count, &
                                   z_count, &
                                   ghosted_neighbors,option)
  ! 
  ! GridGetNeighbors: Returns an array of neighboring cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/11
  ! 

  use Option_module

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscEnum :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: x_count
  PetscInt :: y_count
  PetscInt :: z_count
  PetscInt :: ghosted_neighbors(*)
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridGetGhostedNeighbors(grid%structured_grid, &
                                         ghosted_id,stencil_type, &
                                         stencil_width_i, &
                                         stencil_width_j,stencil_width_k, &
                                         x_count,y_count,z_count, &
                                         ghosted_neighbors,option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridGetNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridGetGhostedNeighbors

! ************************************************************************** !

subroutine GridGetGhostedNeighborsWithCorners(grid,ghosted_id,stencil_type, &
                                   stencil_width_i,stencil_width_j, &
                                   stencil_width_k,icount, &
                                   ghosted_neighbors,option)
  ! 
  ! GridGetNeighborsWithCorners: Returns an array of neighboring cells along with corner
  ! cells
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/19/12
  ! 

  use Option_module

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscEnum :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: icount
  PetscInt :: ghosted_neighbors(*)
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridGetGhostedNeighborsCorners(grid%structured_grid, &
                                         ghosted_id,stencil_type, &
                                         stencil_width_i, &
                                         stencil_width_j, &
                                         stencil_width_k, &
                                         icount, &
                                         ghosted_neighbors,option)
    case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID) 
      option%io_buffer = 'GridGetNeighbors not currently supported for ' // &
        'unstructured grids.'
      call printErrMsg(option)
  end select

end subroutine GridGetGhostedNeighborsWithCorners

! ************************************************************************** !

subroutine GridDestroy(grid)
  ! 
  ! Deallocates a grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(grid_type), pointer :: grid
  PetscErrorCode :: ierr
  PetscInt :: ghosted_id
    
  if (.not.associated(grid)) return
      
  call DeallocateArray(grid%nL2G)
  call DeallocateArray(grid%nG2L)
  call DeallocateArray(grid%nG2A)

  call DeallocateArray(grid%x)
  call DeallocateArray(grid%y)
  call DeallocateArray(grid%z)
  
  if (associated(grid%hash)) call GridDestroyHashTable(grid)
  
  call UGridDestroy(grid%unstructured_grid)    
  call StructGridDestroy(grid%structured_grid)
                                           
  call ConnectionDestroyList(grid%internal_connection_set_list)

  if (associated(grid)) deallocate(grid)
  nullify(grid)

end subroutine GridDestroy

! ************************************************************************** !

function GridIndexToCellID(vec,index,grid,vec_type)
  ! 
  ! Returns the local grid cell id of a Petsc Vec index
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/09
  ! 

  implicit none
  
  Vec :: vec
  PetscInt :: index
  type(grid_type) :: grid
  PetscInt :: vec_type
  
  PetscInt :: GridIndexToCellID
  
  PetscInt :: low
  PetscInt :: high
  PetscInt :: ndof
  PetscInt :: cell_id
  PetscErrorCode :: ierr

  
  cell_id = -1
  call VecGetOwnershipRange(vec,low,high,ierr);CHKERRQ(ierr)
  call VecGetBlockSize(vec,ndof,ierr);CHKERRQ(ierr)
  if (index >= low .and. index < high) then
    cell_id = (index-low)/ndof+1
    if (vec_type == GLOBAL) then
      cell_id = grid%nG2A(grid%nL2G(cell_id))
    else if (vec_type == LOCAL) then
      cell_id = grid%nG2A(cell_id) !nG2A is 1-based
    endif
  endif
  
  call MPI_Allreduce(cell_id,GridIndexToCellID,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,PETSC_COMM_WORLD,ierr)
                     
end function GridIndexToCellID

! ************************************************************************** !

subroutine GridLocalizeRegionFromBlock(grid,region,option)
  ! 
  ! This routine resticts regions to cells local to processor when the region
  ! was defined using a BLOCK from inputfile.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/04/12
  ! 

  use Option_module
  use Region_module

  implicit none
  
  type(region_type), pointer :: region
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr

  if (grid%itype /= STRUCTURED_GRID) then
     option%io_buffer='Region definition using BLOCK is only supported for ' //&
       ' structured grids'
     call printErrMsg(option)
  endif
  
  ! convert indexing from global (entire domain) to local processor
  region%i1 = region%i1 - grid%structured_grid%lxs
  region%i2 = region%i2 - grid%structured_grid%lxs
  region%j1 = region%j1 - grid%structured_grid%lys
  region%j2 = region%j2 - grid%structured_grid%lys
  region%k1 = region%k1 - grid%structured_grid%lzs
  region%k2 = region%k2 - grid%structured_grid%lzs
          
  ! clip region to within local processor domain
  region%i1 = max(region%i1,1)
  region%i2 = min(region%i2,grid%structured_grid%nlx)
  region%j1 = max(region%j1,1)
  region%j2 = min(region%j2,grid%structured_grid%nly)
  region%k1 = max(region%k1,1)
  region%k2 = min(region%k2,grid%structured_grid%nlz)
   
  count = 0  
  if (region%i1 <= region%i2 .and. &
      region%j1 <= region%j2 .and. &
      region%k1 <= region%k2) then
    region%num_cells = (region%i2-region%i1+1)* &
                       (region%j2-region%j1+1)* &
                       (region%k2-region%k1+1)
    allocate(region%cell_ids(region%num_cells))
    if (region%iface /= 0) then
      allocate(region%faces(region%num_cells))
      region%faces = region%iface
    endif
    region%cell_ids = 0
    do k=region%k1,region%k2
      do j=region%j1,region%j2
        do i=region%i1,region%i2
          count = count + 1
          region%cell_ids(count) = &
                 i + (j-1)*grid%structured_grid%nlx + &
                 (k-1)*grid%structured_grid%nlxy
        enddo
      enddo
    enddo
!   if (region%num_cells > 0) then
!     region%coordinates(1)%x = grid%x(region%cell_ids(ONE_INTEGER))
!     region%coordinates(1)%y = grid%y(region%cell_ids(ONE_INTEGER))
!     region%coordinates(1)%z = grid%z(region%cell_ids(ONE_INTEGER))
!   endif
  else
    region%num_cells = 0
  endif

  if (count /= region%num_cells) then
    option%io_buffer = 'Mismatch in number of cells in block region'
    call printErrMsg(option)
  endif

end subroutine GridLocalizeRegionFromBlock

! ************************************************************************** !

subroutine GridLocalizeRegionFromCartBound(grid,region,option)
  ! 
  ! This routine resticts regions to cells local to processor when the region
  ! was defined using a BLOCK from inputfile.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/04/12
  ! 

  use Option_module
  use Region_module

  implicit none
  
  type(region_type), pointer :: region
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr

  if (grid%itype /= STRUCTURED_GRID) then
    option%io_buffer='Region definition using CARTESIAN_BOUNDARY is ' // &
      'only supported for structured grids.'
    call printErrMsg(option)
  endif

  region%i1 = 1
  region%i2 = grid%structured_grid%nx
  region%j1 = 1
  region%j2 = grid%structured_grid%ny
  region%k1 = 1
  region%k2 = grid%structured_grid%nz
  
  select case(region%iface)
    case(WEST_FACE)
      region%i2 = region%i1
    case(EAST_FACE)
      region%i1 = region%i2
    case(SOUTH_FACE)
      region%j2 = region%j1
    case(NORTH_FACE)
      region%j1 = region%j2
    case(BOTTOM_FACE)
      region%k2 = region%k1
    case(TOP_FACE)
      region%k1 = region%k2
  end select

  call GridLocalizeRegionFromBlock(grid,region,option)

end subroutine GridLocalizeRegionFromCartBound

! ************************************************************************** !

subroutine GridLocalizeRegionFromCoordinates(grid,region,option)
  ! 
  ! This routine resticts regions to cells local to processor when the region
  ! was defined using COORDINATES from inputfile.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/04/12
  ! 

  use Option_module
  use Region_module

  implicit none
  
  type(region_type), pointer :: region
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscReal :: del_x, del_y, del_z
  PetscInt :: iflag
  PetscBool :: same_point
  PetscErrorCode :: ierr

  iflag = 0
  if (size(region%coordinates) > TWO_INTEGER) then
    option%io_buffer = 'GridLocalizeRegions: more than 2 coordinates' // &
                       ' not supported in region object'
    call printErrMsg(option)
  endif

  same_point = PETSC_FALSE
  ! if two coordinates, determine whether they are the same point
  if (size(region%coordinates) == TWO_INTEGER) then
    if (dabs(region%coordinates(ONE_INTEGER)%x - &
             region%coordinates(TWO_INTEGER)%x) < tol .and. &
        dabs(region%coordinates(ONE_INTEGER)%y - &
             region%coordinates(TWO_INTEGER)%y) < tol .and. &
        dabs(region%coordinates(ONE_INTEGER)%z - &
             region%coordinates(TWO_INTEGER)%z) < tol) then
      same_point = PETSC_TRUE
    endif
  endif
   
  ! treat two identical coordinates the same as a single coordinate
  if (size(region%coordinates) == ONE_INTEGER .or. same_point) then
    call GridGetLocalIDFromCoordinate(grid,region%coordinates(ONE_INTEGER), &
                                      option,local_id)
    if (INITIALIZED(local_id)) then
      region%num_cells = 1
      allocate(region%cell_ids(region%num_cells))
      region%cell_ids(1) = local_id
      select case(grid%itype)
        case(STRUCTURED_GRID)
          if (region%iface /= 0) then
            allocate(region%faces(region%num_cells))
            region%faces = region%iface
          endif
      end select
    else
      region%num_cells = 0
    endif
    ! the next test as designed will only work on a uniform grid
    call MPI_Allreduce(region%num_cells,count, &
                        ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                        option%mycomm,ierr)   
    if (count == 0) then
      write(option%io_buffer,*) 'Region: (coord)', &
            region%coordinates(ONE_INTEGER)%x, &
            region%coordinates(ONE_INTEGER)%y, &
            region%coordinates(ONE_INTEGER)%z, &
            ' not found in global domain.', count
      call printErrMsg(option)
    else if (count > 1) then
      write(option%io_buffer,*) 'Region: (coord)', &
            region%coordinates(ONE_INTEGER)%x, &
            region%coordinates(ONE_INTEGER)%y, &
            region%coordinates(ONE_INTEGER)%z, &
            ' duplicated across ', count, &
            ' procs in global domain.'
      call printErrMsg(option)
    endif
  else ! 2 coordinates
    x_min = min(region%coordinates(ONE_INTEGER)%x, &
                region%coordinates(TWO_INTEGER)%x)
    x_max = max(region%coordinates(ONE_INTEGER)%x, &
                region%coordinates(TWO_INTEGER)%x)
    y_min = min(region%coordinates(ONE_INTEGER)%y, &
                region%coordinates(TWO_INTEGER)%y)
    y_max = max(region%coordinates(ONE_INTEGER)%y, &
                region%coordinates(TWO_INTEGER)%y)
    z_min = min(region%coordinates(ONE_INTEGER)%z, &
                region%coordinates(TWO_INTEGER)%z)
    z_max = max(region%coordinates(ONE_INTEGER)%z, &
                region%coordinates(TWO_INTEGER)%z)
                
    if (grid%itype == STRUCTURED_GRID) then
      ! shift box slightly inward
      x_shift = 1.d-8*(grid%x_max_global-grid%x_min_global)
      x_min = x_min+x_shift            
      x_max = x_max-x_shift
      y_shift = 1.d-8*(grid%y_max_global-grid%y_min_global)
      y_min = y_min+y_shift            
      y_max = y_max-y_shift
      z_shift = 1.d-8*(grid%z_max_global-grid%z_min_global)
      z_min = z_min+z_shift            
      z_max = z_max-z_shift
         
      ! if plane or line, ensure it is within the grid cells     
      if (x_max-x_min < 1.d-10) then
        x_max = region%coordinates(ONE_INTEGER)%x
        x_shift = 1.d-8*(grid%x_max_global-grid%x_min_global)
        if (region%iface == WEST_FACE) then
          x_max = x_max + x_shift
        else if (region%iface == EAST_FACE) then
          x_max = x_max - x_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (x_max > grid%x_min_global + x_shift) then
            x_max = x_max - x_shift
          else
            x_max = x_max + x_shift
          endif
        endif
        x_min = x_max
      endif
      if (y_max-y_min < 1.d-10) then
        y_max = region%coordinates(ONE_INTEGER)%y
        y_shift = 1.d-8*(grid%y_max_global-grid%y_min_global)
        if (region%iface == SOUTH_FACE) then
          y_max = y_max + y_shift
        else if (region%iface == NORTH_FACE) then
          y_max = y_max - y_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (y_max > grid%y_min_global + y_shift) then
            y_max = y_max - y_shift
          else
            y_max = y_max + y_shift
          endif
        endif
        y_min = y_max
      endif
      if (z_max-z_min < 1.d-10) then
        z_max = region%coordinates(ONE_INTEGER)%z
        z_shift = 1.d-8*(grid%z_max_global-grid%z_min_global)
        if (region%iface == BOTTOM_FACE) then
          z_max = z_max + z_shift
        else if (region%iface == TOP_FACE) then
          z_max = z_max - z_shift
        ! otherwise, shift upwind, unless at upwind physical boundary
        else
          if (z_max > grid%z_min_global + z_shift) then
            z_max = z_max - z_shift
          else
            z_max = z_max + z_shift
          endif
        endif
        z_min = z_max
      endif
    endif   
             
    ! ensure overlap
    if (x_min <= grid%x_max_local .and. &
        x_max >= grid%x_min_local .and. &
        y_min <= grid%y_max_local .and. &
        y_max >= grid%y_min_local .and. &
        z_min <= grid%z_max_local .and. &
        z_max >= grid%z_min_local) then
        
      ! get I,J,K bounds
      select case(grid%itype)
        case(STRUCTURED_GRID)
          ! local, non-ghosted i,j,k's are returned
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                    max(x_min,grid%x_min_local+x_shift), &
                                    max(y_min,grid%y_min_local+y_shift), &
                                    max(z_min,grid%z_min_local+z_shift), &
                                              i_min,j_min,k_min)
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                    min(x_max,grid%x_max_local-x_shift), &
                                    min(y_max,grid%y_max_local-y_shift), &
                                    min(z_max,grid%z_max_local-z_shift), &
                                              i_max,j_max,k_max)
          if (i_min > 0 .and. j_min > 0 .and. k_min > 0 .and. &
              i_max > 0 .and. j_max > 0 .and. k_max > 0) then
            region%num_cells = (i_max-i_min+1)*(j_max-j_min+1)*(k_max-k_min+1)
            allocate(region%cell_ids(region%num_cells))
            if (region%iface /= 0) then
              allocate(region%faces(region%num_cells))
              region%faces = region%iface
            endif
            region%cell_ids = 0
            count = 0
            do k = k_min, k_max
              do j = j_min, j_max
                do i = i_min, i_max
                  count = count+1
                  region%cell_ids(count) = i + (j-1)*grid%structured_grid%nlx + &
                                      (k-1)*grid%structured_grid%nlxy
                enddo
              enddo
            enddo
          else
            iflag = 1
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID,EXPLICIT_UNSTRUCTURED_GRID, &
             POLYHEDRA_UNSTRUCTURED_GRID)
          del_x = x_max-x_min
          del_y = y_max-y_min
          del_z = z_max-z_min
          ! 3D box
          if (del_x > 1.d-10 .and. &
              del_y > 1.d-10 .and. &
              del_z > 1.d-10) then
            ! geh: if the coordinates define a 3D box, add all cell centers
            ! that reside within the box
            count = 0
            do local_id = 1, grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (grid%x(ghosted_id) >= x_min .and. &
                  grid%x(ghosted_id) <= x_max .and. &
                  grid%y(ghosted_id) >= y_min .and. &
                  grid%y(ghosted_id) <= y_max .and. &
                  grid%z(ghosted_id) >= z_min .and. &
                  grid%z(ghosted_id) <= z_max) then
                count = count + 1
              endif
            enddo
            allocate(region%cell_ids(count))
            region%cell_ids = 0
            count = 0
            do local_id = 1, grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (grid%x(ghosted_id) >= x_min .and. &
                  grid%x(ghosted_id) <= x_max .and. &
                  grid%y(ghosted_id) >= y_min .and. &
                  grid%y(ghosted_id) <= y_max .and. &
                  grid%z(ghosted_id) >= z_min .and. &
                  grid%z(ghosted_id) <= z_max) then
                count = count + 1
                region%cell_ids(count) = local_id
              endif
            enddo
            region%num_cells = count
          ! 2D plane
          elseif ((del_x < 1.d-10 .and. del_y > 1.d-10 .and. &
                   del_z > 1.d-10) .or. &
                  (del_x > 1.d-10 .and. del_y < 1.d-10 .and. &
                   del_z > 1.d-10) .or. &
                  (del_x > 1.d-10 .and. del_y > 1.d-10 .and. &
                   del_z < 1.d-10)) then
            if (grid%itype == EXPLICIT_UNSTRUCTURED_GRID) then
              option%io_buffer = 'Regions defined with 2D planes are not ' // &
                'supported with explicit unstructured grids.'
              call printErrMsg(option)
            endif
            if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
              call UGridGetCellsInRectangle(x_min,x_max,y_min,y_max, &
                                            z_min,z_max, &
                                            grid%unstructured_grid,option, &
                                            region%num_cells,region%cell_ids, &
                                            region%faces)
            else
              call UGridPolyhedraGetCellsInRectangle(x_min,x_max,y_min,y_max, &
                                             z_min,z_max, &
                                             grid%unstructured_grid,option, &
                                             region%num_cells,region%cell_ids, &
                                             region%faces)
            endif

          endif
      end select
    endif

    call MPI_Allreduce(iflag,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
    iflag = i
    if (iflag > 0) then
      option%io_buffer = 'GridLocalizeRegions, between two points'
      call printErrMsg(option)
    endif
  endif

end subroutine GridLocalizeRegionFromCoordinates

! ************************************************************************** !

subroutine GridMapCellsInPolVol(grid,polygonal_volume, &
                                region_name,option,cell_ids)
  ! 
  ! Maps all global boundary cells within a polygonal volume to a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/16/15
  ! 
  use Option_module
  use Geometry_module

  implicit none

  type(grid_type) :: grid
  type(polygonal_volume_type) :: polygonal_volume
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)

  PetscInt :: local_id, ghosted_id, icount
  PetscBool :: found
  PetscInt, allocatable :: temp_int(:)
  
  allocate(temp_int(grid%nlmax))
  temp_int = 0
  icount = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    found = GeometryPointInPolygonalVolume(grid%x(ghosted_id), &
                                           grid%y(ghosted_id), &
                                           grid%z(ghosted_id), &
                                           polygonal_volume,option)
    if (found) then
      icount = icount + 1
      temp_int(icount) = local_id  
    endif
  enddo
  allocate(cell_ids(icount))
  cell_ids = temp_int(1:icount)
  deallocate(temp_int)
  
end subroutine GridMapCellsInPolVol

! ************************************************************************** !

subroutine GridGetLocalIDFromCoordinate(grid,coordinate,option,local_id)
  ! 
  ! Returns the local id of the grid cell occupied by a coordinate
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/16/15
  ! 
  use Option_module
  use Geometry_module

  implicit none

  type(grid_type) :: grid
  type(point3d_type) :: coordinate
  type(option_type) :: option
  PetscInt :: local_id
  
  PetscReal, parameter :: pert = 1.d-8, tol = 1.d-20
  PetscReal :: x_shift, y_shift, z_shift
  PetscInt :: i, j, k

  local_id = UNINITIALIZED_INTEGER
  if (coordinate%x >= grid%x_min_global .and. &
      coordinate%x <= grid%x_max_global .and. &
      coordinate%y >= grid%y_min_global .and. &
      coordinate%y <= grid%y_max_global .and. &
      coordinate%z >= grid%z_min_global .and. &
      coordinate%z <= grid%z_max_global) then
    ! If a point is on the corner of 4 or 8 patches in AMR, the region
    ! will be assigned to all 4/8 patches...a problem.  To avoid this, 
    ! we are going to perturb all point coordinates slightly upwind, as
    ! long as they are not on a global boundary (i.e. boundary condition)
    ! -- shift the coorindate slightly upwind
    x_shift = coordinate%x - &
              pert*(grid%x_max_global-grid%x_min_global)
    y_shift = coordinate%y - &
              pert*(grid%y_max_global-grid%y_min_global)
    z_shift = coordinate%z - &
              pert*(grid%z_max_global-grid%z_min_global)
    ! if the coodinate is shifted out of the global domain or 
    ! onto an exterior edge, set it back to the original value
    if (x_shift - grid%x_min_global < tol) &
      x_shift = coordinate%x
    if (y_shift - grid%y_min_global < tol) &
      y_shift = coordinate%y
    if (z_shift - grid%z_min_global < tol) &
      z_shift = coordinate%z
    select case(grid%itype)
      case(STRUCTURED_GRID)
        call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                            x_shift,y_shift,z_shift, &
                                            i,j,k)
        if (i > 0 .and. j > 0 .and. k > 0) then
          local_id = i + (j-1)*grid%structured_grid%nlx + &
                      (k-1)*grid%structured_grid%nlxy
        endif
      case(IMPLICIT_UNSTRUCTURED_GRID)
        !geh: must check each cell individually
        call UGridGetCellFromPoint(coordinate%x, &
                                   coordinate%y, &
                                   coordinate%z, &
                                   grid%unstructured_grid,option,local_id)
      case(EXPLICIT_UNSTRUCTURED_GRID)
        if (grid%itype == EXPLICIT_UNSTRUCTURED_GRID) then
          option%io_buffer = 'Locating a grid cell through a specified &
            &coordinate (GridGetLocalIDFromCoordinate)is not supported for &
            &explicit (primal) unstructured grids.'
          call printErrMsg(option)
        endif
      case(POLYHEDRA_UNSTRUCTURED_GRID)
          option%io_buffer = &
            'add code POLYHDERA in GridGetLocalIDFromCoordinate'
          call printErrMsg(option)
    end select
  endif
  
  ! Several of the subroutines above may return local_id = 0, therefore, we 
  ! need to ensure that local_id is reset back to uninitiazed if 0
  if (local_id <= 0) then
    local_id = UNINITIALIZED_INTEGER
  endif
    
end subroutine GridGetLocalIDFromCoordinate

! ************************************************************************** !

subroutine GridPrintExtents(grid,option)
  ! 
  ! Prints the extents of the gridded domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/18
  ! 
  use Option_module

  implicit none

  type(grid_type) :: grid
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word1, word2

  write(word1,*) grid%x_min_global
  write(word2,*) grid%x_max_global
  write(string,*) 'X: ', trim(adjustl(word1)), ' - ', trim(adjustl(word2))
  if (OptionPrintToScreen(option)) then
    write(*,*) 'Extent of Gridded Domain'
    write(*,*) trim(string)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,*) 'Extent of Gridded Domain'
    write(option%fid_out,*) trim(string)
  endif
  write(word1,*) grid%y_min_global
  write(word2,*) grid%y_max_global
  write(string,*) 'Y: ', trim(adjustl(word1)), ' - ', trim(adjustl(word2))
  if (OptionPrintToScreen(option)) then
    write(*,*) trim(string)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,*) trim(string)
  endif
  write(word1,*) grid%z_min_global
  write(word2,*) grid%z_max_global
  write(string,*) 'Z: ', trim(adjustl(word1)), ' - ', trim(adjustl(word2))
  if (OptionPrintToScreen(option)) then
    write(*,*) trim(string)
    write(*,*)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,*) trim(string)
    write(option%fid_out,*)
  endif

end subroutine GridPrintExtents

end module Grid_module
