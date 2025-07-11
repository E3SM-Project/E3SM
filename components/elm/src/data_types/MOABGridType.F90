#ifdef HAVE_MOAB
! Only build the module if MOAB is enabled
module MOABGridType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gridcell data type allocation
  ! --------------------------------------------------------
  ! gridcell types can have values of
  ! --------------------------------------------------------
  !   1 => default
  !
  ! PET: 9 Feb 2015: Preparing to change the sub-grid hierarchy to include
  !    topographic units between gridcell and landunit.
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_kind_mod     , only : CXX => SHR_KIND_CXX
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdmod  , only: iam, masterproc, mpicom
  use elm_varctl  ,  only : iulog, fatmlndfrc  ! for messages and domain file name
  use abortutils , only : endrun

  use iMOAB  , only: iMOAB_LoadMesh, iMOAB_WriteMesh, iMOAB_RegisterApplication, &
  iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage, iMOAB_SynchronizeTags, &
  iMOAB_UpdateMeshInfo, iMOAB_GetMeshInfo, &
  iMOAB_DetermineGhostEntities, iMOAB_WriteLocalMesh, iMOAB_GetVisibleElementsInfo, &
  iMOAB_GetNeighborElements, iMOAB_GetElementConnectivity, iMOAB_GetDoubleTagStorage

  use mpi
  use iso_c_binding

  !
  ! !PUBLIC TYPES:
  implicit none
  ! save
  private
  !
  public :: elm_moab_initialize
  public :: elm_moab_load_grid_file
  public :: elm_moab_finalize
  !
  integer :: mlndghostid     ! ID of the MOAB ELM application with ghost-cell regions

  ! local variables to fill in data
  ! retrieve everything we need from land domain mct_ldom
  ! number of vertices is the size of land domain
  character(CXX) ::  tagname ! hold all fields
  ! TODO: should size it to the number of actual fields we want to exchange
  ! between ghost layers on the component side
  integer, dimension(100) :: tag_indices
  integer :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
  integer, parameter :: max_num_neighbor = 10

  type, public :: grid_cell
     integer           :: num_owned     ! number of owned "active" grid cells
     integer           :: num_ghost     ! number of ghost "active" grid cells
     integer           :: num_ghosted   ! number of owned + ghost "active" grid cells
     integer           :: num_global    ! number of global "active" grid cells     

     integer, pointer  :: natural_id(:) ! [num_ghosted] ID of cell in the full mesh, while accounting for active and inactive cells
     logical, pointer  :: is_owned(:)   ! [num_ghosted] true if the grid cell locally owned
     integer, pointer  :: owner_rank(:) ! [num_ghosted] MPI rank that owns the cell

     integer, pointer  :: num_neighbor (:) ! [num_ghosted] number of neighboring grid cells
     integer, pointer  :: neighbor_id(:,:) ! [num_ghosted, max_num_neighbor] ghosted ID of neighboring grid cells

     integer, pointer  :: num_vertex(:)    ! [num_ghosted] number of vertices
     integer, pointer  :: vertex_id(:,:)   ! [num_ghosted, max_num_neigbor] IDs of vertices

     real(r8), pointer :: lat(:)        ! [num_ghosted] latitude of the cell
     real(r8), pointer :: lon(:)        ! [num_ghosted] longitude of the cell

  end type grid_cell
  
  type, public :: grid_vertex
     integer           :: num_owned     ! number of owned "active" grid vertices
     integer           :: num_ghost     ! number of ghost "active" grid vertices
     integer           :: num_ghosted   ! number of owned + ghost "active" grid vertices
     integer           :: num_global    ! number of global "active" grid vertices

     integer, pointer  :: natural_id(:) ! [num_ghosted] ID of vertex in the full mesh, while accounting for active and inactiver vertices
     integer, pointer  :: is_owned(:)   ! [num_ghosted] true if the grid cell locally owned
     integer, pointer  :: owner_rank(:) ! [num_ghosted] MPI rank that owns the vertex

     real(r8), pointer :: lat(:)        ! [num_ghosted] latitude of the cell
     real(r8), pointer :: lon(:)        ! [num_ghosted] longitude of the cell

  end type grid_vertex

  type, public :: grid_edge
     integer           :: num             ! number of edges
     integer, pointer  :: cell_ids(:,:)   ! ghosted cell IDs left and right of the edge
     real(r8), pointer :: vertex_ids(:,:) ! ghosted ID vertices forming the edge
  end type grid_edge

  ! topological entity list
  integer :: proc_offset  ! current task offset in the global index space

  integer :: topodim      ! topological dimension of the mesh
  integer :: bridgedim    ! bridge dimension to compute ghost layers
  integer :: nghostlayers ! number of ghost layer regions

  ! topological mapping functionality, local 1d gdc arrays
  integer , pointer :: sblock   (:) => null() ! subgrid block parent of element

  !------------------------------------------------------------------------

  public :: mlndghostid
  public :: proc_offset
  public :: sblock

  type(grid_cell)  , public :: moab_gcell
  type(grid_vertex), public :: moab_gvertex
  type(grid_edge)  , public :: moab_edge_internal

  !------------------------------------------------------------------------

contains

  subroutine elm_moab_initialize()

    integer       :: LNDGHOSTID ! id of the ghosted land app
    character*32  :: appname
    integer       :: ierr

    topodim = 2      ! topological dimension = 2: manifold mesh on the sphere
    bridgedim = 1    ! use vertices = 0 as the bridge (other options: edges = 1)
    nghostlayers = 1 ! initialize to zero (default)

    ! next define MOAB app for the ghosted one
    ! We do this so that coupling does not have to deal with halos
    appname="MOAB_ELM_GHOSTED"//C_NULL_CHAR
    LNDGHOSTID = 55120 ! this is arbitrary (but needs to be unique)
    ierr = iMOAB_RegisterApplication(appname, mpicom, LNDGHOSTID, mlndghostid)
    if (ierr > 0 )  &
       call endrun('Error: cannot register ELM-MOAB halo application')
    if(masterproc) then
       write(iulog,*) " "
       write(iulog,*) "register MOAB application:", trim(appname), ", id=", mlndghostid
       write(iulog,*) " "
    endif

  end subroutine elm_moab_initialize

  subroutine elm_moab_load_grid_file(meshfile)

    use elm_varctl  ,  only : iulog  ! for messages and domain file name
    integer   ::  ierr
    character(1024), intent(in) :: meshfile
    integer tagtype, numco !, mbtype, block_ID
    integer :: num_components !, mbtype, block_ID
    character(1024)  :: outfile, wopts
    character(1024) :: tagname ! hold all fields
    integer :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
    integer, dimension(5) :: entity_type
    real(r8), pointer :: data(:)  ! temporary
    integer :: g, n, num_neighbor, num_vertex, num_internal_edges
    integer :: ghosted_id_left, ghosted_id_right
    integer :: nat_id_left, nat_id_right
    integer, pointer :: neighbor_id(:), vertex_id(:)

    if(masterproc) &
        write(iulog,*) "elm_moab_load_grid_file(): reading mesh file: ", trim(meshfile)

    ! load the mesh file
    ierr = iMOAB_LoadMesh( mlndghostid, trim(meshfile)//C_NULL_CHAR, &
                          "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;REPARTITION"//C_NULL_CHAR, &
                          0 )
    if (ierr > 0 )  &
        call endrun('load_grid_file: fail to load the domain file for land model')

    if (masterproc) &
        write(iulog,*) "elm_moab_load_grid_file(): generating ", nghostlayers, " ghost layers"

    ! After the ghost cell exchange, the halo regions are computed and
    ! mesh info will get updated correctly so that we can query the data
    ierr = iMOAB_DetermineGhostEntities(mlndghostid, topodim, & ! topological dimension
                                        nghostlayers, &     ! number of ghost layers
                                        bridgedim )         ! bridge dimension (vertex=0)
    if (ierr > 0)  &
        call endrun('elm_moab_load_grid_file(): failed to generate the ghost layer entities')

#ifdef MOABDEBUG
      ! write out the full repartitioned mesh file to disk, in parallel
      outfile = 'wholeLndGhost.h5m'//C_NULL_CHAR
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(mlndghostid, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('elm_moab_load_grid_file(): failed to write the land mesh file')
#endif

    ! let us get some information about the partitioned mesh and print
    ierr = iMOAB_GetMeshInfo(mlndghostid, nverts, nelem, nblocks, nsbc, ndbc)
    if (ierr > 0 )  &
      call endrun('elm_moab_load_grid_file(): failed to get mesh info ')

    ! set the local size (owned, ghosted) and total entity list (vertices/elements)
    moab_gvertex%num_owned   = nverts(1) ! owned vertices
    moab_gvertex%num_ghost   = nverts(2) ! ghost vertices
    moab_gvertex%num_ghosted = nverts(3) ! owned + ghosted vertices

    moab_gcell%num_owned     = nelem(1)  ! owned elements
    moab_gcell%num_ghost     = nelem(2)  ! ghost elements
    moab_gcell%num_ghosted   = nelem(3)  ! owned + ghost elements

    proc_offset = 0      ! initialize proc_offset

    ! now consolidate/reduce data to root and print information
    ! not really necessary for actual code -- for verbose info only
    ! TODO: combine the Allreduce calls
    call MPI_Allreduce(moab_gvertex%num_owned, moab_gvertex%num_global, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    call MPI_Allreduce(moab_gcell%num_owned  , moab_gcell%num_global  , 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (masterproc) then
      write(iulog, *)  "elm_moab_load_grid_file()(): Number of gridcell vertices: owned=", moab_gvertex%num_owned, &
                        ", ghost=", moab_gvertex%num_ghost, ", global=", moab_gvertex%num_global
      write(iulog, *)  "elm_moab_load_grid_file()(): Number of gridcell elements: owned=", moab_gcell%num_owned, &
                        ", ghost=", moab_gcell%num_ghost, ", global=", moab_gcell%num_global
    endif

    ! Determine the cell id offset on each processor
    !call MPI_Scan(neoproc, proc_offset, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    call MPI_Scan(moab_gcell%num_owned, proc_offset, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (ierr /= 0) then
       write(iulog,*) 'load_grid_file(): MPI_Scan error failed to get proc_offset'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    entity_type(:) = 1 ! default: Element-based tags

    ! add more domain fields that are missing from domain fields: lat, lon, mask, hgt
    !tagname = 'lat:lon:mask:hgt'//C_NULL_CHAR
    tagtype = 0 ! dense, integer
    numco = 1
    tagname='GLOBAL_ID'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage( mlndghostid, tagname, tagtype, numco, tag_indices(1) )
    if (ierr > 0 )  &
      call endrun('Error: fail to retrieve GLOBAL_ID tag ')

    ! use data array as a data holder
    ! Note that loop bounds are typical for locally owned points
    ! allocate(data(neproc*3))
    ! do i = 1, neproc
    !   n = i-1 + bounds%begg
    !   data(i) = ldomain%frac(n)               ! frac = area fractions
    !   data(i+neproc) = ldomain%area(n)/(re*re)   ! area = element area
    !   data(i+2*neproc) = data(i+neproc)             ! aream = model area
    ! enddo

    ! set the values on the internal mesh, halo values aren't set
    ! ierr = iMOAB_SetDoubleTagStorage( mlndghostid, tagname, neproc*3, entity_type, data )
    ! if (ierr > 0 )  &
    !   call endrun('Error: fail to set frac:area:aream tag ')

    ! ierr = iMOAB_UpdateMeshInfo( mlndghostid )
    ! if (ierr > 0 )  &
    !   call endrun('Error: fail to update mesh info ')

    ! synchronize: GLOBAL_ID on vertices in the mesh with ghost layers
    entity_type(1) = 0 ! Vertex tag for GLOBAL_ID
    ierr = iMOAB_SynchronizeTags(mlndghostid, 1, tag_indices, entity_type)
    if (ierr > 0 )  &
      call endrun('Error: fail to synchronize vertex tags for ELM ')

    ! ! synchronize: GLOBAL_ID tag defined on elements in the ghost layers
    entity_type(1) = 1 ! Element tag for GLOBAL_ID
    ierr = iMOAB_SynchronizeTags(mlndghostid, 1, tag_indices, entity_type)
    if (ierr > 0 )  &
      call endrun('Error: fail to synchronize element tags for ELM ')

    ! deallocate(data)

    allocate(moab_gcell%natural_id(moab_gcell%num_ghosted), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for moab_gcell%natural_id'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(moab_gcell%owner_rank(moab_gcell%num_ghosted), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for moab_gcell%owner_rank'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(sblock(moab_gcell%num_ghosted), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for sblock'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    ierr = iMOAB_GetVisibleElementsInfo(mlndghostid, moab_gcell%num_ghosted, &
                    moab_gcell%natural_id, moab_gcell%owner_rank, sblock )
    if (ierr > 0 )  &
    call endrun('Error: fail to query information for visible elements')

    ! allocate memory for storing data about neigbors and vertices for each grid cell
    allocate(moab_gcell%num_neighbor(moab_gcell%num_ghosted))
    allocate(moab_gcell%num_vertex(moab_gcell%num_ghosted))
    allocate(moab_gcell%is_owned(moab_gcell%num_ghosted))
    allocate(moab_gcell%neighbor_id(moab_gcell%num_ghosted, max_num_neighbor))
    allocate(moab_gcell%vertex_id(moab_gcell%num_ghosted, max_num_neighbor))

    ! initialize
    moab_gcell%num_neighbor(:)  = 0
    moab_gcell%num_vertex(:)    = 0
    moab_gcell%is_owned(:)      = .false.
    moab_gcell%neighbor_id(:,:) = -1
    moab_gcell%vertex_id(:,:)   = -1

    ! allocate memory for temporay variables
    allocate(neighbor_id(max_num_neighbor))
    allocate(vertex_id  (max_num_neighbor))

    num_internal_edges = 0

    do g = 1, moab_gcell%num_ghosted

       ! get data about cell neighbors
       num_neighbor = max_num_neighbor
       ierr = iMOAB_GetNeighborElements(mlndghostid, g - 1, num_neighbor, neighbor_id) ! convert g from 1- to 0-based index

       moab_gcell%num_neighbor(g)                = num_neighbor
       moab_gcell%neighbor_id(g, 1:num_neighbor) = neighbor_id(1:num_neighbor) + 1 ! convert to 1-based index

       do n = 1, moab_gcell%num_neighbor(g)
          ghosted_id_left  = g
          ghosted_id_right = moab_gcell%neighbor_id(g, n)

          nat_id_left  = moab_gcell%natural_id(ghosted_id_left)
          nat_id_right = moab_gcell%natural_id(ghosted_id_right)

          if (nat_id_left < nat_id_right) then
             num_internal_edges = num_internal_edges + 1
          endif
       end do

       ! get data about vertices
       num_vertex = max_num_neighbor
       ierr = iMOAB_GetElementConnectivity(mlndghostid, g - 1, num_vertex, vertex_id) ! convert g from 1- to 0-based index

       moab_gcell%num_vertex(g)              = num_vertex
       moab_gcell%vertex_id(g, 1:num_vertex) = vertex_id(1:num_vertex) + 1 ! convert to 1-based index

       if (moab_gcell%owner_rank(g) == iam) moab_gcell%is_owned(g) = .true.

    enddo

    moab_edge_internal%num = num_internal_edges
    allocate(moab_edge_internal%cell_ids(num_internal_edges,2))
    allocate(moab_edge_internal%vertex_ids(num_internal_edges,2))

    num_internal_edges = 0
    do g = 1, moab_gcell%num_ghosted

       do n = 1, moab_gcell%num_neighbor(g)
          ghosted_id_left  = g
          ghosted_id_right = moab_gcell%neighbor_id(g, n)

          nat_id_left  = moab_gcell%natural_id(ghosted_id_left)
          nat_id_right = moab_gcell%natural_id(ghosted_id_right)

          if (nat_id_left < nat_id_right) then
             num_internal_edges = num_internal_edges + 1
             moab_edge_internal%cell_ids(num_internal_edges, 1) = ghosted_id_left
             moab_edge_internal%cell_ids(num_internal_edges, 2) = ghosted_id_right
          endif
       end do
    end do

    ! free up temporary memory
    deallocate(neighbor_id)
    deallocate(vertex_id)

    !  Define and Set Fraction on each mesh
    tagname='frac:area:aream'//C_NULL_CHAR
    tagtype = 1 ! dense, double
    ierr = iMOAB_DefineTagStorage( mlndghostid, tagname, tagtype, numco, tag_indices(2) )
    if (ierr > 0 )  &
      call endrun('Error: fail to create frac:area:aream tags')

    allocate(data(moab_gcell%num_ghosted * 3))
    data(:) = -1
    do g = 1, moab_gcell%num_ghosted
       if (moab_gcell%is_owned(g)) then
          data((g-1)*3 + 1) = iam
          data((g-1)*3 + 2) = iam + 0.1_r8
          data((g-1)*3 + 3) = iam + 0.4_r8
       end if
    end do
    ierr = iMOAB_SetDoubleTagStorage( mlndghostid, tagname, moab_gcell%num_ghosted*3, entity_type(1), data )
    if (ierr > 0) call endrun('Error: setting values failed')

    ierr = iMOAB_SynchronizeTags(mlndghostid, 3, tag_indices(1:3), entity_type(2))
    if (ierr > 0) call endrun('Error: synchronize failed')

    ! reset the data to some arbitrary value
    data(:) = -100.0

    ierr = iMOAB_GetDoubleTagStorage( mlndghostid, tagname, moab_gcell%num_ghosted*3, entity_type(2), data )
    if (ierr > 0) call endrun('Error: getting values failed')

    do g = 1, moab_gcell%num_ghosted
       write(iulog,*)'g: ',g, moab_gcell%is_owned(g), data((g-1)*3+1:(g-1)*3+3)
    end do
    deallocate(data)

#ifdef MOABDEBUG
      ! write out the local mesh file to disk (np tasks produce np files)
      outfile = 'elm_local_mesh'//CHAR(0)
      ierr = iMOAB_WriteLocalMesh(mlndghostid, trim(outfile))
      if (ierr > 0 )  &
        call endrun('Error: fail to write ELM local meshes in h5m format')
#endif

  end subroutine elm_moab_load_grid_file

  !------------------------------------------------------------------------
  subroutine elm_moab_finalize()
  !------------------------------------------------------------------------

    deallocate(moab_gcell%natural_id   )
    deallocate(moab_gcell%owner_rank   )
    deallocate(moab_gcell%num_neighbor )
    deallocate(moab_gcell%neighbor_id  )
    deallocate(moab_gcell%is_owned     )

    deallocate(sblock)

  end subroutine elm_moab_finalize

end module MOABGridType
#endif
