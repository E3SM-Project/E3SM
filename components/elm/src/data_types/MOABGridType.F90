#ifdef HAVE_MOAB
! Only build the module if MOAB is enabled
module MOABGridType

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_kind_mod , only : CXX => SHR_KIND_CXX
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use spmdmod      , only : iam, masterproc, mpicom
  use elm_varctl  ,  only : iulog, fatmlndfrc  ! for messages and domain file name
  use abortutils   , only : endrun
  use ncdio_pio    , only : PIO_subsystem, io_type
  use pio          , only : file_desc_t, var_desc_t, io_desc_t
  use pio          , only : pio_nowrite, pio_noerr
  use pio          , only : pio_openfile, pio_closefile, pio_inq_varid, pio_inq_vartype
  use pio          , only : pio_inq_dimid, pio_inq_dimlen, pio_initdecomp, pio_freedecomp
  use pio          , only : pio_read_darray

  use iMOAB        , only : iMOAB_LoadMesh, iMOAB_WriteMesh, iMOAB_RegisterApplication, &
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
  integer                 :: mlndghostid     ! ID of the MOAB ELM application with ghost-cell regions

  ! local variables to fill in data
  ! retrieve everything we need from land domain mct_ldom
  ! number of vertices is the size of land domain
  character(CXX)          ::  tagname ! hold all fields
  ! TODO: should size it to the number of actual fields we want to exchange
  ! between ghost layers on the component side
  integer, dimension(100) :: tag_indices
  integer                 :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
  integer, parameter      :: max_num_neighbor = 10

  type, public :: grid_cell
     integer           :: num_owned        ! number of owned "active" grid cells
     integer           :: num_ghost        ! number of ghost "active" grid cells
     integer           :: num_ghosted      ! number of owned + ghost "active" grid cells
     integer           :: num_global       ! number of global "active" grid cells     

     integer, pointer  :: natural_id(:)    ! [num_ghosted] ID of cell in the full mesh, while accounting for active and inactive cells
     logical, pointer  :: is_owned(:)      ! [num_ghosted] true if the grid cell locally owned
     integer, pointer  :: owner_rank(:)    ! [num_ghosted] MPI rank that owns the cell

     integer, pointer  :: num_neighbor (:) ! [num_ghosted] number of neighboring grid cells
     integer, pointer  :: neighbor_id(:,:) ! [num_ghosted, max_num_neighbor] ghosted ID of neighboring grid cells

     integer, pointer  :: num_vertex(:)    ! [num_ghosted] number of vertices
     integer, pointer  :: vertex_id(:,:)   ! [num_ghosted, max_num_neigbor] IDs of vertices

     real(r8), pointer :: lat(:)           ! [num_ghosted] latitude of the cell
     real(r8), pointer :: lon(:)           ! [num_ghosted] longitude of the cell

     integer           :: nv               ! length of 'nv' dimension in the mesh
     real(r8), pointer :: latv(:,:)        ! [num_ghosted, nv] latitude of cell vertices
     real(r8), pointer :: lonv(:,:)        ! [num_ghosted, nv] longitude of cell vertices

     integer, pointer  :: elm2moab(:)      ! [num_ghosted] for ELM's grid cell given by g (where begg <= g <= endg),
                                           !               return the correponding MOAB grid cell 'i' (where 1 <= i <= num_ghosted)
     integer, pointer  :: moab2elm(:)      ! [num_ghosted] vice-versa of elm2moab
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

     real(r8), pointer :: lat_vertex(:,:) ! latitude of vertices forming the edge
     real(r8), pointer :: lon_vertex(:,:) ! longitude of vertices forming the edge

     real(r8), pointer :: dc(:)           ! distance between the cell centers left and right of the edge
     real(r8), pointer :: lv(:)           ! length of the edge shared by the cells
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
    bridgedim = 0    ! use vertices = 0 as the bridge (other options: edges = 1)
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
    !
    use elm_varctl  ,  only : iulog  ! for messages and domain file name
    !
    character(len=*), intent(in)  :: meshfile
    !
    integer                 :: tagtype, numco
    integer                 :: num_components
    character(1024)         :: outfile, wopts
    character(1024)         :: tagname                                            ! hold all fields
    integer                 :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
    integer  , dimension(2) :: entity_type
    real(r8) , allocatable  :: data(:,:)                                          ! for tags to be set in MOAB
    integer                 :: g, n, num_neighbor, num_vertex, num_internal_edges
    integer                 :: ghosted_id_left, ghosted_id_right
    integer                 :: nat_id_left, nat_id_right
    integer  , pointer      :: neighbor_id(:), vertex_id(:)
    integer                 :: begg, endg, count
    integer                 :: ierr

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
    call MPI_Allreduce(moab_gvertex%num_owned, moab_gvertex%num_global, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    call MPI_Allreduce(moab_gcell%num_owned  , moab_gcell%num_global  , 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (masterproc) then
       write(iulog, *)  "elm_moab_load_grid_file()(): Number of gridcell vertices: owned=", moab_gvertex%num_owned, &
            ", ghost=", moab_gvertex%num_ghost, ", global=", moab_gvertex%num_global
       write(iulog, *)  "elm_moab_load_grid_file()(): Number of gridcell elements: owned=", moab_gcell%num_owned, &
            ", ghost=", moab_gcell%num_ghost, ", global=", moab_gcell%num_global
    endif

    ! Determine the cell id offset on each processor
    call MPI_Scan(moab_gcell%num_owned, proc_offset, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (ierr /= 0) then
       write(iulog,*) 'load_grid_file(): MPI_Scan error failed to get proc_offset'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    entity_type(:) = 1 ! default: Element-based tags

    ! add more domain fields that are missing from domain fields: lat, lon, mask, hgt
    tagtype = 0 ! dense, integer
    numco = 1
    tagname='GLOBAL_ID'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage( mlndghostid, tagname, tagtype, numco, tag_indices(1) )
    if (ierr > 0 )  &
         call endrun('Error: fail to retrieve GLOBAL_ID tag ')

    ! synchronize: GLOBAL_ID on vertices in the mesh with ghost layers
    entity_type(1) = 0 ! Vertex tag for GLOBAL_ID
    ierr = iMOAB_SynchronizeTags(mlndghostid, 1, tag_indices(1), entity_type(1))
    if (ierr > 0 )  &
         call endrun('Error: fail to synchronize vertex tags for ELM ')

    ! ! synchronize: GLOBAL_ID tag defined on elements in the ghost layers
    entity_type(2) = 1 ! Element tag for GLOBAL_ID
    ierr = iMOAB_SynchronizeTags(mlndghostid, 1, tag_indices(1), entity_type(2))
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

    ! poppulate the moab_gcell data structure and determine number of unique
    ! internal edges
    do g = 1, moab_gcell%num_ghosted

       ! get data about cell neighbors
       num_neighbor = max_num_neighbor
       ierr = iMOAB_GetNeighborElements(mlndghostid, g - 1, num_neighbor, neighbor_id) ! convert g from 1- to 0-based index

       moab_gcell%num_neighbor(g)                = num_neighbor
       moab_gcell%neighbor_id(g, 1:num_neighbor) = neighbor_id(1:num_neighbor) + 1 ! convert to 1-based index

       ! determine unqiue edges between cells
       do n = 1, moab_gcell%num_neighbor(g)
          ghosted_id_left  = g
          ghosted_id_right = moab_gcell%neighbor_id(g, n)

          nat_id_left  = moab_gcell%natural_id(ghosted_id_left)
          nat_id_right = moab_gcell%natural_id(ghosted_id_right)

          ! check if the edge hasn't been already included
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

    ! Set up ELM-to-MOAB and MOAB-to-ELM mapping
    begg = proc_offset - moab_gcell%num_owned + 1
    endg = proc_offset + moab_gcell%num_ghost

    allocate(moab_gcell%elm2moab(begg:endg))
    allocate(moab_gcell%moab2elm(1:moab_gcell%num_ghosted))

    count = 0
    begg = begg - 1

    ! first, put the owned cells
    do g = 1, moab_gcell%num_ghosted
       if (moab_gcell%is_owned(g)) then

          count = count + 1

          moab_gcell%elm2moab(begg + count) = g
          moab_gcell%moab2elm(g)            = begg + count
       end if
    end do

    ! next, put the ghost cells
    do g = 1, moab_gcell%num_ghosted
       if (.not.moab_gcell%is_owned(g)) then
          count = count + 1

          moab_gcell%elm2moab(begg + count) = g
          moab_gcell%moab2elm(g)            = begg + count
       end if
    end do

    ! populate data structure for internal edges
    moab_edge_internal%num = num_internal_edges
    allocate(moab_edge_internal%cell_ids(num_internal_edges,2))

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

#ifdef MOABDEBUG
    ! write out the local mesh file to disk (np tasks produce np files)
    outfile = 'elm_local_mesh'//CHAR(0)
    ierr = iMOAB_WriteLocalMesh(mlndghostid, trim(outfile))
    if (ierr > 0 )  &
         call endrun('Error: fail to write ELM local meshes in h5m format')
#endif

    ! Read lat/lon for cell centers and cell vertices
    call read_grid_cell_lat_lon()

    ! Now determine the lat/lon of vertices of internal edges
    call set_internal_edge_lat_lon()

    call compute_length_for_internal_edge()

  end subroutine elm_moab_load_grid_file

  !------------------------------------------------------------------------------
  subroutine read_grid_cell_lat_lon()
    !
    ! !DESCRIPTION:
    ! Reads latitude/longitude for:
    ! - Cell centers, and
    ! - Cell vertices
    !
    implicit none
    !
    type(file_desc_t) :: ncid
    type(var_desc_t)  :: varid
    type(io_desc_t)   :: iodescNCells
    character (1024)  :: varname
    integer           :: vartype
    integer           :: dim2d(2), dim3d(3)
    integer, pointer  :: compdof(:)
    integer           :: g, v
    integer           :: ierr
    integer           :: ni, nj, dimid, count
    real(r8), pointer :: data2d(:,:)

    ierr = pio_openfile(pio_subsystem, ncid, io_type, trim(fatmlndfrc), PIO_NOWRITE)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: Unable to open file : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimid(ncid, 'ni', dimid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: ni dimension not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimlen(ncid, dimid, ni)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: Unable to determine length of ni dimension : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimid(ncid, 'nj', dimid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: nj dimension not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimlen(ncid, dimid, nj)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: Unable to determine length of nj dimension : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimid(ncid, 'nv', dimid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: nv dimension not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_dimlen(ncid, dimid, moab_gcell%nv)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: Unable to determine length of nv dimension : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! set dims for reading lat/lon for cell centers
    dim2d(1) = ni
    dim2d(2) = nj

    ! set dims for reading lat/lon for cell vertices
    dim3d(1) = moab_gcell%nv
    dim3d(2) = ni
    dim3d(3) = nj

    !
    ! Read lat/lon for cell centers
    !

    allocate(moab_gcell%lat(moab_gcell%num_ghosted))
    allocate(moab_gcell%lon(moab_gcell%num_ghosted))

    varname = 'xc'
    ierr = pio_inq_varid(ncid, trim(varname), varid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: xc variable not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_vartype(ncid, varid, vartype)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: could not determine the type of xc variable in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Create the I/O decomposition for owned+ghost cells
    call pio_initdecomp(pio_subsystem, vartype, dim2d, moab_gcell%natural_id, iodescNCells)

    ! Read the data
    call pio_read_darray(ncid, varid, iodescNCells, moab_gcell%lon, ierr)

    varname = 'yc'
    ierr = pio_inq_varid(ncid, trim(varname), varid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: yc variable not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Read the data using the previously created I/O decomposition
    call pio_read_darray(ncid, varid, iodescNCells, moab_gcell%lat, ierr)

    ! Free up memory
    call pio_freedecomp(pio_subsystem, iodescNCells)

    !
    ! Read lat/lon for cell vertices
    !

    allocate(moab_gcell%latv(moab_gcell%num_ghosted, moab_gcell%nv))
    allocate(moab_gcell%lonv(moab_gcell%num_ghosted, moab_gcell%nv))
    allocate(data2d(moab_gcell%nv, moab_gcell%num_ghosted))

    varname = 'xv'
    ierr = pio_inq_varid(ncid, trim(varname), varid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: xv variable not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_vartype(ncid, varid, vartype)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: could not determine the type of xv variable in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Create a temporary variable that has the indices for reading the data
    ! for owned+ghost cells
    allocate(compdof(moab_gcell%num_ghosted * moab_gcell%nv))
    count = 0
    do g = 1, moab_gcell%num_ghosted
       do v = 1, moab_gcell%nv
          count = count + 1;
          compdof(count) = (moab_gcell%natural_id(g) - 1 ) * moab_gcell%nv + v
       enddo
    enddo

    ! Create the I/O decomposition
    call pio_initdecomp(pio_subsystem, vartype, dim3d, compdof, iodescNCells)

    ! Read the data
    call pio_read_darray(ncid, varid, iodescNCells, data2d, ierr)

    ! Copy the data
    do g = 1, moab_gcell%num_ghosted
       do v = 1, moab_gcell%nv
          moab_gcell%lonv(g, v) = data2d(v, g)
       end do
    end do

    varname = 'yv'
    ierr = pio_inq_varid(ncid, trim(varname), varid)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: yv variable not present in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ierr = pio_inq_vartype(ncid, varid, vartype)
    if (ierr /= pio_noerr) then
       write(iulog,*) 'Error: could not determine the type of xc variable in : ' // trim(fatmlndfrc)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Read the data
    call pio_read_darray(ncid, varid, iodescNCells, data2d, ierr)

    ! Copy the data
    do g = 1, moab_gcell%num_ghosted
       do v = 1, moab_gcell%nv
          moab_gcell%latv(g, v) = data2d(v, g)
       end do
    end do

    ! Clean up
    call pio_freedecomp(pio_subsystem, iodescNCells)
    call pio_closefile(ncid)
    deallocate(compdof)

  end subroutine read_grid_cell_lat_lon

  !------------------------------------------------------------------------------
  subroutine set_internal_edge_lat_lon()
    !
    ! !DESCRIPTION:
    ! For each internal edge, determine the lat/lon of vertex forming an edge
    !
    implicit none
    !
    integer             :: e, gup, gdn, vertex_count
    integer             :: v1, v2, min_v2
    real(r8)            :: latv_up, lonv_up, latv_dn, lonv_dn, dist, min_dist
    real(r8), parameter :: dist_threshold = 1.e-10

    allocate(moab_edge_internal%lat_vertex(moab_edge_internal%num, 2))
    allocate(moab_edge_internal%lon_vertex(moab_edge_internal%num, 2))

    do e = 1, moab_edge_internal%num
       gup = moab_edge_internal%cell_ids(e, 1)
       gdn = moab_edge_internal%cell_ids(e, 2)

       vertex_count = 0

       do v1 = 1, moab_gcell%nv
          latv_up = moab_gcell%latv(gup, v1)
          lonv_up = moab_gcell%lonv(gup, v1)

          do v2 = 1, moab_gcell%nv
             latv_dn = moab_gcell%latv(gdn, v2)
             lonv_dn = moab_gcell%lonv(gdn, v2)

             dist = (latv_dn - latv_up)**2._r8 + (lonv_dn - lonv_up)**2._r8
             if (v2 == 1) then
                min_dist = dist;
                min_v2   = v2;
             else if (dist < min_dist) then
                min_dist = dist;
                min_v2   = v2;
             end if
          end do ! v2

          if (min_dist < dist_threshold) then
             vertex_count = vertex_count + 1
             if (vertex_count > 2) then
                call endrun('Found more than 2 vertices that are shared between cells')
             end if
             moab_edge_internal%lat_vertex(e, vertex_count) = latv_up
             moab_edge_internal%lon_vertex(e, vertex_count) = lonv_up
          end if
       end do ! v1
    end do ! e

  end subroutine set_internal_edge_lat_lon

  !------------------------------------------------------------------------------
  subroutine haversine_dist(lat1, lon1, lat2, lon2, dist)
    !
    use shr_const_mod   , only : SHR_CONST_PI, SHR_CONST_REARTH
    !
    implicit none
    !
    real(r8) :: lat1, lon1, lat2, lon2, dist
    !
    real(r8) :: a, c, dlat, dlon

    ! Convert degrees to radians
    dlat = (lat2 - lat1) * SHR_CONST_PI / 180._r8
    dlon = (lon2 - lon1) * SHR_CONST_PI / 180._r8

    dist = 0.0_r8
    a = sin(dlat / 2.0_r8)**2._r8 + cos(lat1 * SHR_CONST_PI / 180.0_r8) * cos(lat2 * SHR_CONST_PI / 180._r8) * &
         sin(dlon / 2.0_r8)**2._r8
    c = 2._r8 * atan2(sqrt(a), sqrt(1._r8 - a))

    dist = SHR_CONST_REARTH * c

  end subroutine haversine_dist

  !------------------------------------------------------------------------------
  subroutine compute_length_for_internal_edge()
    !
    ! !DESCRIPTION:
    ! For each internal edge, compute:
    ! - Distance between cells, and
    ! - Length of the edge shared by the cells
    !
    implicit none
    !
    integer  :: e, gup, gdn
    real(r8) :: lat1, lon1, lat2, lon2

    allocate(moab_edge_internal%dc(moab_edge_internal%num))
    allocate(moab_edge_internal%lv(moab_edge_internal%num))

    do e = 1, moab_edge_internal%num
       gup = moab_edge_internal%cell_ids(e, 1)
       gdn = moab_edge_internal%cell_ids(e, 2)

       lat1 = moab_gcell%lat(gup)
       lon1 = moab_gcell%lon(gup)

       lat2 = moab_gcell%lat(gdn)
       lon2 = moab_gcell%lon(gdn)

       call haversine_dist(lat1, lon1, lat2, lon2, moab_edge_internal%dc(e))

       lat1 = moab_edge_internal%lat_vertex(e, 1)
       lon1 = moab_edge_internal%lon_vertex(e, 1)

       lat2 = moab_edge_internal%lat_vertex(e, 2)
       lon2 = moab_edge_internal%lon_vertex(e, 2)

       call haversine_dist(lat1, lon1, lat2, lon2, moab_edge_internal%lv(e))

    end do

  end subroutine compute_length_for_internal_edge

  !------------------------------------------------------------------------
  subroutine elm_moab_finalize()
    !
    ! !DESCRIPTION:
    ! Frees up memory

    deallocate(moab_gcell%natural_id   )
    deallocate(moab_gcell%owner_rank   )
    deallocate(moab_gcell%num_neighbor )
    deallocate(moab_gcell%neighbor_id  )
    deallocate(moab_gcell%is_owned     )

    deallocate(sblock)

  end subroutine elm_moab_finalize

end module MOABGridType
#endif
