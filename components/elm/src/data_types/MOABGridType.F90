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
  use spmdmod  , only: masterproc, mpicom
  use elm_varctl  ,  only : iulog, fatmlndfrc  ! for messages and domain file name
  use abortutils , only : endrun

  use iMOAB  , only: iMOAB_LoadMesh, iMOAB_WriteMesh, iMOAB_RegisterApplication, &
  iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage, iMOAB_SynchronizeTags, &
  iMOAB_UpdateMeshInfo, iMOAB_GetMeshInfo, &
  iMOAB_DetermineGhostEntities, iMOAB_WriteLocalMesh, iMOAB_GetVisibleElementsInfo

  use mpi
  use iso_c_binding

  !
  ! !PUBLIC TYPES:
  implicit none
  ! save
  private
  !
  public :: initialize
  public :: load_grid_file
  public :: finalize
  !
  integer :: mlndghostid     ! ID of the MOAB ELM application with ghost-cell regions

  !integer n
  ! local variables to fill in data
  ! retrieve everything we need from land domain mct_ldom
  ! number of vertices is the size of land domain
  character(CXX) ::  tagname ! hold all fields
  ! TODO: should size it to the number of actual fields we want to exchange
  ! between ghost layers on the component side
  integer, dimension(100) :: tag_indices
  integer :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)

  ! topological entity list
  integer :: neg, neoproc, negproc, neproc   ! number of elements: global, owned, ghosted, total
  integer :: nvg, nvoproc, nvgproc, nvproc   ! number of vertices: global, owned, ghosted, total
  integer :: proc_offset  ! current task offset in the global index space

  integer :: topodim      ! topological dimension of the mesh
  integer :: bridgedim    ! bridge dimension to compute ghost layers
  integer :: nghostlayers ! number of ghost layer regions

  ! topological mapping functionality, local 1d gdc arrays
  integer , pointer :: globid   (:) => null() ! global index of element
  integer , pointer :: eowner   (:) => null() ! task owner of element
  integer , pointer :: sblock   (:) => null() ! subgrid block parent of element

  !------------------------------------------------------------------------

  public :: mlndghostid
  public :: neg, neoproc, negproc, neproc
  public :: nvg, nvoproc, nvgproc, nvproc
  public :: proc_offset
  public :: globid, eowner, sblock

  !------------------------------------------------------------------------

contains

  subroutine initialize()

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

  end subroutine initialize

  subroutine load_grid_file(meshfile)
    use elm_varctl  ,  only : iulog, fatmlndfrc  ! for messages and domain file name

    integer   ::  ierr
    character*100, intent(in) :: meshfile
    integer tagtype, numco !, mbtype, block_ID
    integer :: num_components !, mbtype, block_ID
    character*100  :: outfile, wopts
    character(CXX) :: tagname ! hold all fields
    integer :: nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
    integer, dimension(5) :: entity_type
    real(r8), pointer :: data(:)  ! temporary

    if(masterproc) &
        write(iulog,*) "MOABGridType%load_grid_file: reading mesh file: ", trim(meshfile)

    ! load the mesh file
    ierr = iMOAB_LoadMesh( mlndghostid, trim(meshfile)//C_NULL_CHAR, &
                          "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;REPARTITION"//C_NULL_CHAR, &
                          0 )
    if (ierr > 0 )  &
        call endrun('load_grid_file: fail to load the domain file for land model')

    if (masterproc) &
        write(iulog,*) "MOABGridType%load_grid_file: generating ", nghostlayers, " ghost layers"

    ! After the ghost cell exchange, the halo regions are computed and
    ! mesh info will get updated correctly so that we can query the data
    ierr = iMOAB_DetermineGhostEntities(mlndghostid, topodim, & ! topological dimension
                                        nghostlayers, &     ! number of ghost layers
                                        bridgedim )         ! bridge dimension (vertex=0)
    if (ierr > 0)  &
        call endrun('MOABGridType%load_grid_file: failed to generate the ghost layer entities')

#ifdef MOABDEBUG
      ! write out the full repartitioned mesh file to disk, in parallel
      outfile = 'wholeLndGhost.h5m'//C_NULL_CHAR
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(mlndghostid, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('MOABGridType%load_grid_file: failed to write the land mesh file')
#endif

    ! let us get some information about the partitioned mesh and print
    ierr = iMOAB_GetMeshInfo(mlndghostid, nverts, nelem, nblocks, nsbc, ndbc)
    if (ierr > 0 )  &
      call endrun('MOABGridType%load_grid_file: failed to get mesh info ')

    ! set the local size (owned, ghosted) and total entity list (vertices/elements)
    nvoproc = nverts(1)  ! owned vertices
    nvgproc = nverts(2)  ! ghosted vertices
    nvproc = nverts(3)   ! owned + ghosted vertices
    neoproc = nelem(1)   ! owned elements
    negproc = nelem(2)   ! ghosted elements
    neproc = nelem(3)    ! owned + ghosted elements
    proc_offset = 0      ! initialize proc_offset

    ! now consolidate/reduce data to root and print information
    ! not really necessary for actual code -- for verbose info only
    ! TODO: combine the calls
    call MPI_Allreduce(nvoproc, nvg, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    call MPI_Allreduce(neoproc, neg, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (masterproc) then
      write(iulog, *)  "MOABGridType%load_grid_file(): Number of gridcell vertices: owned=", nvoproc, &
                        ", ghosted=", nvgproc, ", global=", nvg
      write(iulog, *)  "MOABGridType%load_grid_file(): Number of gridcell elements: owned=", neoproc, &
                        ", ghosted=", negproc, ", global=", neg
    endif

    ! Determine the cell id offset on each processor
    call MPI_Scan(neoproc, proc_offset, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
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

    !  Define and Set Fraction on each mesh
    tagname='frac:area:aream'//C_NULL_CHAR
    tagtype = 1 ! dense, double
    ierr = iMOAB_DefineTagStorage( mlndghostid, tagname, tagtype, numco, tag_indices(2) )
    if (ierr > 0 )  &
      call endrun('Error: fail to create frac:area:aream tags')

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

    allocate(globid(neproc), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for globid'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(eowner(neproc), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for eowner'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(sblock(neproc), stat=ierr)
    if (ierr /= 0) then
        write(iulog,*) 'load_grid_file(): allocation error for sblock'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    ierr = iMOAB_GetVisibleElementsInfo(mlndghostid, neproc, &
                    globid, eowner, sblock )
    if (ierr > 0 )  &
    call endrun('Error: fail to query information for visible elements')

#ifdef MOABDEBUG
      ! write out the local mesh file to disk (np tasks produce np files)
      outfile = 'elm_local_mesh'//CHAR(0)
      ierr = iMOAB_WriteLocalMesh(mlndghostid, trim(outfile))
      if (ierr > 0 )  &
        call endrun('Error: fail to write ELM local meshes in h5m format')
#endif

  end subroutine load_grid_file


  !------------------------------------------------------------------------
  subroutine finalize()
    !------------------------------------------------------------------------

    deallocate(globid)
    deallocate(eowner)
    deallocate(sblock)

  end subroutine finalize

end module MOABGridType
#endif