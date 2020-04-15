module dshr_strdata_mod

  ! holds data and methods to advance data models

  use ESMF
  use shr_const_mod , only : SHR_CONST_PI
  use shr_kind_mod  , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod  , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_kind_mod  , only : CXX=>SHR_KIND_CXX
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod   , only : shr_mpi_bcast
  use shr_file_mod  , only : shr_file_getunit, shr_file_freeunit
  use shr_log_mod   , only : loglev  => shr_log_Level
  use shr_log_mod   , only : logunit => shr_log_Unit
  use shr_cal_mod   , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod   , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod   , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod   , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_nl_mod    , only : shr_nl_find_group_name
  use shr_pio_mod   , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
  use shr_mpi_mod   , only : shr_mpi_bcast
  use shr_ncread_mod, only : shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes, shr_ncread_tField
  use shr_ncread_mod, only : shr_ncread_domain, shr_ncread_vardimsizes
  use shr_ncread_mod, only : shr_ncread_varexists, shr_ncread_vardimnum,  shr_ncread_field4dG
  use shr_ncread_mod, only : shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes
  use shr_ncread_mod, only : shr_ncread_tField
  use shr_string_mod
  use dshr_methods_mod , only : chkerr
  use dshr_stream_mod ! stream data type and methods
  use dshr_dmodel_mod ! shr data model stuff
  use dshr_tinterp_mod
  use shr_mct_mod
  use mct_mod        ! mct
  use perf_mod       ! timing
  use pio            ! pio

  implicit none
  private

  ! !PUBLIC TYPES:
  public ::  shr_strdata_type

  public ::  shr_strdata_readnml
  public ::  shr_strdata_restRead
  public ::  shr_strdata_restWrite
  public ::  shr_strdata_setOrbs
  public ::  shr_strdata_print
  public ::  shr_strdata_init_model_domain
  public ::  shr_strdata_init_streams
  public ::  shr_strdata_init_mapping
  public ::  shr_strdata_getfrac_from_stream
  public ::  shr_strdata_advance
  public ::  shr_strdata_clean
  public ::  shr_strdata_pioinit

  private :: shr_strdata_init_stream_domain

  ! !PUBLIC DATA MEMBERS:

  integer                              :: debug    = 0  ! local debug flag
  integer          ,parameter          :: nStrMax = 30
  integer          ,parameter          :: nVecMax = 30
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  real(R8)         ,parameter, private :: dtlimit_default = 1.5_R8

  type shr_strdata_type
     ! --- set by input ---
     character(CL)                  :: dataMode          ! flags physics options wrt input data
     character(CL)                  :: streams (nStrMax) ! stream description file names
     character(CL)                  :: taxMode (nStrMax) ! time axis cycling mode
     real(R8)                       :: dtlimit (nStrMax) ! dt max/min limit
     character(CL)                  :: vectors (nVecMax) ! define vectors to vector map
     character(CL)                  :: fillalgo(nStrMax) ! fill algorithm
     character(CL)                  :: fillmask(nStrMax) ! fill mask
     character(CL)                  :: fillread(nStrMax) ! fill mapping file to read
     character(CL)                  :: fillwrit(nStrMax) ! fill mapping file to write
     character(CL)                  :: mapalgo (nStrMax) ! scalar map algorithm
     character(CL)                  :: mapmask (nStrMax) ! scalar map mask
     character(CL)                  :: mapread (nStrMax) ! regrid mapping file to read
     character(CL)                  :: mapwrit (nStrMax) ! regrid mapping file to write
     character(CL)                  :: tintalgo(nStrMax) ! time interpolation algorithm
     character(CL)                  :: readmode(nStrMax) ! file read mode
     integer                        :: io_type
     integer                        :: io_format

     ! --- mpi info
     integer                        :: mpicom
     integer                        :: ntasks 
     integer                        :: my_task
     integer                        :: master_task

     ! --- standard output info
     integer                        :: logunit

     !--- data required by stream  cosz t-interp method, set by user ---
     real(R8)                       :: eccen
     real(R8)                       :: mvelpp
     real(R8)                       :: lambm0
     real(R8)                       :: obliqr
     integer                        :: modeldt           ! model dt in seconds

     ! --- model domain info, internal, public ---
     character(CL)                  :: domainFile        ! file containing domain info (set my input)
     integer                        :: nxg               ! model grid lon size
     integer                        :: nyg               ! model grid lat size
     integer                        :: nzg               ! model grid vertical size
     integer                        :: lsize             ! model grid local size
     type(mct_gsmap)                :: gsmap             ! model grid global seg map
     type(mct_ggrid)                :: grid              ! model grid ggrid
     type(mct_avect)                :: avs(nStrMax)      ! model grid stream attribute vectors
                                                         ! stream attribute vectors that are time and spatially
                                                         ! interpolated to model grid
     ! --- stream info, internal ---
     type(shr_stream_streamType)    :: stream(nStrMax)
     type(iosystem_desc_t), pointer :: pio_subsystem => null()
     type(io_desc_t)                :: pio_iodesc(nStrMax)
     integer                        :: nstreams          ! number of streams
     integer                        :: strnxg(nStrMax)
     integer                        :: strnyg(nStrMax)
     integer                        :: strnzg(nStrMax)
     logical                        :: dofill(nStrMax)
     logical                        :: domaps(nStrMax)
     integer                        :: lsizeR(nStrMax)
     type(mct_gsmap)                :: gsmapR(nStrMax)
     type(mct_rearr)                :: rearrR(nStrMax)
     type(mct_ggrid)                :: gridR(nStrMax)
     type(mct_avect)                :: avRFile(nStrMax)  ! Read attrvect for multiple time slices
     type(mct_avect)                :: avRLB(nStrMax)    ! Read attrvect
     type(mct_avect)                :: avRUB(nStrMax)    ! Read attrvect
     type(mct_avect)                :: avFUB(nStrMax)    ! Final attrvect
     type(mct_avect)                :: avFLB(nStrMax)    ! Final attrvect
     type(mct_avect)                :: avCoszen(nStrMax) ! data assocaited with coszen time interp
     type(mct_sMatP)                :: sMatPf(nStrMax)
     type(mct_sMatP)                :: sMatPs(nStrMax)
     integer                        :: ymdLB(nStrMax),todLB(nStrMax)
     integer                        :: ymdUB(nStrMax),todUB(nStrMax)
     real(R8)                       :: dtmin(nStrMax)
     real(R8)                       :: dtmax(nStrMax)
     integer                        :: ymd  ,tod
     character(CL)                  :: calendar          ! model calendar for ymd,tod
     integer                        :: nvectors          ! number of vectors
     integer                        :: ustrm (nVecMax)
     integer                        :: vstrm (nVecMax)
     character(CL)                  :: allocstring
  end type shr_strdata_type

  interface shr_strdata_init_model_domain
     module procedure :: shr_strdata_init_model_domain_mesh
     module procedure :: shr_strdata_init_model_domain_scol
  end interface shr_strdata_init_model_domain

  real(R8),parameter,private :: deg2rad = SHR_CONST_PI/180.0_R8
  character(len=*),parameter :: allocstring_value = 'strdata_allocated'

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_strdata_init_model_domain_mesh( mesh, mpicom, compid, sdat, rc)

    ! ---------------------------------------------------------------------
    ! Initialize sdat%lsize, sdat%gsmap and sdat%grid
    ! sdat%nxg, sdat%nyg and sdat%nzg are initialized in shr_strdata_readnl
    ! sdat%avs(:) is initialized in shr_strdata_init_mapping
    ! ---------------------------------------------------------------------

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: compid
    type(shr_strdata_type) , intent(inout) :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    integer              :: n,k          ! generic counters
    integer              :: lsize        ! local size
    integer              :: gsize
    type(ESMF_DistGrid)  :: distGrid
    integer              :: dimCount
    integer              :: tileCount
    integer              :: deCount
    integer, allocatable :: elementCountPTile(:)
    integer, allocatable :: indexCountPDE(:,:)
    integer              :: spatialDim         ! number of dimension in mesh
    integer              :: numOwnedElements   ! size of mesh
    real(r8), pointer    :: ownedElemCoords(:) ! mesh lat and lons
    real(r8), pointer    :: lat(:), lon(:)     ! mesh lats and lons
    integer              :: klon
    integer              :: klat
    integer              :: khgt
    integer              :: kmask
    integer              :: my_task
    integer              :: ierr
    integer              :: nxg,nyg,nzg  ! size of input fields
    integer, pointer     :: idata(:)     ! temporary
    type(ESMF_Array)     :: elemMaskArray
    integer, pointer     :: elemMask(:)
    integer, allocatable, target :: gindex(:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize sdat%lsize and sdat%gsmap (the data model gsmap)
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    sdat%lsize = lsize
    allocate(gindex(sdat%lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, dimCount=dimCount, deCount=deCount, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    gsize = gsize
    deallocate(elementCountPTile)
    call mct_gsMap_init(sdat%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    ! initialize sdat%ggrid
    ! to obtain the mask, create an esmf array from a distgrid and a pointer
    ! the following call will automatically fill in the pointer value for elemMask
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elemMask(lsize)) 
    elemMaskArray = ESMF_ArrayCreate(distGrid, elemMask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create the model ggrid currently
    call mpi_comm_rank(mpicom, my_task, ierr)
    call mct_gGrid_init(GGrid=sdat%grid, CoordChars='lat:lon:hgt', OtherChars='mask', lsize=lsize )
    call mct_gsMap_orderedPoints(sdat%gsMap, my_task, idata)
    call mct_gGrid_importIAttr(sdat%Grid,'GlobGridNum', idata, lsize)
    deallocate(idata)
    klon  = mct_aVect_indexRA(sdat%Grid%data,'lon')
    klat  = mct_aVect_indexRA(sdat%Grid%data,'lat')
    kmask = mct_aVect_indexRA(sdat%Grid%data,'mask')
    khgt  = mct_aVect_indexRA(sdat%Grid%data,'hgt')
    do n = 1, lsize
       sdat%Grid%data%rAttr(klon ,n) = ownedElemCoords(2*n-1)
       sdat%Grid%data%rAttr(klat ,n) = ownedElemCoords(2*n)
       sdat%Grid%data%rAttr(kmask,n) = real(elemMask(n),r8)
       !sdat%Grid%data%rAttr(khgt ,n) = -1 ! currently hard-wired
    end do
    deallocate(elemMask)

  end subroutine shr_strdata_init_model_domain_mesh

  !===============================================================================
  subroutine shr_strdata_init_model_domain_scol(scmlon, scmlat, mpicom, compid, sdat) 

    !----------------------------------------------------------------------------
    ! Create mct ggrid for model grid and set model gsmap if not input
    ! - assumes a very specific netCDF domain file format wrt var names, etc.
    !----------------------------------------------------------------------------

    ! input/output variables
    real(R8)               , intent(in)    :: scmlon     ! single column lon
    real(R8)               , intent(in)    :: scmlat     ! single column lat
    integer                , intent(in)    :: mpicom     ! mpi communicator
    integer                , intent(in)    :: compid
    type(shr_strdata_type) , intent(inout) :: sdat

    !----- local -----
    integer              :: n,k,j,i    ! indices
    integer              :: my_task
    integer              :: ierr       ! error code
    logical              :: fileexists ! true if input domain file exists
    character(CS)        :: lonname    ! name of  lon variable
    character(CS)        :: latname    ! name of  lat variable
    integer              :: nlon
    integer              :: nlat
    integer              :: nmask
    integer              :: nxg, nyg 
    real(R8)             :: dist,mind  ! scmmode point search
    integer              :: ni,nj      ! scmmode point search
    real(R8)             :: lscmlon    ! local copy of scmlon
    real(R8),allocatable :: lon(:,:)   ! temp array for domain lon  info
    real(R8),allocatable :: lat(:,:)   ! temp array for domain lat  info
    integer, pointer     :: idata(:)   ! temporary
    integer              :: gsize 
    integer, allocatable, target :: gindex(:)
    character(*), parameter :: subname = '(shr_strdata_init_model_domain_scol) '
    character(*), parameter :: F00   = "('(shr_strdata_init_model_domain_scol) ',8a)"
    !-------------------------------------------------------------------------------

    ! error check
    call mpi_comm_rank(mpicom,my_task,ierr)
    if (my_task > 0) then
       write(logunit,*) subname,' ERROR: scmmode must be run on one pe'
       call shr_sys_abort(subname//' ERROR: scmmode2 tasks')
    endif

    ! reset sdat%nxg, sdat%nyg and sdat%nzg
    sdat%nxg  =  1
    sdat%nyg  =  1
    sdat%nzg  = -1

    ! initialize sdat%lsize and sdat%gsmap
    sdat%lsize = 1
    gsize = 1
    allocate(gindex(1)); gindex(1) = 1
    call mct_gsMap_init(sdat%gsmap, gindex, mpicom, compid, sdat%lsize, gsize)

    ! sdat%grid is for a single point - but the input model domain from the namelist
    ! is for the total grid-  since will be finding the nearest neighbor
    ! Read full model domain as specified in namelist
    ! Need to allocate the total grid size of the domain file - ane
    ! will find the nearest neighbor fot the single column

    inquire(file=trim(sdat%domainfile), exist=fileExists)
    if (.not. fileExists) then
       write(logunit,F00) "ERROR: file does not exist: ", trim(sdat%domainfile)
       call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(sdat%domainfile))
    end if
    lonname  = "xc" 
    latname  = "yc"
    call shr_ncread_varDimSizes(sdat%domainfile, lonName, n1=nxg)
    call shr_ncread_varDimSizes(sdat%domainfile, latName, n1=nyg)
    allocate(lon(nxg,nyg))
    allocate(lat(nxg,nyg))
    call shr_ncread_domain(sdat%domainfile, lonName, lon, latName, lat)

    ! Find nearest neigbor of input domain for scmlon and scmlat
    ! determine whether dealing with 2D input files (typical of Eulerian
    ! want lon values between 0 and 360, assume 1440 is enough (start with wraparound)

    lscmlon = mod(scmlon+1440.0_r8,360.0_r8)
    lon     = mod(lon   +1440.0_r8,360.0_r8)
    if (nyg .ne. 1) then
       ! lat and lon are on 2D logically rectanular arrays
       ! assumes regular 2d grid for compatability with shr_scam_getCloseLatLon ---
       ni = 1
       mind = abs(lscmlon - (lon(1,1)+360.0_r8))
       do i=1,nxg
          dist = abs(lscmlon - lon(i,1))
          if (dist < mind) then
             mind = dist
             ni = i
          endif
       enddo
       nj = -1
       mind = 1.0e20
       do j=1,nyg
          dist = abs(scmlat - lat(1,j))
          if (dist < mind) then
             mind = dist
             nj = j
          endif
       enddo
       j = nj
    else 
       ! lat and lon are on 1D arrays (e.g. spectral element grids)
       mind = 1.0e20
       do i=1,nxg
          dist=abs(lscmlon - lon(i,1)) + abs(scmlat - lat(i,1))
          if (dist < mind) then
             mind = dist
             ni = i
          endif
       enddo
       j = 1
    endif
    n = 1
    i = ni

    ! Initialize model ggrid for single column
    ! This will be used to interpolate the stream(s) to the single column value

    call mct_gGrid_init(ggrid=sdat%Grid, CoordChars='lat:lon', OtherChars='mask', lsize=1)
    allocate(idata(1)); idata = 1
    call mct_gGrid_importIAttr(sdat%Grid,'GlobGridNum', idata, 1)
    deallocate(idata)
    sdat%Grid%data%rAttr = -9999.0_R8
    nlon  = mct_aVect_indexRA(sdat%Grid%data,'lon')
    nlat  = mct_aVect_indexRA(sdat%Grid%data,'lat')
    nmask = mct_aVect_indexRA(sdat%Grid%data,'mask')

    sdat%Grid%data%rAttr(nlat ,1) = lat(i,j)
    sdat%Grid%data%rAttr(nlon ,1) = lon(i,j)
    sdat%Grid%data%rAttr(nmask,1) = 1
    deallocate(lon)
    deallocate(lat)

  end subroutine shr_strdata_init_model_domain_scol

  !===============================================================================
  subroutine shr_strdata_init_streams(SDAT, compid, mpicom, my_task)

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout) :: SDAT
    integer                ,intent(in)    :: compid
    integer                ,intent(in)    :: mpicom
    integer                ,intent(in)    :: my_task

    ! local variables
    integer          :: n,m,k    ! generic index
    integer          :: nfiles
    integer, pointer :: dof(:)
    integer, parameter :: master_task = 0
    character(len=*), parameter :: subname = "(shr_strdata_init_streams) "
    !-------------------------------------------------------------------------------

    ! Count streams again in case user made changes
    if (my_task == master_task) then
       do n=1,nStrMax

          ! check if a streams string is defined in strdata
          if (trim(SDAT%streams(n)) /= trim(shr_strdata_nullstr)) then
             SDAT%nstreams = max(SDAT%nstreams,n)
          end if

          ! check if a filename is defined in the stream
          call shr_stream_getNFiles(SDAT%stream(n),nfiles)
          if (nfiles > 0) then
             SDAT%nstreams = max(SDAT%nstreams,n)
          end if
          if (trim(SDAT%taxMode(n)) == trim(shr_stream_taxis_extend)) then
             SDAT%dtlimit(n) = 1.0e30
          end if
       end do

       ! Determine vector size for stream n
       SDAT%nvectors = 0
       do n = 1,nVecMax
          if (trim(SDAT%vectors(n)) /= trim(shr_strdata_nullstr)) then
             SDAT%nvectors = n
          end if
       end do
    endif
    call shr_mpi_bcast(SDAT%nstreams  ,mpicom, 'nstreams')
    call shr_mpi_bcast(SDAT%nvectors  ,mpicom, 'nvectors')
    call shr_mpi_bcast(SDAT%dtlimit   ,mpicom, 'dtlimit')

    ! Initialize domain, pio and calendar for each stream
    do n = 1,SDAT%nstreams

       ! Initialize stream domain info for stream n
       call shr_strdata_init_stream_domain(SDAT%stream(n), compid, mpicom, &
            SDAT%gridR(n), SDAT%gsmapR(n), SDAT%strnxg(n), SDAT%strnyg(n), SDAT%strnzg(n), SDAT%lsizeR(n))

       ! Initialize pio settings for stream n
       call mct_gsmap_OrderedPoints(SDAT%gsmapR(n), my_task, dof)
       if (SDAT%strnzg(n) <= 0) then
          call pio_initdecomp(SDAT%pio_subsystem, pio_double, (/SDAT%strnxg(n),SDAT%strnyg(n)/), &
               dof, SDAT%pio_iodesc(n))
       else
          call pio_initdecomp(SDAT%pio_subsystem, pio_double, (/SDAT%strnxg(n),SDAT%strnyg(n),SDAT%strnzg(n)/), &
               dof, SDAT%pio_iodesc(n))
       endif
       deallocate(dof)

       ! Initialize calendar for stream n
       call shr_mpi_bcast(SDAT%stream(n)%calendar, mpicom)
    enddo

  end subroutine shr_strdata_init_streams

  !===============================================================================

  subroutine shr_strdata_getfrac_from_stream(sdat, mpicom, my_task, fracname, fracdata)

    ! The following initializes the fracname from the first stream
    ! This is applicable for data models that read the model domain from the
    ! domain of the first stream; we read the data model's domain fraction from the
    ! first stream file, and this variable provides the name of the frac field on this
    ! file. Also see ESMCI/cime#2515).
    ! Note: this assumes that the first stream is on the same mesh as the model
    ! Note: This also assumes that the field with fracname exists on the domain file

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout) :: sdat
    integer                ,intent(in)    :: mpicom
    integer                ,intent(in)    :: my_task
    character(len=*)       ,intent(in)    :: fracname
    real(r8)               ,pointer       :: fracdata(:)

    ! local variables
    character(CXX)     :: domainfile
    integer            :: lsize
    integer            :: gsize
    integer            :: fid
    integer            :: rcode
    real(r8), pointer  :: data2d(:,:)
    type(mct_avect)    :: avG
    type(mct_avect)    :: avtmp
    integer            :: ierr
    integer            :: nx,ny
    integer            :: i,j
    integer            :: master_task = 0
    character(len=*), parameter :: subname = "(shr_strdata_init_streams) "
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom, my_task, ierr)
    master_task = 0

    if (my_task == master_task) then
       call shr_stream_getDomainFile(sdat%stream(1), domainfile)
       call shr_ncread_varDimSizes(trim(domainfile), trim(fracname), n1=nx, n2=ny)
       gsize = mct_gsmap_gsize(sdat%gsMap)
       if (gsize == nx*ny) then
          allocate(data2d(nx,ny))
       else
          write(logunit,*) "ERROR in nx,ny,gsize data sizes ",nx,ny,gsize
          call shr_sys_abort(subname//"ERROR in data sizes")
       endif
       call shr_ncread_open(trim(domainfile), fid, rCode)
       call shr_ncread_tField(domainfile, 1, fracname, data2d, fidi=fid, rc=rCode)
       call mct_aVect_init(avG, rlist=trim(fracname), lsize=gsize)
       avG%rAttr(1,:) = reshape(data2d, (/gsize/))
       call shr_ncread_close(fid, rCode)
       deallocate(data2d)
    end if
    ! The following call creates avtmp
    call mct_aVect_scatter(avG, avtmp, sdat%gsMap, master_task, mpicom)
    fracdata(:) = avtmp%rattr(1,:)
    call mct_aVect_clean(avtmp)
    if (my_task == master_task) call mct_aVect_clean(avG)

  end subroutine shr_strdata_getfrac_from_stream

  !===============================================================================
  subroutine shr_strdata_init_stream_domain(stream, compid, mpicom, &
       gGrid, gsMap, nxg, nyg, nzg, lsize)

    !----------------------------------------------------------------------------
    ! Create mct ggrid for model grid and set model gsmap if not input
    ! o assumes a very specific netCDF domain file format wrt var names, etc.
    !----------------------------------------------------------------------------

    ! input/output variables
    type(shr_stream_streamType) , intent(in)    :: stream
    integer                     , intent(in)    :: compid
    integer                     , intent(in)    :: mpicom
    type(mct_gGrid)             , intent(inout) :: gGrid
    type(mct_gsMap)             , intent(inout) :: gsMap
    integer                     , intent(out)   :: nxg
    integer                     , intent(out)   :: nyg
    integer                     , intent(out)   :: nzg
    integer                     , intent(out)   :: lsize

    ! local variables
    character(CL)        :: filePath     ! file path of stream domain file
    character(CXX)       :: filename     ! file name of stream domain file
    character(CS)        :: timeName     ! domain file: time variable name
    character(CS)        :: lonname      ! name of  lon variable in file
    character(CS)        :: latname      ! name of  lat variable in file
    character(CS)        :: hgtname      ! name of  hgt variable in file
    character(CS)        :: maskname     ! name of mask variable in file
    integer              :: n,k,j,i      ! indices
    integer              :: gsize        ! gsize
    integer              :: my_task
    integer              :: master_task
    integer              :: ierr         ! error code
    logical              :: fileexists   !
    logical              :: maskexists   ! is mask on dataset
    integer              :: ndims        ! number of dims
    integer              :: nlon
    integer              :: nlat
    integer              :: nmask
    integer              :: nhgt
    integer              :: ni,nj        ! scmmode point search
    real(R8),allocatable :: lon(:,:)     ! temp array for domain lon  info
    real(R8),allocatable :: lat(:,:)     ! temp array for domain lat  info
    integer,allocatable  :: mask(:,:)    ! temp array for domain mask info
    real(R8),allocatable :: hgt(:)       ! temp array for domain height info
    real(R8),allocatable :: a4d(:,:,:,:) ! temp array for reading generic stuff
    integer, pointer     :: idata(:)     ! temporary
    type(mct_ggrid)      :: gGridRoot    ! global mct ggrid
    character(*), parameter :: subname = '(shr_strdata_readgrid_stream) '
    character(*), parameter :: F00   = "('(shr_strdata_readgrid_stream) ',8a)"
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    ! Determine stream filename and namelist of domain variables

    if (my_task == master_task) then
       call shr_stream_getDomainInfo(stream, filePath, fileName,&
            timeName, lonName, latName, hgtName, maskName)
       call shr_stream_getFile(filePath, fileName)
       write(logunit,*) subname,' stream data'
       write(logunit,*) subname,' filePath = ',trim(filePath)
       write(logunit,*) subname,' fileName = ',trim(fileName)
       write(logunit,*) subname,' timeName = ',trim(timeName)
       write(logunit,*) subname,' lonName  = ',trim(lonName)
       write(logunit,*) subname,' latName  = ',trim(latName)
       write(logunit,*) subname,' hgtName  = ',trim(hgtName)
       write(logunit,*) subname,' maskName = ',trim(maskName)
    endif
    call shr_mpi_bcast(fileName ,mpicom)
    call shr_mpi_bcast(lonName  ,mpicom)
    call shr_mpi_bcast(latName  ,mpicom)
    call shr_mpi_bcast(hgtName  ,mpicom)
    call shr_mpi_bcast(maskName ,mpicom)

    if (my_task == master_task) then
       inquire(file=trim(fileName), exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif

    ! Determine stream lat/lon (degrees), area (radians^2), mask ( 1 for ocean, 0 for non-ocean)

    if (my_task == master_task) then
       if (shr_ncread_varexists(fileName, maskname)) then
          maskexists = .true.
          call shr_ncread_varDimSizes(fileName,maskname,nxg,nyg)
       else
          maskexists = .false.
          call shr_ncread_varDimNum(fileName,lonName,ndims)
          if (ndims == 1) then
             call shr_ncread_varDimSizes(fileName,lonName,nxg)
             call shr_ncread_varDimSizes(fileName,latName,nyg)
          else
             call shr_ncread_varDimSizes(fileName,lonName,nxg,nyg)
          endif
       endif
       if (shr_ncread_varexists(fileName,hgtName)) then
          call shr_ncread_varDimSizes(fileName,hgtname,nzg)
       else
          nzg = -1
       endif
    endif
    call shr_mpi_bcast(nxg,mpicom)
    call shr_mpi_bcast(nyg,mpicom)
    call shr_mpi_bcast(nzg,mpicom)
    gsize = abs(nxg*nyg*nzg)

    ! Create stream gsmap using 1d decomp of 2d grid plus 3rd dim if exists

    if (my_task == master_task) then
       write(logunit,*)trim(subname) // ': Creating gsmap for input stream'
    end if
    call shr_strdata_gsmapCreate(gsMap, nxg, nyg, nzg, compid, mpicom)
    lsize = mct_gsmap_lsize(gsmap, mpicom)

    ! Create stream ggrid

    lsize = mct_gsMap_lsize(gsMap, mpicom)
    call mct_gGrid_init(GGrid=Ggrid, CoordChars='lat:lon:hgt', OtherChars='mask', lsize=lsize )
    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(gGrid,'GlobGridNum',idata,lsize)
    deallocate(idata)
    gGrid%data%rAttr = -9999.0_R8

    if (my_task == master_task) then
       allocate(lon(nxg,nyg))
       allocate(lat(nxg,nyg))
       allocate(mask(nxg,nyg))
       allocate(hgt(abs(nzg)))
       if (.not. maskexists) then
          call shr_ncread_domain(fileName, lonname, lon, latname, lat)
          mask = 1
       else
          call shr_ncread_domain(fileName, lonname, lon, latname, lat, maskName=maskname, mask=mask)
       endif
       if (nzg > 1) then
          allocate(a4d(nzg, 1, 1, 1))
          call shr_ncread_field4dG(fileName, hgtName, rfld=a4d)
          hgt(:) = a4d(:, 1, 1, 1)
          deallocate(a4d)
       else
          hgt = 1
       endif
       call mct_gGrid_init(gGridRoot, gGrid, gsize)
       gGridRoot%data%rAttr = -9999.0_R8 !to avoid errors when using strict compiler checks
       nlon  = mct_aVect_indexRA(gGridRoot%data, 'lon')
       nlat  = mct_aVect_indexRA(gGridRoot%data, 'lat')
       nmask = mct_aVect_indexRA(gGridRoot%data, 'mask')
       nhgt  = mct_aVect_indexRA(gGridRoot%data, 'hgt')
       n = 0
       do k = 1,abs(nzg)
          do j = 1,nyg
             do i = 1,nxg
                n = n+1
                gGridRoot%data%rAttr(nlat ,n) = lat(i,j)
                gGridRoot%data%rAttr(nlon ,n) = lon(i,j)
                gGridRoot%data%rAttr(nmask,n) = real(mask(i,j),R8)
                gGridRoot%data%rAttr(nhgt ,n) = hgt(k)
             enddo
          enddo
       enddo
       deallocate(lon)
       deallocate(lat)
       deallocate(mask)
       deallocate(hgt)
    endif
    call mct_gGrid_scatter(gGridRoot, gGrid, gsMap, master_task, mpicom)
    if (my_task == master_task) then
       call mct_gGrid_clean(gGridRoot)
    end if

  contains

    subroutine shr_strdata_gsmapCreate(gsmap, nxg, nyg, nzg, compid, mpicom)

      ! input/output variables
      type(mct_gsMap) , intent(inout) :: gsmap
      integer         , intent(in)    :: nxg,nyg,nzg
      integer         , intent(in)    :: compid
      integer         , intent(in)    :: mpicom

      ! local
      integer :: n,nz,nb,npes,ierr,gsize,dsize,ngseg,lnzg
      integer, pointer :: start(:)     ! for gsmap initialization
      integer, pointer :: length(:)    ! for gsmap initialization
      integer, pointer :: pe_loc(:)    ! for gsmap initialization
      ! ---------------------------------------------

      gsize = abs(nxg*nyg*nzg)
      dsize = nxg*nyg
      lnzg = 1
      if (nzg > 1) lnzg = nzg  ! check for 3d

      call mpi_comm_size(mpicom,npes,ierr)

      !--- 1d decomp of 2d grid plus 3rd dim if exists ---
      ngseg = npes*lnzg
      allocate(start(ngseg),length(ngseg),pe_loc(ngseg))
      start = 0
      length = 0
      pe_loc = 0
      do n = 1,npes
         length(n)  = dsize/npes
         if (n <= mod(dsize,npes)) then
            length(n) = length(n) + 1
         end if
         if (n == 1) then
            start(n) = 1
         else
            start(n) = start(n-1) + length(n-1)
         endif
         pe_loc(n) = n-1
         do nz = 2,lnzg
            nb = (nz-1)*npes + n
            start(nb)  = start(n) + (nz-1)*dsize
            length(nb) = length(n)
            pe_loc(nb) = pe_loc(n)
         enddo
      enddo
      call mct_gsmap_init( gsmap, compid, ngseg, gsize, start, length, pe_loc)
      deallocate(start,length,pe_loc)

    end subroutine shr_strdata_gsmapCreate

  end subroutine shr_strdata_init_stream_domain

  !===============================================================================
  subroutine shr_strdata_init_mapping(SDAT, compid, mpicom, my_task)

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout) :: SDAT
    integer                ,intent(in)    :: compid
    integer                ,intent(in)    :: mpicom
    integer                ,intent(in)    :: my_task

    ! local variables
    type(mct_sMat) :: sMati
    integer        :: n,m,k   ! generic index
    integer        :: nu,nv   ! u,v index
    integer        :: method  ! mapping method
    character(CXX) :: fldList ! list of fields
    character(CS)  :: uname   ! u vector field name
    character(CS)  :: vname   ! v vector field name
    logical        :: compare1
    logical        :: compare2
    integer          ,parameter :: master_task = 0
    character(*)     ,parameter :: F00 = "('(shr_strdata_init_mapping) ',8a)"
    character(len=*) ,parameter :: subname = "(shr_strdata_init_mapping) "
    !-------------------------------------------------------------------------------

    do n = 1,SDAT%nstreams

       method = CompareMaskSubset
       if ( shr_dmodel_gGridCompare(SDAT%gridR(n), SDAT%gsmapR(n), SDAT%grid,SDAT%gsmap, method, mpicom) .or. &
            trim(SDAT%fillalgo(n))=='none') then
          SDAT%dofill(n) = .false.
       else
          SDAT%dofill(n) = .true.
       endif

       if (trim(SDAT%mapmask(n)) == 'dstmask') then
          method = CompareXYabsMask
       else
          method = CompareXYabs
       endif
       if ( shr_dmodel_gGridCompare(SDAT%gridR(n),SDAT%gsmapR(n),SDAT%grid,SDAT%gsmap, method, mpicom, 0.01_r8) .or. &
            trim(SDAT%mapalgo(n))=='none') then
          SDAT%domaps(n) = .false.
       else
          SDAT%domaps(n) = .true.
       endif

       ! Set up fills
       if (SDAT%dofill(n)) then
          if (SDAT%strnzg(n) > 1) then
             write(logunit,*) trim(subname),' do fill called with 3d data, not allowed'
             call shr_sys_abort(subname//': do fill called with 3d data, not allowed')
          endif

          if (trim(SDAT%fillread(n)) == trim(shr_strdata_unset)) then
             if (my_task == master_task) then
                write(logunit,F00) ' calling shr_dmodel_mapSet for fill'
                call shr_sys_flush(logunit)
             endif

             ! Set up fill
             call shr_dmodel_mapSet(SDAT%sMatPf(n), &
                  SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
                  SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
                  name='mapFill', &
                  type='cfill', &
                  algo=trim(SDAT%fillalgo(n)),&
                  mask=trim(SDAT%fillmask(n)),&
                  vect='scalar', &
                  compid=compid,&
                  mpicom=mpicom)

             if (trim(SDAT%fillwrit(n)) /= trim(shr_strdata_unset)) then
                if (my_task == master_task) then
                   write(logunit,F00) ' writing ',trim(SDAT%fillwrit(n))
                   call shr_sys_flush(logunit)
                endif

                call shr_mct_sMatWritednc(&
                     SDAT%sMatPf(n)%Matrix,&
                     SDAT%pio_subsystem,&
                     SDAT%io_type,&
                     SDAT%io_format, &
                     SDAT%fillwrit(n),&
                     compid,&
                     mpicom)
             endif
          else
             if (my_task == master_task) then
                write(logunit,F00) ' reading ',trim(SDAT%fillread(n))
                call shr_sys_flush(logunit)
             endif
             call shr_mct_sMatReaddnc(sMati,SDAT%gsmapR(n),SDAT%gsmapR(n),'src', &
                  filename=trim(SDAT%fillread(n)),mytask=my_task,mpicom=mpicom)

             call mct_sMatP_Init(SDAT%sMatPf(n),sMati,SDAT%gsMapR(n),&
                  SDAT%gsmapR(n),0, mpicom, compid)

             call mct_sMat_Clean(sMati)
          endif
       endif

       ! Set up maps
       if (SDAT%domaps(n)) then
          if (SDAT%strnzg(n) > 1) then
             write(logunit,*) trim(subname),' do maps called with 3d data, not allowed'
             call shr_sys_abort(subname//': do maps called with 3d data, not allowed')
          endif

          if (trim(SDAT%mapread(n)) == trim(shr_strdata_unset)) then
             if (my_task == master_task) then
                write(logunit,F00) ' calling shr_dmodel_mapSet for remap'
                call shr_sys_flush(logunit)
             endif

             call shr_dmodel_mapSet(SDAT%sMatPs(n), &
                  SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
                  SDAT%grid    ,SDAT%gsmap    ,SDAT%nxg      ,SDAT%nyg,       &
                  name='mapScalar', &
                  type='remap',                             &
                  algo=trim(SDAT%mapalgo(n)),&
                  mask=trim(SDAT%mapmask(n)), &
                  vect='scalar', &
                  compid=compid,&
                  mpicom=mpicom)

             if (trim(SDAT%mapwrit(n)) /= trim(shr_strdata_unset)) then
                if (my_task == master_task) then
                   write(logunit,F00) ' writing ',trim(SDAT%mapwrit(n))
                   call shr_sys_flush(logunit)
                endif
                call shr_mct_sMatWritednc(&
                     SDAT%sMatPs(n)%Matrix,&
                     sdat%pio_subsystem,&
                     sdat%io_type,&
                     SDAT%io_format,&
                     SDAT%mapwrit(n),&
                     compid,&
                     mpicom)
             endif
          else
             if (my_task == master_task) then
                write(logunit,F00) ' reading ',trim(SDAT%mapread(n))
                call shr_sys_flush(logunit)
             endif
             call shr_mct_sMatReaddnc(sMati,SDAT%gsmapR(n),SDAT%gsmap,'src', &
                  filename=trim(SDAT%mapread(n)),mytask=my_task,mpicom=mpicom)
             call mct_sMatP_Init(SDAT%sMatPs(n),sMati,SDAT%gsMapR(n),SDAT%gsmap,0, mpicom, compid)
             call mct_sMat_Clean(sMati)
          endif
       else
          call mct_rearr_init(SDAT%gsmapR(n), SDAT%gsmap, mpicom, SDAT%rearrR(n))
       endif
    enddo

    ! --- setup datatypes ---

    do n = 1,SDAT%nstreams
       if (my_task == master_task) then
          call shr_stream_getModelFieldList(SDAT%stream(n),fldList)
       endif
       call shr_mpi_bcast(fldList,mpicom)
       call mct_aVect_init(SDAT%avs(n)  ,rlist=fldList,lsize=SDAT%lsize)
       call mct_aVect_init(SDAT%avFLB(n),rlist=fldList,lsize=SDAT%lsize)
       call mct_aVect_init(SDAT%avFUB(n),rlist=fldList,lsize=SDAT%lsize)
       call mct_aVect_init(SDAT%avRLB(n),rlist=fldList,lsize=SDAT%lsizeR(n))
       call mct_aVect_init(SDAT%avRUB(n),rlist=fldList,lsize=SDAT%lsizeR(n))
       if (trim(SDAT%tintalgo(n)) == 'coszen') then
          call mct_aVect_init(SDAT%avCoszen(n),rlist="tavCosz",lsize=SDAT%lsize)
       endif
    enddo

    ! --- check vectors and compute ustrm,vstrm ---

    do m = 1,SDAT%nvectors
       if (.not. shr_string_listIsValid(SDAT%vectors(m))) then
          write(logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(SDAT%vectors(m))
          call shr_sys_abort(subname//': vec fldlist invalid:'//trim(SDAT%vectors(m)))
       endif
       if (shr_string_listGetNum(SDAT%vectors(m)) /= 2) then
          write(logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(SDAT%vectors(m))
          call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(SDAT%vectors(m)))
       endif
       call shr_string_listGetName(SDAT%vectors(m),1,uname)
       call shr_string_listGetName(SDAT%vectors(m),2,vname)
       nu = 0
       nv = 0
       do n = 1,SDAT%nstreams
          k = mct_aVect_indexRA(SDAT%avRLB(n),trim(uname),perrWith='quiet')
          if (k > 0) nu = n
          k = mct_aVect_indexRA(SDAT%avRLB(n),trim(vname),perrWith='quiet')
          if (k > 0) nv = n
       enddo
       if (nu == 0  .or. nv == 0) then
          write(logunit,*) trim(subname),' vec flds not found  m=',m,trim(SDAT%vectors(m))
          call shr_sys_abort(subname//': vec flds not found:'//trim(SDAT%vectors(m)))
       endif
       if (nu /= nv) then
          compare1 = shr_dmodel_gGridCompare(SDAT%gridR(nu), SDAT%gsmapR(nu), SDAT%gridR(nv),SDAT%gsmapR(nv), &
               CompareXYabs, mpicom, 0.01_r8)
          compare2 = shr_dmodel_gGridCompare(SDAT%gridR(nu), SDAT%gsmapR(nu), SDAT%gridR(nv),SDAT%gsmapR(nv), &
               CompareMaskZeros, mpicom)
          if ((.not. compare1) .or. (.not. compare2)) then
             write(logunit,*) trim(subname),' vec fld doms not same m=',m,trim(SDAT%vectors(m))
             call shr_sys_abort(subname//': vec fld doms not same:'//trim(SDAT%vectors(m)))
          endif
       endif
       SDAT%ustrm(m) = nu
       SDAT%vstrm(m) = nv
    enddo

  end subroutine shr_strdata_init_mapping

  !===============================================================================

  subroutine shr_strdata_advance(SDAT,ymd,tod,mpicom,istr,timers)

    type(shr_strdata_type) ,intent(inout)       :: SDAT
    integer                ,intent(in)          :: ymd    ! current model date
    integer                ,intent(in)          :: tod    ! current model date
    integer                ,intent(in)          :: mpicom
    character(len=*)       ,intent(in),optional :: istr
    logical                ,intent(in),optional :: timers

    ! local variables
    integer                    :: n,m,i,kf  ! generic index
    integer                    :: my_task,npes
    integer,parameter          :: master_task = 0
    logical                    :: mssrmlf
    logical,allocatable        :: newData(:)
    integer                    :: ierr
    integer                    :: nu,nv
    integer                    :: lsize,lsizeR,lsizeF
    integer,allocatable        :: ymdmod(:) ! modified model dates to handle Feb 29
    integer                    :: todmod    ! modified model dates to handle Feb 29
    type(mct_avect)            :: avRtmp
    type(mct_avect)            :: avRV,avFV
    character(len=32)          :: lstr
    logical                    :: ltimers
    real(R8)                   :: flb,fub   ! factor for lb and ub

    !--- for cosz method ---
    real(R8),pointer           :: lonr(:)              ! lon radians
    real(R8),pointer           :: latr(:)              ! lat radians
    real(R8),pointer           :: cosz(:)              ! cosz
    real(R8),pointer           :: tavCosz(:)           ! cosz, time avg over [LB,UB]
    real(R8),pointer           :: xlon(:),ylon(:)
    real(R8),parameter         :: solZenMin = 0.001_R8 ! minimum solar zenith angle
    type(ESMF_Time)            :: timeLB, timeUB       ! lb and ub times
    type(ESMF_TimeInterval)    :: timeint              ! delta time
    integer                    :: dday                 ! delta days
    real(R8)                   :: dtime                ! delta time
    integer                    :: uvar,vvar
    character(CS)              :: uname                ! u vector field name
    character(CS)              :: vname                ! v vector field name
    integer                    :: year,month,day       ! date year month day
    character(len=*),parameter :: timname = "_strd_adv"
    integer,parameter          :: tadj = 2
    character(*),parameter     :: subname = "(shr_strdata_advance) "
    !-------------------------------------------------------------------------------

    if (SDAT%nstreams < 1) return

    lstr = ''
    if (present(istr)) then
       lstr = trim(istr)
    endif

    ltimers = .true.
    if (present(timers)) then
       ltimers = timers
    endif

    if (.not.ltimers) call t_adj_detailf(tadj)

    call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
    call t_startf(trim(lstr)//trim(timname)//'_total')

    call MPI_COMM_SIZE(mpicom,npes,ierr)
    call MPI_COMM_RANK(mpicom,my_task,ierr)

    mssrmlf = .false.

    SDAT%ymd = ymd
    SDAT%tod = tod

    if (SDAT%nstreams > 0) then
       allocate(newData(SDAT%nstreams))
       allocate(ymdmod(SDAT%nstreams))

       do n = 1,SDAT%nstreams
          ! ------------------------------------------------------- !
          ! tcraig, Oct 11 2010.  Mismatching calendars: 4 cases    !
          ! ------------------------------------------------------- !
          ! ymdmod and todmod are the ymd and tod to time           !
          ! interpolate to.  Generally, these are just the model    !
          ! date and time.  Also, always use the stream calendar    !
          ! for time interpolation for reasons described below.     !
          ! When there is a calendar mismatch, support Feb 29 in a  !
          ! special way as needed to get reasonable values.         !
          ! Note that when Feb 29 needs to be treated specially,    !
          ! a discontinuity will be introduced.  The size of that   !
          ! discontinuity will depend on the time series input data.!
          ! ------------------------------------------------------- !
          ! (0) The stream calendar and model calendar are          !
          ! identical.  Proceed in the standard way.                !
          ! ------------------------------------------------------- !
          ! (1) If the stream is a no leap calendar and the model   !
          ! is gregorian, then time interpolate on the noleap       !
          ! calendar.  Then if the model date is Feb 29, compute    !
          ! stream data for Feb 28 by setting ymdmod and todmod to  !
          ! Feb 28.  This results in duplicate stream data on       !
          ! Feb 28 and Feb 29 and a discontinuity at the start of   !
          ! Feb 29.                                                 !
          ! This could be potentially updated by using the gregorian!
          ! calendar for time interpolation when the input data is  !
          ! relatively infrequent (say greater than daily) with the !
          ! following concerns.
          !   - The forcing will not be reproduced identically on   !
          !     the same day with climatological inputs data        !
          !   - Input data with variable input frequency might      !
          !     behave funny
          !   - An arbitrary discontinuity will be introduced in    !
          !     the time interpolation method based upon the        !
          !     logic chosen to transition from reproducing Feb 28  !
          !     on Feb 29 and interpolating to Feb 29.              !
          !   - The time gradient of data will change by adding a   !
          !     day arbitrarily.
          ! ------------------------------------------------------- !
          ! (2) If the stream is a gregorian calendar and the model !
          ! is a noleap calendar, then just time interpolate on the !
          ! gregorian calendar.  The causes Feb 29 stream data      !
          ! to be skipped and will lead to a discontinuity at the   !
          ! start of March 1.                                       !
          ! ------------------------------------------------------- !
          ! (3) If the calendars mismatch and neither of the three  !
          ! cases above are recognized, then abort.                 !
          ! ------------------------------------------------------- !

          ! case(0)
          ymdmod(n) = ymd
          todmod    = tod
          if (trim(SDAT%calendar) /= trim(SDAT%stream(n)%calendar)) then
             if ((trim(SDAT%calendar) == trim(shr_cal_gregorian)) .and. &
                  (trim(SDAT%stream(n)%calendar) == trim(shr_cal_noleap))) then
                ! case (1), set feb 29 = feb 28
                call shr_cal_date2ymd (ymd,year,month,day)
                if (month == 2 .and. day == 29) then
                   call shr_cal_ymd2date(year,2,28,ymdmod(n))
                endif
             else if ((trim(SDAT%calendar) == trim(shr_cal_noleap)) .and. &
                  (trim(SDAT%stream(n)%calendar) == trim(shr_cal_gregorian))) then
                ! case (2), feb 29 input data will be skipped automatically
             else
                ! case (3), abort
                write(logunit,*) trim(subname),' ERROR: mismatch calendar ', &
                     trim(SDAT%calendar),':',trim(SDAT%stream(n)%calendar)
                call shr_sys_abort(trim(subname)//' ERROR: mismatch calendar ')
             endif
          endif

          call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
          call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

          call shr_dmodel_readLBUB(SDAT%stream(n),SDAT%pio_subsystem,SDAT%io_type,SDAT%pio_iodesc(n), &
               ymdmod(n),todmod,mpicom,SDAT%gsmapR(n),&
               SDAT%avRLB(n),SDAT%ymdLB(n),SDAT%todLB(n), &
               SDAT%avRUB(n),SDAT%ymdUB(n),SDAT%todUB(n), &
               SDAT%avRFile(n), trim(SDAT%readmode(n)), newData(n), &
               istr=trim(lstr)//'_readLBUB')

          if (debug > 0) then
             write(logunit,*) trim(subname),' newData flag = ',n,newData(n)
             write(logunit,*) trim(subname),' LB ymd,tod = ',n,SDAT%ymdLB(n),SDAT%todLB(n)
             write(logunit,*) trim(subname),' UB ymd,tod = ',n,SDAT%ymdUB(n),SDAT%todUB(n)
          endif

          if (newData(n)) then
             if (debug > 0) then
                write(logunit,*) trim(subname),' newData RLB = ',n,minval(SDAT%avRLB(n)%rAttr), &
                     maxval(SDAT%avRLB(n)%rAttr),sum(SDAT%avRLB(n)%rAttr)
                write(logunit,*) trim(subname),' newData RUB = ',n,minval(SDAT%avRUB(n)%rAttr), &
                     maxval(SDAT%avRUB(n)%rAttr),sum(SDAT%avRUB(n)%rAttr)
             endif
             call shr_cal_date2ymd(SDAT%ymdLB(n),year,month,day)
             call shr_cal_timeSet(timeLB,SDAT%ymdLB(n),0,SDAT%stream(n)%calendar)
             call shr_cal_timeSet(timeUB,SDAT%ymdUB(n),0,SDAT%stream(n)%calendar)
             timeint = timeUB-timeLB
             call ESMF_TimeIntervalGet(timeint,StartTimeIn=timeLB,d=dday)
             dtime = abs(real(dday,R8) + real(SDAT%todUB(n)-SDAT%todLB(n),R8)/shr_const_cDay)

             SDAT%dtmin(n) = min(SDAT%dtmin(n),dtime)
             SDAT%dtmax(n) = max(SDAT%dtmax(n),dtime)
             if ((SDAT%dtmax(n)/SDAT%dtmin(n)) > SDAT%dtlimit(n)) then
                write(logunit,*) trim(subname),' ERROR: for stream ',n
                write(logunit,*) trim(subName),' ERROR: dt limit1 ',SDAT%dtmax(n),SDAT%dtmin(n),SDAT%dtlimit(n)
                write(logunit,*) trim(subName),' ERROR: dt limit2 ',SDAT%ymdLB(n),SDAT%todLB(n), &
                     SDAT%ymdUB(n),SDAT%todUB(n)
                call shr_sys_abort(trim(subName)//' ERROR dt limit for stream')
             endif
          endif
          call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')
       enddo

       do n = 1,SDAT%nstreams
          if (newData(n)) then

             if (SDAT%doFill(n)) then
                call t_startf(trim(lstr)//trim(timname)//'_fill')
                lsize = mct_aVect_lsize(SDAT%avRLB(n))
                call mct_aVect_init(avRtmp,SDAT%avRLB(n),lsize)
                call mct_aVect_copy(SDAT%avRLB(n),avRtmp)
                call mct_sMat_avMult(avRtmp,SDAT%sMatPf(n),SDAT%avRLB(n))
                call mct_aVect_copy(SDAT%avRUB(n),avRtmp)
                call mct_sMat_avMult(avRtmp,SDAT%sMatPf(n),SDAT%avRUB(n))
                call mct_aVect_clean(avRtmp)
                call t_stopf(trim(lstr)//trim(timname)//'_fill')
             endif

             if (SDAT%domaps(n)) then
                call t_startf(trim(lstr)//trim(timname)//'_map')
                call mct_sMat_avMult(SDAT%avRLB(n),SDAT%sMatPs(n),SDAT%avFLB(n))
                call mct_sMat_avMult(SDAT%avRUB(n),SDAT%sMatPs(n),SDAT%avFUB(n))
                call t_stopf(trim(lstr)//trim(timname)//'_map')
             else
                call t_startf(trim(lstr)//trim(timname)//'_rearr')
                call mct_rearr_rearrange(SDAT%avRLB(n),SDAT%avFLB(n),SDAT%rearrR(n))
                call mct_rearr_rearrange(SDAT%avRUB(n),SDAT%avFUB(n),SDAT%rearrR(n))
                call t_stopf(trim(lstr)//trim(timname)//'_rearr')
             endif

             if (debug > 0) then
                write(logunit,*) trim(subname),' newData FLB = ',n,minval(SDAT%avFLB(n)%rAttr), &
                     maxval(SDAT%avFLB(n)%rAttr),sum(SDAT%avFLB(n)%rAttr)
                write(logunit,*) trim(subname),' newData FUB = ',n,minval(SDAT%avFUB(n)%rAttr), &
                     maxval(SDAT%avFUB(n)%rAttr),sum(SDAT%avFUB(n)%rAttr)
             endif
          endif
       enddo

       do m = 1,SDAT%nvectors
          nu = SDAT%ustrm(m)
          nv = SDAT%vstrm(m)
          if ((SDAT%domaps(nu) .or. SDAT%domaps(nv)) .and. &
               (newdata(nu) .or. newdata(nv))) then

             call t_startf(trim(lstr)//trim(timname)//'_vect')
             call shr_string_listGetName(SDAT%vectors(m),1,uname)
             call shr_string_listGetName(SDAT%vectors(m),2,vname)
             lsizeR = mct_aVect_lsize(SDAT%avRLB(nu))
             lsizeF = mct_aVect_lsize(SDAT%avFLB(nu))
             call mct_aVect_init(avRV,rlist=SDAT%vectors(m),lsize=lsizeR)
             call mct_aVect_init(avFV,rlist=SDAT%vectors(m),lsize=lsizeF)
             allocate(xlon(lsizeR))
             allocate(ylon(lsizeF))
             call mct_aVect_exportRattr(SDAT%gridR(nu)%data,'lon',xlon)
             call mct_aVect_exportRattr(SDAT%grid     %data,'lon',ylon)
             xlon = xlon * deg2rad
             ylon = ylon * deg2rad

             !--- map LB ---

             uvar = mct_aVect_indexRA(SDAT%avRLB(nu),trim(uname))
             vvar = mct_aVect_indexRA(SDAT%avRLB(nv),trim(vname))
             do i = 1,lsizeR
                avRV%rAttr(1,i) =  SDAT%avRLB(nu)%rAttr(uvar,i) * cos(xlon(i))  &
                     -SDAT%avRLB(nv)%rAttr(vvar,i) * sin(xlon(i))
                avRV%rAttr(2,i) =  SDAT%avRLB(nu)%rAttr(uvar,i) * sin(xlon(i))  &
                     +SDAT%avRLB(nv)%rAttr(vvar,i) * cos(xlon(i))
             enddo
             call mct_sMat_avMult(avRV,SDAT%sMatPs(nu),avFV)
             ! ---   don't need to recompute uvar and vvar, should be the same
             !       uvar = mct_aVect_indexRA(SDAT%avFLB(nu),trim(uname))
             !       vvar = mct_aVect_indexRA(SDAT%avFLB(nv),trim(vname))
             do i = 1,lsizeF
                SDAT%avFLB(nu)%rAttr(uvar,i) =  avFV%rAttr(1,i) * cos(ylon(i))  &
                     +avFV%rAttr(2,i) * sin(ylon(i))
                SDAT%avFLB(nv)%rAttr(vvar,i) = -avFV%rAttr(1,i) * sin(ylon(i))  &
                     +avFV%rAttr(2,i) * cos(ylon(i))
             enddo

             !--- map UB ---

             uvar = mct_aVect_indexRA(SDAT%avRUB(nu),trim(uname))
             vvar = mct_aVect_indexRA(SDAT%avRUB(nv),trim(vname))
             do i = 1,lsizeR
                avRV%rAttr(1,i) =  SDAT%avRUB(nu)%rAttr(uvar,i) * cos(xlon(i))  &
                     -SDAT%avRUB(nv)%rAttr(vvar,i) * sin(xlon(i))
                avRV%rAttr(2,i) =  SDAT%avRUB(nu)%rAttr(uvar,i) * sin(xlon(i))  &
                     +SDAT%avRUB(nv)%rAttr(vvar,i) * cos(xlon(i))
             enddo
             call mct_sMat_avMult(avRV,SDAT%sMatPs(nu),avFV)
             ! ---   don't need to recompute uvar and vvar, should be the same
             !       uvar = mct_aVect_indexRA(SDAT%avFUB(nu),trim(uname))
             !       vvar = mct_aVect_indexRA(SDAT%avFUB(nv),trim(vname))
             do i = 1,lsizeF
                SDAT%avFUB(nu)%rAttr(uvar,i) =  avFV%rAttr(1,i) * cos(ylon(i))  &
                     +avFV%rAttr(2,i) * sin(ylon(i))
                SDAT%avFUB(nv)%rAttr(vvar,i) = -avFV%rAttr(1,i) * sin(ylon(i))  &
                     +avFV%rAttr(2,i) * cos(ylon(i))
             enddo

             call mct_aVect_clean(avRV)
             call mct_aVect_clean(avFV)
             deallocate(xlon,ylon)

             call t_stopf(trim(lstr)//trim(timname)//'_vect')
          endif
       enddo

       do n = 1,SDAT%nstreams

          !--- method: coszen -------------------------------------------------------
          if (trim(SDAT%tintalgo(n)) == 'coszen') then
             call t_startf(trim(lstr)//trim(timname)//'_coszen')

             !--- make sure orb info has been set ---
             if (SDAT%eccen == SHR_ORB_UNDEF_REAL) then
                call shr_sys_abort(subname//' ERROR in orb params for coszen tinterp')
             else if (SDAT%modeldt < 1) then
                call shr_sys_abort(subname//' ERROR: model dt < 1 for coszen tinterp')
             endif

             !--- allocate avg cosz array ---
             lsizeF = mct_aVect_lsize(SDAT%avFLB(n))
             allocate(tavCosz(lsizeF),cosz(lsizeF),lonr(lsizeF),latr(lsizeF))

             !--- get lat/lon data ---
             kf = mct_aVect_indexRA(SDAT%grid%data,'lat')
             latr(:) = SDAT%grid%data%rAttr(kf,:) * deg2rad
             kf = mct_aVect_indexRA(SDAT%grid%data,'lon')
             lonr(:) = SDAT%grid%data%rAttr(kf,:) * deg2rad

             call t_startf(trim(lstr)//trim(timname)//'_coszenC')
             cosz = 0.0_r8
             call shr_tInterp_getCosz(cosz,lonr,latr,ymdmod(n),todmod, &
                  SDAT%eccen,SDAT%mvelpp,SDAT%lambm0,SDAT%obliqr,SDAT%stream(n)%calendar)
             call t_stopf(trim(lstr)//trim(timname)//'_coszenC')

             if (newdata(n)) then
                !--- compute a new avg cosz ---
                call t_startf(trim(lstr)//trim(timname)//'_coszenN')
                call shr_tInterp_getAvgCosz(tavCosz,lonr,latr, &
                     SDAT%ymdLB(n),SDAT%todLB(n), SDAT%ymdUB(n),SDAT%todUB(n), &
                     SDAT%eccen,SDAT%mvelpp,SDAT%lambm0,SDAT%obliqr,SDAT%modeldt,&
                     SDAT%stream(n)%calendar)
                call mct_avect_importRAttr(SDAT%avCoszen(n),'tavCosz',tavCosz,lsizeF)
                call t_stopf(trim(lstr)//trim(timname)//'_coszenN')
             else
                !--- reuse existing avg cosz ---
                call mct_avect_exportRAttr(SDAT%avCoszen(n),'tavCosz',tavCosz)
             endif

             !--- t-interp is LB data normalized with this factor: cosz/tavCosz ---
             do i = 1,lsizeF
                if (cosz(i) > solZenMin) then
                   SDAT%avs(n)%rAttr(:,i) = SDAT%avFLB(n)%rAttr(:,i)*cosz(i)/tavCosz(i)
                else
                   SDAT%avs(n)%rAttr(:,i) =  0._r8
                endif
             enddo
             deallocate(tavCosz,cosz,lonr,latr)
             call t_stopf(trim(lstr)//trim(timname)//'_coszen')

             !--- method: not coszen ---------------------------------------------------
          elseif (trim(SDAT%tintalgo(n)) /= trim(shr_strdata_nullstr)) then

             call t_startf(trim(lstr)//trim(timname)//'_tint')
             call shr_tInterp_getFactors(SDAT%ymdlb(n),SDAT%todlb(n),SDAT%ymdub(n),SDAT%todub(n), &
                  ymdmod(n),todmod,flb,fub, &
                  calendar=SDAT%stream(n)%calendar,algo=trim(SDAT%tintalgo(n)))
             if (debug > 0) then
                write(logunit,*) trim(subname),' interp = ',n,flb,fub
             endif
             SDAT%avs(n)%rAttr(:,:) = SDAT%avFLB(n)%rAttr(:,:)*flb + SDAT%avFUB(n)%rAttr(:,:)*fub
             call t_stopf(trim(lstr)//trim(timname)//'_tint')

          else
             call t_startf(trim(lstr)//trim(timname)//'_zero')
             call mct_avect_zero(SDAT%avs(n))
             call t_stopf(trim(lstr)//trim(timname)//'_zero')
          endif
          if (debug > 0) then
             write(logunit,*) trim(subname),' SDAT av = ',n,minval(SDAT%avs(n)%rAttr),&
                  maxval(SDAT%avs(n)%rAttr),sum(SDAT%avs(n)%rAttr)
          endif

       enddo

       deallocate(newData)
       deallocate(ymdmod)

    endif    ! nstreams > 0

    call t_stopf(trim(lstr)//trim(timname)//'_total')
    if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine shr_strdata_advance

  !===============================================================================

  subroutine shr_strdata_clean(SDAT)

    type(shr_strdata_type),intent(inout) :: SDAT

    integer :: n
    character(len=*),parameter :: subname = "(shr_strdata_clean) "
    !-------------------------------------------------------------------------------

    if (SDAT%nxg * SDAT%nyg == 0) then
       return
    endif

    ! Free MCT and PIO data first, while we still know which objects were
    ! allocated for which streams.
    call mct_ggrid_clean(SDAT%grid)
    call mct_gsmap_clean(SDAT%gsmap)

    do n = 1, SDAT%nstreams
       call pio_freedecomp(SDAT%pio_subsystem, SDAT%pio_iodesc(n))
       call mct_avect_clean(SDAT%avs(n))
       call mct_avect_clean(SDAT%avRLB(n))
       call mct_avect_clean(SDAT%avRUB(n))
       call mct_avect_clean(SDAT%avFLB(n))
       call mct_avect_clean(SDAT%avFUB(n))
       call mct_ggrid_clean(SDAT%gridR(n))
       if (SDAT%dofill(n)) call mct_sMatP_clean(SDAT%sMatPf(n))
       if (SDAT%domaps(n)) call mct_sMatP_clean(SDAT%sMatPs(n))
       call mct_gsmap_clean(SDAT%gsmapR(n))
    end do

    ! Now that all sub-objects are freed, clear components of the strdata
    ! object itself.
    SDAT%nxg = 0
    SDAT%nyg = 0
    SDAT%nzg = 0
    SDAT%strnxg = 0
    SDAT%strnyg = 0
    SDAT%strnzg = 0

    SDAT%nstreams = 0
    SDAT%nvectors = 0
    SDAT%ustrm = 0
    SDAT%vstrm = 0

    SDAT%dofill = .false.
    SDAT%domaps = .false.

  end subroutine shr_strdata_clean

  !===============================================================================

  subroutine shr_strdata_restWrite(filename,SDAT,mpicom,str1,str2)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: SDAT
    integer           ,intent(in)    :: mpicom
    character(len=*)      ,intent(in)    :: str1
    character(len=*)      ,intent(in)    :: str2

    !--- local ----
    integer :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restWrite) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)

    if (my_task == 0) then
       call shr_stream_restWrite(SDAT%stream,trim(filename),trim(str1),trim(str2),SDAT%nstreams)
    endif

  end subroutine shr_strdata_restWrite

  !===============================================================================

  subroutine shr_strdata_restRead(filename,SDAT,mpicom)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: SDAT
    integer           ,intent(in)    :: mpicom

    !--- local ----
    integer :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restRead) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)

    if (my_task == 0) then
       call shr_stream_restRead(SDAT%stream,trim(filename),SDAT%nstreams)
    endif

  end subroutine shr_strdata_restRead

  !===============================================================================

  subroutine shr_strdata_setOrbs(SDAT,eccen,mvelpp,lambm0,obliqr,modeldt)

    type(shr_strdata_type),intent(inout) :: SDAT
    real(R8),intent(in) :: eccen
    real(R8),intent(in) :: mvelpp
    real(R8),intent(in) :: lambm0
    real(R8),intent(in) :: obliqr
    integer,intent(in) :: modeldt

    ! local variables
    character(len=*),parameter :: subname = "(shr_strdata_setOrbs) "
    !-------------------------------------------------------------------------------

    SDAT%eccen   = eccen
    SDAT%mvelpp  = mvelpp
    SDAT%lambm0  = lambm0
    SDAT%obliqr  = obliqr
    SDAT%modeldt = modeldt

  end subroutine shr_strdata_setOrbs

  !===============================================================================
  subroutine shr_strdata_readnml(SDAT, file, rc, mpicom)

    ! Reads shr_strdata_nml namelist input

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: SDAT   ! strdata data data-type
    character(*) ,optional ,intent(in)   :: file   ! file to read strdata from
    integer      ,optional ,intent(out)  :: rc     ! return code
    integer      ,optional ,intent(in)   :: mpicom ! mpi comm

    ! local variables
    integer    :: rCode         ! return code
    integer    :: nUnit         ! fortran i/o unit number
    integer    :: n             ! generic loop index
    integer    :: my_task       ! my task number, 0 is default
    integer    :: master_task   ! master task number, 0 is default
    integer    :: ntasks        ! total number of tasks

    ! shr_strdata_nml namelist variables
    character(CL) :: dataMode           ! flags physics options wrt input data
    character(CL) :: domainFile         ! file   containing domain info
    integer       :: nx_global          ! global size of nx
    integer       :: ny_global          ! global size of ny
    character(CL) :: streams(nStrMax)   ! stream description file names
    character(CL) :: taxMode(nStrMax)   ! time axis cycling mode
    real(R8)      :: dtlimit(nStrMax)   ! delta time limiter
    character(CL) :: vectors(nVecMax)   ! define vectors to vector map
    character(CL) :: fillalgo(nStrMax)  ! fill algorithm
    character(CL) :: fillmask(nStrMax)  ! fill mask
    character(CL) :: fillread(nStrMax)  ! fill mapping file to read
    character(CL) :: fillwrite(nStrMax) ! fill mapping file to write
    character(CL) :: mapalgo(nStrMax)   ! scalar map algorithm
    character(CL) :: mapmask(nStrMax)   ! scalar map mask
    character(CL) :: mapread(nStrMax)   ! regrid mapping file to read
    character(CL) :: mapwrite(nStrMax)  ! regrid mapping file to write
    character(CL) :: tintalgo(nStrMax)  ! time interpolation algorithm
    character(CL) :: readmode(nStrMax)  ! file read mode
    character(CL) :: fileName           ! generic file name
    integer       :: yearFirst          ! first year to use in data stream
    integer       :: yearLast           ! last  year to use in data stream
    integer       :: yearAlign          ! data year that aligns with yearFirst

    !----- define namelist -----
    namelist / shr_strdata_nml / &
           dataMode        &
         , domainFile      &
         , nx_global       &
         , ny_global       &
         , streams         &
         , taxMode         &
         , dtlimit         &
         , vectors         &
         , fillalgo        &
         , fillmask        &
         , fillread        &
         , fillwrite       &
         , mapalgo         &
         , mapmask         &
         , mapread         &
         , mapwrite        &
         , tintalgo        &
         , readmode

    !----- formats -----
    character(*),parameter :: subName = "(shr_strdata_readnml) "
    character(*),parameter ::   F00 = "('(shr_strdata_readnml) ',8a)"
    character(*),parameter ::   F01 = "('(shr_strdata_readnml) ',a,i6,a)"
    character(*),parameter ::   F02 = "('(shr_strdata_readnml) ',a,es13.6)"
    character(*),parameter ::   F03 = "('(shr_strdata_readnml) ',a,l6)"
    character(*),parameter ::   F04 = "('(shr_strdata_readnml) ',a,i2,a,a)"
    character(*),parameter ::   F20 = "('(shr_strdata_readnml) ',a,i6,a)"
    character(*),parameter ::   F90 = "('(shr_strdata_readnml) ',58('-'))"
    !-------------------------------------------------------------------------------

    if (present(rc)) rc = 0

    my_task = 0
    master_task = 0
    ntasks = 1
    if (present(mpicom)) then
       call mpi_comm_rank(mpicom, my_task, rCode)
       call mpi_comm_size(mpicom, ntasks, rCode)
    endif

    sdat%my_task = my_task
    sdat%ntasks = ntasks
    sdat%master_task = master_task

    !--master--task--
    if (my_task == master_task) then

       !----------------------------------------------------------------------------
       ! set default values for namelist vars
       !----------------------------------------------------------------------------
       dataMode    = 'NULL'
       domainFile  = trim(shr_strdata_nullstr)
       streams(:)  = trim(shr_strdata_nullstr)
       taxMode(:)  = trim(shr_stream_taxis_cycle)
       dtlimit(:)  = dtlimit_default
       vectors(:)  = trim(shr_strdata_nullstr)
       fillalgo(:) = 'nn'
       fillmask(:) = 'nomask'
       fillread(:) = trim(shr_strdata_unset)
       fillwrite(:)= trim(shr_strdata_unset)
       mapalgo(:)  = 'bilinear'
       mapmask(:)  = 'dstmask'
       mapread(:)  = trim(shr_strdata_unset)
       mapwrite(:) = trim(shr_strdata_unset)
       tintalgo(:) = 'linear'
       readmode(:) = 'single'

       !----------------------------------------------------------------------------
       ! read input namelist
       !----------------------------------------------------------------------------
       if (present(file)) then
          write(logunit,F00) 'reading input namelist file: ',trim(file)
          open (newunit=nUnit,file=trim(file),status="old",action="read")
          call shr_nl_find_group_name(nUnit, 'shr_strdata_nml', status=rCode)
          if (rCode == 0) then
             read (nUnit, nml=shr_strdata_nml, iostat=rCode)
             if (rCode /= 0) then
                write(logunit,F01) 'ERROR: reading input namelist shr_strdata_input from file, &
                     &'//trim(file)//' iostat=',rCode
                call shr_sys_abort(subName//": namelist read error "//trim(file))
             end if
          end if
          close(nUnit)
       endif

       !----------------------------------------------------------------------------
       ! copy temporary/local namelist vars into data structure
       !----------------------------------------------------------------------------
       SDAT%nstreams    = 0
       do n=1,nStrMax
          call shr_stream_default(SDAT%stream(n))
       enddo
       SDAT%dataMode    = dataMode
       SDAT%domainFile  = domainFile
       SDAT%nxg         = nx_global
       SDAT%nyg         = ny_global
       SDAT%streams(:)  = streams(:)
       SDAT%taxMode(:)  = taxMode(:)
       SDAT%dtlimit(:)  = dtlimit(:)
       SDAT%vectors(:)  = vectors(:)
       SDAT%fillalgo(:) = fillalgo(:)
       SDAT%fillmask(:) = fillmask(:)
       SDAT%fillread(:) = fillread(:)
       SDAT%fillwrit(:) = fillwrite(:)
       SDAT%mapalgo(:)  = mapalgo(:)
       SDAT%mapmask(:)  = mapmask(:)
       SDAT%mapread(:)  = mapread(:)
       SDAT%mapwrit(:)  = mapwrite(:)
       SDAT%tintalgo(:) = tintalgo(:)
       SDAT%readmode(:) = readmode(:)
       do n=1,nStrMax
          if (trim(streams(n)) /= trim(shr_strdata_nullstr)) then
             SDAT%nstreams = max(SDAT%nstreams,n)
          end if
          if (trim(SDAT%taxMode(n)) == trim(shr_stream_taxis_extend)) then
             SDAT%dtlimit(n) = 1.0e30
          end if
       end do
       SDAT%nvectors = 0
       do n=1,nVecMax
          if (trim(vectors(n)) /= trim(shr_strdata_nullstr)) then
             SDAT%nvectors = n
          end if
       end do

       do n = 1,SDAT%nstreams
          if (trim(SDAT%streams(n)) /= shr_strdata_nullstr) then
             ! extract fileName (stream description text file), yearAlign, yearFirst, yearLast from SDAT%streams(n)
             call shr_stream_parseInput(SDAT%streams(n), fileName, yearAlign, yearFirst, yearLast)

             ! initialize stream datatype, read description text file
             call shr_stream_init(SDAT%stream(n), fileName, yearFirst, yearLast, yearAlign, trim(SDAT%taxMode(n)))
          end if
       enddo
    endif   ! master_task

    if (present(mpicom)) then
       call shr_mpi_bcast(SDAT%dataMode  ,mpicom ,'dataMode')
       call shr_mpi_bcast(SDAT%domainFile,mpicom ,'domainFile')
       call shr_mpi_bcast(SDAT%nxg       ,mpicom ,'nxg')
       call shr_mpi_bcast(SDAT%nyg       ,mpicom ,'nyg')
       call shr_mpi_bcast(SDAT%calendar  ,mpicom ,'calendar')
       call shr_mpi_bcast(SDAT%nstreams  ,mpicom ,'nstreams')
       call shr_mpi_bcast(SDAT%nvectors  ,mpicom ,'nvectors')
       call shr_mpi_bcast(SDAT%streams   ,mpicom ,'streams')
       call shr_mpi_bcast(SDAT%taxMode   ,mpicom ,'taxMode')
       call shr_mpi_bcast(SDAT%dtlimit   ,mpicom ,'dtlimit')
       call shr_mpi_bcast(SDAT%vectors   ,mpicom ,'vectors')
       call shr_mpi_bcast(SDAT%fillalgo  ,mpicom ,'fillalgo')
       call shr_mpi_bcast(SDAT%fillmask  ,mpicom ,'fillmask')
       call shr_mpi_bcast(SDAT%fillread  ,mpicom ,'fillread')
       call shr_mpi_bcast(SDAT%fillwrit  ,mpicom ,'fillwrit')
       call shr_mpi_bcast(SDAT%mapalgo   ,mpicom ,'mapalgo')
       call shr_mpi_bcast(SDAT%mapmask   ,mpicom ,'mapmask')
       call shr_mpi_bcast(SDAT%mapread   ,mpicom ,'mapread')
       call shr_mpi_bcast(SDAT%mapwrit   ,mpicom ,'mapwrit')
       call shr_mpi_bcast(SDAT%tintalgo  ,mpicom ,'tintalgo')
       call shr_mpi_bcast(SDAT%readmode  ,mpicom ,'readmode')
    endif

    SDAT%ymdLB = -1
    SDAT%todLB = -1
    SDAT%ymdUB = -1
    SDAT%todUB = -1
    SDAT%dtmin = 1.0e30
    SDAT%dtmax = 0.0
    SDAT%nzg   = 0
    SDAT%eccen  = SHR_ORB_UNDEF_REAL
    SDAT%mvelpp = SHR_ORB_UNDEF_REAL
    SDAT%lambm0 = SHR_ORB_UNDEF_REAL
    SDAT%obliqr = SHR_ORB_UNDEF_REAL
    SDAT%modeldt = 0
    SDAT%calendar = shr_cal_noleap

  end subroutine shr_strdata_readnml

  !===============================================================================

  subroutine shr_strdata_pioinit(SDAT, compid )

    ! !DESCRIPTION: Initialize PIO for a component model

    ! input/output arguments
    type(shr_strdata_type),intent(inout) :: SDAT  ! strdata data data-type
    type(iosystem_desc_t) , pointer      :: io_subsystem
    integer               , intent(in)   :: compid
    !-------------------------------------------------------------------------------

    SDAT%pio_subsystem => shr_pio_getiosys(compid)
    SDAT%io_type       =  shr_pio_getiotype(compid)
    SDAT%io_format     =  shr_pio_getioformat(compid)

  end subroutine shr_strdata_pioinit

  !===============================================================================

  subroutine shr_strdata_print(SDAT,name)

    ! !DESCRIPTION:
    !  Print strdata common to all data models

    ! !INPUT/OUTPUT PARAMETERS:
    type(shr_strdata_type)  ,intent(in) :: SDAT  ! strdata data data-type
    character(len=*),optional,intent(in) :: name  ! just a name for tracking

    integer   :: n
    character(CL) :: lname

    !----- formats -----
    character(*),parameter :: subName = "(shr_strdata_print) "
    character(*),parameter ::   F00 = "('(shr_strdata_print) ',8a)"
    character(*),parameter ::   F01 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter ::   F02 = "('(shr_strdata_print) ',a,es13.6)"
    character(*),parameter ::   F03 = "('(shr_strdata_print) ',a,l6)"
    character(*),parameter ::   F04 = "('(shr_strdata_print) ',a,i2,a,a)"
    character(*),parameter ::   F05 = "('(shr_strdata_print) ',a,i2,a,i6)"
    character(*),parameter ::   F06 = "('(shr_strdata_print) ',a,i2,a,l2)"
    character(*),parameter ::   F07 = "('(shr_strdata_print) ',a,i2,a,es13.6)"
    character(*),parameter ::   F20 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter ::   F90 = "('(shr_strdata_print) ',58('-'))"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lname = 'unknown'
    if (present(name)) then
       lname = trim(name)
    endif
    !----------------------------------------------------------------------------
    ! document datatype settings
    !----------------------------------------------------------------------------
    write(logunit,F90)
    write(logunit,F00) "name        = ",trim(lname)
    write(logunit,F00) "dataMode    = ",trim(SDAT%dataMode)
    write(logunit,F00) "domainFile  = ",trim(SDAT%domainFile)
    write(logunit,F01) "nxg         = ",SDAT%nxg
    write(logunit,F01) "nyg         = ",SDAT%nyg
    write(logunit,F01) "nzg         = ",SDAT%nzg
    write(logunit,F00) "calendar    = ",trim(SDAT%calendar)
    write(logunit,F01) "io_type     = ",SDAT%io_type
    write(logunit,F02) "eccen       = ",SDAT%eccen
    write(logunit,F02) "mvelpp      = ",SDAT%mvelpp
    write(logunit,F02) "lambm0      = ",SDAT%lambm0
    write(logunit,F02) "obliqr      = ",SDAT%obliqr
    write(logunit,F01) "nstreams    = ",SDAT%nstreams
    write(logunit,F01) "pio_iotype  = ",sdat%io_type

    do n=1, SDAT%nstreams
       write(logunit,F04) "  streams (",n,") = ",trim(SDAT%streams(n))
       write(logunit,F04) "  taxMode (",n,") = ",trim(SDAT%taxMode(n))
       write(logunit,F07) "  dtlimit (",n,") = ",SDAT%dtlimit(n)
       write(logunit,F05) "  strnxg  (",n,") = ",SDAT%strnxg(n)
       write(logunit,F05) "  strnyg  (",n,") = ",SDAT%strnyg(n)
       write(logunit,F05) "  strnzg  (",n,") = ",SDAT%strnzg(n)
       write(logunit,F06) "  dofill  (",n,") = ",SDAT%dofill(n)
       write(logunit,F04) "  fillalgo(",n,") = ",trim(SDAT%fillalgo(n))
       write(logunit,F04) "  fillmask(",n,") = ",trim(SDAT%fillmask(n))
       write(logunit,F04) "  fillread(",n,") = ",trim(SDAT%fillread(n))
       write(logunit,F04) "  fillwrit(",n,") = ",trim(SDAT%fillwrit(n))
       write(logunit,F06) "  domaps  (",n,") = ",SDAT%domaps(n)
       write(logunit,F04) "  mapalgo (",n,") = ",trim(SDAT%mapalgo(n))
       write(logunit,F04) "  mapmask (",n,") = ",trim(SDAT%mapmask(n))
       write(logunit,F04) "  mapread (",n,") = ",trim(SDAT%mapread(n))
       write(logunit,F04) "  mapwrit (",n,") = ",trim(SDAT%mapwrit(n))
       write(logunit,F04) "  tintalgo(",n,") = ",trim(SDAT%tintalgo(n))
       write(logunit,F04) "  readmode(",n,") = ",trim(SDAT%readmode(n))
       write(logunit,F01) " "
    end do
    write(logunit,F01) "nvectors    = ",SDAT%nvectors
    do n=1, SDAT%nvectors
       write(logunit,F04) "  vectors (",n,") = ",trim(SDAT%vectors(n))
    end do
    write(logunit,F90)
    call shr_sys_flush(logunit)

  end subroutine shr_strdata_print

end module dshr_strdata_mod
