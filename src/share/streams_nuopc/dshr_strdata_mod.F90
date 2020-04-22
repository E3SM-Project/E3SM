module dshr_strdata_mod

  ! holds data and methods to advance data models
  ! Obtain the model domain and the stream domain for each stream
  ! For the model domain 
  !  - if single column - read in the data model domain file - and find the nearest neighbor
  !  - if not single column - will obtain it directly from the mesh input but will still need 
  !    to read in the domain file to obtain the global nx and ny that need to be passed as scalar
  !    data back to the mediator

  use ESMF
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod      , only : shr_sys_abort
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_file_mod     , only : shr_file_getunit, shr_file_freeunit
  use shr_const_mod    , only : shr_const_pi, shr_const_cDay
  use shr_log_mod      , only : logunit => shr_log_Unit
  use shr_cal_mod      , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod      , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod      , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_nl_mod       , only : shr_nl_find_group_name
  use shr_pio_mod      , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
  use shr_string_mod   , only : shr_string_listgetname, shr_string_listisvalid, shr_string_listgetnum
  use dshr_stream_mod  , only : shr_stream_streamtype, shr_stream_getModelFieldList
  use dshr_stream_mod  , only : shr_stream_taxis_cycle, shr_stream_taxis_extend, shr_stream_default
  use dshr_stream_mod  , only : shr_stream_getNFiles, shr_stream_getFile
  use dshr_stream_mod  , only : shr_stream_getCurrFile, shr_stream_setCurrFile
  use dshr_stream_mod  , only : shr_stream_getDomainFile, shr_stream_getDomainInfo
  use dshr_stream_mod  , only : shr_stream_init_from_infiles, shr_stream_init_from_fortran
  use dshr_stream_mod  , only : shr_stream_restWrite, shr_stream_restRead
  use dshr_stream_mod  , only : shr_stream_getFilePath, shr_stream_findBounds
  use dshr_stream_mod  , only : shr_stream_getnextfilename, shr_stream_getprevfilename 
  use dshr_stream_mod  , only : shr_stream_getfilefieldname 
  use dshr_tinterp_mod , only : shr_tInterp_getCosz, shr_tInterp_getAvgCosz, shr_tInterp_getFactors
  use dshr_methods_mod , only : dshr_fldbun_getfldptr, chkerr
  use pio              , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
  use pio              , only : pio_openfile, pio_closefile, pio_nowrite
  use pio              , only : pio_seterrorhandling, pio_initdecomp, pio_freedecomp
  use pio              , only : pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
  use pio              , only : pio_inq_dimlen, pio_double, pio_offset_kind
  use pio              , only : pio_read_darray, pio_get_var, pio_setframe
  use pio              , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
  use shr_map_mod      , only : shr_map_maptype, shr_map_mapset, shr_map_get, shr_map_clean
  use shr_map_mod      , only : shr_map_fs_ndst, shr_map_fs_nsrc, shr_map_fs_nwts
  use shr_mct_mod
  use mct_mod
  use perf_mod

  implicit none
  private

  public  :: shr_strdata_type
  public  :: shr_strdata_init_from_infiles
  public  :: shr_strdata_init_from_fortran
  public  :: shr_strdata_restRead
  public  :: shr_strdata_restWrite
  public  :: shr_strdata_setOrbs
  public  :: shr_strdata_print
  public  :: shr_strdata_advance
  public  :: shr_strdata_get_stream_domain
  public  :: shr_strdata_get_griddata
  public  :: shr_strdata_set_griddata
  public  :: shr_strdata_mapSet

  private :: shr_strdata_readnml
  private :: shr_strdata_init_model_domain
  private :: shr_strdata_init_streams_from_infiles
  private :: shr_strdata_init_streams_from_fortran
  private :: shr_strdata_init_mapping
  private :: shr_strdata_init_stream_domain
  private :: shr_strdata_gGridCompare
  private :: shr_strdata_readLBUB

  ! public data members:
  integer                              :: debug    = 0  ! local debug flag
  integer          ,parameter          :: nStrMax = 30
  integer          ,parameter          :: nVecMax = 30
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  real(R8)         ,parameter, private :: dtlimit_default = 1.5_R8
  integer          ,parameter          :: master_task = 0

  type shr_strdata_type
     ! stream info set by shr_strdata_nml
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
     integer                        :: nstreams          ! number of streams set in shr_strdata_readnml
     integer                        :: nvectors          ! number of vectors set in shr_strdata_readnml

     ! pio info
     integer                        :: io_type
     integer                        :: io_format
     type(iosystem_desc_t), pointer :: pio_subsystem => null()
     type(io_desc_t)                :: pio_iodesc(nStrMax)

     ! data required by stream  cosz t-interp method, set by user
     real(R8)                       :: eccen
     real(R8)                       :: mvelpp
     real(R8)                       :: lambm0
     real(R8)                       :: obliqr
     integer                        :: modeldt           ! model dt in seconds

     ! model info
     integer                        :: nxg                   ! model domain lon size
     integer                        :: nyg                   ! model domain lat size
     integer                        :: nzg                   ! model domain vertical size
     integer                        :: lsize                 ! model domain local size
     integer                        :: gsize                 ! model domain global
     integer, pointer               :: gindex(:)             ! model domain global index spzce
     type(mct_gsmap)                :: gsmap                 ! model domain mct global seg map
     type(mct_ggrid)                :: grid                  ! model domain mct ggrid
     type(ESMF_Mesh)                :: mesh_model            ! model mesh
     type(ESMF_FieldBundle)         :: fldbun_model(nStrMax) ! stream data time interpolated and mapped to model domain
     type(mct_avect)                :: avFUB(nStrMax)        ! stream data LB mapped to model domain
     type(mct_avect)                :: avFLB(nStrMax)        ! stream data UB mapped to model domain
     type(mct_avect)                :: avCoszen(nStrMax)     ! data assocaited with coszen time interp

     ! stream info
     type(shr_stream_streamType)    :: stream(nStrMax)
     integer                        :: strnxg(nStrMax)
     integer                        :: strnyg(nStrMax)
     integer                        :: strnzg(nStrMax)
     logical                        :: dofill(nStrMax)
     logical                        :: domaps(nStrMax)
     integer                        :: lsizeR(nStrMax)
     type(mct_gsmap)                :: gsmapR(nStrMax)
     type(mct_rearr)                :: rearrR(nStrMax)
     type(mct_ggrid)                :: gridR(nStrMax)    ! stream grid
     type(mct_avect)                :: avRFile(nStrMax)  ! Read attrvect for multiple time slices - stream grid
     type(mct_avect)                :: avRLB(nStrMax)    ! Read attrvect - stream grid
     type(mct_avect)                :: avRUB(nStrMax)    ! Read attrvect - stream grid
     type(mct_sMatP)                :: sMatPf(nStrMax)
     type(mct_sMatP)                :: sMatPs(nStrMax)
     integer                        :: ymdLB(nStrMax)
     integer                        :: todLB(nStrMax)
     integer                        :: ymdUB(nStrMax)
     integer                        :: todUB(nStrMax)
     integer                        :: ustrm(nVecMax)
     integer                        :: vstrm(nVecMax)

     ! time invo
     real(R8)                       :: dtmin(nStrMax)
     real(R8)                       :: dtmax(nStrMax)
     integer                        :: ymd  ,tod
     character(CL)                  :: calendar          ! model calendar for ymd,tod
     character(CL)                  :: allocstring
  end type shr_strdata_type

  integer          ,parameter :: CompareXYabs      = 1   ! X,Y  relative error
  integer          ,parameter :: CompareXYrel      = 2   ! X,Y  absolute error
  integer          ,parameter :: CompareAreaAbs    = 3   ! area relative error
  integer          ,parameter :: CompareAreaRel    = 4   ! area absolute error
  integer          ,parameter :: CompareMaskIdent  = 5   ! masks are identical
  integer          ,parameter :: CompareMaskZeros  = 6   ! masks have same zeros
  integer          ,parameter :: CompareMaskSubset = 7   ! mask is subset of other
  integer          ,parameter :: CompareXYabsMask  = 101 ! X,Y  relative error
  integer          ,parameter :: iotype_std_netcdf = -99 ! non pio option
  real(R8)         ,parameter :: deg2rad = SHR_CONST_PI/180.0_R8
  character(len=*) ,parameter :: allocstring_value = 'strdata_allocated'
  character(*)     ,parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_strdata_init_from_infiles(sdat, nlfilename, mesh, clock, &
        mpicom, compid, logunit, reset_domain_mask, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    character(len=*)       , intent(in)    :: nlfilename ! for shr_strdata_nml namelist
    type(ESMF_Mesh)        , intent(inout) :: mesh 
    type(ESMF_Clock)       , intent(in)    :: clock
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: compid
    integer                , intent(in)    :: logunit
    logical, optional      , intent(in)    :: reset_domain_mask
    integer                , intent(out)   :: rc

    ! local varaibles
    type(ESMF_Calendar)          :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag)      :: esmf_caltype  ! esmf calendar type
    character(CS)                :: calendar      ! calendar name
    integer                      :: my_task
    integer                      :: ierr
    integer                      :: master_task = 0
    character(len=*), parameter  :: subname='(dshr_mod:dshr_sdat_init)'
    character(*)    , parameter  :: F01="('(dshr_init_strdata) ',a,2f10.4)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call mpi_comm_rank(mpicom, my_task, ierr)

    ! Read shr_strdata_nml from nlfilename
    ! Read sdat namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdat, trim(nlfilename), mpicom=mpicom)

    ! Initialize sdat  pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! Set sdat mesh_model
    sdat%mesh_model = mesh

    ! Initialize the sdat model domain info
    call shr_strdata_init_model_domain(mesh,  mpicom, compid, sdat, reset_domain_mask=reset_domain_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize sdat stream domains
    call shr_strdata_init_streams_from_infiles(sdat, compid, mpicom, my_task)

    ! initialize sdat attributes mapping of streams to model domain
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

    ! initialize sdat calendar
    call ESMF_ClockGet(clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(calendar))
    end if
    sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))

    ! print sdat output
    if (my_task == master_task) then
       call shr_strdata_print(sdat,'sdat data ')
       write(logunit,*) ' successfully initialized sdat'
    endif

  end subroutine shr_strdata_init_from_infiles

  !===============================================================================
  subroutine shr_strdata_init_from_fortran(                    &
       sdat, mpicom, compid, mesh, nxg, nyg, clock,            &
       !--- streams stuff required ---
       yearFirst, yearLast, yearAlign, offset,                 &
       DomFilePath, DomFileName,                               &
       DomTvarName, DomXvarName, DomYvarName, domMaskName,     &
       strmFldFilePath, strmFldFilenames,                      &
       strmFldNamesInFile, strmfldNamesInModel,                &
       !--- strdata optional ---
       nzg, domZvarName,                                       &
       taxMode, dtlimit, tintalgo, readmode,                   &
       fillalgo, fillmask, mapalgo, mapmask) 

    ! Set strdata and stream info from fortran interface.
    ! Note: When this is called, previous settings are reset to defaults
    ! and then the values passed are used.

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: sdat                ! strdata data data-type
    integer                ,intent(in)   :: mpicom              ! mpi comm
    integer                ,intent(in)   :: compid
    type(ESMF_Mesh)        ,intent(in)   :: mesh
    integer                ,intent(in)   :: nxg
    integer                ,intent(in)   :: nyg
    type(ESMF_Clock)       ,intent(in)   :: clock
    integer                ,intent(in)   :: yearFirst           ! first year to use
    integer                ,intent(in)   :: yearLast            ! last  year to use
    integer                ,intent(in)   :: yearAlign           ! align yearFirst with this model year
    integer                ,intent(in)   :: offset              ! offset in seconds of stream data
    character(*)           ,intent(in)   :: domFilePath         ! domain file path
    character(*)           ,intent(in)   :: domFileName         ! domain file name
    character(*)           ,intent(in)   :: domTvarName         ! domain time dim name
    character(*)           ,intent(in)   :: domXvarName         ! domain x dim name
    character(*)           ,intent(in)   :: domYvarName         ! domain y dim name
    character(*)           ,intent(in)   :: domMaskName         ! domain mask name
    character(*)           ,intent(in)   :: strmFldFilePath     ! path to stream data 
    character(*)           ,intent(in)   :: strmFldFileNames(:) ! filenames for stream data
    character(*)           ,intent(in)   :: strmFldNamesInFile  ! file field names, colon delim list
    character(*)           ,intent(in)   :: strmFldNamesInModel ! model field names, colon delim list
    integer                ,intent(in)   :: nzg
    character(*)           ,intent(in)   :: domZvarName         ! domain z dim name
    character(*)           ,intent(in)   :: taxMode
    real(R8)               ,intent(in)   :: dtlimit
    character(*)           ,intent(in)   :: fillalgo            ! fill algorithm
    character(*)           ,intent(in)   :: fillmask            ! fill mask
    character(*)           ,intent(in)   :: mapalgo             ! scalar map algorithm
    character(*)           ,intent(in)   :: mapmask             ! scalar map mask
    character(*)           ,intent(in)   :: tintalgo            ! time interpolation algorithm
    character(*)           ,intent(in)   :: readmode            ! file read mode

    ! local variables
    type(ESMF_Calendar)     :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype  ! esmf calendar type
    character(CS)           :: calendar      ! calendar name
    character(CS)           :: zname
    integer                 :: my_task
    integer                 :: ierr
    integer                 :: rc
    character(*),parameter  :: subName = "(shr_strdata_create) "
    character(*),parameter  :: F00 = "('(shr_strdata_create) ',8a)"
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom, my_task, ierr)

    ! Assume only 1 stream
    sdat%nstreams = 1

    ! Initialize pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! set defaults - but calling shr_strdata_readnml - but not reading int the namelist
    call shr_strdata_readnml(sdat) 

    ! Reset sdat values if they appear as optional arguments
    ! TODO: pass these values to shr_strdata routine rather than reading them here
    sdat%taxMode(1)  = taxMode
    sdat%dtlimit(1)  = dtlimit
    sdat%fillalgo(1) = fillalgo
    sdat%fillmask(1) = fillmask
    sdat%mapalgo(1)  = mapalgo
    sdat%mapmask(1)  = mapmask
    sdat%tintalgo(1) = tintalgo
    sdat%readmode(1) = readmode ! single or full_file
    if (trim(sdat%taxMode(1)) == trim(shr_stream_taxis_extend)) then
       ! reset dtlimit if necessary
       sdat%dtlimit(1) = 1.0e30
    end if

    ! Initialize sdat mesh
    sdat%mesh_model = mesh

    ! Initialize sdat model domain info
    call shr_strdata_init_model_domain(mesh, mpicom, compid, sdat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize sdat stream domain info
    call shr_strdata_init_streams_from_fortran(sdat, &
       yearFirst, yearLast, yearAlign, offset, sdat%taxmode(1), &
       domFilePath, domFileName, &
       domTvarName, domXvarName, domYvarName, domMaskName, &
       domZvarName, nzg, &
       strmFldNamesInFile, strmFldNamesInModel, &
       strmFldFilePath, strmFldFileNames)

    ! Initialize sdat mapping of stream to model domain
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

    ! Initialize sdat calendar
    call ESMF_ClockGet(clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(calendar))
    end if
    sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))

  end subroutine shr_strdata_init_from_fortran

  !===============================================================================
  subroutine shr_strdata_init_model_domain( mesh, mpicom, compid, sdat, reset_domain_mask, rc)

    ! ----------------------------------------------
    ! Initialize sdat%lsize, sdat%gsmap and sdat%grid
    ! sdat%nxg, sdat%nyg and sdat%nzg are initialized in shr_strdata_readnl
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: compid
    type(shr_strdata_type) , intent(inout) :: sdat
    logical, optional      , intent(in)    :: reset_domain_mask
    integer                , intent(out)   :: rc

    ! local variables
    integer              :: n,k          ! generic counters
    integer              :: lsize        ! local size
    integer              :: gsize
    type(ESMF_DistGrid)  :: distGrid
    integer              :: dimCount
    integer              :: tileCount
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
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize sdat%lsize, sdat%gsize and sdat%gindex
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    sdat%lsize = lsize
    allocate(sdat%gindex(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=sdat%gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    sdat%gsize = gsize
    deallocate(elementCountPTile)

    ! initialize sdat%gsmap
    call mct_gsMap_init(sdat%gsmap, sdat%gindex, mpicom, compid, sdat%lsize, sdat%gsize)

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

    if (present(reset_domain_mask)) then
       if (reset_domain_mask) then
          if (my_task == master_task) then
             write(logunit,*) ' Resetting the component domain mask to 1'
          end if
          sdat%grid%data%rattr(kmask,:) = 1
       end if
    end if

  end subroutine shr_strdata_init_model_domain

  !===============================================================================
  subroutine shr_strdata_init_streams_from_infiles(sdat, compid, mpicom, my_task)

    ! Initialize streams for input stream txt files

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout) :: sdat
    integer                ,intent(in)    :: compid
    integer                ,intent(in)    :: mpicom
    integer                ,intent(in)    :: my_task

    ! local variables
    integer       :: n         ! generic index
    character(CL) :: str       ! namelist input for "streams" 
    integer       :: i         ! index in string
    integer       :: yearFirst ! first year to use in data stream
    integer       :: yearLast  ! last  year to use in data stream
    integer       :: yearAlign ! data year that aligns with yearFirst
    character(CL) :: filename  ! stream txt filename
    character(len=*), parameter :: subname = "(shr_strdata_init_streams) "
    !-------------------------------------------------------------------------------

    ! Loop over all active streams and initialize domain and calendar for each stream
    do n = 1,sdat%nstreams
       if (trim(sdat%streams(n)) /= shr_strdata_nullstr) then

          if (my_task == master_task) then
             ! Set yearAlign, yearFirst, yearLast for sdat%streams(n)
             str = adjustL(sdat%streams(n))  ! stream info from namelist
             i = index(str," ")
             fileName = str(:i)
             read(str(i:),*) yearAlign, yearFirst, yearLast

             ! Initialize stream datatype by reading stream txt files
             call shr_stream_init_from_infiles(sdat%stream(n), &
                  fileName, yearFirst, yearLast, yearAlign, trim(sdat%taxMode(n)))
          end if

          ! Initialize stream domain info for stream n
          call shr_strdata_init_stream_domain( sdat%stream(n), &
               sdat%pio_subsystem, sdat%io_type, compid, mpicom, &
               sdat%gridR(n), sdat%gsmapR(n), sdat%pio_iodesc(n), &
               sdat%strnxg(n), sdat%strnyg(n), sdat%strnzg(n), sdat%lsizeR(n))

          ! Initialize calendar for stream n
          call shr_mpi_bcast(sdat%stream(n)%calendar, mpicom)
       end if
    enddo

  end subroutine shr_strdata_init_streams_from_infiles

  !===============================================================================
  subroutine shr_strdata_init_streams_from_fortran(sdat,   &
       yearFirst, yearLast, yearAlign, offset, taxmode,    &
       domFilePath, domFileName,                           &
       domTvarName, domXvarName, domYvarName, domMaskName, &
       domZvarName, nzg, &
       strmFldNamesInFile, strmFldNamesInModel, &
       strmFldfilePath, strmFldFileNames)

    ! -----------------------------------------------
    ! initialize stream info from fortran interfaces 
    ! rather than from stream text files
    ! -----------------------------------------------

    ! input/output variables
    type(shr_strdata_type) ,intent(inout):: sdat                ! strdata data data-type
    integer                ,intent(in)   :: yearFirst           ! first year to use
    integer                ,intent(in)   :: yearLast            ! last  year to use
    integer                ,intent(in)   :: yearAlign           ! align yearFirst with this model year
    integer                ,intent(in)   :: offset              ! offset in seconds of stream data
    character(*)           ,intent(in)   :: taxMode             ! time axis mode
    character(*)           ,intent(in)   :: domFilePath         ! domain file path
    character(*)           ,intent(in)   :: domFileName         ! domain file name
    character(*)           ,intent(in)   :: domTvarName         ! domain time dim name
    character(*)           ,intent(in)   :: domXvarName         ! domain x dim name
    character(*)           ,intent(in)   :: domYvarName         ! domain y dim name
    character(*)           ,intent(in)   :: domMaskName         ! domain mask name
    character(*)           ,intent(in)   :: domZvarName         ! domain z dim name
    integer                ,intent(in)   :: nzg                 ! domain z dim size
    character(*)           ,intent(in)   :: strmFldNamesInFile  ! file field names, colon delim list
    character(*)           ,intent(in)   :: strmFldNamesInModel ! model field names, colon delim list
    character(*)           ,intent(in)   :: strmFldfilePath     ! path to filenames
    character(*)           ,intent(in)   :: strmFldFileNames(:) ! filename for index filenumber
    ! --------------------------------------------------------

    call shr_stream_init_from_fortran(sdat%stream(1), &
       yearFirst, yearLast, yearAlign, offset, taxmode,    &
       domFilePath, domFileName, &
       domTvarName, domXvarName, domYvarName, domZvarName, nzg, domMaskName, &
       strmFldNamesInFile, strmFldNamesInModel, &
       strmFldfilePath, strmFldFileNames)

  end subroutine shr_strdata_init_streams_from_fortran

  !===============================================================================
  subroutine shr_strdata_init_stream_domain(stream, &
       pio_subsystem, pio_type, compid, mpicom, &
       gGrid, gsMap, pio_iodesc, &
       nxg, nyg, nzg, lsize)

    !----------------------------------------------------------------------------
    ! Create mct ggrid for model grid and set model gsmap if not input
    ! o assumes a very specific netCDF domain file format wrt var names, etc.
    !----------------------------------------------------------------------------

    ! input/output variables
    type(shr_stream_streamType) , intent(in)    :: stream
    type(iosystem_desc_t)       , pointer       :: pio_subsystem
    integer                     , intent(in)    :: pio_type
    integer                     , intent(in)    :: compid
    integer                     , intent(in)    :: mpicom
    type(mct_gGrid)             , intent(inout) :: gGrid
    type(mct_gsMap)             , intent(inout) :: gsMap
    type(io_desc_t)             , intent(inout) :: pio_iodesc
    integer                     , intent(out)   :: nxg
    integer                     , intent(out)   :: nyg
    integer                     , intent(out)   :: nzg
    integer                     , intent(out)   :: lsize

    ! local variables
    character(CL)        :: filePath   ! file path of stream domain file
    character(CXX)       :: filename   ! file name of stream domain file
    character(CS)        :: timeName   ! domain file: time variable name
    character(CS)        :: lonname    ! name of  lon variable in file
    character(CS)        :: latname    ! name of  lat variable in file
    character(CS)        :: hgtname    ! name of  hgt variable in file
    character(CS)        :: maskname   ! name of mask variable in file
    integer              :: n,k,j,i    ! indices
    integer              :: gsize      ! gsize
    integer              :: my_task
    integer              :: npes
    integer              :: ierr       ! error code
    logical              :: fileexists !
    logical              :: maskexists ! is mask on dataset
    integer              :: klon
    integer              :: klat
    integer              :: kmask
    integer              :: khgt
    real(R8),allocatable :: lon1d(:)   ! temp array for domain lon  info
    real(R8),allocatable :: lat1d(:)   ! temp array for domain lat  info
    real(R8),allocatable :: lon(:,:)   ! temp array for domain lon  info
    real(R8),allocatable :: lat(:,:)   ! temp array for domain lat  info
    integer ,allocatable :: mask(:,:)  ! temp array for domain mask info
    real(R8),allocatable :: hgt(:)     ! temp array for domain height info
    integer, allocatable :: dimid(:)
    type(var_desc_t)     :: varid
    type(file_desc_t)    :: pioid
    integer              :: rcode
    integer              :: ndims      ! number of dims
    integer, pointer     :: gindex_stream(:)
    integer              :: my_start, my_end
    integer              :: npt
    type(mct_ggrid)      :: gGridRoot    ! global mct ggrid
    character(*), parameter :: subname = '(shr_strdata_readgrid_stream) '
    character(*), parameter :: F00   = "('(shr_strdata_readgrid_stream) ',8a)"
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom, my_task, ierr)
    call mpi_comm_size(mpicom, npes, ierr)

    ! Determine stream filename and namelist of domain variables
    if (my_task == master_task) then

       ! Obtain filename and variables names for stream domain info
       call shr_stream_getDomainInfo(stream, filePath, fileName, &
            timeName, lonName, latName, hgtName, maskName)

       ! Determine if stream domain file exists and if not exit
       call shr_stream_getFile(filePath, fileName)
       inquire(file=trim(filename), exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if

       ! Write out stream domain info to logunit
       write(logunit,*) subname,' stream data'
       write(logunit,*) subname,' filePath = ',trim(filePath)
       write(logunit,*) subname,' fileName = ',trim(fileName)
       write(logunit,*) subname,' timeName = ',trim(timeName)
       write(logunit,*) subname,' lonName  = ',trim(lonName)
       write(logunit,*) subname,' latName  = ',trim(latName)
       write(logunit,*) subname,' hgtName  = ',trim(hgtName)
       write(logunit,*) subname,' maskName = ',trim(maskName)
    endif
    call shr_mpi_bcast(filePath ,mpicom)
    call shr_mpi_bcast(fileName ,mpicom)
    call shr_mpi_bcast(lonName  ,mpicom)
    call shr_mpi_bcast(latName  ,mpicom)
    call shr_mpi_bcast(hgtName  ,mpicom)
    call shr_mpi_bcast(maskName ,mpicom)

    ! Open stream domain file
    rcode = pio_openfile(pio_subsystem, pioid, pio_type, trim(filename), pio_nowrite)
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)

    ! Obtain stream domain lon/lat arrays (in degrees)
    call pio_seterrorhandling(pioid, PIO_RETURN_ERROR)
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, trim(lonname), varid)
    rcode = pio_inq_varndims(pioid, varid, ndims)
    allocate(dimid(ndims))
    if (ndims == 1) then
       ! lon and lat in stream domain file are 1d arrays
       rcode = pio_inq_varid(pioid, trim(lonname), varid)
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), nxg)
       rcode = pio_inq_varid(pioid, trim(latname), varid)
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), nyg)
       allocate(lon1d(nxg)); allocate(lat1d(nyg))
       rcode = pio_inq_varid(pioid, trim(lonname), varid)
       rcode = pio_get_var(pioid, varid, (/1/), (/nxg/), lon1d)
       rcode = pio_inq_varid(pioid, trim(latname), varid)
       rcode = pio_get_var(pioid, varid, (/1/), (/nyg/), lat1d)
       allocate(lon(nxg,nyg)); allocate(lat(nxg,nyg))
       do j = 1,nyg
          do i = 1,nxg
             lon(i,j) = lon1d(i)
             lat(i,j) = lat1d(j)
          end do
       end do
       deallocate(lon1d); deallocate(lat1d)
    else
       ! lon and lat in stream domain file are 2d arrays 
       ! obtain both nxg and nyg from lon info
       rcode = pio_inq_varid(pioid, trim(lonname), varid)
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), nxg)
       rcode = pio_inq_dimlen(pioid, dimid(2), nyg)
       allocate(lon(nxg,nyg)); allocate(lat(nxg,nyg))
       rcode = pio_get_var(pioid, varid, lon)
       rcode = pio_inq_varid(pioid, trim(latname), varid)
       rcode = pio_get_var(pioid, varid, lat)
    endif

    ! Create mask array (  for ocean, 0 for non-ocean)
    allocate(mask(nxg,nyg))

    call pio_seterrorhandling(pioid, PIO_RETURN_ERROR)
    rcode = pio_inq_varid(pioid, trim(maskname), varid)
    if (rcode == PIO_NOERR) then
       ! mask exists on file
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid, trim(maskname), varid)
       rcode = pio_get_var(pioid, varid, mask)
    else
       ! mask does not exist on file - set to 1
       mask(:,:) = 1
    end if

    ! Create hgt (vertical dimension array)
    call pio_seterrorhandling(pioid, PIO_RETURN_ERROR)
    rcode = pio_inq_varid(pioid, trim(hgtName), varid)
    if (rcode == PIO_NOERR) then
       ! vertical coordinate exists in file
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), nzg)
       allocate(hgt(nzg))
       rCode = pio_inq_varid(pioid, trim(hgtName), varid)
       rcode = pio_get_var(pioid, varid, hgt)
    else
       nzg = -1
       allocate(hgt(1))
       hgt(1) = 1
    endif

    ! Create 1d decomp of grid (gindex_stream)
    gsize = abs(nxg*nyg*nzg)
    lsize = gsize/npes
    my_start = lsize*my_task + min(my_task, mod(gsize, npes)) + 1
    if (my_task < mod(gsize, npes)) then
       lsize = lsize + 1
    end if
    my_end = my_start + lsize - 1
    allocate(gindex_stream(lsize))
    npt = 0
    do n = 1,gsize
       npt = npt + 1
       if (npt >= my_start .and. npt <= my_end) then
          gindex_stream(npt - my_start + 1) = n
       end if
    end do

    ! Create stream pio_iodesc using gindex_stream
    if (nzg <= 0) then
       call pio_initdecomp(pio_subsystem, pio_double, (/nxg,nyg/), gindex_stream, pio_iodesc)
    else
       call pio_initdecomp(pio_subsystem, pio_double, (/nxg,nyg,nzg/), gindex_stream, pio_iodesc)
    endif
    
    ! Create mct stream gsmap using gindex_stream
    if (my_task == master_task) write(logunit,*)trim(subname) // ': Creating gsmap for input stream'
    call mct_gsMap_init(gsmap, gindex_stream, mpicom, compid, lsize, gsize)

    ! Create mct stream ggrid
    if (my_task == master_task) write(logunit,*)trim(subname) // ': Creating ggrid for input stream'
    call mct_gGrid_init(gGrid=gGridRoot, CoordChars='lat:lon:hgt', OtherChars='mask', lsize=gsize)
    call mct_gGrid_init(gGrid=Ggrid    , CoordChars='lat:lon:hgt', OtherChars='mask', lsize=lsize)
    call mct_gGrid_importIAttr(gGrid, 'GlobGridNum', gindex_stream, lsize)
    gGridRoot%data%rAttr = -9999.0_R8 !to avoid errors when using strict compiler checks
    gGrid%data%rAttr     = -9999.0_R8
    klon  = mct_aVect_indexRA(gGrid%data, 'lon')
    klat  = mct_aVect_indexRA(gGrid%data, 'lat')
    kmask = mct_aVect_indexRA(gGrid%data, 'mask')
    khgt  = mct_aVect_indexRA(gGrid%data, 'hgt')
    n = 0
    do k = 1,abs(nzg)
       do j = 1,nyg
          do i = 1,nxg
             n = n+1
             gGridRoot%data%rAttr(klat ,n) = lat(i,j)
             gGridRoot%data%rAttr(klon ,n) = lon(i,j)
             gGridRoot%data%rAttr(kmask,n) = real(mask(i,j),R8)
             gGridRoot%data%rAttr(khgt ,n) = hgt(k)
          enddo
       enddo
    enddo
    call mct_gGrid_scatter(gGridRoot, gGrid, gsMap, master_task, mpicom)

    ! Deallocate memory
    if (my_task == master_task) then
       call mct_gGrid_clean(gGridRoot)
    end if
    deallocate(lon)
    deallocate(lat)
    deallocate(mask)
    deallocate(hgt)
    deallocate(gindex_stream)
    deallocate(dimid)

    if (my_task == master_task) write(logunit,*)trim(subname) // ': Successfully initialized stream domain'

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
    character(CS)  :: fldname
    logical        :: compare1
    logical        :: compare2
    integer        :: nfld
    integer        :: rc
    type(ESMF_Field) :: lfield
    character(*)     ,parameter :: F00 = "('(shr_strdata_init_mapping) ',8a)"
    character(len=*) ,parameter :: subname = "(shr_strdata_init_mapping) "
    !-------------------------------------------------------------------------------

    do n = 1,SDAT%nstreams

       method = CompareMaskSubset
       if ( shr_strdata_gGridCompare(SDAT%gridR(n), SDAT%gsmapR(n), &
                                    SDAT%grid    , SDAT%gsmap    , method, mpicom) .or. &
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
       if ( shr_strdata_gGridCompare(SDAT%gridR(n),SDAT%gsmapR(n),SDAT%grid,SDAT%gsmap, method, mpicom, 0.01_r8) .or. &
            trim(SDAT%mapalgo(n))=='none') then
          SDAT%domaps(n) = .false.
       else
          SDAT%domaps(n) = .true.
       endif

       ! ---------------------
       ! Set up fills
       ! ---------------------

       if (SDAT%dofill(n)) then
          if (SDAT%strnzg(n) > 1) then
             write(logunit,*) trim(subname),' do fill called with 3d data, not allowed'
             call shr_sys_abort(subname//': do fill called with 3d data, not allowed')
          endif

          if (trim(SDAT%fillread(n)) == trim(shr_strdata_unset)) then
             if (my_task == master_task) then
                write(logunit,F00) ' calling shr_strdata_mapSet for fill'
             endif

             call shr_strdata_mapSet(SDAT%sMatPf(n), &
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
             endif
             call shr_mct_sMatReaddnc(sMati, SDAT%gsmapR(n), SDAT%gsmapR(n), 'src',  &
                  filename=trim(SDAT%fillread(n)), mytask=my_task, mpicom=mpicom)

             call mct_sMatP_Init(SDAT%sMatPf(n), sMati, SDAT%gsMapR(n), SDAT%gsmapR(n), 0, mpicom, compid)

             call mct_sMat_Clean(sMati)
          endif
       endif

       ! ---------------------
       ! Set up maps
       ! ---------------------

       if (SDAT%domaps(n)) then
          if (SDAT%strnzg(n) > 1) then
             write(logunit,*) trim(subname),' do maps called with 3d data, not allowed'
             call shr_sys_abort(subname//': do maps called with 3d data, not allowed')
          endif

          if (trim(SDAT%mapread(n)) == trim(shr_strdata_unset)) then
             if (my_task == master_task) then
                write(logunit,F00) ' calling shr_strdata_mapSet for remap'
             endif

             call shr_strdata_mapSet(SDAT%sMatPs(n), &
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

    ! Create model datatypes in sdat

    do n = 1,SDAT%nstreams
       ! Determine colon delimited field name string for stream fields
       if (my_task == master_task) then
          call shr_stream_getModelFieldList(SDAT%stream(n),fldList)
       endif
       call shr_mpi_bcast(fldList,mpicom)

       ! Create attribute vectors on model mesh
       call mct_aVect_init(SDAT%avFLB(n),rlist=fldList,lsize=SDAT%lsize)
       call mct_aVect_init(SDAT%avFUB(n),rlist=fldList,lsize=SDAT%lsize)
       call mct_aVect_init(SDAT%avRLB(n),rlist=fldList,lsize=SDAT%lsizeR(n))
       call mct_aVect_init(SDAT%avRUB(n),rlist=fldList,lsize=SDAT%lsizeR(n))

       if (trim(SDAT%tintalgo(n)) == 'coszen') then
          call mct_aVect_init(SDAT%avCoszen(n), rlist="tavCosz", lsize=SDAT%lsize)
       endif

       ! Create field bundle on model mesh for spatially and time
       ! interpolated streams fields to model mesh
       sdat%fldbun_model(n) = ESMF_FieldBundleCreate(rc=rc)
       do nfld = 1,shr_string_listGetNum(fldList)
          ! get nth fldname in colon delimited string
          call shr_string_listGetName(fldlist, nfld, fldname)

          ! create temporary field with name fldname on model mesh
          ! add the field to the lbound and ubound field bundle as well as the time interpolated field bundle
          lfield = ESMF_FieldCreate(sdat%mesh_model, ESMF_TYPEKIND_R8, name=trim(fldname), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          call ESMF_FieldBundleAdd(sdat%fldbun_model(n), (/lfield/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end do
    end do

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
          compare1 = shr_strdata_gGridCompare(SDAT%gridR(nu), SDAT%gsmapR(nu), SDAT%gridR(nv),SDAT%gsmapR(nv), &
               CompareXYabs, mpicom, 0.01_r8)
          compare2 = shr_strdata_gGridCompare(SDAT%gridR(nu), SDAT%gsmapR(nu), SDAT%gridR(nv),SDAT%gsmapR(nv), &
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
  subroutine shr_strdata_get_stream_domain(sdat, stream_index, mpicom, my_task, fldname, flddata)

    ! Obtain the data for fldname from the stream domain data

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    integer                , intent(in)    :: stream_index
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: my_task
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , pointer       :: flddata(:)

    ! local variables
    type(var_desc_t)  :: varid
    type(file_desc_t) :: pioid
    integer           :: rcode
    character(CL)     :: filename
    ! ----------------------------------------------

    if (my_task == master_task) then
       call shr_stream_getDomainFile(sdat%stream(stream_index), filename)
    end if
    call shr_mpi_bcast(filename, mpicom, 'streamfile')
    rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(filename), pio_nowrite)
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, trim(fldname), varid)
    call pio_read_darray(pioid, varid, sdat%pio_iodesc(stream_index), flddata, rcode)
    call pio_closefile(pioid)

  end subroutine shr_strdata_get_stream_domain

  !===============================================================================
  subroutine shr_strdata_advance(SDAT,ymd,tod,mpicom,istr,timers)

    type(shr_strdata_type) ,intent(inout)       :: SDAT
    integer                ,intent(in)          :: ymd    ! current model date
    integer                ,intent(in)          :: tod    ! current model date
    integer                ,intent(in)          :: mpicom
    character(len=*)       ,intent(in),optional :: istr
    logical                ,intent(in),optional :: timers

    ! local variables
    integer                    :: n,m,i,kf,nf          ! generic index
    integer                    :: my_task,npes
    logical                    :: mssrmlf
    logical,allocatable        :: newData(:)
    integer                    :: ierr
    integer                    :: nu,nv
    integer                    :: lsize,lsizeR,lsizeF
    integer,allocatable        :: ymdmod(:)            ! modified model dates to handle Feb 29
    integer                    :: todmod               ! modified model dates to handle Feb 29
    type(mct_avect)            :: avRtmp
    type(mct_avect)            :: avRV,avFV
    character(len=32)          :: lstr
    logical                    :: ltimers
    real(R8)                   :: flb,fub              ! factor for lb and ub
    real(R8),pointer           :: lonr(:)              ! lon radians, cosz method
    real(R8),pointer           :: latr(:)              ! lat radians, cosz method
    real(R8),pointer           :: cosz(:)              ! cosz, cosz method
    real(R8),pointer           :: tavCosz(:)           ! cosz, time avg over [LB,UB], cosz method
    real(R8),pointer           :: xlon(:),ylon(:)
    type(ESMF_Time)            :: timeLB, timeUB       ! lb and ub times
    type(ESMF_TimeInterval)    :: timeint              ! delta time
    integer                    :: dday                 ! delta days
    real(R8)                   :: dtime                ! delta time
    integer                    :: uvar,vvar
    character(CS)              :: uname                ! u vector field name
    character(CS)              :: vname                ! v vector field name
    integer                    :: year,month,day       ! date year month day
    real(r8), pointer          :: dataptr(:) 
    integer                    :: rc
    integer                    :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: fieldnamelist(:)
    real(R8)         ,parameter :: solZenMin = 0.001_R8 ! minimum solar zenith angle
    character(len=*) ,parameter :: timname = "_strd_adv"
    integer          ,parameter :: tadj = 2
    character(*)     ,parameter :: subname = "(shr_strdata_advance) "
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

          call shr_strdata_readLBUB(SDAT%stream(n),&
               SDAT%pio_subsystem, SDAT%io_type, SDAT%pio_iodesc(n),  &
               ymdmod(n), todmod, mpicom, SDAT%gsmapR(n), &
               SDAT%avRLB(n), SDAT%ymdLB(n), SDAT%todLB(n),  &
               SDAT%avRUB(n), SDAT%ymdUB(n), SDAT%todUB(n),  &
               SDAT%avRFile(n),  trim(SDAT%readmode(n)),  newData(n),  &
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
          if ((SDAT%domaps(nu) .or. SDAT%domaps(nv)) .and. (newdata(nu) .or. newdata(nv))) then

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

             ! do not need to recompute uvar and vvar, should be the same
             do i = 1,lsizeF
                SDAT%avFUB(nu)%rAttr(uvar,i) =  avFV%rAttr(1,i) * cos(ylon(i)) + avFV%rAttr(2,i) * sin(ylon(i))
                SDAT%avFUB(nv)%rAttr(vvar,i) = -avFV%rAttr(1,i) * sin(ylon(i)) + avFV%rAttr(2,i) * cos(ylon(i))
             enddo

             call mct_aVect_clean(avRV)
             call mct_aVect_clean(avFV)
             deallocate(xlon,ylon)

             call t_stopf(trim(lstr)//trim(timname)//'_vect')
          endif
       enddo

       ! ---------------------------------------------------------
       ! Do time interpolation to create FB_model
       ! ---------------------------------------------------------

       do n = 1,SDAT%nstreams

          ! Get field namelist
          call ESMF_FieldBundleGet(sdat%fldbun_model(n), fieldCount=fieldCount, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
          allocate(fieldnamelist(fieldCount))
          call ESMF_FieldBundleGet(sdat%fldbun_model(n), fieldNameList=fieldnamelist, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
               call ESMF_Finalize(endflag=ESMF_END_ABORT)

          if (trim(SDAT%tintalgo(n)) == 'coszen') then

             ! ------------------------------------------
             ! time interpolation method is coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_coszen')

             ! make sure orb info has been set
             if (SDAT%eccen == SHR_ORB_UNDEF_REAL) then
                call shr_sys_abort(subname//' ERROR in orb params for coszen tinterp')
             else if (SDAT%modeldt < 1) then
                call shr_sys_abort(subname//' ERROR: model dt < 1 for coszen tinterp')
             endif

             !--- allocate avg cosz array ---
             lsizeF = sdat%lsize
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

             ! t-interp is LB data normalized with this factor: cosz/tavCosz

             do nf = 1,fieldcount
                call dshr_fldbun_getfldptr(sdat%fldbun_model(n), fieldnamelist(nf), dataptr, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
                     call ESMF_Finalize(endflag=ESMF_END_ABORT)
                kf = mct_aVect_indexRA(sdat%avFLB(n), trim(fieldnamelist(nf)))
                do i = 1,lsize
                   if (cosz(i) > solZenMin) then
                      dataptr(i) = SDAT%avFLB(n)%rAttr(kf,i)*cosz(i)/tavCosz(i)
                   else
                      dataptr(i) = 0._r8
                   endif
                end do
             end do

             deallocate(tavCosz,cosz,lonr,latr)
             call t_stopf(trim(lstr)//trim(timname)//'_coszen')

          elseif (trim(SDAT%tintalgo(n)) /= trim(shr_strdata_nullstr)) then

             ! ------------------------------------------
             ! time interpolation method is not coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_tint')
             call shr_tInterp_getFactors(SDAT%ymdlb(n),SDAT%todlb(n),SDAT%ymdub(n),SDAT%todub(n), &
                  ymdmod(n),todmod,flb,fub, calendar=SDAT%stream(n)%calendar,algo=trim(SDAT%tintalgo(n)))
             if (debug > 0) then
                write(logunit,*) trim(subname),' interp = ',n,flb,fub
             endif

             do nf = 1,fieldcount
                call dshr_fldbun_getfldptr(sdat%fldbun_model(n), fieldnamelist(nf), dataptr, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
                     call ESMF_Finalize(endflag=ESMF_END_ABORT)
                kf = mct_aVect_indexRA(sdat%avFLB(n), trim(fieldnamelist(nf)))
                dataptr(:) = sdat%avFLB(n)%rAttr(kf,:)*flb + sdat%avFUB(n)%rAttr(kf,:)*fub
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_tint')

          else

             ! ------------------------------------------
             ! zero out stream data for this field
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_zero')
             do nf = 1,fieldcount
                call dshr_fldbun_getfldptr(sdat%fldbun_model(n), fieldnamelist(nf), dataptr, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
                     call ESMF_Finalize(endflag=ESMF_END_ABORT)
                dataptr(:) = 0._r8
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_zero')
          endif

       enddo

       deallocate(newData)
       deallocate(ymdmod)

    endif    ! nstreams > 0

    call t_stopf(trim(lstr)//trim(timname)//'_total')
    if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine shr_strdata_advance

  !===============================================================================

  subroutine shr_strdata_restWrite(filename,SDAT,mpicom,str1,str2)

    character(len=*)       ,intent(in)    :: filename
    type(shr_strdata_type) ,intent(inout) :: SDAT
    integer                ,intent(in)    :: mpicom
    character(len=*)       ,intent(in)    :: str1
    character(len=*)       ,intent(in)    :: str2

    !--- local ----
    integer :: my_task,ier
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom,my_task,ier)
    if (my_task == 0) then
       call shr_stream_restWrite(SDAT%stream,trim(filename),trim(str1),trim(str2),SDAT%nstreams)
    endif

  end subroutine shr_strdata_restWrite

  !===============================================================================

  subroutine shr_strdata_restRead(filename,SDAT,mpicom)

    character(len=*)       ,intent(in)    :: filename
    type(shr_strdata_type) ,intent(inout) :: SDAT
    integer                ,intent(in)    :: mpicom

    !--- local ----
    integer :: my_task,ier
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom,my_task,ier)
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
  subroutine shr_strdata_readnml(SDAT, file, mpicom)

    ! Reads shr_strdata_nml namelist input

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: SDAT   ! strdata data data-type
    character(*) ,optional ,intent(in)   :: file   ! file to read strdata from
    integer      ,optional ,intent(in)   :: mpicom ! mpi comm

    ! local variables
    integer    :: rCode         ! return code
    integer    :: nUnit         ! fortran i/o unit number
    integer    :: n             ! generic loop index
    integer    :: my_task       ! my task number, 0 is default
    integer    :: ntasks        ! total number of tasks

    ! shr_strdata_nml namelist variables
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
           nx_global ,           &
           ny_global ,           &
           streams   ,           &
           taxMode   ,           &
           dtlimit   ,           &
           vectors   ,           &
           fillalgo  ,           &
           fillmask  ,           &
           fillread  ,           &
           fillwrite ,           &
           mapalgo   ,           &
           mapmask   ,           &
           mapread   ,           &
           mapwrite  ,           &
           tintalgo  ,           &
           readmode   

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

    my_task = 0
    ntasks = 1

    if (present(mpicom)) then
       call mpi_comm_rank(mpicom, my_task, rCode)
       call mpi_comm_size(mpicom, ntasks, rCode)
    endif

    !--master--task--
    if (my_task == master_task) then

       ! set default values for namelist vars
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

       ! read input namelist
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

       ! copy temporary/local namelist vars into data structure
       SDAT%nstreams    = 0
       do n=1,nStrMax
          call shr_stream_default(SDAT%stream(n))
       enddo
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
    endif   ! master_task

    if (present(mpicom)) then
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

  end subroutine shr_strdata_print

  subroutine shr_strdata_readLBUB(stream,pio_subsystem,pio_iotype,pio_iodesc,&
       mDate,mSec,mpicom,gsMap, &
       avLB,mDateLB,mSecLB,avUB,mDateUB,mSecUB,avFile,readMode, &
       newData,rmOldFile,istr)

    !-------------------------------------------------------------------------
    ! Read LB and UB of stream data
    !-------------------------------------------------------------------------

    !----- arguments -----
    type(shr_stream_streamType) ,intent(inout) :: stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)    :: pio_iotype
    type(io_desc_t)             ,intent(inout) :: pio_iodesc
    integer                     ,intent(in)    :: mDate  ,mSec
    integer                     ,intent(in)    :: mpicom
    type(mct_gsMap)             ,intent(in)    :: gsMap
    type(mct_aVect)             ,intent(inout) :: avLB
    integer                     ,intent(inout) :: mDateLB,mSecLB
    type(mct_aVect)             ,intent(inout) :: avUB
    integer                     ,intent(inout) :: mDateUB,mSecUB
    type(mct_aVect)             ,intent(inout) :: avFile
    character(len=*)            ,intent(in)    :: readMode
    logical                     ,intent(out)   :: newData
    logical          ,optional  ,intent(in)    :: rmOldFile
    character(len=*) ,optional  ,intent(in)    :: istr

    !----- local -----
    integer           :: my_task
    integer           :: ierr       ! error code
    integer           :: rCode      ! return code
    logical           :: localCopy,fileexists
    integer           :: ivals(6)   ! bcast buffer
    integer           :: oDateLB,oSecLB,dDateLB,oDateUB,oSecUB,dDateUB
    real(R8)          :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer           :: n_lb, n_ub
    character(CL)     :: fn_lb,fn_ub,fn_next,fn_prev
    character(CL)     :: path
    character(len=32) :: lstr
    real(R8)          :: spd
    character(*), parameter :: subname = '(shr_strdata_readLBUB) '
    character(*), parameter :: F00   = "('(shr_strdata_readLBUB) ',8a)"
    character(*), parameter :: F01   = "('(shr_strdata_readLBUB) ',a,5i8)"
    !-------------------------------------------------------------------------

    lstr = 'shr_strdata_readLBUB'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    call t_startf(trim(lstr)//'_setup')
    call MPI_COMM_RANK(mpicom,my_task,ierr)
    spd = shr_const_cday

    newData = .false.
    n_lb = -1
    n_ub = -1
    fn_lb = 'undefinedlb'
    fn_ub = 'undefinedub'

    oDateLB = mDateLB
    oSecLB  = mSecLB
    oDateUB = mDateUB
    oSecUB  = mSecUB

    rDateM  = real(mDate  ,R8) + real(mSec  ,R8)/spd
    rDateLB = real(mDateLB,R8) + real(mSecLB,R8)/spd
    rDateUB = real(mDateUB,R8) + real(mSecUB,R8)/spd
    call t_stopf(trim(lstr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(lstr)//'_fbound')
       if (my_task == master_task) then
          call shr_stream_findBounds(stream,mDate,mSec,                 &
               ivals(1),dDateLB,ivals(2),ivals(5),fn_lb, &
               ivals(3),dDateUB,ivals(4),ivals(6),fn_ub  )
          call shr_stream_getFilePath(stream,path)
       endif
       call t_stopf(trim(lstr)//'_fbound')
       call t_startf(trim(lstr)//'_bcast')
       call shr_mpi_bcast(stream%calendar,mpicom)
       call shr_mpi_bcast(ivals,mpicom)
       mDateLB = ivals(1)
       mSecLB  = ivals(2)
       mDateUB = ivals(3)
       mSecUB  = ivals(4)
       n_lb    = ivals(5)
       n_ub    = ivals(6)
       call t_stopf(trim(lstr)//'_bcast')
    endif

    if (mDateLB /= oDateLB .or. mSecLB /= oSecLB) then
       newdata = .true.
       if (mDateLB == oDateUB .and. mSecLB == oSecUB) then
          call t_startf(trim(lstr)//'_LB_copy')
          avLB%rAttr(:,:) = avUB%rAttr(:,:)
          call t_stopf(trim(lstr)//'_LB_copy')
       else
          select case(readMode)
          case ('single')
             call shr_strdata_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, &
                  gsMap, avLB, mpicom, &
                  path, fn_lb, n_lb,istr=trim(lstr)//'_LB', boundstr = 'lb')
          case ('full_file')
             call shr_strdata_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
                  gsMap, avLB, avFile, mpicom, &
                  path, fn_lb, n_lb,istr=trim(lstr)//'_LB', boundstr = 'lb')
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
          end select
       endif
    endif

    if (mDateUB /= oDateUB .or. mSecUB /= oSecUB) then
       newdata = .true.
       select case(readMode)
       case ('single')
          call shr_strdata_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, avUB, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case ('full_file')
          call shr_strdata_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
               gsMap, avUB, avFile, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case default
          write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
          call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
       end select
    endif

    call t_startf(trim(lstr)//'_filemgt')
    !--- determine previous & next data files in list of files ---
    if (my_task == master_task .and. newdata) then
       call shr_stream_getFilePath(stream,path)
       call shr_stream_getPrevFileName(stream,fn_lb,fn_prev,path)
       call shr_stream_getNextFileName(stream,fn_ub,fn_next,path)
       inquire(file=trim(fn_next),exist=fileExists)
       if ( trim(fn_next) == "unknown" .or. fileExists) then
          ! do nothing
       end if
    endif
    call t_stopf(trim(lstr)//'_filemgt')

  end subroutine shr_strdata_readLBUB

  !===============================================================================

  subroutine shr_strdata_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, av, mpicom, &
       path, fn, nt, istr, boundstr)

    !----- arguments -----
    type(shr_stream_streamType) ,intent(inout)         :: stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)            :: pio_iotype
    type(io_desc_t)             ,intent(inout)         :: pio_iodesc
    type(mct_gsMap)             ,intent(in)            :: gsMap
    type(mct_aVect)             ,intent(inout)         :: av
    integer                     ,intent(in)            :: mpicom
    character(len=*)            ,intent(in)            :: path
    character(len=*)            ,intent(in)            :: fn
    integer                     ,intent(in)            :: nt
    character(len=*),optional   ,intent(in)            :: istr
    character(len=*),optional   ,intent(in)            :: boundstr

    !----- local -----
    integer              :: my_task
    integer              :: ierr
    logical              :: fileexists
    integer              :: gsize,nx,ny,nz
    integer              :: k
    integer              :: rCode      ! return code
    character(CL)        :: fileName
    character(CL)        :: sfldName
    type(mct_avect)      :: avtmp
    character(len=32)    :: lstr
    character(len=32)    :: bstr
    logical              :: fileopen
    character(CL)        :: currfile
    integer              :: ndims
    integer,pointer      :: dimid(:)
    type(file_desc_t)    :: pioid
    type(var_desc_t)     :: varid
    integer(kind=pio_offset_kind) :: frame
    character(*), parameter :: subname = '(shr_strdata_readstrm) '
    character(*), parameter :: F00   = "('(shr_strdata_readstrm) ',8a)"
    character(*), parameter :: F01   = "('(shr_strdata_readstrm) ',a,5i8)"
    character(*), parameter :: F02   = "('(shr_strdata_readstrm) ',2a,i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_strdata_readstrm'
    if (present(istr)) then
       lstr = trim(istr)
    endif
    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    call mpi_comm_rank(mpicom,my_task,ierr)

    if (my_task == master_task) then
       fileName = trim(path)//trim(fn)
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif

    if (my_task == master_task) then
       call shr_stream_getFileFieldName(stream,1,sfldName)
    endif

    call t_stopf(trim(lstr)//'_setup')

    call t_startf(trim(lstr)//'_readpio')
    call shr_mpi_bcast(sfldName,mpicom,'sfldName')
    call shr_mpi_bcast(filename,mpicom,'filename')

    call shr_stream_getCurrFile(stream,fileopen=fileopen,currfile=currfile,currpioid=pioid)

    if (fileopen .and. currfile==filename) then
       ! don't reopen file, all good
    else
       ! otherwise close the old file if open and open new file
       if (fileopen) then
          if (my_task == master_task) write(logunit,F00) 'close  : ',trim(currfile)
          call pio_closefile(pioid)
       endif
       if (my_task == master_task) write(logunit,F00) 'open   : ',trim(filename)
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call shr_stream_setCurrFile(stream, fileopen=.true., currfile=trim(filename), currpioid=pioid)
    endif
    if (my_task == master_task) write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename),nt

    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

    rcode = pio_inq_varid(pioid,trim(sfldName),varid)
    rcode = pio_inq_varndims(pioid, varid, ndims)
    allocate(dimid(ndims))
    rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
    if (ndims >= 1) rcode = pio_inq_dimlen(pioid, dimid(1), nx)
    if (ndims >= 2) rcode = pio_inq_dimlen(pioid, dimid(2), ny)
    if (ndims >= 3) rcode = pio_inq_dimlen(pioid, dimid(3), nz)
    deallocate(dimid)
    gsize = mct_gsmap_gsize(gsMap)
    if (gsize /= nx*ny .and. gsize /= nx*ny*nz) then
       write(logunit,F01) "ERROR in data sizes ",nx,ny,nz,gsize
       call shr_sys_abort(subname//"ERROR in data sizes")
    endif

    do k = 1,mct_aVect_nRAttr(av)
       if (my_task == master_task) then
          call shr_stream_getFileFieldName(stream,k,sfldName)
       endif
       call shr_mpi_bcast(sfldName,mpicom,'sfldName')
       rcode = pio_inq_varid(pioid,trim(sfldName),varid)
       frame = nt
       call pio_setframe(pioid,varid,frame)
       call pio_read_darray(pioid,varid,pio_iodesc,av%rattr(k,:),rcode)
    enddo

    call t_stopf(trim(lstr)//'_readpio')

  end subroutine shr_strdata_readstrm

  !===============================================================================

  subroutine shr_strdata_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
       gsMap, av, avFile, mpicom, &
       path, fn, nt, istr, boundstr)

    !----- arguments -----
    type(shr_stream_streamType) ,intent(inout)         :: stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)            :: pio_iotype
    type(mct_gsMap)             ,intent(in)            :: gsMap
    type(mct_aVect)             ,intent(inout)         :: av
    type(mct_aVect)             ,intent(inout)         :: avFile
    integer                     ,intent(in)            :: mpicom
    character(len=*)            ,intent(in)            :: path
    character(len=*)            ,intent(in)            :: fn
    integer                     ,intent(in)            :: nt
    character(len=*)            ,intent(in) ,optional  :: istr
    character(len=*)            ,intent(in) ,optional  :: boundstr

    !----- local -----
    integer                       :: my_task
    integer                       :: ierr
    logical                       :: localCopy,fileexists
    integer                       :: gsize,nx,ny,nz
    integer                       :: k
    integer                       :: rCode   ! return code
    character(CL)                 :: fileName
    character(CL)                 :: sfldName
    character(len=32)             :: lstr
    character(len=32)             :: bstr
    logical                       :: fileopen
    character(CL)                 :: currfile
    character(CXX)                :: fldList ! list of fields
    integer                       :: ndims
    integer,pointer               :: dimid(:)
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    type(io_desc_t)               :: pio_iodesc_local
    integer                       :: avFile_beg, avFile_end
    integer                       :: lsize, cnt,m,n
    integer, allocatable          :: count(:), compDOF(:)
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    character(*), parameter :: subname = ' (shr_strdata_readstrm_fullfile) '
    character(*), parameter :: F00   = "(' (shr_strdata_readstrm_fullfile) ',8a)"
    character(*), parameter :: F01   = "(' (shr_strdata_readstrm_fullfile) ',a,5i8)"
    character(*), parameter :: F02   = "(' (shr_strdata_readstrm_fullfile) ',2a,2i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_strdata_readstrm_fullfile'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    call mpi_comm_rank(mpicom,my_task,ierr)

    if (my_task == master_task) then
       fileName = trim(path) // trim(fn)
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif
    if (my_task == master_task) then
       call shr_stream_getFileFieldName(stream,1,sfldName)
    endif
    call t_stopf(trim(lstr)//'_setup')

    call t_startf(trim(lstr)//'_readpio')
    call shr_mpi_bcast(sfldName,mpicom,'sfldName')
    call shr_mpi_bcast(filename,mpicom,'filename')
    call shr_stream_getCurrFile(stream,fileopen=fileopen,currfile=currfile,currpioid=pioid)

    if (fileopen .and. currfile==filename) then
       ! don't reopen file, all good
    else
       ! otherwise close the old file if open, open the new file,
       ! and read all time slices of a temporal dataset within the new file.
       if (fileopen) then
          if (my_task == master_task) write(logunit,F00) 'close  : ',trim(currfile)
          call pio_closefile(pioid)
       endif
       if (my_task == master_task) write(logunit,F00) 'open   : ',trim(filename)
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call shr_stream_setCurrFile(stream,fileopen=.true.,currfile=trim(filename),currpioid=pioid)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
       rcode = pio_inq_varid(pioid,trim(sfldName),varid)
       rcode = pio_inq_varndims(pioid, varid, ndims)

       allocate(dimid(ndims))
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       nx = 1
       ny = 1
       nz = 1
       if (ndims >= 1) rcode = pio_inq_dimlen(pioid, dimid(1), nx)
       if (ndims >= 2) rcode = pio_inq_dimlen(pioid, dimid(2), ny)
       if (ndims >= 3) rcode = pio_inq_dimlen(pioid, dimid(3), nz)
       deallocate(dimid)

       gsize = mct_gsmap_gsize(gsMap)
       if (gsize /= nx*ny) then
          write(logunit,F01) "ERROR in data sizes ",nx,ny,nz,gsize
          call shr_sys_abort(subname//"ERROR in data sizes")
       endif

       if (my_task == master_task) then
          call shr_stream_getModelFieldList(stream,fldList)
       endif
       call shr_mpi_bcast(fldList,mpicom)
       call mct_avect_clean(avFile)
       call mct_aVect_init(avFile,rlist=fldList,lsize=nx*ny*nz)

       lsize = mct_gsmap_lsize(gsMap,mpicom)
       allocate(compDOF(lsize*nz),stat=rcode)
       if (rcode /= 0) call shr_sys_abort(subname//"ERROR insufficient memory")
       if (my_task == master_task) then
          write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename),1,nz
          call shr_sys_flush(logunit)
       endif

       ! Create a 3D MCT component DOF corresponding to "2D(=gsmOP) x nz"
       call mct_gsmap_orderedPoints(gsMap,my_task,gsmOP)
       cnt = 0
       do n = 1,nz
          do m = 1,lsize
             cnt = cnt + 1
             compDOF(cnt) = (n-1)*gsize + gsmOP(m)
          enddo
       enddo
       
       ! Initialize the decomposition
       allocate(count(3))
       count(1) = nx
       count(2) = ny
       count(3) = nz
       call pio_initdecomp(pio_subsystem, pio_double, count, compDOF, pio_iodesc_local)
       deallocate(count)
       deallocate(compDOF)

       ! For each attribute, read all frames in one go
       frame = 1
       do k = 1, mct_aVect_nRAttr(avFile)
          if (my_task == master_task) then
             call shr_stream_getFileFieldName(stream,k,sfldName)
          endif
          call shr_mpi_bcast(sfldName,mpicom,'sfldName')
          rcode = pio_inq_varid(pioid,trim(sfldName),varid)
          call pio_setframe(pioid,varid,frame)
          call pio_read_darray(pioid, varid, pio_iodesc_local, avFile%rattr(k,:), rcode)
       enddo
       call pio_freedecomp(pio_subsystem, pio_iodesc_local)

    endif

    ! Copy the `nt` time slice data from avFile into av
    avFile_beg = lsize*(nt-1) + 1
    avFile_end = lsize*nt
    do k = 1, mct_aVect_nRAttr(avFile)
       av%rattr(k,1:lsize) = avFile%rattr(k,avFile_beg:avFile_end)
    enddo

    call t_stopf(trim(lstr)//'_readpio')

  end subroutine shr_strdata_readstrm_fullfile

  !===============================================================================

  logical function shr_strdata_gGridCompare(ggrid1,gsmap1,ggrid2,gsmap2,method,mpicom,eps)

    ! Returns TRUE if two domains are the the same (within tolerance).

    ! input/output variables
    type(mct_gGrid)   ,intent(in)  :: ggrid1   ! 1st ggrid
    type(mct_gsmap)   ,intent(in)  :: gsmap1   ! 1st gsmap
    type(mct_gGrid)   ,intent(in)  :: ggrid2   ! 2nd ggrid
    type(mct_gsmap)   ,intent(in)  :: gsmap2   ! 2nd gsmap
    integer           ,intent(in)  :: method   ! selects what to compare
    integer           ,intent(in)  :: mpicom   ! mpicom
    real(R8),optional ,intent(in)  :: eps      ! epsilon compare value

    ! local variables
    real(R8)        :: leps         ! local epsilon
    integer         :: n            ! counters
    integer         :: my_task
    integer         :: gsize
    integer         :: ierr
    integer         :: nlon1, nlon2, nlat1, nlat2, nmask1, nmask2  ! av field indices
    logical         :: compare               ! local compare logical
    real(R8)        :: lon1,lon2             ! longitudes to compare
    real(R8)        :: lat1,lat2             ! latitudes to compare
    real(R8)        :: msk1,msk2             ! masks to compare
    integer         :: nx,ni1,ni2            ! i grid size, i offset for 1 vs 2 and 2 vs 1
    integer         :: n1,n2,i,j             ! local indices
    type(mct_aVect) :: avG1                  ! global av
    type(mct_aVect) :: avG2                  ! global av
    integer         :: lmethod               ! local method
    logical         :: maskmethod, maskpoint ! masking on method
    character(*),parameter :: subName = '(shr_strdata_gGridCompare) '
    character(*),parameter :: F01     = "('(shr_strdata_gGridCompare) ',4a)"
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)

    leps = 1.0e-6_R8
    if (present(eps)) leps = eps

    lmethod = mod(method,100)
    if (method > 100) then
       maskmethod=.true.
    else
       maskmethod=.false.
    endif

    call mct_aVect_gather(gGrid1%data,avG1,gsmap1,master_task,mpicom)
    call mct_aVect_gather(gGrid2%data,avG2,gsmap2,master_task,mpicom)

    if (my_task == master_task) then

       compare = .true.
       gsize = mct_aVect_lsize(avG1)
       if (gsize /= mct_aVect_lsize(avG2)) then
          compare = .false.
       endif

       if (.not. compare ) then
          !--- already failed the comparison test, check no futher ---
       else
          nlon1 = mct_aVect_indexRA(avG1,'lon')
          nlat1 = mct_aVect_indexRA(avG1,'lat')
          nlon2 = mct_aVect_indexRA(avG2,'lon')
          nlat2 = mct_aVect_indexRA(avG2,'lat')
          nmask1 = mct_aVect_indexRA(avG1,'mask')
          nmask2 = mct_aVect_indexRA(avG2,'mask')

          ! To compare, want to be able to treat longitude wraparound generally.
          ! So we need to compute i index offset and we need to compute the size of the nx dimension
          ! First adjust the lon so it's in the range [0,360), add 1440 to lon to take into
          ! accounts lons less than 1440.  if any lon is less than -1440, abort.  1440 is arbitrary
          ! Next, comute ni1 and ni2.  These are the offsets of grid1 relative to grid2 and
          ! grid2 relative to grid1.  The sum of those offsets is nx.  Use ni1 to offset grid2
          ! in comparison and compute new grid2 index from ni1 and nx.  If ni1 is zero, then
          ! there is no offset, don't need to compute ni2, and nx can be anything > 0.

          !--- compute offset of grid2 compared to pt 1 of grid 1
          lon1 = minval(avG1%rAttr(nlon1,:))
          lon2 = minval(avG2%rAttr(nlon2,:))
          if ((lon1 < -1440.0_R8) .or. (lon2 < -1440.0_R8)) then
             write(logunit,*) subname,' ERROR: lon1 lon2 lt -1440 ',lon1,lon2
             call shr_sys_abort(subname//' ERROR: lon1 lon2 lt -1440')
          endif

          lon1 = mod(avG1%rAttr(nlon1,1)+1440.0_R8,360.0_R8)
          lat1 = avG1%rAttr(nlat1,1)
          ni1 = -1
          do n = 1,gsize
             lon2 = mod(avG2%rAttr(nlon2,n)+1440.0_R8,360.0_R8)
             lat2 = avG2%rAttr(nlat2,n)
             if ((ni1 < 0) .and. abs(lon1-lon2) <= leps .and. abs(lat1-lat2) <= leps) then
                ni1 = n - 1  ! offset, compare to first gridcell in grid 1
             endif
          enddo

          if (ni1 < 0) then        ! no match for grid point 1, so fails without going further
             compare = .false.
          elseif (ni1 == 0) then   ! no offset, set nx to anything > 0
             nx = 1
          else                     ! now compute ni2
             ! compute offset of grid1 compared to pt 1 of grid 2
             lon2 = mod(avG2%rAttr(nlon2,1)+1440.0_R8,360.0_R8)
             lat2 = avG2%rAttr(nlat2,1)
             ni2 = -1
             do n = 1,gsize
                lon1 = mod(avG1%rAttr(nlon1,n)+1440.0_R8,360.0_R8)
                lat1 = avG1%rAttr(nlat1,n)
                if ((ni2 < 0) .and. abs(lon1-lon2) <= leps .and. abs(lat1-lat2) <= leps) then
                   ni2 = n - 1  ! offset, compare to first gridcell in grid 1
                endif
             enddo
             if (ni2 < 0) then
                write(logunit,*) subname,' ERROR in ni2 ',ni1,ni2
                call shr_sys_abort(subname//' ERROR in ni2')
             endif
             nx = ni1 + ni2
          endif

          if (compare) then
             do n = 1,gsize
                j = ((n-1)/nx) + 1
                i = mod(n-1,nx) + 1
                n1 = (j-1)*nx + mod(n-1,nx) + 1
                n2 = (j-1)*nx + mod(n-1+ni1,nx) + 1
                if (n1 /= n) then    ! sanity check, could be commented out
                   write(logunit,*) subname,' ERROR in n1 n2 ',n,i,j,n1,n2
                   call shr_sys_abort(subname//' ERROR in n1 n2')
                endif
                lon1 = mod(avG1%rAttr(nlon1,n1)+1440.0_R8,360.0_R8)
                lat1 = avG1%rAttr(nlat1,n1)
                lon2 = mod(avG2%rAttr(nlon2,n2)+1440.0_R8,360.0_R8)
                lat2 = avG2%rAttr(nlat2,n2)
                msk1 = avG1%rAttr(nmask1,n1)
                msk2 = avG2%rAttr(nmask2,n2)

                maskpoint = .true.
                if (maskmethod .and. (msk1 == 0 .or. msk2 == 0)) then
                   maskpoint = .false.
                endif

                if (maskpoint) then
                   if (lmethod == CompareXYabs      ) then
                      if (abs(lon1 - lon2) > leps) compare = .false.
                      if (abs(lat1 - lat2) > leps) compare = .false.
                   else if (lmethod == CompareXYrel      ) then
                      if (rdiff(lon1,lon2) > leps) compare = .false.
                      if (rdiff(lat1,lat2) > leps) compare = .false.
                   else if (lmethod == CompareMaskIdent  ) then
                      if (msk1 /= msk2)compare = .false.
                   else if (lmethod == CompareMaskZeros  ) then
                      if (msk1 == 0 .and. msk2 /= 0) compare = .false.
                      if (msk1 /= 0 .and. msk2 == 0) compare = .false.
                   else if (lmethod == CompareMaskSubset ) then
                      if (msk1 /= 0 .and. msk2 == 0) compare = .false.
                   else
                      write(logunit,F01) "ERROR: compare method not recognized, method = ",method
                      call shr_sys_abort(subName//"ERROR: compare method not recognized")
                   endif  ! lmethod
                endif  ! maskpoint
             enddo ! gsize
          endif  ! compare
       endif   ! compare
    endif   ! master_task

    if (my_task == master_task) call mct_avect_clean(avG1)
    if (my_task == master_task) call mct_avect_clean(avG2)

    call shr_mpi_bcast(compare,mpicom)
    shr_strdata_gGridCompare = compare

  contains   ! internal subprogram

    real(R8) function rdiff(v1,v2) ! internal function
      !------------------------------------------
      real(R8),intent(in) :: v1,v2                 ! two values to compare
      real(R8),parameter  :: c0           = 0.0_R8 ! zero
      real(R8),parameter  :: large_number = 1.0e20_R8 ! infinity
      !------------------------------------------
      if (v1 == v2) then
         rdiff = c0
      elseif (v1 == c0 .and. v2 /= c0) then
         rdiff = large_number
      elseif (v2 == c0 .and. v1 /= c0) then
         rdiff = large_number
      else
         !        rdiff = abs((v2-v1)/v1)   ! old version, but rdiff(v1,v2) /= vdiff(v2,v1)
         rdiff = abs(2.0_R8*(v2-v1)/(v1+v2))
      endif
      !------------------------------------------
    end function rdiff

  end function shr_strdata_gGridCompare

  !===============================================================================

  subroutine shr_strdata_mapSet(smatp,&
       ggridS,gsmapS,nxgS,nygS,&
       ggridD,gsmapD,nxgD,nygD, &
       name,type,algo,mask,vect,&
       compid,mpicom,strategy)

    !-------------------------------------------------------------------------------
    ! Initialize sparse mapping routine handle, sMatP, from source and destination
    ! mct gGridS and gGridD
    !-------------------------------------------------------------------------------

    ! input/output variables
    type(mct_sMatP)  , intent(inout)       :: smatp
    type(mct_gGrid)  , intent(in)          :: ggridS
    type(mct_gsmap)  , intent(in)          :: gsmapS
    integer          , intent(in)          :: nxgS
    integer          , intent(in)          :: nygS
    type(mct_gGrid)  , intent(in)          :: ggridD
    type(mct_gsmap)  , intent(in)          :: gsmapD
    integer          , intent(in)          :: nxgD
    integer          , intent(in)          :: nygD
    character(len=*) , intent(in)          :: name
    character(len=*) , intent(in)          :: type
    character(len=*) , intent(in)          :: algo
    character(len=*) , intent(in)          :: mask
    character(len=*) , intent(in)          :: vect
    integer          , intent(in)          :: compid
    integer          , intent(in)          :: mpicom
    character(len=*) , intent(in),optional :: strategy

    ! local variables
    integer               :: n,i,j
    integer               :: lsizeS,gsizeS,lsizeD,gsizeD
    integer               :: nlon,nlat,nmsk
    integer               :: my_task,ierr
    real(R8) , pointer    :: Xsrc(:,:)
    real(R8) , pointer    :: Ysrc(:,:)
    integer  , pointer    :: Msrc(:,:)
    real(R8) , pointer    :: Xdst(:,:)
    real(R8) , pointer    :: Ydst(:,:)
    integer  , pointer    :: Mdst(:,:)
    type(shr_map_mapType) :: shrmap
    type(mct_aVect)       :: AVl
    type(mct_aVect)       :: AVg
    character(len=32)     :: lstrategy
    integer               :: nsrc,ndst,nwts
    integer  , pointer    :: isrc(:)
    integer  , pointer    :: idst(:)
    real(R8) , pointer    :: wgts(:)
    type(mct_sMat)        :: sMat0
    character(*), parameter :: subname = '(shr_strdata_mapSet) '
    character(*), parameter :: F00   = "('(shr_strdata_mapSet) ',8a)"
    character(*), parameter :: F01   = "('(shr_strdata_mapSet) ',a,5i8)"
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)

    !--- get sizes and allocate for SRC ---

    lsizeS = mct_aVect_lsize(ggridS%data)
    call mct_avect_init(AVl,rList='lon:lat:mask',lsize=lsizeS)
    call mct_avect_copy(ggridS%data,AVl,rList='lon:lat:mask')
    call mct_avect_gather(AVl,AVg,gsmapS,master_task,mpicom)

    if (my_task == master_task) then
       gsizeS = mct_aVect_lsize(AVg)
       if (gsizeS /= nxgS*nygS) then
          write(logunit,F01) ' ERROR: gsizeS ',gsizeS,nxgS,nygS
          call shr_sys_abort(subname//' ERROR gsizeS')
       endif
       allocate(Xsrc(nxgS,nygS),Ysrc(nxgS,nygS),Msrc(nxgS,nygS))

       nlon = mct_avect_indexRA(AVg,'lon')
       nlat = mct_avect_indexRA(AVg,'lat')
       nmsk = mct_avect_indexRA(AVg,'mask')

       n = 0
       Msrc = 1
       do j = 1,nygS
          do i = 1,nxgS
             n = n + 1
             Xsrc(i,j) = AVg%rAttr(nlon,n)
             Ysrc(i,j) = AVg%rAttr(nlat,n)
             if (abs(AVg%rAttr(nmsk,n)) < 0.5_R8) Msrc(i,j) = 0
          enddo
       enddo
    endif

    if (my_task == master_task) call mct_aVect_clean(AVg)
    call mct_aVect_clean(AVl)

    !--- get sizes and allocate for DST ---

    lsizeD = mct_aVect_lsize(ggridD%data)
    call mct_avect_init(AVl,rList='lon:lat:mask',lsize=lsizeD)
    call mct_avect_copy(ggridD%data,AVl,rList='lon:lat:mask')
    call mct_avect_gather(AVl,AVg,gsmapD,master_task,mpicom)

    if (my_task == master_task) then
       gsizeD = mct_aVect_lsize(AVg)
       if (gsizeD /= nxgD*nygD) then
          write(logunit,F01) ' ERROR: gsizeD ',gsizeD,nxgD,nygD
          call shr_sys_abort(subname//' ERROR gsizeD')
       endif
       allocate(Xdst(nxgD,nygD),Ydst(nxgD,nygD),Mdst(nxgD,nygD))

       nlon = mct_avect_indexRA(AVg,'lon')
       nlat = mct_avect_indexRA(AVg,'lat')
       nmsk = mct_avect_indexRA(AVg,'mask')

       n = 0
       Mdst = 1
       do j = 1,nygD
          do i = 1,nxgD
             n = n + 1
             Xdst(i,j) = AVg%rAttr(nlon,n)
             Ydst(i,j) = AVg%rAttr(nlat,n)
             if (abs(AVg%rAttr(nmsk,n)) < 0.5_R8) Mdst(i,j) = 0
          enddo
       enddo
    endif

    if (my_task == master_task) call mct_aVect_clean(AVg)
    call mct_aVect_clean(AVl)

    !--- set map ---

    if (my_task == master_task) then
       call shr_map_mapSet(shrmap,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst, &
            trim(name),trim(type),trim(algo),trim(mask),trim(vect))
       deallocate(Xsrc,Ysrc,Msrc)
       deallocate(Xdst,Ydst,Mdst)
    endif

    !--- convert map to sMatP ---

    lstrategy = 'Xonly'
    if (present(strategy)) then
       lstrategy = trim(strategy)
    endif

    if (my_task == master_task) then
       call shr_map_get(shrmap,shr_map_fs_nsrc,nsrc)
       call shr_map_get(shrmap,shr_map_fs_ndst,ndst)
       call shr_map_get(shrmap,shr_map_fs_nwts,nwts)
       allocate(isrc(nwts),idst(nwts),wgts(nwts))
       call shr_map_get(shrmap,isrc,idst,wgts)
       call shr_map_clean(shrmap)

       call mct_sMat_init(sMat0,ndst,nsrc,nwts)

       call mct_sMat_ImpGColI (sMat0,isrc,nwts)
       call mct_sMat_ImpGRowI (sMat0,idst,nwts)
       call mct_sMat_ImpMatrix(sMat0,wgts,nwts)
       deallocate(isrc,idst,wgts)
    endif

    call mct_sMatP_Init(sMatP,sMat0,gsmapS,gsmapD,lstrategy,master_task,mpicom,compid)

    if (my_task == master_task) then
       call mct_sMat_clean(sMat0)
    endif

  end subroutine shr_strdata_mapSet

  !===============================================================================
  subroutine shr_strdata_set_griddata(sdat, fldname, rvalue) 
    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , intent(in)    :: rvalue

    ! local variables
    integer :: kf

    kf = mct_aVect_indexRA(sdat%grid%data, trim(fldname))
    sdat%grid%data%rAttr(kf,:) = rvalue
  end subroutine shr_strdata_set_griddata

  !===============================================================================
  subroutine shr_strdata_get_griddata(sdat, fldname, data) 
    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , intent(out)   :: data(:)

    ! local variables
    integer :: kf

    kf = mct_aVect_indexRA(sdat%grid%data, trim(fldname))
    data(:) = sdat%grid%data%rAttr(kf,:)
  end subroutine shr_strdata_get_griddata

end module dshr_strdata_mod
