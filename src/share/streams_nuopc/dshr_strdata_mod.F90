module dshr_strdata_mod

  ! holds data and methods to advance data models
  ! Obtain the model domain and the stream domain for each stream

  use ESMF

  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use shr_const_mod    , only : shr_const_pi, shr_const_cDay, shr_const_spval
  use shr_cal_mod      , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod      , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod      , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_pio_mod      , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
  use shr_string_mod   , only : shr_string_listgetname, shr_string_listisvalid, shr_string_listgetnum

  use dshr_stream_mod  , only : shr_stream_streamtype, shr_stream_getModelFieldList, shr_stream_getStreamFieldList
  use dshr_stream_mod  , only : shr_stream_taxis_cycle, shr_stream_taxis_extend, shr_stream_findBounds
  use dshr_stream_mod  , only : shr_stream_getCurrFile, shr_stream_setCurrFile, shr_stream_getMeshFilename
  use dshr_stream_mod  , only : shr_stream_init_from_xml, shr_stream_init_from_inline
!  use dshr_stream_mod  , only : shr_stream_restWrite, shr_stream_restRead
  use dshr_stream_mod  , only : shr_stream_getnextfilename, shr_stream_getprevfilename, shr_stream_getData
  use dshr_tinterp_mod , only : shr_tInterp_getCosz, shr_tInterp_getAvgCosz, shr_tInterp_getFactors
  use dshr_methods_mod , only : dshr_fldbun_getfldptr, dshr_fldbun_getfieldN, dshr_fldbun_fldchk, chkerr
  use dshr_methods_mod , only : dshr_fldbun_diagnose, dshr_fldbun_regrid, dshr_field_getfldptr

  use pio              , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
  use pio              , only : pio_openfile, pio_closefile, pio_nowrite
  use pio              , only : pio_seterrorhandling, pio_initdecomp, pio_freedecomp
  use pio              , only : pio_inquire, pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
  use pio              , only : pio_inq_dimlen, pio_inq_vartype, pio_inq_dimname
  use pio              , only : pio_double, pio_real, pio_int, pio_offset_kind
  use pio              , only : pio_read_darray, pio_setframe, pio_fill_double, pio_get_att
  use pio              , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
  use perf_mod         , only : t_startf, t_stopf, t_adj_detailf

  implicit none
  private

  public  :: shr_strdata_type
  public  :: shr_strdata_init_from_xml
  public  :: shr_strdata_init_from_inline
!  public  :: shr_strdata_restRead
!  public  :: shr_strdata_restWrite
  public  :: shr_strdata_setOrbs
  public  :: shr_strdata_advance
  public  :: shr_strdata_get_stream_domain  ! public since needed by dshr_mod
  public  :: shr_strdata_get_stream_pointer ! get a pointer into a stream's fldbun_model field bundle
  public  :: shr_strdata_print
  public  :: shr_strdata_get_stream_count
  public  :: shr_strdata_get_stream_fieldbundle
  private :: shr_strdata_init_model_domain
  private :: shr_strdata_readLBUB

  ! public data members:
  integer                              :: debug    = 0  ! local debug flag
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  integer          ,parameter          :: master_task = 0

  ! note that the fields in fldbun_stream_lb and fldbun_stream_ub contain the the names fldlist_model

  type shr_strdata_perstream
     character(CL)                       :: stream_meshfile                 ! stream mesh file from stream txt file
     type(ESMF_Mesh)                     :: stream_mesh                     ! stream mesh created from stream mesh file
     type(io_desc_t)                     :: stream_pio_iodesc               ! stream pio descriptor
     logical                             :: stream_pio_iodesc_set =.false.  ! true=>pio iodesc has been set
     type(ESMF_RouteHandle)              :: routehandle                     ! stream n -> model mesh mapping
     character(CL), allocatable          :: fldlist_stream(:)               ! names of stream file fields
     character(CL), allocatable          :: fldlist_model(:)                ! names of stream model fields
     integer                             :: stream_lb                       ! index of the Lowerbound (LB) in fldlist_stream
     integer                             :: stream_ub                       ! index of the Upperbound (UB) in fldlist_stream
     type(ESMF_Field)                    :: field_stream                    ! a field on the stream data domain
     type(ESMF_Field)                    :: stream_vector                   ! a vector field on the stream data domain
     type(ESMF_FieldBundle), allocatable :: fldbun_data(:)                  ! stream field bundle interpolated to model grid
     type(ESMF_FieldBundle)              :: fldbun_model                    ! stream n field bundle interpolated to model grid and time
     integer                             :: ucomp = -1                      ! index of vector u in stream
     integer                             :: vcomp = -1                      ! index of vector v in stream
     integer                             :: ymdLB = -1                      ! stream ymd lower bound
     integer                             :: todLB = -1                      ! stream tod lower bound
     integer                             :: ymdUB = -1                      ! stream ymd upper bound
     integer                             :: todUB = -1                      ! stream tod upper bound
     real(r8)                            :: dtmin = 1.0e30_r8
     real(r8)                            :: dtmax = 0.0_r8
     type(ESMF_Field)                    :: field_coszen                    ! needed for coszen time interp
  end type shr_strdata_perstream

  type shr_strdata_type
     type(shr_strdata_perstream), allocatable :: pstrm(:)              ! stream info
     type(shr_stream_streamType), pointer :: stream(:)=> null()        ! stream datatype
     integer                        :: nvectors                        ! number of vectors
     logical                        :: masterproc
     integer                        :: logunit                         ! stdout unit
     integer                        :: io_type                         ! pio info
     integer                        :: io_format                       ! pio info
     integer                        :: modeldt = 0                     ! model dt in seconds
     type(ESMF_Mesh)                :: model_mesh                      ! model mesh
     real(r8), pointer              :: model_lon(:) => null()          ! model longitudes
     real(r8), pointer              :: model_lat(:) => null()          ! model latitudes
     integer                        :: model_nxg                       ! model global domain lon size
     integer                        :: model_nyg                       ! model global domain lat size
     integer                        :: model_nzg                       ! model global domain vertical size
     integer                        :: model_lsize                     ! model local domain size
     integer, pointer               :: model_gindex(:)                 ! model global index spzce
     integer                        :: model_gsize                     ! model global domain size
     type(ESMF_CLock)               :: model_clock                     ! model clock
     character(CL)                  :: model_calendar = shr_cal_noleap ! model calendar for ymd,tod
     integer                        :: ymd, tod                        ! model time
     type(iosystem_desc_t), pointer :: pio_subsystem => null()         ! pio info
     real(r8)                       :: eccen  = SHR_ORB_UNDEF_REAL     ! cosz t-interp info
     real(r8)                       :: mvelpp = SHR_ORB_UNDEF_REAL     ! cosz t-interp info
     real(r8)                       :: lambm0 = SHR_ORB_UNDEF_REAL     ! cosz t-interp info
     real(r8)                       :: obliqr = SHR_ORB_UNDEF_REAL     ! cosz t-interp info
     real(r8), allocatable          :: tavCoszen(:)                    ! cosz t-interp data
  end type shr_strdata_type

  integer          ,parameter :: iotype_std_netcdf = -99 ! non pio option
  real(r8)         ,parameter :: deg2rad = SHR_CONST_PI/180.0_r8
  character(*)     ,parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  integer function shr_strdata_get_stream_count(sdat)
    type(shr_strdata_type)     , intent(in) :: sdat
    shr_strdata_get_stream_count = size(sdat%stream)
  end function shr_strdata_get_stream_count

  !===============================================================================
  type(ESMF_FieldBundle) function shr_strdata_get_stream_fieldbundle(sdat, ns, name)

    ! input/output variables
    type(shr_strdata_type)     , intent(in) :: sdat
    integer                    , intent(in) :: ns ! stream number
    character(len=*)           , intent(in) :: name

    if(trim(name) .eq. 'model') then
       shr_strdata_get_stream_fieldbundle = sdat%pstrm(ns)%fldbun_model
    else if (trim(name) .eq. 'model_lb') then
       shr_strdata_get_stream_fieldbundle = sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_lb)
    else if (trim(name) .eq. 'model_ub') then
       shr_strdata_get_stream_fieldbundle = sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_ub)
    else if (trim(name) .eq. 'stream_lb') then
       call shr_sys_abort("should not be here")
!       shr_strdata_get_stream_fieldbundle = sdat%pstrm(ns)%fldbun_stream(sdat%pstrm(ns)%stream_lb)
    else if (trim(name) .eq. 'stream_ub') then
       call shr_sys_abort("should not be here")
!       shr_strdata_get_stream_fieldbundle = sdat%pstrm(ns)%fldbun_stream(sdat%pstrm(ns)%stream_ub)
    else
       call shr_sys_abort(trim(name)//' is not a recognized stream bundle name')
    endif

  end function shr_strdata_get_stream_fieldbundle

  !===============================================================================
  subroutine shr_strdata_init_from_xml(sdat, xmlfilename, model_mesh, clock, compid, logunit, rc)

    ! input/output variables
    type(shr_strdata_type)     , intent(inout) :: sdat
    character(len=*)           , intent(in)    :: xmlfilename
    type(ESMF_Mesh)            , intent(in)    :: model_mesh
    type(ESMF_Clock)           , intent(in)    :: clock
    integer                    , intent(in)    :: compid
    integer                    , intent(in)    :: logunit
    integer                    , intent(out)   :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: i
    integer       :: localPet
    integer       :: ierr
    character(len=*), parameter  :: subname='(shr_strdata_init_from_xml)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize log unit
    sdat%logunit = logunit

    ! Initialize sdat  pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize sdat streams (read xml file for streams)
    sdat%masterproc = (localPet == master_task)

    call shr_stream_init_from_xml(xmlfilename, sdat%stream, sdat%masterproc, &
         sdat%logunit, compid, rc=rc)

    allocate(sdat%pstrm(shr_strdata_get_stream_count(sdat)))

    ! Initialize sdat model domain
    sdat%model_mesh = model_mesh
    call shr_strdata_init_model_domain(sdat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now finish initializing sdat
    call shr_strdata_init(sdat, clock, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_strdata_init_from_xml

  !===============================================================================
  subroutine shr_strdata_init_from_inline(sdat, my_task, logunit, compid, model_clock, model_mesh,&
       stream_meshfile, stream_filenames, stream_fldlistFile, stream_fldListModel, &
       stream_yearFirst, stream_yearLast, stream_yearAlign, stream_offset, stream_taxmode, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat                   ! stream data type
    integer                , intent(in)    :: my_task                ! my mpi task
    integer                , intent(in)    :: logunit                ! stdout logunit
    integer                , intent(in)    :: compid                 ! component id (needed by pio)
    type(ESMF_Clock)       , intent(in)    :: model_clock            ! model clock
    type(ESMF_Mesh)        , intent(in)    :: model_mesh             ! model mesh
    character(*)           , intent(in)    :: stream_meshFile        ! full pathname to stream mesh file
    character(*)           , intent(in)    :: stream_filenames(:)    ! stream data filenames (full pathnamesa)
    character(*)           , intent(in)    :: stream_fldListFile(:)  ! file field names, colon delim list
    character(*)           , intent(in)    :: stream_fldListModel(:) ! model field names, colon delim list
    integer                , intent(in)    :: stream_yearFirst       ! first year to use
    integer                , intent(in)    :: stream_yearLast        ! last  year to use
    integer                , intent(in)    :: stream_yearAlign       ! align yearFirst with this model year
    integer                , intent(in)    :: stream_offset          ! offset in seconds of stream data
    character(*)           , intent(in)    :: stream_taxMode         ! time axis mode
    integer                , intent(out)   :: rc                     ! error code
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize sdat%logunit and sdat%masterproc
    sdat%logunit = logunit
    sdat%masterproc = (my_task == master_task)

    ! Initialize sdat pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! Initialize sdat%pstrm - ASSUME only 1 stream
    allocate(sdat%pstrm(1))

    ! Initialize sdat model domain
    sdat%model_mesh = model_mesh
    call shr_strdata_init_model_domain(sdat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize sdat stream - ASSUME only 1 stream
    write(6,*)'DEBUG: stream_meshfile = ',stream_meshfile
    call shr_stream_init_from_inline(sdat%stream, stream_meshfile, &
       stream_yearFirst, stream_yearLast, stream_yearAlign, stream_offset, stream_taxmode, &
       stream_fldlistFile, stream_fldListModel, stream_fileNames, logunit)

    ! Now finish initializing sdat
    call shr_strdata_init(sdat, model_clock, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_strdata_init_from_inline

  !===============================================================================
  subroutine shr_strdata_init_model_domain( sdat, rc)

    ! ----------------------------------------------
    ! Initialize sdat model domain info
    ! ----------------------------------------------

    ! input/output variables
    type(shr_strdata_type)     , intent(inout) :: sdat
    integer                    , intent(out)   :: rc

    ! local variables
    integer               :: n,k          ! generic counters
    type(ESMF_DistGrid)   :: distGrid
    integer               :: dimCount
    integer               :: tileCount
    integer, allocatable  :: elementCountPTile(:)
    integer, allocatable  :: indexCountPDE(:,:)
    integer               :: spatialDim         ! number of dimension in mesh
    integer               :: numOwnedElements   ! local size of mesh
    real(r8), allocatable :: ownedElemCoords(:) ! mesh lat and lons
    integer               :: my_task
    integer               :: ierr
    integer               :: rcode
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize sdat%lsize
    call ESMF_MeshGet(sdat%model_mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=sdat%model_lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initialize sdat%model_gindex
    allocate(sdat%model_gindex(sdat%model_lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=sdat%model_gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initialize sdat%model_gsize
    call ESMF_DistGridGet(distGrid, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    sdat%model_gsize = 0
    do n = 1,size(elementCountPTile)
       sdat%model_gsize = sdat%model_gsize + elementCountPTile(n)
    end do
    deallocate(elementCountPTile)

    ! determine sdat%model_lon, sdat%model_lat
    call ESMF_MeshGet(sdat%model_mesh, spatialDim=spatialDim, &
         numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(sdat%model_mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(sdat%model_lon(numOwnedElements))
    allocate(sdat%model_lat(numOwnedElements))
    do n = 1, numOwnedElements
       sdat%model_lon(n) = ownedElemCoords(2*n-1)
       sdat%model_lat(n) = ownedElemCoords(2*n)
    end do

  end subroutine shr_strdata_init_model_domain

  !===============================================================================
  subroutine shr_strdata_init(sdat, model_clock, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(inout), target :: sdat
    type(ESMF_Clock)       , intent(in)    :: model_clock
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Calendar)          :: esmf_calendar   ! esmf calendar
    type(ESMF_CalKind_Flag)      :: esmf_caltype    ! esmf calendar type
    type(ESMF_DistGrid)          :: distgrid
    type(ESMF_RegridMethod_Flag) :: regridmethod
    type(ESMF_PoleMethod_Flag)   :: polemethod
    character(CS)                :: calendar        ! calendar name
    integer                      :: dimcount
    integer, allocatable         :: minIndexPTile(:,:)
    integer, allocatable         :: maxIndexPTile(:,:)
    integer                      :: lnx, lny        ! global mesh dimensions
    integer                      :: ne              ! number of local mesh elements
    integer                      :: ns              ! stream index
    integer                      :: n,m,k           ! generic index
    character(CL)                :: fileName        ! generic file name
    integer                      :: nfiles          ! number of data files for a given stream
    character(CS)                :: uname           ! u vector field name
    character(CS)                :: vname           ! v vector field name
    integer                      :: nu, nv          ! vector indices
    integer                      :: nstream         ! loop stream index
    integer                      :: nvector         ! loop vector index
    integer                      :: nfld            ! loop stream field index
    integer                      :: nflds           ! total number of fields in a given stream
    type(ESMF_Field)             :: lfield          ! temporary
    type(ESMF_Field)             :: lfield_dst      ! temporary
    integer                      :: srcTermProcessing_Value = 0 ! should this be a module variable?
    integer , pointer            :: stream_gindex(:)
    integer                      :: stream_lsize
    character(CS)                :: tmpstr
    integer                      :: ierr
    integer                      :: localpet
    logical                      :: fileExists
    type(ESMF_VM)                :: vm
    logical                      :: masterproc
    integer                      :: nvars
    integer                      :: i
    character(len=*), parameter  :: subname='(shr_strdata_mod:shr_sdat_init)'
    character(*)    , parameter  :: F00 = "('(shr_sdat_init) ',a)"
    character(*)    , parameter  :: F01  = "('(shr_sdat) ',a,2x,i8)"
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_VmGetCurrent(vm, rc=rc)
    call ESMF_VMGet(vm, localpet=localPet, rc=rc)
    masterproc= localPet==master_task

    do ns = 1,shr_strdata_get_stream_count(sdat)
       ! Initialize calendar for stream n
       call ESMF_VMBroadCast(vm, sdat%stream(ns)%calendar, CS, 0, rc=rc)

       ! Create the target stream mesh from the stream mesh file
       ! TODO: add functionality if the stream mesh needs to be created from a grid
       call shr_stream_getMeshFileName (sdat%stream(ns), filename)

       if (masterproc) then
          inquire(file=trim(filename),exist=fileExists)
          if (.not. fileExists) then
             write(sdat%logunit,F00) "ERROR: file does not exist: ", trim(fileName)
             call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
          end if
       endif

       sdat%pstrm(ns)%stream_mesh = ESMF_MeshCreate(trim(filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine field names for stream fields with both stream file names and data model names
       nvars = sdat%stream(ns)%nvars

       allocate(sdat%pstrm(ns)%fldList_model(nvars))
       call shr_stream_getModelFieldList(sdat%stream(ns), sdat%pstrm(ns)%fldlist_model)

       allocate(sdat%pstrm(ns)%fldlist_stream(nvars))
       call shr_stream_getStreamFieldList(sdat%stream(ns), sdat%pstrm(ns)%fldlist_stream)

       ! Create field bundles on model mesh
       if(sdat%stream(ns)%readmode=='single') then
          sdat%pstrm(ns)%stream_lb = 1
          sdat%pstrm(ns)%stream_ub = 2
          allocate(sdat%pstrm(ns)%fldbun_data(2))
       else if(sdat%stream(ns)%readmode=='full_file') then


       endif
       do i=1,size(sdat%pstrm(ns)%fldbun_data)
          sdat%pstrm(ns)%fldbun_data(i) = ESMF_FieldBundleCreate(rc=rc) ! stream mesh
       enddo
       sdat%pstrm(ns)%fldbun_model    = ESMF_FieldBundleCreate(rc=rc) ! time interpolation on model mesh
       do nfld = 1, nvars
          ! create temporary fields on model mesh and add the fields to the field bundle
          do i=1,size(sdat%pstrm(ns)%fldbun_data)
             lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
                  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_data(i), (/lfield/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          enddo
          lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_model   , (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do


       ! Create a field on the model mesh for coszen time interpolation for this stream if needed
       if (trim(sdat%stream(ns)%tinterpalgo) == 'coszen') then
          sdat%pstrm(ns)%field_coszen = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name='tavCosz', &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       endif

       ! ------------------------------------
       ! Create the sdat route handles for mapping the stream -> model
       ! ------------------------------------

       ! create the source and destination fields needed to create the route handles
       ! assume that all fields in a stream share the same mesh and there is only a unique model mesh
       ! can do this outside of a stream loop by just using the first stream index
       sdat%pstrm(ns)%field_stream = ESMF_FieldCreate(sdat%pstrm(ns)%stream_mesh, &
            ESMF_TYPEKIND_r8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_fldbun_getFieldN(sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_lb), 1, lfield_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       !sdat%stream(ns)%mapalgo = "redist"
       if (trim(sdat%stream(ns)%mapalgo) == "bilinear") then
          call ESMF_FieldRegridStore(sdat%pstrm(ns)%field_stream, lfield_dst, &
               routehandle=sdat%pstrm(ns)%routehandle, &
               regridmethod=ESMF_REGRIDMETHOD_BILINEAR,  &
               polemethod=ESMF_POLEMETHOD_ALLAVG, &
               extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
               dstMaskValues = (/0/), &  ! ignore destination points where the mask is 0
               srcMaskValues = (/0/), &  ! ignore source points where the mask is 0
               srcTermProcessing=srcTermProcessing_Value, ignoreDegenerate=.true., rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else if (trim(sdat%stream(ns)%mapalgo) == 'redist') then
          call ESMF_FieldRedistStore(sdat%pstrm(ns)%field_stream, lfield_dst, &
               routehandle=sdat%pstrm(ns)%routehandle, &
               ignoreUnmatchedIndices = .true., rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort('ERROR: only bilinear regrid or redist is supported for now')
       end if

    end do ! end of loop over streams
    !
    ! Check for vector pairs in the stream - both ucomp and vcomp must be in the same stream
    !
    do m = 1,shr_strdata_get_stream_count(sdat)
       ! check that vector field list is a valid colon delimited string
       if(trim(sdat%stream(m)%stream_vectors).eq.'null') cycle

       if (.not. shr_string_listIsValid(sdat%stream(m)%stream_vectors)) then
          write(sdat%logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(sdat%stream(m)%stream_vectors)
          call shr_sys_abort(subname//': vec fldlist invalid:'//trim(sdat%stream(m)%stream_vectors))
       endif

       ! check that only 2 fields are contained for any vector pairing
       if (shr_string_listGetNum(sdat%stream(m)%stream_vectors) /= 2) then
          write(sdat%logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(sdat%stream(m)%stream_vectors)
          call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(sdat%stream(m)%stream_vectors))
       endif

       sdat%pstrm(m)%stream_vector = ESMF_FieldCreate(sdat%pstrm(m)%stream_mesh, ESMF_TYPEKIND_r8, name='stream_vector', &
            ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! initialize sdat model clock and calendar
    sdat%model_clock = model_clock
    call ESMF_ClockGet(sdat%model_clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       sdat%model_calendar = trim(shr_cal_noleap)
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       sdat%model_calendar = trim(shr_cal_gregorian)
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(calendar))
    end if

    ! print sdat output
    if (masterproc) then
       call shr_strdata_print(sdat,'sdat data ')
       write(sdat%logunit,*) ' successfully initialized sdat'
    endif

  end subroutine shr_strdata_init

  !===============================================================================
  subroutine shr_strdata_get_stream_domain(sdat, stream_index, fldname, flddata, rc)

    ! Obtain the data for fldname from the stream domain data

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    integer                , intent(in)    :: stream_index
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , pointer       :: flddata(:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    type(var_desc_t)        :: varid
    type(file_desc_t)       :: pioid
    integer                 :: rcode
    character(CL)           :: filename
    type(io_desc_t)         :: pio_iodesc
    real(r4), allocatable   :: data_real(:)
    real(r8), allocatable   :: data_double(:)
    integer                 :: pio_iovartype
    integer                 :: lsize
    character(*), parameter :: subname = '(shr_strdata_set_stream_domain) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine the file to open
    if (sdat%masterproc) then
       call shr_stream_getData(sdat%stream(stream_index), 1, filename)
    end if
    call ESMF_VMBroadCast(vm, filename, CL, 0, rc=rc)


    ! Open the file
    rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(filename), pio_nowrite)

    ! Create the pio iodesc for fldname
    call shr_strdata_set_stream_iodesc(sdat, sdat%pio_subsystem, pioid, &
         trim(fldname), sdat%pstrm(stream_index)%stream_mesh, pio_iodesc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now read in the data for fldname
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    lsize = size(flddata)
    rcode = pio_inq_varid(pioid, trim(fldname), varid)
    rcode = pio_inq_vartype(pioid, varid, pio_iovartype)
    if (pio_iovartype == PIO_REAL) then
       allocate(data_real(lsize))
       call pio_read_darray(pioid, varid, pio_iodesc, data_real, rcode)
       flddata(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    else if (pio_iovartype == PIO_DOUBLE) then
       allocate(data_double(lsize))
       call pio_read_darray(pioid, varid, pio_iodesc, data_double, rcode)
       flddata(:) = data_double(:)
       deallocate(data_double)
    else
       call shr_sys_abort(subName//"ERROR: only real and double types are supported for stream domain read")
    end if

    ! Free the memory associate with the iodesc and close the file
    call pio_freedecomp(pioid, pio_iodesc)
    call pio_closefile(pioid)

  end subroutine shr_strdata_get_stream_domain

  !===============================================================================
  subroutine shr_strdata_advance(sdat, ymd, tod, logunit, istr, timers, rc)

    ! -------------------------------------------------------
    ! Mismatching calendars: 4 cases
    ! (0) The stream calendar and model calendar are identical
    ! (1) The stream is a no leap calendar and the model is gregorian
    ! (2) The stream is a gregorian calendar and the model is a noleap calendar
    ! (3) The calendars mismatch and none of the above
    ! -------------------------------------------------------
    !
    ! ymdmod and todmod are the ymd and tod to time interpolate to.
    ! Generally, these are just the model date and time.  Also, always
    ! use the stream calendar for time interpolation for reasons
    ! described below.  When there is a calendar mismatch, support Feb
    ! 29 in a special way as needed to get reasonable values.  Note
    ! that when Feb 29 needs to be treated specially, a discontinuity
    ! will be introduced.  The size of that discontinuity will depend
    ! on the time series input data.
    !
    ! (0) The stream calendar and model calendar are identical:
    ! Proceed in the standard way.
    !
    ! (1) The stream is a no leap calendar and the model is gregorian:
    ! Time interpolate on the noleap calendar.  If the model date is Feb 29,
    ! compute stream data for Feb 28 by setting ymdmod and todmod to Feb 28.
    ! This results in duplicate stream data on Feb 28 and Feb 29 and a
    ! discontinuity at the start of Feb 29.  This could potentially be fixed
    ! by using the gregorian calendar for time interpolation when the input data
    ! is relatively infrequent (say greater than daily) with the following concerns.
    !   - The forcing will not be reproduced identically on the same day with
    !     with climatological inputs data
    !   - Input data with variable input frequency might behave funny
    !   - An arbitrary discontinuity will be introduced in the time
    !     interpolation method based upon the logic chosen to transition
    !     from reproducing Feb 28 on Feb 29 and interpolating to Feb 29.
    !   - The time gradient of data will change by adding a day arbitrarily.
    !
    ! (2) The stream is a gregorian calendar and the model is a noleap calendar:
    ! Then just time interpolate on the gregorian calendar. This causes Feb 29
    ! stream data to be skipped and lead to a discontinuity at the start of March 1.
    !
    ! (3) If the calendars mismatch and neither of the three cases above
    ! are recognized, then abort.
    ! -------------------------------------------------------

    ! input/output variables
    type(shr_strdata_type) ,intent(inout)       :: sdat
    integer                ,intent(in)          :: ymd    ! current model date
    integer                ,intent(in)          :: tod    ! current model date
    integer                ,intent(in)          :: logunit
    character(len=*)       ,intent(in)          :: istr
    logical                ,intent(in),optional :: timers
    integer                ,intent(out)         :: rc

    ! local variables
    integer                             :: ns               !stream index
    integer                             :: nf               ! field index
    integer                             :: m                ! vector index
    integer                             :: n,i              ! generic indices
    integer                             :: ierr
    integer                             :: nu,nv
    integer                             :: lsize
    logical , allocatable               :: newData(:)
    integer , allocatable               :: ymdmod(:)        ! modified model dates to handle Feb 29
    real(r8), allocatable               :: coszen(:)        ! cosine of zenith angle
    integer                             :: todmod           ! modified model dates to handle Feb 29
    character(len=32)                   :: lstr             ! local string
    logical                             :: ltimers          ! local logical for timers
    real(r8)                            :: flb,fub          ! factor for lb and ub
    real(r8) ,pointer                   :: dataptr(:)       ! pointer into field bundle
    real(r8) ,pointer                   :: dataptr_lb(:)    ! pointer into field bundle
    real(r8) ,pointer                   :: dataptr_ub(:)    ! pointer into field bundle
    real(r8), pointer                   :: nu_coords(:)     ! allocatable local element mesh lat and lons
    real(r8), pointer                   :: nv_coords(:)     ! allocatable local element mesh lat and lons
    real(r8), pointer                   :: data2d_src(:,:)  ! pointer into field bundle
    real(r8), pointer                   :: data2d_dst(:,:)  ! pointer into field bundle
    real(r8), pointer                   :: data_u_src(:)    ! pointer into field bundle
    real(r8), pointer                   :: data_v_src(:)    ! pointer into field bundle
    real(r8), pointer                   :: data_u_dst(:)    ! pointer into field bundle
    real(r8), pointer                   :: data_v_dst(:)    ! pointer into field bundle
    type(ESMF_Time)                     :: timeLB, timeUB   ! lb and ub times
    type(ESMF_TimeInterval)             :: timeint          ! delta time
    integer                             :: dday             ! delta days
    real(r8)                            :: dtime            ! delta time
    integer                             :: year,month,day   ! date year month day
    integer                             :: spatialDim       ! spatial dimension of mesh
    integer                             :: numOwnedElements ! local size of mesh
    character(CS)                       :: uname            ! u vector field name
    character(CS)                       :: vname            ! v vector field name
    type(ESMF_Field)                    :: field_src
    type(ESMF_Field)                    :: field_dst
    real(r8)                            :: lon, lat
    real(r8)                            :: sinlon, sinlat
    real(r8)                            :: coslon, coslat
    real(r8)                            :: ux, uy
    logical                             :: checkflag = .false.
    integer                             :: npes
    integer                             :: my_task
    real(r8)         ,parameter         :: solZenMin = 0.001_r8 ! minimum solar zenith angle
    integer          ,parameter         :: tadj = 2
    character(len=*) ,parameter         :: timname = "_strd_adv"
    character(*)     ,parameter         :: subname = "(shr_strdata_advance) "
    character(*)     ,parameter         :: F00  = "('(shr_strdata_advance) ',a)"
    character(*)     ,parameter         :: F01  = "('(shr_strdata_advance) ',a,a,i4,2(f10.5,2x))"
    real(r8), pointer :: dataptr_temp1(:)
    real(r8), pointer :: dataptr_temp2(:)
    integer :: nstreams
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(dataptr)
    nullify(dataptr_ub)
    nullify(dataptr_lb)
    nullify(nu_coords)
    nullify(nv_coords)
    nullify(data2d_src)
    nullify(data2d_dst)
    nullify(data_u_src)
    nullify(data_v_src)
    nullify(data_u_dst)
    nullify(data_v_dst)
    nstreams = shr_strdata_get_stream_count(sdat)
    if (nstreams < 1) return ! TODO: is this needed

    lstr = trim(istr)

    ltimers = .true.
    if (present(timers)) then
       ltimers = timers
    endif

    if (.not.ltimers) call t_adj_detailf(tadj)

!    call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
    call t_startf(trim(lstr)//trim(timname)//'_total')


    sdat%ymd = ymd
    sdat%tod = tod

    if (nstreams > 0) then
       allocate(newData(nstreams))
       allocate(ymdmod(nstreams))

       do ns = 1,nstreams
          ! ---------------------------------------------------------
          ! Consistency checks
          ! ---------------------------------------------------------

          ! case(0)
          ymdmod(ns) = ymd
          todmod    = tod
          if (trim(sdat%model_calendar) /= trim(sdat%stream(ns)%calendar)) then
             if (( trim(sdat%model_calendar) == trim(shr_cal_gregorian)) .and. &
                  (trim(sdat%stream(ns)%calendar) == trim(shr_cal_noleap))) then
                ! case (1), set feb 29 = feb 28
                call shr_cal_date2ymd (ymd,year,month,day)
                if (month == 2 .and. day == 29) then
                   call shr_cal_ymd2date(year,2,28,ymdmod(ns))
                endif
             else if ((trim(sdat%model_calendar) == trim(shr_cal_noleap)) .and. &
                      (trim(sdat%stream(ns)%calendar) == trim(shr_cal_gregorian))) then
                ! case (2), feb 29 input data will be skipped automatically
             else
                ! case (3), abort
                write(logunit,*) trim(subname),' ERROR: mismatch calendar ', &
                     trim(sdat%model_calendar),':',trim(sdat%stream(ns)%calendar)
                call shr_sys_abort(trim(subname)//' ERROR: mismatch calendar ')
             endif
          endif

          ! ---------------------------------------------------------
          ! Determine if new data is read in - if so then copy
          ! fldbun_stream_ub to fldbun_stream_lb and read in new fldbun_stream_ub data
          ! ---------------------------------------------------------

!          call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
          call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

          select case(sdat%stream(ns)%readmode)
          case ('single')
             call shr_strdata_readLBUB(sdat, ns, &
                  ymdmod(ns), todmod, &
                  newData(ns), trim(lstr)//'_readLBUB', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          case ('full_file')

             ! TODO: need to put in capability to read all stream data at once
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(sdat%stream(ns)%readmode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(sdat%stream(ns)%readmode))
          end select

          if (debug > 0 .and. sdat%masterproc) then
             write(sdat%logunit,*) trim(subname),' newData flag = ',ns,newData(ns)
             write(sdat%logunit,*) trim(subname),' LB ymd,tod = ',ns,sdat%pstrm(ns)%ymdLB,sdat%pstrm(ns)%todLB
             write(sdat%logunit,*) trim(subname),' UB ymd,tod = ',ns,sdat%pstrm(ns)%ymdUB,sdat%pstrm(ns)%todUB
          endif

          ! ---------------------------------------------------------
          ! If new data is read in:
          ! ---------------------------------------------------------

          if (newData(ns)) then
             if(sdat%pstrm(ns)%ymdLB <= 0 ) then
                call shr_sys_abort('time out of bounds')
             endif
             ! Reset time bounds if newdata read in
             call shr_cal_date2ymd(sdat%pstrm(ns)%ymdLB,year,month,day)
             call shr_cal_timeSet(timeLB,sdat%pstrm(ns)%ymdLB,0,sdat%stream(ns)%calendar,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call shr_cal_timeSet(timeUB,sdat%pstrm(ns)%ymdUB,0,sdat%stream(ns)%calendar,rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             timeint = timeUB-timeLB
             call ESMF_TimeIntervalGet(timeint,StartTimeIn=timeLB,d=dday)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             dtime = abs(real(dday,r8) + real(sdat%pstrm(ns)%todUB-sdat%pstrm(ns)%todLB,r8)/shr_const_cDay)

             sdat%pstrm(ns)%dtmin = min(sdat%pstrm(ns)%dtmin,dtime)
             sdat%pstrm(ns)%dtmax = max(sdat%pstrm(ns)%dtmax,dtime)
             if ((sdat%pstrm(ns)%dtmax/sdat%pstrm(ns)%dtmin) > sdat%stream(ns)%dtlimit) then
                write(logunit,*) trim(subname),' ERROR: for stream ',n
                write(logunit,*) trim(subName),' ERROR: dt limit1 ',&
                     sdat%pstrm(ns)%dtmax, sdat%pstrm(ns)%dtmin, sdat%stream(ns)%dtlimit
                write(logunit,*) trim(subName),' ERROR: dt limit2 ',&
                     sdat%pstrm(ns)%ymdLB, sdat%pstrm(ns)%todLB, sdat%pstrm(ns)%ymdUB, sdat%pstrm(ns)%todUB
                call shr_sys_abort(trim(subName)//' ERROR dt limit for stream')
             endif
          endif

          call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')

       enddo ! end of loop over streams

       ! ---------------------------------------------------------
       ! Do time interpolation to create fldbun_model
       ! ---------------------------------------------------------

       do ns = 1,nstreams

          if (trim(sdat%stream(ns)%tinterpalgo) == 'coszen') then

             ! ------------------------------------------
             ! time interpolation method is coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_coszen')
             allocate(coszen(sdat%model_lsize))
             write(6,*)'DEBUG: model_size = ',sdat%model_lsize

             ! get coszen
             call t_startf(trim(lstr)//trim(timname)//'_coszenC')
             call shr_tInterp_getCosz(coszen, sdat%model_lon, sdat%model_lat, ymdmod(ns), todmod, &
                  sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%stream(ns)%calendar)
             call t_stopf(trim(lstr)//trim(timname)//'_coszenC')

             ! get avg cosz factor
             if (newdata(ns)) then
                ! compute a new avg cosz
                call t_startf(trim(lstr)//trim(timname)//'_coszenN')
                if (.not. allocated(sdat%tavCoszen)) then
                   allocate(sdat%tavCoszen(sdat%model_lsize))
                end if
                call shr_tInterp_getAvgCosz(sdat%tavCoszen, sdat%model_lon, sdat%model_lat,  &
                     sdat%pstrm(ns)%ymdLB, sdat%pstrm(ns)%todLB,  sdat%pstrm(ns)%ymdUB, sdat%pstrm(ns)%todUB,  &
                     sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%modeldt, &
                     sdat%stream(ns)%calendar, rc=rc)
                call t_stopf(trim(lstr)//trim(timname)//'_coszenN')
             endif

             ! compute time interperpolate value - LB data normalized with this factor: cosz/tavCosz
             do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model   , sdat%pstrm(ns)%fldlist_model(nf), dataptr   , rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_lb), sdat%pstrm(ns)%fldlist_model(nf), dataptr_lb, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                do i = 1,size(dataptr)
                   if (coszen(i) > solZenMin) then
                      dataptr(i) = dataptr_lb(i)*coszen(i)/sdat%tavCoszen(i)
                   else
                      dataptr(i) = 0._r8
                   endif
                end do
             end do

             deallocate(coszen)
             call t_stopf(trim(lstr)//trim(timname)//'_coszen')

          elseif (trim(sdat%stream(ns)%tinterpalgo) /= trim(shr_strdata_nullstr)) then

             ! ------------------------------------------
             ! time interpolation method is not coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_tint')
             call shr_tInterp_getFactors(sdat%pstrm(ns)%ymdlb, sdat%pstrm(ns)%todlb, sdat%pstrm(ns)%ymdub, sdat%pstrm(ns)%todub, &
                  ymdmod(ns), todmod, flb, fub, calendar=sdat%stream(ns)%calendar, logunit=sdat%logunit, &
                  algo=trim(sdat%stream(ns)%tinterpalgo), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (debug > 0 .and. sdat%masterproc) then
                write(sdat%logunit,F01) trim(subname),' interp = ',ns,flb,fub
             endif

             do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model   , sdat%pstrm(ns)%fldlist_model(nf), dataptr   , rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_lb), sdat%pstrm(ns)%fldlist_model(nf), dataptr_lb, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_ub), sdat%pstrm(ns)%fldlist_model(nf), dataptr_ub, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                dataptr(:) = dataptr_lb(:) * flb + dataptr_ub(:) * fub

             end do
             call t_stopf(trim(lstr)//trim(timname)//'_tint')

          else

             ! ------------------------------------------
             ! zero out stream data for this field
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_zero')
             do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model, sdat%pstrm(ns)%fldlist_model(nf), dataptr, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                dataptr(:) = 0._r8
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_zero')

          endif

       end do  ! loop over ns (number of streams)

       deallocate(newData)
       deallocate(ymdmod)

    endif  ! nstreams > 0

    call t_stopf(trim(lstr)//trim(timname)//'_total')
    if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine shr_strdata_advance

  !===============================================================================

  subroutine shr_strdata_setOrbs(sdat,eccen,mvelpp,lambm0,obliqr,modeldt)

    type(shr_strdata_type),intent(inout) :: sdat
    real(r8),intent(in) :: eccen
    real(r8),intent(in) :: mvelpp
    real(r8),intent(in) :: lambm0
    real(r8),intent(in) :: obliqr
    integer,intent(in) :: modeldt

    ! local variables
    character(len=*),parameter :: subname = "(shr_strdata_setOrbs) "
    !-------------------------------------------------------------------------------

    sdat%eccen   = eccen
    sdat%mvelpp  = mvelpp
    sdat%lambm0  = lambm0
    sdat%obliqr  = obliqr
    sdat%modeldt = modeldt

  end subroutine shr_strdata_setOrbs

  !===============================================================================
  subroutine shr_strdata_print(sdat, name)

    !  Print strdata common to all data models

    ! input/output parameters:
    type(shr_strdata_type) , intent(in) :: sdat  ! strdata data data-type
    character(len=*)       , intent(in) :: name

    ! local variables
    integer  :: ns,n
    character(*),parameter :: subName = "(shr_strdata_print) "
    character(*),parameter ::   F00 = "('(shr_strdata_print) ',8a)"
    character(*),parameter ::   F01 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter ::   F02 = "('(shr_strdata_print) ',a,es13.6)"
    character(*),parameter ::   F04 = "('(shr_strdata_print) ',a,i2,a,a)"
    character(*),parameter ::   F07 = "('(shr_strdata_print) ',a,i2,a,es13.6)"
    character(*),parameter ::   F90 = "('(shr_strdata_print) ',58('-'))"
    !-------------------------------------------------------------------------------

    write(sdat%logunit,*)
    write(sdat%logunit,F90)
    write(sdat%logunit,F00) "name        = ",trim(name)
    write(sdat%logunit,F00) "calendar    = ",trim(sdat%model_calendar)
    write(sdat%logunit,F02) "eccen       = ",sdat%eccen
    write(sdat%logunit,F02) "mvelpp      = ",sdat%mvelpp
    write(sdat%logunit,F02) "lambm0      = ",sdat%lambm0
    write(sdat%logunit,F02) "obliqr      = ",sdat%obliqr
    write(sdat%logunit,F01) "pio_iotype  = ",sdat%io_type

    write(sdat%logunit,F01) "nstreams    = ",shr_strdata_get_stream_count(sdat)
    do ns = 1, shr_strdata_get_stream_count(sdat)
       write(sdat%logunit,F04) "  taxMode (",ns,") = ",trim(sdat%stream(ns)%taxmode)
       write(sdat%logunit,F07) "  dtlimit (",ns,") = ",sdat%stream(ns)%dtlimit
       write(sdat%logunit,F04) "  mapalgo (",ns,") = ",trim(sdat%stream(ns)%mapalgo)
       write(sdat%logunit,F04) "  tintalgo(",ns,") = ",trim(sdat%stream(ns)%tinterpalgo)
       write(sdat%logunit,F04) "  readmode(",ns,") = ",trim(sdat%stream(ns)%readmode)
       write(sdat%logunit,F01) " "
    end do

    write(sdat%logunit,F01) "nvectors    = ",sdat%nvectors
    do n=1, sdat%nvectors
       write(sdat%logunit,F04) "  vectors (",n,") = ",trim(sdat%stream(n)%stream_vectors)
    end do
    write(sdat%logunit,F90)

  end subroutine shr_strdata_print

  !===============================================================================

  subroutine shr_strdata_readLBUB(sdat, ns, mDate, mSec, newData, istr, rc)

    !-------------------------------------------------------------------------
    ! Read LB and UB of stream data
    !-------------------------------------------------------------------------

    ! input/output variables
    type(shr_strdata_type) , intent(inout), target :: sdat  ! strdata data data-type
    integer                , intent(in) :: ns    ! stream number
    integer                ,intent(in)    :: mDate  ,mSec
    logical                ,intent(out)   :: newData
    character(len=*)       ,intent(in)    :: istr
    integer                ,intent(out)   :: rc

    ! local variables
    type(shr_stream_streamType),pointer :: stream
    type(ESMF_Mesh), pointer            :: stream_mesh
    type(ESMF_FieldBundle)        ,pointer :: fldbun_stream_lb
    type(ESMF_FieldBundle)        ,pointer :: fldbun_stream_ub
    type(ESMF_VM)                       :: vm
    integer                             :: nf
    integer                             :: rCode      ! return code
    logical                             :: fileexists
    integer                             :: ivals(6)   ! bcast buffer
    integer                             :: oDateLB,oSecLB,dDateLB
    integer                             :: oDateUB,oSecUB,dDateUB
    real(r8)                            :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer                             :: n_lb, n_ub
    integer                             :: i
    character(CL)                       :: filename_lb
    character(CL)                       :: filename_ub
    character(CL)                       :: filename_next
    character(CL)                       :: filename_prev
    real(r8), pointer                   :: dataptr_lb(:)
    real(r8), pointer                   :: dataptr_ub(:)
    character(*), parameter             :: subname = '(shr_strdata_readLBUB) '
    character(*), parameter             :: F00   = "('(shr_strdata_readLBUB) ',8a)"
    character(*), parameter             :: F01   = "('(shr_strdata_readLBUB) ',a,5i8)"
    !-------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    stream => sdat%stream(ns)
    stream_mesh => sdat%pstrm(ns)%stream_mesh

    call t_startf(trim(istr)//'_setup')
    ! allocate streamdat instance on all tasks
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    newData = .false.
    n_lb = -1
    n_ub = -1
    filename_lb = 'undefinedlb'
    filename_ub = 'undefinedub'

    oDateLB = sdat%pstrm(ns)%ymdLB
    oSecLB  = sdat%pstrm(ns)%todLB
    oDateUB = sdat%pstrm(ns)%ymdUB
    oSecUB  = sdat%pstrm(ns)%todUB

    rDateM  = real(mDate  ,r8) + real(mSec  ,r8)/shr_const_cday
    rDateLB = real(sdat%pstrm(ns)%ymdLB,r8) + real(sdat%pstrm(ns)%todLB,r8)/shr_const_cday
    rDateUB = real(sdat%pstrm(ns)%ymdUB,r8) + real(sdat%pstrm(ns)%todUB,r8)/shr_const_cday
    call t_stopf(trim(istr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(istr)//'_fbound')
       call shr_stream_findBounds(stream, mDate, mSec,  &
            sdat%pstrm(ns)%ymdLB, dDateLB, sdat%pstrm(ns)%todLB, n_lb, filename_lb, &
            sdat%pstrm(ns)%ymdUB, dDateUB, sdat%pstrm(ns)%todUB, n_ub, filename_ub)
       call t_stopf(trim(istr)//'_fbound')
    endif
    if (sdat%pstrm(ns)%ymdLB /= oDateLB .or. sdat%pstrm(ns)%todLB /= oSecLB) then
       newdata = .true.
       if (sdat%pstrm(ns)%ymdLB == oDateUB .and. sdat%pstrm(ns)%todLB == oSecUB) then
          ! copy fldbun_stream_ub to fldbun_stream_lb
          i = sdat%pstrm(ns)%stream_ub
          sdat%pstrm(ns)%stream_ub = sdat%pstrm(ns)%stream_lb
          sdat%pstrm(ns)%stream_lb = i
!          i = sdat%pstrm(ns)%n_ub
!          sdat%pstrm(ns)%n_lb = sdat%pstrm(ns)%n_ub
!          sdat%pstrm(ns)%n_ub = i
       else
          ! read lower bound of data
          call shr_strdata_readstrm(sdat, ns, stream, stream_mesh, &
               sdat%pstrm(ns)%fldlist_stream, sdat%pstrm(ns)%fldlist_model, &
               sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_lb), &
               sdat%pio_subsystem, sdat%io_type, sdat%pstrm(ns)%stream_pio_iodesc_set, sdat%pstrm(ns)%stream_pio_iodesc, &
               filename_lb, n_lb, istr=trim(istr)//'_LB', boundstr='lb', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    if (sdat%pstrm(ns)%ymdUB /= oDateUB .or. sdat%pstrm(ns)%todUB /= oSecUB) then
       newdata = .true.
       call shr_strdata_readstrm(sdat, ns, stream, stream_mesh, &
            sdat%pstrm(ns)%fldlist_stream, sdat%pstrm(ns)%fldlist_model, &
            sdat%pstrm(ns)%fldbun_data(sdat%pstrm(ns)%stream_ub), &
            sdat%pio_subsystem, sdat%io_type, sdat%pstrm(ns)%stream_pio_iodesc_set, sdat%pstrm(ns)%stream_pio_iodesc, &
            filename_ub, n_ub, istr=trim(istr)//'_UB', boundstr='ub', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    ! determine previous & next data files in list of files

    call t_startf(trim(istr)//'_filemgt')
    if (sdat%masterproc .and. newdata) then
       call shr_stream_getPrevFileName(stream, filename_lb, filename_prev)
       call shr_stream_getNextFileName(stream, filename_ub, filename_next)
       inquire(file=trim(filename_next),exist=fileExists)
       if ( trim(filename_next) == "unknown" .or. fileExists) then
          ! do nothing
       end if
    endif
    call t_stopf(trim(istr)//'_filemgt')

  end subroutine shr_strdata_readLBUB

  !===============================================================================
  subroutine shr_strdata_readstrm(sdat, ns, stream, stream_mesh, &
       fldlist_stream, fldlist_model, fldbun_model, &
       pio_subsystem, pio_iotype, pio_iodesc_set, pio_iodesc, &
       filename, nt, istr, boundstr, rc)
    use shr_const_mod         , only : r8fill => SHR_CONST_SPVAL
    use shr_infnan_mod        , only : shr_infnan_isnan

    ! Read the stream data and initialize the strea pio_iodesc the first time
    ! the stream is read

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat  ! strdata data data-type
    integer,                 intent(in) :: ns    ! current stream index
    type(shr_stream_streamType) , intent(inout)         :: stream
    type(ESMF_Mesh)             , intent(in)            :: stream_mesh
    character(len=*)            , intent(in)            :: fldlist_stream(:)
    character(len=*)            , intent(in)            :: fldlist_model(:)
    type(ESMF_FieldBundle)      , intent(inout)         :: fldbun_model
    type(iosystem_desc_t)       , intent(inout), target :: pio_subsystem
    integer                     , intent(in)            :: pio_iotype
    logical                     , intent(inout)         :: pio_iodesc_set
    type(io_desc_t)             , intent(inout)         :: pio_iodesc
    character(len=*)            , intent(in)            :: filename
    integer                     , intent(in)            :: nt
    character(len=*)            , intent(in)            :: istr
    character(len=*)            , intent(in)            :: boundstr
    integer                     , intent(out)           :: rc

    ! local variables
    type(ESMF_Field)              :: field_dst, vector_dst
    character(CL)                 :: currfile
    logical                       :: fileexists
    logical                       :: fileopen
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer                       :: nf
    integer                       :: rCode
    real(r4), allocatable         :: data_real(:)
    real(r4)                      :: fillvalue_r4
    real(r8)                      :: fillvalue_r8
    logical                       :: handlefill = .false.
    integer                       :: old_error_handle
    real(r8), pointer             :: dataptr(:)=>null()
    real(r8), pointer             :: dataptr2d_src(:,:) => null(), dataptr2d_dst(:,:) => null()
    integer                       :: lsize, n
    integer                       :: spatialDim, numOwnedElements
    integer                       :: pio_iovartype
    real(r8), pointer             :: nv_coords(:), nu_coords(:), data_u_dst(:), data_v_dst(:)
    real(r8)                      :: lat, lon, sinlat, sinlon, coslat, coslon
    character(CS)                 :: uname, vname
    character(*), parameter       :: subname = '(shr_strdata_readstrm) '
    character(*), parameter       :: F00   = "('(shr_strdata_readstrm) ',8a)"
    character(*), parameter       :: F02   = "('(shr_strdata_readstrm) ',2a,i8)"
    integer :: i
    logical :: checkflag = .false.
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set up file to read from
!    call t_barrierf(trim(istr)//'_BARRIER', mpicom)
    if (sdat%masterproc) then
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(sdat%logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif

    ! Get current file and determine if it is open
    call shr_stream_getCurrFile(stream, fileopen=fileopen, currfile=currfile, currpioid=pioid)

    if (fileopen .and. currfile==filename) then
       ! don't reopen file, all good
    else
       ! otherwise close the old file if open and open new file
       if (fileopen) then
          if (sdat%masterproc) write(sdat%logunit,F00) 'close  : ',trim(currfile)
          call pio_closefile(pioid)
       endif
       if (sdat%masterproc) write(sdat%logunit,F00) 'opening   : ',trim(filename)
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call shr_stream_setCurrFile(stream, fileopen=.true., currfile=trim(filename), currpioid=pioid)
    endif

    ! ******************************************************************************
    ! Determine the pio io descriptor for the stream from the first data field in the stream
    ! ******************************************************************************

    if (.not. pio_iodesc_set) then
       call shr_strdata_set_stream_iodesc(sdat, pio_subsystem, pioid, &
            trim(fldlist_stream(1)), stream_mesh, pio_iodesc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       pio_iodesc_set = .true.
    end if

    ! ******************************************************************************
    ! Read in the stream data for field names in fldname_stream_input - but fill in
    ! the data for fldbun_stream with the field names fldname_stream_model
    ! ******************************************************************************

    call t_startf(trim(istr)//'_readpio')
    if (sdat%masterproc) then
       write(sdat%logunit,F02) 'file ' // trim(boundstr) //': ',trim(filename), nt
    endif

    call dshr_field_getfldptr(sdat%pstrm(ns)%field_stream, fldptr1=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if(sdat%pstrm(ns)%ucomp > 0 .and. sdat%pstrm(ns)%vcomp > 0) then
       call shr_string_listGetName(stream%stream_vectors,1,uname)
       call shr_string_listGetName(stream%stream_vectors,2,vname)
       call dshr_field_getfldptr(sdat%pstrm(ns)%stream_vector, fldptr2=dataptr2d_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    lsize = size(dataptr)
    do nf = 1,size(fldlist_stream)
       call dshr_fldbun_getfieldN(fldbun_model, nf, field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rcode = pio_inq_varid(pioid, trim(fldlist_stream(nf)), varid)
       ! determine type of the variable
       rcode = pio_inq_vartype(pioid, varid, pio_iovartype)

       if (pio_iovartype == PIO_REAL .and. .not. allocated(data_real)) then
          allocate(data_real(lsize))
       endif

       handlefill = .false.
       call PIO_seterrorhandling(pioid, PIO_BCAST_ERROR, old_error_handle)
       if (pio_iovartype == PIO_REAL) then
          rcode = pio_get_att(pioid, varid, "_FillValue", fillvalue_r4)
       else if (pio_iovartype == PIO_DOUBLE) then
          rcode = pio_get_att(pioid, varid, "_FillValue", fillvalue_r8)
       endif
       if(rcode == PIO_NOERR) handlefill=.true.
       call PIO_seterrorhandling(pioid, old_error_handle)

       if (debug>0 .and. sdat%masterproc)  then
          write(sdat%logunit,F02)' reading '//trim(fldlist_stream(nf))//' into '//trim(fldlist_model(nf)),' at time index: ',nt
       end if

       call pio_setframe(pioid, varid, int(nt,kind=Pio_Offset_Kind))

       if (pio_iovartype == PIO_REAL) then
          call pio_read_darray(pioid, varid, pio_iodesc, data_real, rcode)
          if ( rcode /= PIO_NOERR ) then
             call shr_sys_abort(' ERROR: reading in variable: '// trim(fldlist_stream(nf)))
          end if
          if(handlefill) then
             do n=1,size(dataptr)
                if(.not. shr_infnan_isnan(data_real(n)) .and. data_real(n) .ne. fillvalue_r4) then
                   dataptr(n) = real(data_real(n), kind=r8)
                else
                   dataptr(n) = r8fill
                endif
             enddo
          else
             dataptr(:) = real(data_real(:),kind=r8)
          endif
       else if (pio_iovartype == PIO_DOUBLE) then
          call pio_read_darray(pioid, varid, pio_iodesc, dataptr, rcode)
          if ( rcode /= PIO_NOERR ) then
             call shr_sys_abort(' ERROR: reading in variable: '// trim(fldlist_stream(nf)))
          end if
          if(handlefill) then
             do n=1,size(dataptr)
                if(.not. shr_infnan_isnan(dataptr(n)) .and. dataptr(n).eq.fillvalue_r8) then
                   dataptr(n) = r8fill
                endif
             enddo
          endif
       else
          call shr_sys_abort(subName//"ERROR: only real and double types are subborted for stream read")
       end if
       if(associated(dataptr2d_src) .and. trim(fldlist_model(nf)) .eq. uname) then
          ! save in dataptr2d_src
          dataptr2d_src(1,:) = dataptr(:)
       elseif(associated(dataptr2d_src) .and. trim(fldlist_model(nf)) .eq. vname) then
          dataptr2d_src(2,:) = dataptr(:)
       else
          call ESMF_FieldRegrid(sdat%pstrm(ns)%field_stream, field_dst, routehandle=sdat%pstrm(ns)%routehandle, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=.false., zeroregion=ESMF_REGION_TOTAL, rc=rc)
       endif
    enddo

    ! Both components of a vector stream must be in the same input stream file
    if(associated(dataptr2d_src)) then
       ! get lon and lat of stream u and v fields
       allocate(dataptr(lsize))

       call ESMF_MeshGet(sdat%pstrm(ns)%stream_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(nu_coords(spatialDim*numOwnedElements))
       call ESMF_MeshGet(sdat%pstrm(ns)%stream_mesh, ownedElemCoords=nu_coords)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(nv_coords(spatialDim*numOwnedElements))
       call ESMF_MeshGet(sdat%pstrm(ns)%stream_mesh, ownedElemCoords=nv_coords)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       do i=1,lsize
          dataptr(i) = dataptr2d_src(1,i)
          lon = nu_coords(2*i-1)
          lat = nu_coords(2*i)
          sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
          sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
          dataptr2d_src(1,i) = coslon * dataptr(i) - sinlon * dataptr2d_src(2,i)
          dataptr2d_src(2,i) = sinlon * dataptr(i) + coslon * dataptr2d_src(2,i)
       enddo
       vector_dst = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name='vector_dst', &
            ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldRegrid(sdat%pstrm(ns)%stream_vector, vector_dst, sdat%pstrm(ns)%routehandle, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(vector_dst, farrayPtr=dataptr2d_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_fldbun_getFldPtr(fldbun_model, trim(uname), data_u_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_fldbun_getFldPtr(fldbun_model, trim(vname), data_v_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       do i = 1,size(data_u_dst)
          lon = sdat%model_lon(i)
          lat = sdat%model_lat(i)
          sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
          sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
          data_u_dst(i) =  coslon * dataptr2d_dst(1,i) + sinlon * dataptr2d_dst(2,i)
          data_v_dst(i) = -sinlon * dataptr2d_dst(1,i) + coslon * dataptr2d_dst(2,i)
       enddo


    endif
    if (pio_iovartype == PIO_REAL) then
       deallocate(data_real)
    endif
    call t_stopf(trim(istr)//'_readpio')

  end subroutine shr_strdata_readstrm

  !===============================================================================
  subroutine shr_strdata_set_stream_iodesc(sdat, &
       pio_subsystem, pioid, fldname, stream_mesh, pio_iodesc, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(in)    :: sdat
    type(iosystem_desc_t) , intent(inout), target :: pio_subsystem
    type(file_desc_t)     , intent(inout)         :: pioid
    character(len=*)      , intent(in)            :: fldname
    type(ESMF_Mesh)       , intent(in)            :: stream_mesh
    type(io_desc_t)       , intent(inout)         :: pio_iodesc
    integer               , intent(out)           :: rc

    ! local variables
    integer                       :: pio_iovartype
    integer                       :: n
    type(var_desc_t)              :: varid
    integer                       :: ndims
    integer, allocatable          :: dimids(:)
    integer, allocatable          :: dimlens(:)
    integer                       :: unlimdid
    type(ESMF_DistGrid)           :: distGrid
    integer                       :: lsize
    integer, pointer              :: compdof(:)
    character(CS)                 :: dimname
    integer                       :: rCode      ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
    character(*), parameter       :: subname = '(shr_strdata_set_stream_iodesc) '
    character(*), parameter       :: F00  = "('(shr_strdata_set_stream_iodesc) ',a,i8,2x,i8,2x,a)"
    character(*), parameter       :: F01  = "('(shr_strdata_set_stream_iodesc) ',a,i8,2x,i8,2x,a)"
    character(*), parameter       :: F02  = "('(shr_strdata_set_stream_iodesc) ',a,i8,2x,i8,2x,i8,2x,a)"
    !-------------------------------------------------------------------------------
    integer :: old_error_handle

    rc = ESMF_SUCCESS

    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR, old_error_handle)

    ! query the first field in the stream dataset
    rcode = pio_inq_varid(pioid, trim(fldname), varid)
    rcode = pio_inq_varndims(pioid, varid, ndims)

    allocate(dimids(ndims))
    allocate(dimlens(ndims))

    rcode = pio_inq_vardimid(pioid, varid, dimids(1:ndims))

    do n = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
    end do

    ! determine compdof for stream
    call ESMF_MeshGet(stream_mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(compdof(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine type of the variable
    rcode = pio_inq_vartype(pioid, varid, pio_iovartype)

    if (ndims == 2) then
       if (sdat%masterproc) then
          write(sdat%logunit,F00) 'setting iodesc for : '//trim(fldname)// &
               ' with dimlens(1), dimlens2 = ',dimlens(1),dimlens(2),&
               ' variable had no time dimension '
       end if
       call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
    else if (ndims == 3) then
       rcode = pio_inq_dimname(pioid, dimids(ndims), dimname)
       if (trim(dimname) == 'time' .or. trim(dimname) == 'nt') then
          if (sdat%masterproc) then
             write(sdat%logunit,F01) 'setting iodesc for : '//trim(fldname)// &
                  ' with dimlens(1), dimlens2 = ',dimlens(1),dimlens(2),&
                  ' variable had time dimension '//trim(dimname)
          end if
          call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       else
          if (sdat%masterproc) then
             write(sdat%logunit,F02) 'setting iodesc for : '//trim(fldname)// &
                  ' with dimlens(1), dimlens2 = ',dimlens(1),dimlens(2),dimlens(3),&
                  ' variable had no time dimension '
          end if
          call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof, pio_iodesc)
       end if
    else
       write(6,*)'ERROR: dimlens= ',dimlens
       call shr_sys_abort(trim(subname)//' only dimlen of 2 and 3 are currently supported')
    end if
    call pio_seterrorhandling(pioid, old_error_handle)

    deallocate(compdof)
    deallocate(dimids)
    deallocate(dimlens)

  end subroutine shr_strdata_set_stream_iodesc

  !===============================================================================
  subroutine shr_strdata_get_stream_pointer(sdat, strm_fld, strm_ptr, rc)

    ! Set a pointer, strm_ptr, for field, strm_fld, into sdat fldbun_model field bundle

    ! input/output variables
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: strm_fld
    real(r8)               , pointer       :: strm_ptr(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer :: ns, nf
    logical :: found
    character(len=*), parameter :: subname='(shr_strdata_get_stream_pointer)'
    character(*)    , parameter :: F00 = "('(shr_strdata_get_stream_pointer) ',8a)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! loop over all input streams and determine if the strm_fld is in the field bundle of the target stream
    do ns = 1, shr_strdata_get_stream_count(sdat)
       found = .false.
       ! Check if requested stream field is read in - and if it is then point into the stream field bundle
       do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
          if (trim(strm_fld) == trim(sdat%pstrm(ns)%fldlist_model(nf))) then
             call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model, trim(sdat%pstrm(ns)%fldlist_model(nf)), &
                  strm_ptr, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (sdat%masterproc) then
                write(sdat%logunit,F00)' strm_ptr is allocated for stream field strm_'//trim(strm_fld)
             end if
             found = .true.
             exit
          end if
       end do
       if (found) exit
    end do

  end subroutine shr_strdata_get_stream_pointer

  !===============================================================================
  subroutine shr_strdata_handle_error(ierr, errorstr)
    use pio, only: pio_noerr

    ! input/output variables
    integer,          intent(in)  :: ierr
    character(len=*), intent(in)  :: errorstr

    ! local variables
    character(len=256) :: errormsg
    !-------------------------------------------------------------------------------

    if (ierr /= PIO_NOERR) then
      write(errormsg, '(a,i6,2a)') '(PIO:', ierr, ') ', trim(errorstr)
      call shr_sys_abort(errormsg)
    end if
  end subroutine shr_strdata_handle_error

end module dshr_strdata_mod
