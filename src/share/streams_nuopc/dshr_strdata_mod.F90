module dshr_strdata_mod

  ! holds data and methods to advance data models
  ! Obtain the model domain and the stream domain for each stream

  use ESMF

  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_file_mod     , only : shr_file_getunit, shr_file_freeunit
  use shr_const_mod    , only : shr_const_pi, shr_const_cDay, shr_const_spval
  use shr_log_mod      , only : logunit => shr_log_Unit
  use shr_cal_mod      , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod      , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod      , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_nl_mod       , only : shr_nl_find_group_name
  use shr_pio_mod      , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
  use shr_string_mod   , only : shr_string_listgetname, shr_string_listisvalid, shr_string_listgetnum

  use dshr_stream_mod  , only : shr_stream_streamtype, shr_stream_getModelFieldList, shr_stream_getStreamFieldList
  use dshr_stream_mod  , only : shr_stream_taxis_cycle, shr_stream_taxis_extend, shr_stream_default
  use dshr_stream_mod  , only : shr_stream_getNFiles
  use dshr_stream_mod  , only : shr_stream_getCurrFile, shr_stream_setCurrFile
  use dshr_stream_mod  , only : shr_stream_getMeshFilename
  use dshr_stream_mod  , only : shr_stream_init_from_xml, shr_stream_init_from_fortran
  use dshr_stream_mod  , only : shr_stream_restWrite, shr_stream_restRead
  use dshr_stream_mod  , only : shr_stream_findBounds
  use dshr_stream_mod  , only : shr_stream_getnextfilename, shr_stream_getprevfilename
  use dshr_stream_mod  , only : shr_stream_getfilefieldname, shr_stream_getData
  use dshr_tinterp_mod , only : shr_tInterp_getCosz, shr_tInterp_getAvgCosz, shr_tInterp_getFactors
  use dshr_methods_mod , only : dshr_fldbun_getfldptr, dshr_fldbun_getfieldN, dshr_fldbun_fldchk, chkerr
  use dshr_methods_mod , only : dshr_fldbun_diagnose, dshr_fldbun_regrid

  use pio              , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
  use pio              , only : pio_openfile, pio_closefile, pio_nowrite
  use pio              , only : pio_seterrorhandling, pio_initdecomp, pio_freedecomp
  use pio              , only : pio_inquire, pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
  use pio              , only : pio_inq_dimlen, pio_inq_vartype, pio_inq_dimname
  use pio              , only : pio_double, pio_real, pio_int, pio_offset_kind
  use pio              , only : pio_read_darray, pio_get_var, pio_setframe
  use pio              , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
  use perf_mod

  implicit none
  private

  public  :: shr_strdata_type
  public  :: shr_strdata_init_from_infiles
  public  :: shr_strdata_init_from_fortran
  public  :: shr_strdata_restRead
  public  :: shr_strdata_restWrite
  public  :: shr_strdata_setOrbs
  public  :: shr_strdata_advance
  public  :: shr_strdata_get_stream_domain  ! public since needed by dshr_mod
  public  :: shr_strdata_get_stream_pointer ! get a pointer into a stream's fldbun_model field bundle
  public  :: shr_strdata_print

  private :: shr_strdata_readnml_from_infiles
  private :: shr_strdata_init_model_domain
  private :: shr_strdata_readLBUB

  ! public data members:
  integer                              :: debug    = 1  ! local debug flag
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  real(r8)         ,parameter, private :: dtlimit_default = 1.5_r8
  integer          ,parameter          :: master_task = 0

  type shr_strdata_perstream
     ! character(CS)                     :: taxMode                  ! stream time axis cycling mode
     ! real(r8)                          :: dtlimit                  ! stream dt max/min limit
     ! character(CS)                     :: tintalgo                 ! stream time interpolation algorithm
     ! character(CS)                     :: readmode                 ! stream file(s) read mode

     character(CL)                       :: stream_vectors            ! stream vectors names from shr_strdata_nml
     character(CL)                       :: stream_meshfile           ! stream mesh file from stream txt file
     type(ESMF_Mesh)                     :: stream_mesh               ! stream mesh created from stream mesh file
     type(io_desc_t)                     :: stream_pio_iodesc         ! stream pio descriptor
     logical                             :: stream_pio_iodesc_set =.false.  ! true=>pio iodesc has been set

                                                                      ! stream-> model mapping info
     logical                             :: domaps
     ! character(CL)                     :: mapalgo                   ! scalar map algorithm
     ! character(CL)                     :: mapmask                   ! scalar map mask
     type(ESMF_RouteHandle)              :: routehandle               ! stream n -> model mesh mapping

     ! field bundles
     ! note that the fields in fldbun_stream_lb and fldbun_stream_ub contain the the names fldlist_model
     character(CS), allocatable          :: fldlist_stream(:)         ! names of fields read in from stream
     character(CS), allocatable          :: fldlist_model(:)          ! names of fields in model (1/1 correspondence with fldlist_stream)
     type(ESMF_FieldBundle)              :: fldbun_stream_lb          ! stream n field bundle for lb of time period (stream grid)
     type(ESMF_FieldBundle)              :: fldbun_stream_ub          ! stream n field bundle for ub of time period (stream grid)
     type(ESMF_FieldBundle)              :: fldbun_model_lb           ! stream n field bundle for lb of time period (model grid)
     type(ESMF_FieldBundle)              :: fldbun_model_ub           ! stream n field bundle for ub of time period (model grid)
     type(ESMF_FieldBundle)              :: fldbun_model              ! stream n field bundle for model time (model grid)
     type(ESMF_FieldBundle), allocatable :: fldbun_stream_alltimes(:) ! field bundle for stream n for all time slices for stream
     integer                             :: ymdLB                     ! stream n ymd lower bound
     integer                             :: todLB
     integer                             :: ymdUB
     integer                             :: todUB
     integer                             :: ustrm
     integer                             :: vstrm
     real(r8)                            :: dtmin
     real(r8)                            :: dtmax
     type(ESMF_Field)                    :: field_coszen              ! needed for coszen time interp
  end type shr_strdata_perstream


  type shr_strdata_type
     type(shr_strdata_perstream), allocatable :: pstrm(:)

     ! stream domain info
     type(shr_stream_streamType), pointer :: stream(:)=> null()                    ! stream datatype
     integer                        :: nstreams                          ! number of streams set in shr_strdata_readnml
     integer                        :: nvectors                          ! number of vectors set in shr_strdata_readnml

     ! mpi info
     integer                        :: mpicom
     integer                        :: my_task

     ! pio info
     integer                        :: io_type
     integer                        :: io_format
     type(iosystem_desc_t), pointer :: pio_subsystem => null()

     ! data required by stream  cosz t-interp method, set by user
     real(r8)                       :: eccen
     real(r8)                       :: mvelpp
     real(r8)                       :: lambm0
     real(r8)                       :: obliqr
     integer                        :: modeldt                           ! model dt in seconds

     ! model domain info
     type(ESMF_Mesh)                :: model_mesh                        ! model mesh
     real(r8), pointer              :: model_lon(:) => null()            ! model longitudes
     real(r8), pointer              :: model_lat(:) => null()            ! model latitudes
     real(r8), pointer              :: model_lev(:) => null()            ! model levels (if needed)
     integer                        :: model_nxg                         ! model global domain lon size
     integer                        :: model_nyg                         ! model global domain lat size
     integer                        :: model_nzg                         ! model global domain vertical size
     integer                        :: model_lsize                       ! model local domain size
     integer, pointer               :: model_gindex(:)                   ! model global index spzce
     integer                        :: model_gsize                       ! model global domain size

     ! time info
     type(ESMF_CLock)               :: model_clock
     integer                        :: ymd, tod                          ! model time
     character(CL)                  :: model_calendar                    ! model calendar for ymd,tod
     real(r8), allocatable          :: tavCoszen(:)
  end type shr_strdata_type

  integer          ,parameter :: iotype_std_netcdf = -99 ! non pio option
  real(r8)         ,parameter :: deg2rad = SHR_CONST_PI/180.0_r8
  character(*)     ,parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_strdata_init_from_infiles(sdat, xmlfilename, mesh, clock, mpicom, &
       compid, logunit, reset_mask, model_maskfile, rc)

    ! input/output variables
    type(shr_strdata_type)     , intent(inout) :: sdat
    character(len=*)           , intent(in)    :: xmlfilename ! for shr_strdata_nml namelist
    type(ESMF_Mesh)            , intent(inout) :: mesh
    type(ESMF_Clock)           , intent(in)    :: clock
    integer                    , intent(in)    :: mpicom
    integer                    , intent(in)    :: compid
    integer                    , intent(in)    :: logunit
    logical         , optional , intent(in)    :: reset_mask
    character(len=*), optional , intent(in)    :: model_maskfile
    integer                    , intent(out)   :: rc

    ! local variables
    logical       :: masterproc
    integer       :: ns ! stream index
    integer       :: i  ! generic index
    character(CL) :: str
    integer       :: yearFirst     ! first year to use
    integer       :: yearLast      ! last  year to use
    integer       :: yearAlign     ! align yearFirst with this model year
    character(CL) :: filename
    integer       :: ierr
    character(len=*), parameter  :: subname='(shr_strdata_mod:dshr_sdat_init_from_infiles)'
    character(*)    , parameter  :: F01="('(shr_init_strdata) ',a,2f10.4)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize sdat mpi info
    sdat%mpicom = mpicom
    call mpi_comm_rank(sdat%mpicom, sdat%my_task, ierr)
    masterproc = (sdat%my_task == master_task)

    ! Initialize sdat  pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! Read xml file
    call shr_stream_init_from_xml(xmlfilename, sdat%stream, masterproc, rc=rc)
    sdat%nstreams = size(sdat%stream)
    allocate(sdat%pstrm(sdat%nstreams))

    ! Initialize the sdat model domain info
    sdat%pstrm(:)%ymdLB = -1
    sdat%pstrm(:)%todLB = -1
    sdat%pstrm(:)%ymdUB = -1
    sdat%pstrm(:)%todUB = -1
    sdat%pstrm(:)%dtmin = 1.0e30
    sdat%pstrm(:)%dtmax = 0.0
    sdat%eccen          = SHR_ORB_UNDEF_REAL
    sdat%mvelpp         = SHR_ORB_UNDEF_REAL
    sdat%lambm0         = SHR_ORB_UNDEF_REAL
    sdat%obliqr         = SHR_ORB_UNDEF_REAL
    sdat%modeldt        = 0
    sdat%model_calendar = shr_cal_noleap
    sdat%pstrm(:)%dtmin = 1.0e30
    sdat%pstrm(:)%dtmax = 0.0e0

    call shr_strdata_init_model_domain(sdat, mesh, compid, &
         reset_mask=reset_mask, model_maskfile=model_maskfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now finish initializing sdat
    call shr_strdata_init(sdat, clock, compid, logunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_strdata_init_from_infiles

  !===============================================================================
  subroutine shr_strdata_init_from_fortran( sdat, mpicom, compid,                  &
       model_mesh, model_nxg, model_nyg, model_clock,                              &
       stream_meshfile, stream_filenames, stream_FldListFile, stream_fldListModel, &
       stream_yearFirst, stream_yearLast, stream_yearAlign, stream_offset,         &
       stream_taxMode, stream_dtlimit, stream_tintalgo, stream_readmode,           &
       stream_mapalgo, stream_mapmask, rc)

    ! Set strdata and stream info from fortran interface.
    ! Note: When this is called, previous settings are reset to defaults
    ! and then the values passed are used.

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: sdat                ! strdata data data-type
    integer                ,intent(in)   :: mpicom              ! mpi comm
    integer                ,intent(in)   :: compid
    type(ESMF_Mesh)        ,intent(in)   :: model_mesh
    integer                ,intent(in)   :: model_nxg
    integer                ,intent(in)   :: model_nyg
    type(ESMF_Clock)       ,intent(in)   :: model_clock
    character(len=*)       ,intent(in)   :: stream_meshfile
    character(*)           ,intent(in)   :: stream_FileNames(:) ! filenames for stream data
    character(*)           ,intent(in)   :: stream_fldListFile  ! file field names, colon delim list
    character(*)           ,intent(in)   :: stream_fldListModel ! model field names, colon delim list
    integer                ,intent(in)   :: stream_yearFirst    ! first year to use
    integer                ,intent(in)   :: stream_yearLast     ! last  year to use
    integer                ,intent(in)   :: stream_yearAlign    ! align yearFirst with this model year
    integer                ,intent(in)   :: stream_offset       ! offset in seconds of stream data
    character(*)           ,intent(in)   :: stream_taxMode
    real(r8)               ,intent(in)   :: stream_dtlimit
    character(*)           ,intent(in)   :: stream_mapalgo      ! scalar map algorithm
    character(*)           ,intent(in)   :: stream_mapmask      ! scalar map mask
    character(*)           ,intent(in)   :: stream_tintalgo     ! time interpolation algorithm
    character(*)           ,intent(in)   :: stream_readmode     ! file read mode
    integer                ,intent(out)  :: rc                  ! error status

    ! local variables
    character(*),parameter  :: subName = "(shr_strdata_create) "
    character(*),parameter  :: F00 = "('(shr_strdata_create) ',8a)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Assume only 1 stream
    ! TODO: set all default values - dont read shr_strdata_readnml to set default values

    sdat%nstreams    = 1
    sdat%model_mesh  = model_mesh
    sdat%stream(1)%taxmode  = stream_taxMode
    sdat%stream(1)%dtlimit  = stream_dtlimit
    sdat%stream(1)%mapalgo  = stream_mapalgo
    sdat%stream(1)%mapmask  = stream_mapmask
    sdat%stream(1)%tinterpalgo = stream_tintalgo
    sdat%stream(1)%readmode = stream_readmode ! single or full_file
    if (trim(sdat%stream(1)%taxmode) == trim(shr_stream_taxis_extend)) then
       ! reset dtlimit if necessary
       sdat%stream(1)%dtlimit = 1.0e30
    end if

    ! Initialize sdat  pio
    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! Initialize sdat stream info
    call shr_stream_init_from_fortran(sdat%stream(1), stream_meshfile, &
         stream_yearFirst, stream_yearLast, stream_yearAlign, stream_offset, stream_taxmode, &
         stream_fldlistFile, stream_fldlistModel, stream_fileNames)

    ! Initialize sdat model domain info
    call shr_strdata_init_model_domain(sdat, model_mesh, compid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize the rest of sdat
    call shr_strdata_init(sdat, model_clock, compid, logunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_strdata_init_from_fortran

  !===============================================================================
  subroutine shr_strdata_init_model_domain( sdat, model_mesh, compid, reset_mask, model_maskfile, rc)

    ! ----------------------------------------------
    ! Initialize sdat model domain info
    ! ----------------------------------------------

    ! input/output variables
    type(shr_strdata_type)     , intent(inout) :: sdat
    type(ESMF_Mesh)            , intent(in)    :: model_mesh
    integer                    , intent(in)    :: compid
    logical         , optional , intent(in)    :: reset_mask
    character(len=*), optional , intent(in)    :: model_maskfile
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
    type(ESMF_Array)      :: elemMaskArray
    integer, allocatable  :: elemMask(:)
    type(file_desc_t)     :: pioid
    type(var_desc_t)      :: varid
    type(io_desc_t)       :: pio_iodesc
    integer               :: rcode
    logical               :: lreset_mask
    character(CL)         :: lmodel_maskfile
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(reset_mask)) then
       lreset_mask = reset_mask
    else
       lreset_mask = .false.
    end if
    if (present(model_maskfile)) then
       lmodel_maskfile = trim(model_maskfile)
    else
       lmodel_maskfile = ''
    end if

    ! initialize sdat%model_mesh
    sdat%model_mesh = model_mesh

    ! initialize sdat%lsize
    call ESMF_MeshGet(model_mesh, elementdistGrid=distGrid, rc=rc)
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

    ! determine sdat%model_lon, sdat%model_lat and sdat%model_levs
    ! TODO: add levs - assume this is in the mesh - for now not reference here
    call ESMF_MeshGet(model_mesh, spatialDim=spatialDim, &
         numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(model_mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(sdat%model_lon(numOwnedElements))
    allocate(sdat%model_lat(numOwnedElements))
    allocate(sdat%model_lev(numOwnedElements))
    do n = 1, numOwnedElements
       sdat%model_lon(n) = ownedElemCoords(2*n-1)
       sdat%model_lat(n) = ownedElemCoords(2*n)
       sdat%model_lev(n) = shr_const_spval
    end do

    ! initialize the model mask if appropriate
    if (present(model_maskfile)) then
       allocate(elemMask(sdat%model_lsize))
       ! obtain model mask from separate model_maskfile
       rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(model_maskfile), pio_nowrite)
       call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
       ! TODO: check that mask name is on file and if not abort
       rcode = pio_inq_varid(pioid, 'mask', varid) ! assume mask name on domain file
       call pio_initdecomp(sdat%pio_subsystem, pio_int, (/sdat%model_nxg, sdat%model_nyg/), sdat%model_gindex, pio_iodesc)
       allocate(elemMask(sdat%model_lsize))
       call pio_read_darray(pioid, varid, pio_iodesc, elemMask, rcode)
       call pio_closefile(pioid)
       call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
       call ESMF_MeshSet(model_mesh, elementMask=elemMask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(elemMask)
    else if (lreset_mask) then
       allocate(elemMask(sdat%model_lsize))
       elemMask(:) = 1._r8
       call ESMF_MeshSet(model_mesh, elementMask=elemMask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(elemMask)
    end if

  end subroutine shr_strdata_init_model_domain

  !===============================================================================
  subroutine shr_strdata_init(sdat, model_clock, compid, logunit, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_Clock)       , intent(in)    :: model_clock
    integer                , intent(in)    :: compid
    integer                , intent(in)    :: logunit
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Calendar)     :: esmf_calendar   ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype    ! esmf calendar type
    type(ESMF_DistGrid)     :: distgrid
    character(CS)           :: calendar        ! calendar name
    integer                 :: dimcount
    integer, allocatable    :: minIndexPTile(:,:)
    integer, allocatable    :: maxIndexPTile(:,:)
    integer                 :: lnx, lny        ! global mesh dimensions
    integer                 :: ne              ! number of local mesh elements
    integer                 :: ns              ! stream index
    integer                 :: n,m,k           ! generic index
    character(CL)           :: fileName        ! generic file name
    integer                 :: nfiles          ! number of data files for a given stream
    character(CS)           :: uname           ! u vector field name
    character(CS)           :: vname           ! v vector field name
    integer                 :: nu, nv          ! vector indices
    integer                 :: nstream         ! loop stream index
    integer                 :: nvector         ! loop vector index
    integer                 :: nfld            ! loop stream field index
    integer                 :: nflds           ! total number of fields in a given stream
    type(ESMF_Field)        :: lfield          ! temporary
    type(ESMF_Field)        :: lfield_src      ! temporary
    type(ESMF_Field)        :: lfield_dst      ! temporary
    integer                 :: srcTermProcessing_Value = 0 ! should this be a module variable?
    integer , pointer       :: stream_gindex(:)
    integer                 :: stream_lsize
    character(CS)           :: tmpstr
    integer                 :: ierr
    integer                 :: localpet
    logical                 :: fileExists
    type(ESMF_VM)           :: vm
    logical                 :: masterproc
    integer                 :: nvars
    character(len=*), parameter :: subname='(shr_strdata_mod:shr_sdat_init)'
    character(*)    , parameter :: F00 = "('(shr_sdat_init) ',a)"
    character(*)    , parameter :: F01 = "('(shr_sdat) ',a,2x,i8)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_VmGetCurrent(vm, rc=rc)
    call ESMF_VMGet(vm, localpet=localPet, rc=rc)
    masterproc= localPet==master_task

    do ns = 1,sdat%nstreams
       ! Initialize calendar for stream n
       call ESMF_VMBroadCast(vm, sdat%stream(ns)%calendar, CS, 0, rc=rc)

       ! Create the target stream mesh from the stream mesh file
       ! TODO: add functionality if the stream mesh needs to be created from a grid
       call shr_stream_getMeshFileName (sdat%stream(ns), filename)

       if (masterproc) then
          inquire(file=trim(filename),exist=fileExists)
          if (.not. fileExists) then
             write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
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
       sdat%pstrm(ns)%fldbun_model_lb = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%pstrm(ns)%fldbun_model_ub = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%pstrm(ns)%fldbun_model    = ESMF_FieldBundleCreate(rc=rc) ! time interpolation on model mesh
       do nfld = 1, nvars
          ! create temporary fields on model mesh and add the fields to the field bundle
          lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_model_lb, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_model_ub, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          lfield = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_model   , (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do

       ! Create field bundles on stream mesh but with fldlist_model names
       sdat%pstrm(ns)%fldbun_stream_lb = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at lower time bound
       sdat%pstrm(ns)%fldbun_stream_ub = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at upper time bound
       if (masterproc) then
          write(logunit,F00)' initializing fldbun_stream_model_lb and fldbun_stream_model_up on stream mesh'
       end if
       do nfld = 1, nvars
          ! create temporary fields on stream mesh and add the fields to the field bundle
          lfield = ESMF_FieldCreate(sdat%pstrm(ns)%stream_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_stream_lb, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          lfield = ESMF_FieldCreate(sdat%pstrm(ns)%stream_mesh, ESMF_TYPEKIND_r8, name=trim(sdat%pstrm(ns)%fldlist_model(nfld)), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%pstrm(ns)%fldbun_stream_ub, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (masterproc) then
             write(logunit,F00)' adding field '//trim(sdat%pstrm(ns)%fldlist_model(nfld))//&
                  ' to fldbun_stream_lb and fldbun_stream_ub'
          end if
       end do ! end of loop over streams (ns)

       ! Create a field for coszen time interpolation for this stream if needed
       if (trim(sdat%stream(ns)%tinterpalgo) == 'coszen') then
          sdat%pstrm(ns)%field_coszen = ESMF_FieldCreate(sdat%pstrm(ns)%stream_mesh, ESMF_TYPEKIND_r8, name='tavCosz', &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       endif

       ! ------------------------------------
       ! Create the sdat route handles for mapping the stream -> model
       ! For now only  create bilinear route handle for stream(n) -> model mapping
       ! ------------------------------------

       ! create the source and destination fields needed for the route handles
       ! these fields will be used to create the route handles
       ! since all fields in a stream share the same mesh and there is only a unique model mesh
       ! can do this outside of a stream loop by just using the first stream index
       ! TODO: Determine if mask was reset - and if so if the masks are the same - if they are not do the bilinear
       ! interpolation below -

       call dshr_fldbun_getFieldN(sdat%pstrm(ns)%fldbun_stream_lb, 1, lfield_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_fldbun_getFieldN(sdat%pstrm(ns)%fldbun_model_lb, 1, lfield_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldRegridStore(lfield_src, lfield_dst, routehandle=sdat%pstrm(ns)%routehandle, &
           regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
           polemethod=ESMF_POLEMETHOD_ALLAVG, &
         ! extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
         ! regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
         ! polemethod=ESMF_POLEMETHOD_NONE, &
           dstMaskValues = (/0/), &  ! ignore destination points where the mask is 0
           srcTermProcessing=srcTermProcessing_Value, ignoreDegenerate=.true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end do ! end of loop over streams

    ! ------------------------------------
    ! determine the streams for ustrm and vstrm
    ! ------------------------------------
    ! check vectors and compute ustrm,vstrm
    ! determine if vector field names match any field names in the input stream

    do m = 1,sdat%nvectors
       ! check that vector field list is a valid colon delimited string
       if (.not. shr_string_listIsValid(sdat%pstrm(m)%stream_vectors)) then
          write(logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(sdat%pstrm(m)%stream_vectors)
          call shr_sys_abort(subname//': vec fldlist invalid:'//trim(sdat%pstrm(m)%stream_vectors))
       endif

       ! check that only 2 fields are contained for any vector pairing
       if (shr_string_listGetNum(sdat%pstrm(m)%stream_vectors) /= 2) then
          write(logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(sdat%pstrm(m)%stream_vectors)
          call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(sdat%pstrm(m)%stream_vectors))
       endif

       ! get name of the first and second field in the colon delimited string
       call shr_string_listGetName(sdat%pstrm(m)%stream_vectors,1,uname)
       call shr_string_listGetName(sdat%pstrm(m)%stream_vectors,2,vname)

       ! loop through the streams and find which stream(s) contain the vector field pair
       ! normally both fields in the pair will be on one stream
       nu = 0
       nv = 0
       do ns = 1,sdat%nstreams
          if (dshr_fldbun_fldchk(sdat%pstrm(ns)%fldbun_stream_lb, trim(uname), rc=rc)) nu = n
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (dshr_fldbun_fldchk(sdat%pstrm(ns)%fldbun_stream_lb, trim(vname), rc=rc)) nv = n
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo

       ! error checks
       if (nu == 0  .or. nv == 0) then
          ! if the input fields are not contained - then abort
          write(logunit,*) trim(subname),' vec flds not found  m=',m,trim(sdat%pstrm(m)%stream_vectors)
          call shr_sys_abort(subname//': vec flds not found:'//trim(sdat%pstrm(m)%stream_vectors))
       else if (nu /= nv) then
          ! mesh files for nu and nv stream must be the same
          if (trim(sdat%pstrm(nu)%stream_meshfile) /= trim(sdat%pstrm(nv)%stream_meshfile)) then
             write(logunit,*) trim(subname),' vec fld mesh files are not same m=',m,trim(sdat%pstrm(m)%stream_vectors)
             write(logunit,*) trim(subname),' vec mesh file file for nu = ',nu,' is ',trim(sdat%pstrm(nu)%stream_meshfile)
             write(logunit,*) trim(subname),' vec mesh file file for nv = ',nv,' is ',trim(sdat%pstrm(nv)%stream_meshfile)
             call shr_sys_abort(subname//': vec fld mesh files must be the same:'//trim(sdat%pstrm(m)%stream_vectors))
          end if
       else
          ! now set the stream indices for ustrm and vstrm
          sdat%pstrm(m)%ustrm = nu
          sdat%pstrm(m)%vstrm = nv
       end if
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
       write(logunit,*) ' successfully initialized sdat'
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
    type(var_desc_t)      :: varid
    type(file_desc_t)     :: pioid
    integer               :: rcode
    character(CL)         :: filename
    type(io_desc_t)       :: pio_iodesc
    real(r4), allocatable :: data_real(:)
    integer               :: lsize
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine the file to open
    if (sdat%my_task == master_task) then
       call shr_stream_getData(sdat%stream(stream_index), 1, filename)
    end if
    call shr_mpi_bcast(filename, sdat%mpicom, 'streamfile')

    ! Open the file
    rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(filename), pio_nowrite)

    ! Create the pio iodesc for fldname
    call shr_strdata_set_stream_iodesc(sdat%pio_subsystem, pioid, &
         trim(fldname), sdat%pstrm(stream_index)%stream_mesh, pio_iodesc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now read in the data for fldname
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    lsize = size(flddata)
    allocate(data_real(lsize))
    rcode = pio_inq_varid(pioid, trim(fldname), varid)
    call pio_read_darray(pioid, varid, pio_iodesc, data_real, rcode)
    flddata(:) = real(data_real(:), kind=r8)
    deallocate(data_real)

    ! Free the memory associate with the iodesc
    call pio_freedecomp(pioid, pio_iodesc)

    ! Close the file
    call pio_closefile(pioid)

  end subroutine shr_strdata_get_stream_domain

  !===============================================================================
  subroutine shr_strdata_advance(sdat, ymd, tod, mpicom, logunit, istr, timers, rc)

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
    integer                ,intent(in)          :: mpicom
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

    if (sdat%nstreams < 1) return ! TODO: is this needed

    lstr = trim(istr)

    ltimers = .true.
    if (present(timers)) then
       ltimers = timers
    endif

    if (.not.ltimers) call t_adj_detailf(tadj)

    call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
    call t_startf(trim(lstr)//trim(timname)//'_total')

    call mpi_comm_size(mpicom, npes, ierr)
    call mpi_comm_rank(mpicom, my_task,ierr)

    sdat%ymd = ymd
    sdat%tod = tod

    if (sdat%nstreams > 0) then
       allocate(newData(sdat%nstreams))
       allocate(ymdmod(sdat%nstreams))

       do ns = 1,sdat%nstreams
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

          call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
          call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

          select case(sdat%stream(ns)%readmode)
          case ('single')
             call shr_strdata_readLBUB(sdat%mpicom, sdat%my_task, sdat%stream(ns), sdat%pstrm(ns)%stream_mesh, &
                  sdat%pstrm(ns)%fldlist_stream, sdat%pstrm(ns)%fldlist_model, &
                  sdat%pstrm(ns)%fldbun_stream_lb, sdat%pstrm(ns)%fldbun_stream_ub, &
                  sdat%pio_subsystem, sdat%io_type, sdat%pstrm(ns)%stream_pio_iodesc_set, sdat%pstrm(ns)%stream_pio_iodesc, &
                  ymdmod(ns), todmod, sdat%pstrm(ns)%ymdLB, sdat%pstrm(ns)%todLB, &
                  sdat%pstrm(ns)%ymdUB, sdat%pstrm(ns)%todUB, &
                  newData(ns), trim(lstr)//'_readLBUB', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          case ('full_file')
             ! TODO: need to put in capability to read all stream data at once
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(sdat%stream(ns)%readmode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(sdat%stream(ns)%readmode))
          end select

          if (debug > 0) then
             write(logunit,*) trim(subname),' newData flag = ',ns,newData(ns)
             write(logunit,*) trim(subname),' LB ymd,tod = ',ns,sdat%pstrm(ns)%ymdLB,sdat%pstrm(ns)%todLB
             write(logunit,*) trim(subname),' UB ymd,tod = ',ns,sdat%pstrm(ns)%ymdUB,sdat%pstrm(ns)%todUB
          endif

          ! ---------------------------------------------------------
          ! If new data is read in:
          ! ---------------------------------------------------------

          if (newData(ns)) then

             ! Reset time bounds if newdata read in
             call shr_cal_date2ymd(sdat%pstrm(ns)%ymdLB,year,month,day)
             call shr_cal_timeSet(timeLB,sdat%pstrm(ns)%ymdLB,0,sdat%stream(ns)%calendar)
             call shr_cal_timeSet(timeUB,sdat%pstrm(ns)%ymdUB,0,sdat%stream(ns)%calendar)
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

             ! If new data was read in, spatially interpolate the lower and upper bound data to the model grid
             call t_startf(trim(lstr)//trim(timname)//'_map')
             if (debug > 0) then
                call dshr_fldbun_diagnose(sdat%pstrm(ns)%fldbun_stream_lb, subname//':fldbun_stream_lb',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_diagnose(sdat%pstrm(ns)%fldbun_stream_ub, subname//':fldbun_stream_ub',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif

             ! regrid
             call dshr_fldbun_regrid(sdat%pstrm(ns)%fldbun_stream_lb, sdat%pstrm(ns)%fldbun_model_lb, &
                  sdat%pstrm(ns)%routehandle, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_regrid(sdat%pstrm(ns)%fldbun_stream_ub, sdat%pstrm(ns)%fldbun_model_ub, &
                  sdat%pstrm(ns)%routehandle, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             if (debug > 0) then
                call dshr_fldbun_diagnose(sdat%pstrm(ns)%fldbun_model_lb, subname//':fldbun_model_lb',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_diagnose(sdat%pstrm(ns)%fldbun_model_ub, subname//':fldbun_model_ub',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
             call t_stopf(trim(lstr)//trim(timname)//'_map')
          endif

          call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')

       enddo ! end of loop over streams

       ! ---------------------------------------------------------
       ! remap with vectors if needed
       ! ---------------------------------------------------------

       do m = 1,sdat%nvectors
          nu = sdat%pstrm(m)%ustrm ! nu is the stream index that contains the u vector
          nv = sdat%pstrm(m)%vstrm ! nv is the stream index that contains the v vector

          ! TODO: this is not correct logic - need to change it
          if ((nu > 0 .or. nv > 0) .and. (newdata(nu) .or. newdata(nv))) then

             call t_startf(trim(lstr)//trim(timname)//'_vect')

             ! get lon and lat of stream nu
             call ESMF_MeshGet(sdat%pstrm(nu)%stream_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(nu_coords(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%pstrm(nu)%stream_mesh, ownedElemCoords=nu_coords)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get lon and lat of stream nv
             call ESMF_MeshGet(sdat%pstrm(nv)%stream_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(nv_coords(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%pstrm(nv)%stream_mesh, ownedElemCoords=nv_coords)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! create source field and destination fields
             ! TODO: assume that the two meshes are idential - but need to confirm this
             field_src = ESMF_FieldCreate(sdat%pstrm(nu)%stream_mesh, ESMF_TYPEKIND_r8, name='field_src', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             field_dst = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_r8, name='field_dst', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! set pointers to source and destination data that will be filled in with rotation to cart3d
             call ESMF_FieldGet(field_src, farrayPtr=data2d_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(field_dst, farrayPtr=data2d_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get names of vector field pairs
             call shr_string_listGetName(sdat%pstrm(m)%stream_vectors, 1, uname)
             call shr_string_listGetName(sdat%pstrm(m)%stream_vectors, 2, vname)

             ! map lower bs: rotate source data, regrid, then rotate back
             call dshr_fldbun_getFldPtr(sdat%pstrm(nu)%fldbun_stream_lb, trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nu)%fldbun_stream_lb, trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nv)%fldbun_stream_lb, trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nv)%fldbun_stream_lb, trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do i = 1,size(data_u_src)
                lon = nu_coords(2*i-1)
                lat = nu_coords(2*i)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                data2d_src(1,i) = coslon * data_u_src(i) - sinlon * data_v_src(i)
                data2d_src(2,i) = sinlon * data_u_src(i) + coslon * data_v_src(i)
             enddo
             call ESMF_FieldRegrid(field_src, field_dst, sdat%pstrm(nu)%routehandle, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do i = 1,size(data_u_dst)
                lon = sdat%model_lon(i)
                lat = sdat%model_lat(i)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,i)
                uy = data2d_dst(2,i)
                data_u_dst(i) =  coslon * data_u_dst(i) + sinlon * data_v_dst(i)
                data_v_dst(i) = -sinlon * data_u_dst(i) + coslon * data_v_dst(i)
             enddo

             ! map upper bs: rotate source data, regrid, then rotate back
             call dshr_fldbun_getFldPtr(sdat%pstrm(nu)%fldbun_stream_ub, trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nu)%fldbun_stream_ub, trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nv)%fldbun_stream_ub, trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getFldPtr(sdat%pstrm(nv)%fldbun_stream_ub, trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do i = 1,size(data_u_src)
                lon = nu_coords(2*i-1)
                lat = nu_coords(2*i)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                data2d_src(1,i) = coslon * data_u_src(i) - sinlon * data_v_src(i)
                data2d_src(2,i) = sinlon * data_u_src(i) + coslon * data_v_src(i)
             enddo
             call ESMF_FieldRegrid(field_src, field_dst, sdat%pstrm(nu)%routehandle, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do i = 1,size(data_u_dst)
                lon = sdat%model_lon(i)
                lat = sdat%model_lat(i)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,i)
                uy = data2d_dst(2,i)
                data_u_dst(i) =  coslon * data_u_dst(i) + sinlon * data_v_dst(i)
                data_v_dst(i) = -sinlon * data_u_dst(i) + coslon * data_v_dst(i)
             enddo

             deallocate(nu_coords)
             deallocate(nv_coords)

             call t_stopf(trim(lstr)//trim(timname)//'_vect')

          endif ! end block nu > 0 or nv>0
       enddo

       ! ---------------------------------------------------------
       ! Do time interpolation to create fldbun_model
       ! ---------------------------------------------------------

       do ns = 1,sdat%nstreams

          if (trim(sdat%stream(ns)%tinterpalgo) == 'coszen') then

             ! ------------------------------------------
             ! time interpolation method is coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_coszen')
             allocate(coszen(sdat%model_lsize))

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
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model_lb, sdat%pstrm(ns)%fldlist_model(nf), dataptr_lb, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                do i = 1,lsize
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
                  ymdmod(ns), todmod, flb, fub, calendar=sdat%stream(ns)%calendar, algo=trim(sdat%stream(ns)%tinterpalgo), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (debug > 0) then
                write(logunit,F01) trim(subname),' interp = ',ns,flb,fub
             endif

             do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model   , sdat%pstrm(ns)%fldlist_model(nf), dataptr   , rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model_lb, sdat%pstrm(ns)%fldlist_model(nf), dataptr_lb, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model_ub, sdat%pstrm(ns)%fldlist_model(nf), dataptr_ub, rc=rc)
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

  subroutine shr_strdata_restWrite(sdat, filename, mpicom, str1, str2)

    type(shr_strdata_type) ,intent(inout) :: sdat
    character(len=*)       ,intent(in)    :: filename
    integer                ,intent(in)    :: mpicom
    character(len=*)       ,intent(in)    :: str1
    character(len=*)       ,intent(in)    :: str2

    !--- local ----
    integer :: my_task,ier
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom,my_task,ier)
    if (my_task == 0) then
       call shr_stream_restWrite(sdat%stream,trim(filename),trim(str1),trim(str2),sdat%nstreams)
    endif

  end subroutine shr_strdata_restWrite

  !===============================================================================
  subroutine shr_strdata_restRead(sdat, filename, mpicom)

    type(shr_strdata_type) ,intent(inout) :: sdat
    character(len=*)       ,intent(in)    :: filename
    integer                ,intent(in)    :: mpicom

    !--- local ----
    integer :: my_task,ier
    !-------------------------------------------------------------------------------

    call mpi_comm_rank(mpicom,my_task,ier)
    if (my_task == 0) then
       call shr_stream_restRead(sdat%stream, trim(filename), sdat%nstreams)
    endif

  end subroutine shr_strdata_restRead

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
  subroutine shr_strdata_readnml_from_infiles(sdat, file, mpicom)

    ! Reads shr_strdata_nml namelist input

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: sdat   ! strdata data data-type
    character(*) ,optional ,intent(in)   :: file   ! file to read strdata from
    integer      ,optional ,intent(in)   :: mpicom ! mpi comm

    ! local variables
    integer, parameter :: nstrmax=1
    integer, parameter :: nvecmax=1
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
    real(r8)      :: dtlimit(nStrMax)   ! delta time limiter
    character(CL) :: vectors(nVecMax)   ! define vectors to vector map
    character(CL) :: mapalgo(nStrMax)   ! scalar map algorithm
    character(CL) :: mapmask(nStrMax)   ! scalar map mask
    character(CL) :: tintalgo(nStrMax)  ! time interpolation algorithm
    character(CL) :: readmode(nStrMax)  ! file read mode
    character(CL) :: fileName           ! generic file name

    !----- define namelist -----
    namelist / shr_strdata_nml / &
           nx_global ,           &
           ny_global ,           &
           streams   ,           &
           taxMode   ,           &
           dtlimit   ,           &
           vectors   ,           &
           mapalgo   ,           &
           mapmask   ,           &
           tintalgo  ,           &
           readmode

    !----- formats -----
    character(*),parameter :: subName = "(shr_strdata_readnml_from_infiles) "
    character(*),parameter ::   F00 = "('(shr_strdata_readnml_from_infiles) ',8a)"
    character(*),parameter ::   F01 = "('(shr_strdata_readnml_from_infiles) ',a,i6,a)"
    character(*),parameter ::   F02 = "('(shr_strdata_readnml_from_infiles) ',a,es13.6)"
    character(*),parameter ::   F03 = "('(shr_strdata_readnml_from_infiles) ',a,l6)"
    character(*),parameter ::   F04 = "('(shr_strdata_readnml_from_infiles) ',a,i2,a,a)"
    character(*),parameter ::   F20 = "('(shr_strdata_readnml_from_infiles) ',a,i6,a)"
    character(*),parameter ::   F90 = "('(shr_strdata_readnml_from_infiles) ',58('-'))"
    !-------------------------------------------------------------------------------

    my_task = 0
    ntasks = 1

    if (present(mpicom)) then
       call mpi_comm_rank(mpicom, my_task, rCode)
       call mpi_comm_size(mpicom, ntasks, rCode)
    endif
    allocate(sdat%stream(nStrMax))
    allocate(sdat%pstrm(nStrMax))
    !--master--task--
    if (my_task == master_task) then

       ! set default values for namelist vars
       streams(:)  = trim(shr_strdata_nullstr)
       taxMode(:)  = trim(shr_stream_taxis_cycle)
       dtlimit(:)  = dtlimit_default
       vectors(:)  = trim(shr_strdata_nullstr)
       mapalgo(:)  = 'bilinear'
       mapmask(:)  = 'dstmask'
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
       do n=1,nStrMax
          call shr_stream_default(sdat%stream(n))
       enddo
       sdat%stream(:)%taxmode        = taxMode(:)
       sdat%stream(:)%dtlimit        = dtlimit(:)
       sdat%pstrm(:)%stream_vectors = vectors(:)
       sdat%stream(:)%mapalgo        = mapalgo(:)
       sdat%stream(:)%mapmask        = mapmask(:)
       sdat%stream(:)%tinterpalgo       = tintalgo(:)
       sdat%stream(:)%readmode       = readmode(:)

       sdat%nstreams = 0
       do n=1,nStrMax
          if (trim(sdat%stream(n)%taxmode) == trim(shr_stream_taxis_extend)) then
             sdat%stream(n)%dtlimit = 1.0e30
          end if
       end do

       sdat%nvectors = 0
       do n=1,nVecMax
          if (trim(vectors(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nvectors = n
          end if
       end do
    endif   ! master_task

    if (present(mpicom)) then
       call shr_mpi_bcast(sdat%model_calendar ,mpicom ,'calendar')
       call shr_mpi_bcast(sdat%nstreams       ,mpicom ,'nstreams')
       call shr_mpi_bcast(sdat%nvectors       ,mpicom ,'nvectors')
       call shr_mpi_bcast(sdat%stream(1)%taxmode        ,mpicom ,'taxMode')
       call shr_mpi_bcast(sdat%stream(1)%dtlimit        ,mpicom ,'dtlimit')
       call shr_mpi_bcast(sdat%pstrm(1)%stream_vectors ,mpicom ,'vectors')
       call shr_mpi_bcast(sdat%stream(1)%mapalgo        ,mpicom ,'mapalgo')
       call shr_mpi_bcast(sdat%stream(1)%mapmask        ,mpicom ,'mapmask')
       call shr_mpi_bcast(sdat%stream(1)%tinterpalgo       ,mpicom ,'tintalgo')
       call shr_mpi_bcast(sdat%stream(1)%readmode       ,mpicom ,'readmode')
    endif

    sdat%pstrm(:)%ymdLB          = -1
    sdat%pstrm(:)%todLB          = -1
    sdat%pstrm(:)%ymdUB          = -1
    sdat%pstrm(:)%todUB          = -1
    sdat%pstrm(:)%dtmin          = 1.0e30
    sdat%pstrm(:)%dtmax          = 0.0
    sdat%eccen                   = SHR_ORB_UNDEF_REAL
    sdat%mvelpp                  = SHR_ORB_UNDEF_REAL
    sdat%lambm0                  = SHR_ORB_UNDEF_REAL
    sdat%obliqr                  = SHR_ORB_UNDEF_REAL
    sdat%modeldt                 = 0
    sdat%model_calendar          = shr_cal_noleap

  end subroutine shr_strdata_readnml_from_infiles

  !===============================================================================
  subroutine shr_strdata_print(sdat,name)

    !  Print strdata common to all data models

    ! input/output parameters:
    type(shr_strdata_type)  ,intent(in) :: sdat  ! strdata data data-type
    character(len=*),optional,intent(in) :: name  ! just a name for tracking

    ! local variables
    integer   :: ns,n
    character(CL) :: lname
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

    lname = 'unknown'
    if (present(name)) then
       lname = trim(name)
    endif

    write(logunit,*)
    write(logunit,F90)
    write(logunit,F00) "name        = ",trim(lname)
    write(logunit,F00) "calendar    = ",trim(sdat%model_calendar)
    write(logunit,F01) "io_type     = ",sdat%io_type
    write(logunit,F02) "eccen       = ",sdat%eccen
    write(logunit,F02) "mvelpp      = ",sdat%mvelpp
    write(logunit,F02) "lambm0      = ",sdat%lambm0
    write(logunit,F02) "obliqr      = ",sdat%obliqr
    write(logunit,F01) "pio_iotype  = ",sdat%io_type

    write(logunit,F01) "nstreams    = ",sdat%nstreams
    do ns = 1, sdat%nstreams
       write(logunit,F04) "  taxMode (",ns,") = ",trim(sdat%stream(ns)%taxmode)
       write(logunit,F07) "  dtlimit (",ns,") = ",sdat%stream(ns)%dtlimit
       write(logunit,F06) "  domaps  (",ns,") = ",sdat%pstrm(ns)%domaps
       write(logunit,F04) "  mapalgo (",ns,") = ",trim(sdat%stream(ns)%mapalgo)
       write(logunit,F04) "  mapmask (",ns,") = ",trim(sdat%stream(ns)%mapmask)
       write(logunit,F04) "  tintalgo(",ns,") = ",trim(sdat%stream(ns)%tinterpalgo)
       write(logunit,F04) "  readmode(",ns,") = ",trim(sdat%stream(ns)%readmode)
       write(logunit,F01) " "
    end do
    write(logunit,F01) "nvectors    = ",sdat%nvectors

    do n=1, sdat%nvectors
       write(logunit,F04) "  vectors (",n,") = ",trim(sdat%pstrm(n)%stream_vectors)
    end do
    write(logunit,F90)

  end subroutine shr_strdata_print

  !===============================================================================

  subroutine shr_strdata_readLBUB(mpicom, my_task, stream, stream_mesh, &
       fldlist_stream, fldlist_model, fldbun_stream_lb, fldbun_stream_ub, &
       pio_subsystem, pio_iotype, pio_iodesc_set, pio_iodesc, &
       mDate, mSec, mDateLB, mSecLB, mDateUB, mSecUB, newData, istr, rc)

    !-------------------------------------------------------------------------
    ! Read LB and UB of stream data
    !-------------------------------------------------------------------------

    ! input/output variables
    integer                       ,intent(in)    :: mpicom
    integer                       ,intent(in)    :: my_task
    type(shr_stream_streamType)   ,intent(inout) :: stream
    type(ESMF_Mesh)               ,intent(in)    :: stream_mesh
    character(len=*)              ,intent(in)    :: fldlist_stream(:)  
    character(len=*)              ,intent(in)    :: fldlist_model(:)  
    type(ESMF_FieldBundle)        ,intent(inout) :: fldbun_stream_lb
    type(ESMF_FieldBundle)        ,intent(inout) :: fldbun_stream_ub
    type(iosystem_desc_t), target ,intent(inout) :: pio_subsystem
    integer                       ,intent(in)    :: pio_iotype
    logical                       ,intent(inout) :: pio_iodesc_set
    type(io_desc_t)               ,intent(inout) :: pio_iodesc
    integer                       ,intent(in)    :: mDate  ,mSec
    integer                       ,intent(inout) :: mDateLB,mSecLB
    integer                       ,intent(inout) :: mDateUB,mSecUB
    logical                       ,intent(out)   :: newData
    character(len=*)              ,intent(in)    :: istr
    integer                       ,intent(out)   :: rc

    ! local variables
    integer                             :: nf
    integer                             :: rCode      ! return code
    logical                             :: fileexists
    integer                             :: ivals(6)   ! bcast buffer
    integer                             :: oDateLB,oSecLB,dDateLB
    integer                             :: oDateUB,oSecUB,dDateUB
    real(r8)                            :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer                             :: n_lb, n_ub
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

    call t_startf(trim(istr)//'_setup')

    newData = .false.
    n_lb = -1
    n_ub = -1
    filename_lb = 'undefinedlb'
    filename_ub = 'undefinedub'

    oDateLB = mDateLB
    oSecLB  = mSecLB
    oDateUB = mDateUB
    oSecUB  = mSecUB

    rDateM  = real(mDate  ,r8) + real(mSec  ,r8)/shr_const_cday
    rDateLB = real(mDateLB,r8) + real(mSecLB,r8)/shr_const_cday
    rDateUB = real(mDateUB,r8) + real(mSecUB,r8)/shr_const_cday
    call t_stopf(trim(istr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(istr)//'_fbound')
       if (my_task == master_task) then ! Note that the stream bounds is only done on the master task
          call shr_stream_findBounds(stream, mDate, mSec,  &
               ivals(1), dDateLB, ivals(2), ivals(5), filename_lb, &
               ivals(3), dDateUB, ivals(4), ivals(6), filename_ub)
       endif
       call t_stopf(trim(istr)//'_fbound')

       call t_startf(trim(istr)//'_bcast')
       call shr_mpi_bcast(stream%calendar, mpicom)
       call shr_mpi_bcast(ivals, mpicom)
       call shr_mpi_bcast(filename_lb, mpicom)
       call shr_mpi_bcast(filename_ub, mpicom)
       mDateLB = ivals(1) ! Now all processors have the bounds
       mSecLB  = ivals(2)
       mDateUB = ivals(3)
       mSecUB  = ivals(4)
       n_lb    = ivals(5)
       n_ub    = ivals(6)
       call t_stopf(trim(istr)//'_bcast')
    endif

    if (mDateLB /= oDateLB .or. mSecLB /= oSecLB) then
       newdata = .true.
       if (mDateLB == oDateUB .and. mSecLB == oSecUB) then
          ! copy fldbun_stream_ub to fldbun_stream_lb
          call t_startf(trim(istr)//'_LB_copy')
          do nf = 1,size(fldlist_stream)
             call dshr_fldbun_getfldptr(fldbun_stream_ub, trim(fldlist_stream(nf)), fldptr1=dataptr_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_fldbun_getfldptr(fldbun_stream_lb, trim(fldlist_stream(nf)), fldptr1=dataptr_lb, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr_lb(:) = dataptr_ub(:)
          end do
          call t_stopf(trim(istr)//'_LB_copy')
       else
          ! read lower bound of data
          call shr_strdata_readstrm(mpicom, my_task, stream, stream_mesh, &
               fldlist_stream, fldlist_model, fldbun_stream_lb, &
               pio_subsystem, pio_iotype, pio_iodesc_set, pio_iodesc, filename_lb, n_lb, &
               istr=trim(istr)//'_LB', boundstr='lb', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    if (mDateUB /= oDateUB .or. mSecUB /= oSecUB) then
       newdata = .true.
       call shr_strdata_readstrm(mpicom, my_task, stream, stream_mesh, &
            fldlist_stream, fldlist_model, fldbun_stream_ub, &
            pio_subsystem, pio_iotype, pio_iodesc_set, pio_iodesc, filename_ub, n_ub, &
            istr=trim(istr)//'_UB', boundstr='ub', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    ! determine previous & next data files in list of files

    call t_startf(trim(istr)//'_filemgt')
    if (my_task == master_task .and. newdata) then
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
  subroutine shr_strdata_readstrm(mpicom, my_task, stream, stream_mesh, &
       fldlist_stream, fldlist_model, fldbun_stream, &
       pio_subsystem, pio_iotype, pio_iodesc_set, pio_iodesc, &
       filename, nt, istr, boundstr, rc)

    ! Read the stream data and initialize the strea pio_iodesc the first time
    ! the stream is read

    ! input/output variables
    integer                     , intent(in)            :: mpicom
    integer                     , intent(in)            :: my_task
    type(shr_stream_streamType) , intent(inout)         :: stream
    type(ESMF_Mesh)             , intent(in)            :: stream_mesh
    character(len=*)            , intent(in)            :: fldlist_stream(:)
    character(len=*)            , intent(in)            :: fldlist_model(:)
    type(ESMF_FieldBundle)      , intent(inout)         :: fldbun_stream 
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
    character(CL)                 :: currfile
    logical                       :: fileexists
    logical                       :: fileopen
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    integer                       :: nf
    integer                       :: rCode
    real(r4), allocatable         :: data_real(:)
    real(r8), pointer             :: dataptr(:)
    integer                       :: lsize
    character(*), parameter       :: subname = '(shr_strdata_readstrm) '
    character(*), parameter       :: F00   = "('(shr_strdata_readstrm) ',8a)"
    character(*), parameter       :: F02   = "('(shr_strdata_readstrm) ',2a,i8)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set up file to read from
    call t_barrierf(trim(istr)//'_BARRIER', mpicom)
    if (my_task == master_task) then
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
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
          if (my_task == master_task) write(logunit,F00) 'close  : ',trim(currfile)
          call pio_closefile(pioid)
       endif
       if (my_task == master_task) write(logunit,F00) 'opening   : ',trim(filename)
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call shr_stream_setCurrFile(stream, fileopen=.true., currfile=trim(filename), currpioid=pioid)
    endif

    ! ******************************************************************************
    ! Determine the pio io descriptor for the stream from the first data field in the stream
    ! ******************************************************************************

    if (.not. pio_iodesc_set) then
       call shr_strdata_set_stream_iodesc(pio_subsystem, pioid, trim(fldlist_stream(1)), stream_mesh, pio_iodesc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       pio_iodesc_set = .true.
    end if

    ! ******************************************************************************
    ! Read in the stream data for field names in fldname_stream_input - but fill in 
    ! the data for fldbun_stream with the field names fldname_stream_model
    ! ******************************************************************************

    call t_startf(trim(istr)//'_readpio')
    if (my_task == master_task) then
       write(logunit,F02) 'file ' // trim(boundstr) //': ',trim(filename), nt
    endif
    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    do nf = 1,size(fldlist_stream)
       rcode = pio_inq_varid(pioid, trim(fldlist_stream(nf)), varid)
       frame = nt ! set frame to time index
       call pio_setframe(pioid, varid, int(nt,kind=Pio_Offset_Kind))
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       call dshr_fldbun_getfldptr(fldbun_stream, trim(fldlist_model(nf)), fldptr1=dataptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       lsize = size(dataptr)
       allocate(data_real(lsize))
       call pio_read_darray(pioid, varid, pio_iodesc, data_real, rcode)
       if ( rcode /= PIO_NOERR ) then
          call shr_sys_abort(' ERROR: reading in variable: '// trim(fldlist_stream(nf)))
       end if
       dataptr(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    enddo
    call t_stopf(trim(istr)//'_readpio')

  end subroutine shr_strdata_readstrm

  !===============================================================================
  subroutine shr_strdata_set_stream_iodesc(pio_subsystem, pioid, fldname, stream_mesh, pio_iodesc, rc)

    ! input/output variables
    type(iosystem_desc_t) , intent(inout), target :: pio_subsystem
    type(file_desc_t)     , intent(inout)         :: pioid
    character(len=*)      , intent(in)            :: fldname
    type(ESMF_Mesh)       , intent(in)            :: stream_mesh
    type(io_desc_t)       , intent(inout)         :: pio_iodesc
    integer               , intent(out)           :: rc

    ! local variables
    integer                       :: n
    type(var_desc_t)              :: varid
    integer                       :: ndims
    integer, allocatable          :: dimids(:)
    integer, allocatable          :: dimlens(:)
    integer                       :: unlimdid
    integer                       :: itype
    type(ESMF_DistGrid)           :: distGrid
    integer                       :: lsize
    integer, pointer              :: compdof(:)
    character(CS)                 :: dimname
    integer                       :: rCode      ! pio return code
    character(*), parameter       :: subname = '(shr_strdata_set_stream_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the first field in the stream dataset
    call PIO_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, trim(fldname), varid)
    call shr_strdata_handle_error(rcode, 'SHR_STRDATA_PIO_VAR_INFO: Error inquiring varid for '//trim(fldname))
    rcode = pio_inq_varndims(pioid, varid, ndims)
    call shr_strdata_handle_error(rcode, 'SHR_STRDATA_PIO_VAR_INFO: Error inquiring ndims for '//trim(fldname))
    allocate(dimids(ndims))
    allocate(dimlens(ndims))
    rcode = pio_inq_vardimid(pioid, varid, dimids(1:ndims))
    call shr_strdata_handle_error(rcode, 'SHR_STRDATA_PIO_VAR_INFO: Error inquiring dimids for '//trim(fldname))
    do n = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
       call shr_strdata_handle_error(rcode, 'SHR_STRDATA_PIO_VAR_INFO: Error inquiring dimlens for '//trim(fldname))
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
    rcode = pio_inq_vartype(pioid, varid, itype)
    call shr_strdata_handle_error(rcode, 'SHR_STRDATA_PIO_VAR_INFO: Error inquiring type for '//trim(fldname))

    ! now create the io descriptor
    rcode = pio_inquire(pioid, unlimitedDimid=unlimdid)
    if (rcode == PIO_NOERR) then
       if (dimids(ndims) == unlimdid) then
          ! remove unlimited time dimension from the pio_init decomp
          ndims = ndims - 1
       end if
    end if
    if (ndims == 2) then
      !call pio_initdecomp(pio_subsystem, itype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       call pio_initdecomp(pio_subsystem, pio_real, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
    else if (ndims == 3) then
       rcode = pio_inq_dimname(pioid, dimids(ndims), dimname)
       if (trim(dimname) == 'time') then
          call pio_initdecomp(pio_subsystem, pio_real, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       else
          call pio_initdecomp(pio_subsystem, pio_real, (/dimlens(1),dimlens(2),dimlens(3)/), compdof, pio_iodesc)
       end if
    else
       write(6,*)'ERROR: dimlens= ',dimlens
       call shr_sys_abort(trim(subname)//' only dimlen of 2 and 3 are currently supported')
    end if

    deallocate(compdof)
    deallocate(dimids)
    deallocate(dimlens)

  end subroutine shr_strdata_set_stream_iodesc

  !===============================================================================
  subroutine shr_strdata_get_stream_pointer(sdat, strm_fld, strm_ptr, logunit, masterproc, rc)
    
    ! Set a pointer, strm_ptr, for field, strm_fld, into sdat fldbun_model field bundle

    ! input/output variables
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: strm_fld
    real(r8)               , pointer       :: strm_ptr(:)
    integer                , intent(in)    :: logunit
    logical                , intent(in)    :: masterproc
    integer                , intent(out)   :: rc 

    ! local variables
    integer :: ns, nf
    logical :: found
    character(len=*), parameter :: subname='(shr_strdata_get_stream_pointer)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! loop over all input streams and determine if the strm_fld is in the field bundle of the target stream
    do ns = 1, sdat%nstreams
       found = .false.
       ! Check if requested stream field is read in - and if it is then point into the stream field bundle
       do nf = 1,size(sdat%pstrm(ns)%fldlist_model)
          if (trim(strm_fld) == trim(sdat%pstrm(ns)%fldlist_model(nf))) then
             call dshr_fldbun_getfldptr(sdat%pstrm(ns)%fldbun_model, trim(sdat%pstrm(ns)%fldlist_model(nf)), &
                  strm_ptr, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (masterproc) then
                write(logunit,*)'(dshr_addfield_add) strm_ptr is allocated for stream field strm_'//trim(strm_fld)
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
