1module dshr_stream_mod

  ! -------------------------------------------------------------------------------
  ! Data type and methods to manage input data streams.
  ! A "data stream" is a sequence of input files where each file contains the
  ! same set of data fields and all the data fields are on the same grid.
  ! The sequence of input data files provides an uninterupted time series of
  ! data.
  !
  ! A stream data type stores information about one data stream, including the
  ! range of data date years to use and how data dates align with model dates.
  !
  ! Given a model date, this module can return data dates that are upper and
  ! lower time bounds around the given model date and the names of the files
  ! containing those dates.
  ! -------------------------------------------------------------------------------

  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod      , only : shr_sys_abort
  use shr_const_mod    , only : shr_const_cday
  use shr_string_mod   , only : shr_string_leftalign_and_convert_tabs, shr_string_parseCFtunit
  use shr_cal_mod      , only : shr_cal_noleap
  use shr_cal_mod      , only : shr_cal_date2ymd
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_cal_mod      , only : shr_cal_calendarName
  use shr_cal_mod      , only : shr_cal_advDate
  use shr_cal_mod      , only : shr_cal_advdateint
  use dshr_methods_mod , only : chkerr
  use pio              , only : file_desc_t
  use netcdf

  implicit none
  private ! default private

  ! !PUBLIC TYPES:
  public :: shr_stream_streamType        ! stream data type with private components

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: shr_stream_init_from_xml
  public :: shr_stream_init_from_inline  ! initial stream type
  public :: shr_stream_findBounds        ! return lower/upper bounding date info
  public :: shr_stream_getMeshFileName   ! return stream filename
  public :: shr_stream_getModelFieldList ! return model field name list
  public :: shr_stream_getStreamFieldList! return stream file field name list
  public :: shr_stream_getPrevFileName   ! return previous file in sequence
  public :: shr_stream_getNextFileName   ! return next file in sequence
  public :: shr_stream_getNFiles         ! get the number of files in a stream
  public :: shr_stream_getCalendar       ! get the stream calendar
  public :: shr_stream_getCurrFile       ! get the currfile, fileopen, and currpioid
  public :: shr_stream_getData           ! get stream data from target file
  public :: shr_stream_setCurrFile       ! set the currfile, fileopen, and currpioid
  public :: shr_stream_dataDump          ! internal stream data for debugging
  public :: shr_stream_restWrite         ! write a streams restart file
  public :: shr_stream_restRead          ! read  a streams restart file

  character(CS),parameter,public :: shr_stream_taxis_cycle  = 'cycle'
  character(CS),parameter,public :: shr_stream_taxis_extend = 'extend'
  character(CS),parameter,public :: shr_stream_taxis_limit  = 'limit'
  character(CS),parameter,public :: shr_stream_file_null    = 'not_set'

  ! a useful derived type to use inside shr_streamType ---
  type shr_stream_file_type
     character(CL)         :: name = shr_stream_file_null ! the file name (full pathname)
     logical               :: haveData = .false.          ! has t-coord data been read in?
     integer               :: nt = 0                      ! size of time dimension
     integer  ,allocatable :: date(:)                     ! t-coord date: yyyymmdd
     integer  ,allocatable :: secs(:)                     ! t-coord secs: elapsed on date
  end type shr_stream_file_type

  type shr_stream_data_variable
     character(CS) :: nameinfile
     character(CS) :: nameinmodel
  end type shr_stream_data_variable

  type shr_stream_streamType
     !private ! no public access to internal components
     integer           :: logunit                               ! stdout log unit
     logical           :: init         = .false.                ! has stream been initialized
     integer           :: nFiles       = 0                      ! number of data files
     integer           :: yearFirst    = -1                     ! first year to use in t-axis (yyyymmdd)
     integer           :: yearLast     = -1                     ! last  year to use in t-axis (yyyymmdd)
     integer           :: yearAlign    = -1                     ! align yearFirst with this model year
     character(CS)     :: taxMode      = shr_stream_taxis_cycle ! cycling option for time axis
     character(CS)     :: tInterpAlgo  = 'linear'               ! algorithm to use for time interpolation
     character(CS)     :: mapalgo      = 'bilinear'             ! type of mapping - default is 'bilinear'
     character(CS)     :: readMode     = 'single'               ! stream read model - 'single' or 'full_file'
     real(r8)          :: dtlimit      = 1.5_r8                 ! delta time ratio limits for time interpolation
     integer           :: offset       = 0                      ! offset in seconds of stream data
     character(CS)     :: calendar     = shr_cal_noleap         ! stream calendar (obtained from first stream data file)
     character(CL)     :: meshFile     = ' '                    ! filename for mesh for all fields on stream (full pathname)
     integer           :: k_lvd        = -1                     ! file/sample of least valid date
     integer           :: n_lvd        = -1                     ! file/sample of least valid date
     logical           :: found_lvd    = .false.                ! T <=> k_lvd,n_lvd have been set
     integer           :: k_gvd        = -1                     ! file/sample of greatest valid date
     integer           :: n_gvd        = -1                     ! file/sample of greatest valid date
     logical           :: found_gvd    = .false.                ! T <=> k_gvd,n_gvd have been set
     logical           :: fileopen     = .false.                ! is current file open
     character(CL)     :: currfile     = ' '                    ! current filename
     integer           :: nvars                                 ! number of stream variables
     type(file_desc_t) :: currpioid                             ! current pio file desc
     type(shr_stream_file_type)    , allocatable :: file(:)     ! filenames of stream data files (full pathname)
     type(shr_stream_data_variable), allocatable :: varlist(:)  ! stream variable names (on file and in model)
  end type shr_stream_streamType

  !----- parameters -----
  integer          , save      :: debug = 1            ! edit/turn-on for debug write statements
  real(R8)         , parameter :: spd = shr_const_cday ! seconds per day
  character(len=*) , parameter :: sourcefile = &
       __FILE__
  character(*)     ,parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_stream_init_from_xml(xmlfilename, streamdat, mastertask, logunit, rc)

    use FoX_DOM
    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMBroadCast, ESMF_SUCCESS

    ! ---------------------------------------------------------------------
    ! The xml format of a stream txt file will look like the following
    ! <?xml version="1.0"?>
    ! <file id="stream" version="1.0">
    !   <stream_info>
    !    <meshfile>
    !      mesh_filename
    !    </meshfile>
    !    <data_files>
    !       /glade/p/cesmdata/cseg/inputdata/atm/datm7/NYF/nyf.ncep.T62.050923.nc
    !       .....
    !    <data_files>
    !    <data_variables>
    !       u_10  u
    !    </data_variables>
    !    <stream_offset>
    !       0
    !    </stream_offset>
    !  </stream_info>
    ! </file>
    ! ---------------------------------------------------------------------

    ! input/output variables
    type(shr_stream_streamType) , intent(inout), pointer :: streamdat(:)
    character(len=*)            , intent(in)             :: xmlfilename
    logical                     , intent(in)             :: mastertask
    integer                     , intent(in)             :: logunit
    integer                     , intent(out)            :: rc

    ! local variables
    type(ESMF_VM)            :: vm
    type(Node)     , pointer :: Sdoc, p, streamnode
    type(NodeList) , pointer :: streamlist, filelist, varlist
    character(len=CL)        :: tmpstr
    integer                  :: i, n, nstrms
    integer                  :: status
    integer                  :: tmp(6)
    real(r8)                 :: rtmp(1)
    ! --------------------------------------------------------

    rc = ESMF_SUCCESS

    nstrms = 0

    if (mastertask) then

       Sdoc => parseFile(xmlfilename, iostat=status)
       if (status /= 0) then
          call shr_sys_abort("Could not open file "//trim(xmlfilename))
       endif
       streamlist => getElementsByTagname(Sdoc, "stream_info")
       nstrms = getLength(streamlist)

       ! allocate an array of shr_stream_streamtype objects on just mastertask
       allocate(streamdat(nstrms))

       ! fill in non-default values for the streamdat attributes
       do i= 1, nstrms
          streamnode => item(streamlist, i-1)

          p => item(getElementsByTagname(streamnode, "taxmode"), 0)
          if (associated(p)) then
             call extractDataContent(p, streamdat(i)%taxmode)
          endif

          p => item(getElementsByTagname(streamnode, "mapalgo"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%mapalgo)
          endif
          p => item(getElementsByTagname(streamnode, "tInterpAlgo"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%tInterpAlgo)
          endif
          p => item(getElementsByTagname(streamnode, "readMode"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%readMode)
          endif
          p=> item(getElementsByTagname(streamnode, "yearFirst"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearFirst)
          else
             call shr_sys_abort("yearFirst must be provided")
          endif
          p=> item(getElementsByTagname(streamnode, "yearLast"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearLast)
          else
             call shr_sys_abort("yearLast must be provided")
          endif
          p=> item(getElementsByTagname(streamnode, "yearAlign"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearAlign)
          else
             call shr_sys_abort("yearAlign must be provided")
          endif
          p=> item(getElementsByTagname(streamnode, "dtlimit"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%dtlimit)
          endif
          p=> item(getElementsByTagname(streamnode, "stream_offset"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%offset)
          endif
          p=> item(getElementsByTagname(streamnode, "stream_mesh_file"), 0)
          if (associated(p)) then
             call extractDataContent(p, streamdat(i)%meshfile)
          else
             call shr_sys_abort("mesh file name must be provided")
          endif
          p => item(getElementsByTagname(streamnode, "stream_data_files"), 0)
          if (.not. associated(p)) then
             call shr_sys_abort("stream data files must be provided")
          endif

          filelist => getElementsByTagname(p,"file")
          streamdat(i)%nfiles = getLength(filelist)

          allocate(streamdat(i)%file( streamdat(i)%nfiles))
          do n=1, streamdat(i)%nfiles
             p => item(filelist, n-1)
             call extractDataContent(p, streamdat(i)%file(n)%name)
          enddo

          ! Determine name of stream variables in file and model
          p => item(getElementsByTagname(streamnode, "stream_data_variables"), 0)
          varlist => getElementsByTagname(p, "var")
          streamdat(i)%nvars = getLength(varlist)
          allocate(streamdat(i)%varlist(streamdat(i)%nvars))
          do n = 1, streamdat(i)%nvars
             p => item(varlist, n-1)
             call extractDataContent(p, tmpstr)
             streamdat(i)%varlist(n)%nameinfile = tmpstr(1:index(tmpstr, " "))
             streamdat(i)%varlist(n)%nameinmodel = tmpstr(index(trim(tmpstr), " ", .true.)+1:)
          enddo

          call shr_stream_getCalendar(streamdat(i), 1, streamdat(i)%calendar)
       enddo
       call destroy(Sdoc)
    endif

    ! allocate streamdat instance on all tasks
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    tmp(1) = nstrms
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nstrms = tmp(1)
    if (.not. mastertask) then
       allocate(streamdat(nstrms))
    endif

    ! broadcast the contents of streamdat from master to all tasks
    do i=1,nstrms
       tmp(1) = streamdat(i)%nfiles
       tmp(2) = streamdat(i)%nvars
       tmp(3) = streamdat(i)%yearFirst
       tmp(4) = streamdat(i)%yearLast
       tmp(5) = streamdat(i)%yearAlign
       tmp(6) = streamdat(i)%offset
       call ESMF_VMBroadCast(vm, tmp, 6, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       streamdat(i)%nfiles    = tmp(1)
       streamdat(i)%nvars     = tmp(2)
       streamdat(i)%yearFirst = tmp(3)
       streamdat(i)%yearLast  = tmp(4)
       streamdat(i)%yearAlign = tmp(5)
       streamdat(i)%offset    = tmp(6)
       if(.not. mastertask) then
          allocate(streamdat(i)%file(streamdat(i)%nfiles))
          allocate(streamdat(i)%varlist(streamdat(i)%nvars))
       endif
       do n=1,streamdat(i)%nfiles
          call ESMF_VMBroadCast(vm, streamdat(i)%file(n)%name, CL, 0, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo
       do n=1,streamdat(i)%nvars
          call ESMF_VMBroadCast(vm, streamdat(i)%varlist(n)%nameinfile, CS, 0, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_VMBroadCast(vm, streamdat(i)%varlist(n)%nameinmodel, CS, 0, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo
       call ESMF_VMBroadCast(vm, streamdat(i)%meshfile,     CL, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMBroadCast(vm, streamdat(i)%taxmode,      CS, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMBroadCast(vm, streamdat(i)%readmode,     CS, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMBroadCast(vm, streamdat(i)%tinterpAlgo,  CS, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMBroadCast(vm, streamdat(i)%mapalgo,      CS, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rtmp(1) = streamdat(i)%dtlimit
       call ESMF_VMBroadCast(vm, rtmp, 1, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       streamdat(i)%dtlimit = rtmp(1)
    enddo

    ! Set logunit
    streamdat(:)%logunit = logunit

    ! initialize flag that stream has been set
    streamdat(:)%init = .true.

  end subroutine shr_stream_init_from_xml

  !===============================================================================

  subroutine shr_stream_init_from_inline(streamdat, stream_meshfile, &
       stream_yearFirst, stream_yearLast, stream_yearAlign, stream_offset, stream_taxmode, &
       stream_fldlistFile, stream_fldListModel, stream_fileNames, logunit)

    ! --------------------------------------------------------
    ! set values of stream datatype independent of a reading in a stream text file
    ! this is used to initialize a stream directly from fortran interface
    ! --------------------------------------------------------

    ! input/output variables
    type(shr_stream_streamType) , pointer, intent(inout) :: streamdat(:)           ! data streams (assume 1 below)
    character(*)                ,intent(in)              :: stream_meshFile        ! full pathname to stream mesh file
    integer                     ,intent(in)              :: stream_yearFirst       ! first year to use
    integer                     ,intent(in)              :: stream_yearLast        ! last  year to use
    integer                     ,intent(in)              :: stream_yearAlign       ! align yearFirst with this model year
    integer                     ,intent(in)              :: stream_offset          ! offset in seconds of stream data
    character(*)                ,intent(in)              :: stream_taxMode         ! time axis mode
    character(*)                ,intent(in)              :: stream_fldListFile(:)  ! file field names, colon delim list
    character(*)                ,intent(in)              :: stream_fldListModel(:) ! model field names, colon delim list
    character(*)                ,intent(in)              :: stream_filenames(:)    ! stream data filenames (full pathnamesa)
    integer                     ,intent(in)              :: logunit                ! stdout unit

    ! local variables
    integer                :: n
    integer                :: nfiles
    integer                :: nvars
    character(CS)          :: calendar ! stream calendar
    character(*),parameter :: subName = '(shr_stream_init_from_inline) '
    ! --------------------------------------------------------

    ! Assume only 1 stream
    allocate(streamdat(1))

    ! overwrite default values
    streamdat(1)%yearFirst    = stream_yearFirst
    streamdat(1)%yearLast     = stream_yearLast
    streamdat(1)%yearAlign    = stream_yearAlign
    streamdat(1)%offset       = stream_offset
    streamdat(1)%taxMode      = trim(stream_taxMode)
    streamdat(1)%meshFile     = trim(stream_meshFile)

    ! initialize stream filenames
    if (allocated(streamdat(1)%file)) then
       deallocate(streamdat(1)%file)
    end if
    nfiles = size(stream_filenames)
    streamdat(1)%nfiles = nfiles
    allocate(streamdat(1)%file(nfiles))
    do n = 1, nfiles
       streamdat(1)%file(n)%name = trim(stream_filenames(n))
    enddo

    ! Determine name of stream variables in file and model
    nvars = size(stream_fldlistFile)
    streamdat(1)%nvars = nvars
    allocate(streamdat(1)%varlist(nvars))
    do n = 1, nvars
       streamdat(1)%varlist(n)%nameinfile  = trim(stream_fldlistFile(n))
       streamdat(1)%varlist(n)%nameinmodel = trim(stream_fldlistModel(n))
    end do

    ! Get initial calendar value
    call shr_stream_getCalendar(streamdat(1), 1, calendar)
    streamdat(1)%calendar = trim(calendar)

    ! Initialize logunit
    streamdat(1)%logunit = logunit

    ! Initialize flag that stream has been set
    streamdat(1)%init = .true.

  end subroutine shr_stream_init_from_inline

  !===============================================================================
  subroutine shr_stream_findBounds(strm,mDateIn, secIn, &
       mDateLB, dDateLB, secLB, n_lb, fileLB,  mDateUB, dDateUB, secUB, n_ub, fileUB)

    ! Given a stream and a model date, find time coordinates of the upper and
    ! lower time bounds surrounding the models date.  Returns the model date,
    ! data date, elasped seconds, time index, and file names associated with
    ! these upper and lower time bounds.

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(inout):: strm    ! data stream to query
    integer                     ,intent(in)   :: mDateIn ! model date (yyyymmdd)
    integer                     ,intent(in)   ::   secIn ! elapsed sec on model date
    integer                     ,intent(out)  :: mDateLB ! model date    of LB
    integer                     ,intent(out)  :: dDateLB ! data  date    of LB
    integer                     ,intent(out)  ::   secLB ! elap sec      of LB
    integer                     ,intent(out)  ::    n_lb ! t-coord index of LB
    character(*)                ,intent(out)  ::  fileLB ! file containing  LB
    integer                     ,intent(out)  :: mDateUB ! model date    of UB
    integer                     ,intent(out)  :: dDateUB ! data  date    of UB
    integer                     ,intent(out)  ::   secUB ! elap sec      of UB
    integer                     ,intent(out)  ::    n_ub ! t-coord index of UB
    character(*)                ,intent(out)  ::  fileUB ! file containing  UB

    ! local variables
    integer  :: dDateIn       ! model date mapped onto a data date
    integer  :: dDateF        ! first date
    integer  :: dDateL        ! last date
    integer  :: n,nf          ! loop index wrt t-coord array within one file
    integer  :: k,kf          ! loop index wrt list of files
    integer  :: k_ub,k_lb     ! file index of U/L bounds
    integer  :: rCode         ! return code
    integer  :: mYear         ! year of model date
    integer  :: yrFirst       ! first year of data loop
    integer  :: yrLast        ! last year of data loop
    integer  :: yrAlign       ! model year that aligns with yearFirst
    integer  :: nYears        ! number of years in data loop
    integer  :: dYear         ! data year corresponding to model year
    integer  :: yy,mm,dd      ! year,month,day
    real(R8) :: rDateIn       ! model dDateIn + secs/(secs per day)
    real(R8) :: rDate1        ! stream dDateIn + secs/(secs per day)
    real(R8) :: rDate2        ! stream dDateIn + secs/(secs per day)
    real(R8) :: rDatelvd      ! lvd dDate + secs/(secs per day)
    real(R8) :: rDategvd      ! gvd dDate + secs/(secs per day)
    logical  :: cycle         ! is cycling on or off
    logical  :: limit         ! is limiting on or off
    character(*),parameter :: subName = '(shr_stream_findBounds) '
    character(*),parameter :: F00   = "('(shr_stream_findBounds) ',8a)"
    character(*),parameter :: F01   = "('(shr_stream_findBounds) ',a,i9.8,a)"
    character(*),parameter :: F02   = "('(shr_stream_findBounds) ',a,2i9.8,i6,i5,1x,a)"
    character(*),parameter :: F03   = "('(shr_stream_findBounds) ',a,i4)"
    character(*),parameter :: F04   = "('(shr_stream_findBounds) ',2a,i4)"
    !-------------------------------------------------------------------------------
    ! Purpose:
    !   1) take the model date, map it into the data date range
    !   2) find the upper and lower bounding data dates
    !   3) return the bounding data and model dates, file names, & t-coord indicies
    !-------------------------------------------------------------------------------

    if (debug>0) write(strm%logunit,F02) "DEBUG: ---------- enter ------------------"

    if ( .not. strm%init ) then
       call shr_sys_abort(trim(subName)//" ERROR: trying to find bounds of uninitialized stream")
    end if

    if (trim(strm%taxMode) == trim(shr_stream_taxis_cycle)) then
       cycle = .true.
       limit = .false.
    elseif (trim(strm%taxMode) == trim(shr_stream_taxis_extend)) then
       cycle = .false.
       limit = .false.
    elseif (trim(strm%taxMode) == trim(shr_stream_taxis_limit)) then
       cycle = .false.
       limit = .true.
    else
       write(strm%logunit,*) trim(subName),' ERROR: illegal taxMode = ',trim(strm%taxMode)
       call shr_sys_abort(trim(subName)//' ERROR: illegal taxMode = '//trim(strm%taxMode))
    endif

    !----------------------------------------------------------------------------
    ! convert/map the model year/date into a data year/date
    ! note: these values will be needed later to convert data year to model year
    !----------------------------------------------------------------------------
    mYear   = mDateIn/10000                      ! assumes/require F90 truncation
    yrFirst = strm%yearFirst                     ! first year in data sequence
    yrLast  = strm%yearLast                      ! last year in data sequence
    yrAlign = strm%yearAlign                     ! model year corresponding to yearFirst
    nYears  = yrLast - yrFirst + 1               ! number of years in data sequence
    dDateF  = yrFirst * 10000 + 101 ! first date in valid range
    dDateL  = (yrLast+1)  * 10000 + 101 ! last date in valid range

    if (cycle) then
       dYear  = yrFirst + modulo(mYear-yrAlign+(2*nYears),nYears)   ! current data year
    else
       dYear  = yrFirst + mYear - yrAlign
    endif

    if (dYear < 0) then
       write(strm%logunit,*) trim(subName),' ERROR: dyear lt zero = ',dYear
       call shr_sys_abort(trim(subName)//' ERROR: dyear lt one')
    endif

    dDateIn = dYear*10000 + modulo(mDateIn,10000) ! mDateIn mapped to range of data years
    rDateIn = dDateIn + secIn/spd                 ! dDateIn + fraction of a day
    if(debug>0) then
       write(strm%logunit,*) 'fbd1 ',mYear,dYear,dDateIn,rDateIn
       write(strm%logunit,*) 'fbd2 ',yrFirst,yrLast,yrAlign,nYears
    endif

    !----------------------------------------------------------------------------
    ! find least valid date (lvd)
    !----------------------------------------------------------------------------

    if (.not. strm%found_lvd) then
       A:    do k=1,strm%nFiles
          if (.not. strm%file(k)%haveData) then
             call shr_stream_readtCoord(strm, k, rCode)
             if ( rCode /= 0 )then
                call shr_sys_abort(trim(subName)//" ERROR: readtCoord1")
             end if
          end if
          do n=1,strm%file(k)%nt
             if ( dDateF <= strm%file(k)%date(n) ) then
                !--- found a date in or beyond yearFirst ---
                strm%k_lvd = k
                strm%n_lvd = n
                strm%found_lvd = .true.
                exit A
             end if
          end do
       end do A
       if (.not. strm%found_lvd) then
          write(strm%logunit,F00)  "ERROR: LVD not found, all data is before yearFirst"
          call shr_sys_abort(trim(subName)//" ERROR: LVD not found, all data is before yearFirst")
       else
          !--- LVD is in or beyond yearFirst, verify it is not beyond yearLast ---
          if ( dDateL <= strm%file(strm%k_lvd)%date(strm%n_lvd) ) then
             write(strm%logunit,F00)  "ERROR: LVD not found, all data is after yearLast"
             call shr_sys_abort(trim(subName)//" ERROR: LVD not found, all data is after yearLast")
          end if
       end if
       if (debug>1 ) then
          if (strm%found_lvd) write(strm%logunit,F01) "DEBUG: found LVD = ",strm%file(k)%date(n)
       end if
    end if

    if (strm%found_lvd) then
       k = strm%k_lvd
       n = strm%n_lvd
       rDatelvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! LVD date + frac day
    else
       write(strm%logunit,F00)  "ERROR: LVD not found yet"
       call shr_sys_abort(trim(subName)//" ERROR: LVD not found yet")
    endif

    if (strm%found_gvd) then
       k = strm%k_gvd
       n = strm%n_gvd
       rDategvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! GVD date + frac day
    else
       rDategvd = 99991231.0
    endif
    if(debug>0) then
       write(strm%logunit,*) 'fbd3 ',rDateIn,rDatelvd,rDategvd
    endif

    !-----------------------------------------------------------
    ! dateIn < rDatelvd
    !   limit -> abort
    !   extend -> use lvd value, set LB to 00000101
    !   cycle -> lvd is UB, gvd is LB, shift mDateLB by -nYears
    !-----------------------------------------------------------

    if (rDateIn < rDatelvd) then
       if (limit) then
          write(strm%logunit,*)  trim(subName)," ERROR: limit on and rDateIn lt rDatelvd",rDateIn,rDatelvd
          call shr_sys_abort(trim(subName)//" ERROR: rDateIn lt rDatelvd limit true")
       endif

       if (.not.cycle) then
          k_lb = strm%k_lvd
          n_lb = strm%n_lvd
          dDateLB = 00000101
          mDateLB = 00000101
          secLB   = 0
          fileLB  = strm%file(k_lb)%name

          k_ub = strm%k_lvd
          n_ub = strm%n_lvd
          dDateUB = strm%file(k_ub)%date(n_ub)
          call shr_cal_date2ymd(dDateUB,yy,mm,dd)
          yy = yy + (mYear-dYear)
          call shr_cal_ymd2date(yy,mm,dd,mDateUB)
          secUB = strm%file(k_ub)%secs(n_ub)
          fileUB = strm%file(k_ub)%name
          return
       endif

       if (cycle) then
          !--- find greatest valid date (GVD) ---
          if (.not. strm%found_gvd) then
             !--- start search at last file & move toward first file ---
             B:          do k=strm%nFiles,1,-1
                !--- read data for file number k ---
                if (.not. strm%file(k)%haveData) then
                   call shr_stream_readtCoord(strm, k, rCode)
                   if ( rCode /= 0 )then
                      call shr_sys_abort(trim(subName)//" ERROR: readtCoord2")
                   end if
                end if
                !--- start search at greatest date & move toward least date ---
                do n=strm%file(k)%nt,1,-1
                   if ( strm%file(k)%date(n) < dDateL ) then
                      strm%k_gvd = k
                      strm%n_gvd = n
                      strm%found_gvd = .true.
                      rDategvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! GVD date + frac day
                      if (debug>1 ) write(strm%logunit,F01) "DEBUG: found GVD ",strm%file(k)%date(n)
                      exit B
                   end if
                end do
             end do B
          end if

          if (.not. strm%found_gvd) then
             write(strm%logunit,F00)  "ERROR: GVD not found1"
             call shr_sys_abort(trim(subName)//" ERROR: GVD not found1")
          endif

          k_lb = strm%k_gvd
          n_lb = strm%n_gvd
          dDateLB = strm%file(k_lb)%date(n_lb)
          call shr_cal_date2ymd(dDateLB,yy,mm,dd)
          yy = yy + (mYear-dYear-nYears)
          call shr_cal_ymd2date(yy,mm,dd,mDateLB)
          secLB   = strm%file(k_lb)%secs(n_lb)
          fileLB  = strm%file(k_lb)%name

          k_ub = strm%k_lvd
          n_ub = strm%n_lvd
          dDateUB = strm%file(k_ub)%date(n_ub)
          call shr_cal_date2ymd(dDateUB,yy,mm,dd)
          yy = yy + (mYear-dYear)
          call shr_cal_ymd2date(yy,mm,dd,mDateUB)
          secUB   = strm%file(k_ub)%secs(n_ub)
          fileUB  = strm%file(k_ub)%name
          return
       endif

       !-----------------------------------------------------------
       ! dateIn > rDategvd
       !   limit -> abort
       !   extend -> use gvd value, set UB to 99991231
       !   cycle -> lvd is UB, gvd is LB, shift mDateLB by +nYears
       !-----------------------------------------------------------

    else if (strm%found_gvd .and. rDateIn >= rDategvd) then
       if (limit) then
          write(strm%logunit,*) trim(subName)," ERROR: limit on and rDateIn gt rDategvd",rDateIn,rDategvd
          call shr_sys_abort(trim(subName)//" ERROR: rDateIn gt rDategvd limit true")
       endif

       if (.not.cycle) then
          k_lb = strm%k_gvd
          n_lb = strm%n_gvd
          dDateLB = strm%file(k_lb)%date(n_lb)
          call shr_cal_date2ymd(dDateLB,yy,mm,dd)
          yy = yy + (mYear-dYear)
          call shr_cal_ymd2date(yy,mm,dd,mDateLB)
          secLB   = strm%file(k_lb)%secs(n_lb)
          fileLB  = strm%file(k_lb)%name

          k_ub = strm%k_gvd
          n_ub = strm%n_gvd
          dDateUB = 99991231
          mDateUB = 99991231
          secUB   = 0
          fileUB  = strm%file(k_ub)%name
          return
       endif

       if (cycle) then
          k_lb = strm%k_gvd
          n_lb = strm%n_gvd
          dDateLB = strm%file(k_lb)%date(n_lb)
          call shr_cal_date2ymd(dDateLB,yy,mm,dd)
          yy = yy + (mYear-dYear)
          call shr_cal_ymd2date(yy,mm,dd,mDateLB)
          secLB   = strm%file(k_lb)%secs(n_lb)
          fileLB  = strm%file(k_lb)%name

          k_ub = strm%k_lvd
          n_ub = strm%n_lvd
          dDateUB = strm%file(k_ub)%date(n_ub)
          call shr_cal_date2ymd(dDateUB,yy,mm,dd)
          yy = yy + (mYear-dYear+nYears)
          call shr_cal_ymd2date(yy,mm,dd,mDateUB)
          secUB   = strm%file(k_ub)%secs(n_ub)
          fileUB  = strm%file(k_ub)%name
          return
       endif

    else

       !-----------------------------------------------------------
       ! dateIn > rDatelvd
       !-----------------------------------------------------------
       k_lb = strm%k_lvd
       n_lb = strm%n_lvd
       C:    do k=strm%k_lvd,strm%nFiles
          !--- read data for file number k ---
          if (.not. strm%file(k)%haveData) then
             call shr_stream_readtCoord(strm, k, rCode)
             if ( rCode /= 0 )then
                call shr_sys_abort(trim(subName)//" ERROR: readtCoord3")
             end if
          end if
          !--- examine t-coords for file k ---
          n      = strm%file(k)%nt                                 ! last t-index in file
          rDate1 = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! last date + frac day

          if (.not. strm%found_gvd) then
             n = strm%file(k)%nt
             if (dDateL <= strm%file(k)%date(n)) then
                !--- set gvd to last timestep in previous file then advance through current file ---
                if (k > 1) then
                   strm%k_gvd = k-1
                   strm%n_gvd = strm%file(k-1)%nt
                   strm%found_gvd = .true.
                endif
                do n=1,strm%file(k)%nt
                   if ( strm%file(k)%date(n) < dDateL ) then
                      strm%k_gvd = k
                      strm%n_gvd = n
                      strm%found_gvd = .true.
                   endif
                enddo
             elseif (k == strm%nFiles) then
                strm%k_gvd = k
                strm%n_gvd = strm%file(k)%nt
                strm%found_gvd = .true.
             end if
             if (strm%found_gvd) then
                kf = strm%k_gvd
                nf = strm%n_gvd
                rDategvd = strm%file(kf)%date(nf) + strm%file(kf)%secs(nf)/spd ! GVD date + frac day
             endif
          end if

          !-----------------------------------------------------------
          ! dateIn > rDategvd
          !   limit -> abort
          !   extend -> use gvd value, set UB to 99991231
          !   cycle -> lvd is UB, gvd is LB, shift mDateLB by nYears
          !-----------------------------------------------------------

          if (strm%found_gvd .and. rDateIn >= rDategvd) then
             if (limit) then
                write(strm%logunit,*) trim(subName)," ERROR: limit on and rDateIn gt rDategvd",rDateIn,rDategvd
                call shr_sys_abort(trim(subName)//" ERROR: rDateIn gt rDategvd limit true")
             endif

             if (.not.cycle) then
                k_lb = strm%k_gvd
                n_lb = strm%n_gvd
                dDateLB = strm%file(k_lb)%date(n_lb)
                call shr_cal_date2ymd(dDateLB,yy,mm,dd)
                yy = yy + (mYear-dYear)
                call shr_cal_ymd2date(yy,mm,dd,mDateLB)
                secLB   = strm%file(k_lb)%secs(n_lb)
                fileLB  = strm%file(k_lb)%name

                k_ub = strm%k_gvd
                n_ub = strm%n_gvd
                dDateUB = 99991231
                mDateUB = 99991231
                secUB   = 0
                fileUB  = strm%file(k_ub)%name
                return
             endif

             if (cycle) then
                k_lb = strm%k_gvd
                n_lb = strm%n_gvd
                dDateLB = strm%file(k_lb)%date(n_lb)
                call shr_cal_date2ymd(dDateLB,yy,mm,dd)
                yy = yy + (mYear-dYear)
                call shr_cal_ymd2date(yy,mm,dd,mDateLB)
                secLB   = strm%file(k_lb)%secs(n_lb)
                fileLB  = strm%file(k_lb)%name

                k_ub = strm%k_lvd
                n_ub = strm%n_lvd
                dDateUB = strm%file(k_ub)%date(n_ub)
                call shr_cal_date2ymd(dDateUB,yy,mm,dd)
                yy = yy + (mYear-dYear+nYears)
                call shr_cal_ymd2date(yy,mm,dd,mDateUB)
                secUB   = strm%file(k_ub)%secs(n_ub)
                fileUB  = strm%file(k_ub)%name
                return
             endif

          endif

          if ( rDate1 < rDateIn ) then
             !--- increment lb and continue to search ---
             k_lb = k
             n_lb = strm%file(k)%nt
          else
             !--- the greatest lower-bound is in file k, find it ---
             do n=1,strm%file(k)%nt
                rDate2 = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! date + frac day
                if ( rDate2 <= rDateIn ) then
                   !--- found another/greater lower-bound ---
                   k_lb = k
                   n_lb = n
                else
                   !--- found the least upper-bound ---
                   k_ub = k
                   n_ub = n

                   dDateLB = strm%file(k_lb)%date(n_lb)
                   call shr_cal_date2ymd(dDateLB,yy,mm,dd)
                   yy = yy + (mYear-dYear)
                   call shr_cal_ymd2date(yy,mm,dd,mDateLB)
                   secLB = strm%file(k_lb)%secs(n_lb)
                   fileLB = strm%file(k_lb)%name

                   dDateUB = strm%file(k_ub)%date(n_ub)
                   call shr_cal_date2ymd(dDateUB,yy,mm,dd)
                   yy = yy + (mYear-dYear)
                   call shr_cal_ymd2date(yy,mm,dd,mDateUB)
                   secUB = strm%file(k_ub)%secs(n_ub)
                   fileUB = strm%file(k_ub)%name
                   return
                endif
             enddo
          endif
       end do C
    endif

    call shr_sys_abort(trim(subName)//' ERROR: findBounds failed')

  end subroutine shr_stream_findBounds

  !===============================================================================
  subroutine shr_stream_readTCoord(strm,k,rc)

    ! Read in time coordinates with possible offset (require that time coordinate is 'time')

    ! input/output parameters:
    type(shr_stream_streamType)  ,intent(inout) :: strm ! data stream to query
    integer         ,intent(in)    :: k    ! stream file index
    integer,optional,intent(out)   :: rc   ! return code

    ! local variables
    character(CL)          :: fileName    ! filename to read
    integer                :: nt
    integer                :: num,n
    integer                :: din,dout
    integer                :: sin,sout,offin
    integer                :: lrc
    integer                :: fid,vid,ndims,rcode
    integer,allocatable    :: dids(:)
    character(CS)          :: units,calendar
    character(CS)          :: bunits        ! time units (days,secs,...)
    integer                :: bdate         ! base date: calendar date
    real(R8)               :: bsec          ! base date: elapsed secs
    integer                :: ndate         ! calendar date of time value
    real(R8)               :: nsec          ! elapsed secs on calendar date
    real(R8),allocatable   :: tvar(:)
    character(*),parameter :: subname = '(shr_stream_readTCoord) '
    character(*),parameter :: F01   = "('(shr_stream_readTCoord) ',a,2i7)"
    !-------------------------------------------------------------------------------

    lrc = 0

    fileName  = trim(strm%file(k)%name)
    rCode = nf90_open(fileName, nf90_nowrite, fid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_open file '//trim(filename))
    rCode = nf90_inq_varid(fid, 'time', vid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inq_varid')
    rCode = nf90_inquire_variable(fid, vid, ndims=ndims)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_variable1')
    allocate(dids(ndims))
    rCode = nf90_inquire_variable(fid, vid, dimids=dids)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_variable2')

    ! determine number of times in file
    rCode = nf90_inquire_dimension(fid, dids(1), len=nt)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_dimension')
    deallocate(dids)

    ! allocate memory for date and secs
    allocate(strm%file(k)%date(nt), strm%file(k)%secs(nt))
    strm%file(k)%nt = nt

    ! get time units
    units = ' '
    rCode = nf90_get_att(fid, vid, 'units', units)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_att units')
    n = len_trim(units)
    if (ichar(units(n:n)) == 0 ) units(n:n) = ' '
    call shr_string_leftalign_and_convert_tabs(units)

    ! get strm%calendar
    calendar = ' '
    rCode = nf90_inquire_attribute(fid, vid, 'calendar')
    if (rCode == nf90_noerr) then
       rCode = nf90_get_att(fid, vid, 'calendar', calendar)
       if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_att calendar')
    else
       calendar = trim(shr_cal_noleap)
    endif
    n = len_trim(calendar)
    if (ichar(calendar(n:n)) == 0 ) calendar(n:n) = ' '
    call shr_string_leftalign_and_convert_tabs(calendar)
    call shr_string_parseCFtunit(units, bunits, bdate, bsec)
    strm%calendar = trim(shr_cal_calendarName(trim(calendar)))

    ! read in time coordinate values
    allocate(tvar(nt))
    rcode = nf90_get_var(fid,vid,tvar)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_var')
    rCode = nf90_close(fid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_close')

    ! determine strm%file(k)%date(n) and strm%file(k)%secs(n)
    do n = 1,nt
       call shr_cal_advDate(tvar(n), bunits, bdate, bsec, ndate, nsec, calendar)
       strm%file(k)%date(n) = ndate
       strm%file(k)%secs(n) = nint(nsec)
    enddo
    deallocate(tvar)

    ! if offset is not zero, adjust strm%file(k)%date(n) and strm%file(k)%secs(n)
    if (strm%offset /= 0) then
       if (size(strm%file(k)%date) /= size(strm%file(k)%secs)) then
          write(strm%logunit,F01) "Incompatable date and secs sizes",size(strm%file(k)%date),size(strm%file(k)%secs)
          call shr_sys_abort()
       endif
       num = size(strm%file(k)%date)
       offin = strm%offset
       do n = 1,num
          din = strm%file(k)%date(n)
          sin = strm%file(k)%secs(n)
          call shr_cal_advDateInt(offin,'seconds',din,sin,dout,sout,calendar)
          strm%file(k)%date(n) = dout
          strm%file(k)%secs(n) = sout
          ! write(strm%logunit,*) 'debug ',n,strm%offset,din,sin,dout,sout
       enddo
    endif

    ! Verify that time coordinate is valid
    strm%file(k)%haveData = .true.
    call verifyTCoord(strm, k, lrc) ! check new t-coord data

    if (present(rc)) then
       rc = lrc
    endif

  contains

    subroutine verifyTCoord(strm,k,rc)
      ! verify time coordinate data is OK

      ! !input/output parameters:
      type(shr_stream_streamType),intent(in) :: strm  ! data stream
      integer                   :: k     ! index of file to check
      integer                   :: rc    ! return code

      !----- local -----
      integer :: n           ! generic loop index
      integer :: nt          ! size of t-dimension
      integer :: date1,secs1 ! date and seconds for a    time coord
      integer :: date2,secs2 ! date and seconds for next time coord
      logical :: checkIt     ! have data / do comparison
      character(*),parameter :: subName = '(shr_stream_verifyTCoord) '
      character(*),parameter :: F00   = "('(shr_stream_verifyTCoord) ',8a)"
      character(*),parameter :: F01   = "('(shr_stream_verifyTCoord) ',a,2i7)"
      character(*),parameter :: F02   = "('(shr_stream_verifyTCoord) ',a,2i9.8)"
      !-------------------------------------------------------------------------------
      ! Notes:
      !   o checks that dates are increasing (must not decrease)
      !   o does not check for valid dates (eg. day=0 & month = 13 are "OK")
      !   o checks that secs are strictly increasing within any one day
      !   o checks that 0 <= secs <= spd (seconds per day)
      !   o checks all dates from one file plus last date of previous file and
      !     first date of next file
      !-------------------------------------------------------------------------------

      rc = 0
      if (debug>1 ) write(strm%logunit,F01) "checking t-coordinate data   for file k =",k

      if ( .not. strm%file(k)%haveData) then
         rc = 1
         write(strm%logunit,F01) "Don't have data for file ",k
         call shr_sys_abort(subName//"ERROR: can't check -- file not read.")
      end if

      do n=1,strm%file(k)%nt+1
         checkIt = .false.

         !--- do we have data for two consecutive dates? ---
         if (n==1) then
            !--- compare with previous file? ---
            if (k>1) then
               if ( strm%file(k-1)%haveData ) then
                  nt    = strm%file(k-1)%nt
                  date1 = strm%file(k-1)%date(nt)
                  secs1 = strm%file(k-1)%secs(nt)
                  date2 = strm%file(k  )%date(n)
                  secs2 = strm%file(k  )%secs(n)
                  checkIt = .true.
                  if (debug>1 ) write(strm%logunit,F01) "comparing with previous file for file k =",k
               end if
            end if
         else if (n==strm%file(k)%nt+1) then
            !--- compare with next file? ---
            if (k<strm%nFiles) then
               if ( strm%file(k+1)%haveData ) then
                  nt    = strm%file(k  )%nt
                  date1 = strm%file(k  )%date(nt)
                  secs1 = strm%file(k  )%secs(nt)
                  date2 = strm%file(k+1)%date(1)
                  secs2 = strm%file(k+1)%secs(1)
                  checkIt = .true.
                  if (debug>1 ) write(strm%logunit,F01) "comparing with next     file for file k =",k
               end if
            end if
         else
            !--- compare within this file ---
            date1 = strm%file(k)%date(n-1)
            secs1 = strm%file(k)%secs(n-1)
            date2 = strm%file(k)%date(n  )
            secs2 = strm%file(k)%secs(n  )
            checkIt = .true.
         end if

         !--- compare two consecutive dates ---
         if (checkIt) then
            if ( date1 > date2 ) then
               rc = 1
               write(strm%logunit,F01) "ERROR: calendar dates must be increasing"
               write(strm%logunit,F02) "date(n), date(n+1) = ",date1,date2
               call shr_sys_abort(subName//"ERROR: calendar dates must be increasing")
            else if ( date1 == date2 ) then
               if ( secs1 >= secs2 ) then
                  rc = 1
                  write(strm%logunit,F01) "ERROR: elapsed seconds on a date must be strickly increasing"
                  write(strm%logunit,F02) "secs(n), secs(n+1) = ",secs1,secs2
                  call shr_sys_abort(subName//"ERROR: elapsed seconds must be increasing")
               end if
            end if
            if ( secs1 < 0 .or. spd < secs1 ) then
               rc = 1
               write(strm%logunit,F01) "ERROR: elapsed seconds out of valid range [0,spd]"
               write(strm%logunit,F02) "secs(n) = ",secs1
               call shr_sys_abort(subName//"ERROR: elapsed seconds out of range")
            end if
         end if
      end do

      if (debug>0) write(strm%logunit,F01) "data is OK (non-decreasing)  for file k =",k
    end subroutine verifyTCoord

  end subroutine shr_stream_readTCoord

  !===============================================================================
  subroutine shr_stream_getMeshFileName(stream, filename)

    ! Return stream mesh filename

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    character(len=*)            ,intent(out) :: filename
    !-------------------------------------------------------------------------------

    filename = stream%meshfile

  end subroutine shr_stream_getMeshFileName

  !===============================================================================
  subroutine shr_stream_getModelFieldList(stream, list)

    ! Get list of file fields

    !input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    character(*)                ,intent(out) :: list(:)    ! field list

    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    do i=1,stream%nvars
       list(i) = stream%varlist(i)%nameinmodel
    enddo

  end subroutine shr_stream_getModelFieldList

  !===============================================================================
  subroutine shr_stream_getStreamFieldList(stream, list)

    ! Get list of file fields

    !input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    character(*)                ,intent(out) :: list(:)    ! field list
    !-------------------------------------------------------------------------------
    integer :: i

    do i=1,stream%nvars
       list(i) = stream%varlist(i)%nameinfile
    enddo

  end subroutine shr_stream_getStreamFieldList

  !===============================================================================
  subroutine shr_stream_getCalendar(strm, k, calendar)

    ! Returns calendar name

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: strm     ! data stream
    integer                     ,intent(in)  :: k        ! file to query
    character(*)                ,intent(out) :: calendar ! calendar name

    ! local
    integer                :: fid, vid, n
    character(CL)          :: fileName,lcal
    integer                :: rCode
    character(*),parameter :: subName = '(shr_stream_getCalendar) '
    !-------------------------------------------------------------------------------

    lcal = ' '
    calendar = ' '
    if (k > strm%nfiles) call shr_sys_abort(subname//' ERROR: k gt nfiles')
    fileName = trim(strm%file(k)%name)
    rCode = nf90_open(fileName,nf90_nowrite,fid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_open file '//trim(filename))
    rCode = nf90_inq_varid(fid, 'time', vid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inq_varid')
    rCode = nf90_inquire_attribute(fid, vid, 'calendar')
    if (rCode == nf90_noerr) then
       rCode = nf90_get_att(fid, vid, 'calendar', lcal)
       if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_att calendar')
    else
       lcal = trim(shr_cal_noleap)
    endif
    rCode = nf90_close(fid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_close')

    n = len_trim(lcal)
    if (ichar(lcal(n:n)) == 0 ) lcal(n:n) = ' '
    call shr_string_leftalign_and_convert_tabs(lcal)
    calendar = trim(shr_cal_calendarName(trim(lcal)))

  end subroutine shr_stream_getCalendar

  !===============================================================================
  subroutine shr_stream_getCurrFile(strm, fileopen, currfile, currpioid)

    ! returns current file information

    ! input/output parameters:
    type(shr_stream_streamType),intent(in)  :: strm      ! data stream
    logical           ,optional,intent(out) :: fileopen  ! file open flag
    character(*)      ,optional,intent(out) :: currfile  ! current filename
    type(file_desc_t) ,optional,intent(out) :: currpioid ! current pioid
    !-------------------------------------------------------------------------------

    if (present(fileopen  )) fileopen = strm%fileopen
    if (present(currfile  )) currfile = strm%currfile
    if (present(currpioid )) currpioid = strm%currpioid

  end subroutine shr_stream_getCurrFile

  !===============================================================================
  subroutine shr_stream_setCurrFile(strm, fileopen, currfile, currpioid)

    ! set current file information

    ! input/output parameters:
    type(shr_stream_streamType),intent(inout) :: strm      ! data stream
    logical           ,optional,intent(in)    :: fileopen  ! file open flag
    character(*)      ,optional,intent(in)    :: currfile  ! current filename
    type(file_desc_t) ,optional,intent(in)    :: currpioid ! current pioid
    !-------------------------------------------------------------------------------

    if (present(fileopen  )) strm%fileopen = fileopen
    if (present(currfile  )) strm%currfile = currfile
    if (present(currpioid )) strm%currpioid = currpioid

  end subroutine shr_stream_setCurrFile

  !===============================================================================
  subroutine shr_stream_getNextFileName(strm, fn, fnNext,rc)

    ! Returns next file name in sequence
    ! Note: will wrap-around data loop if lvd & gvd are known
    ! otherwise may return file name = "unknown"

    ! !input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: strm   ! data stream
    character(*)                ,intent(in)  :: fn     ! file name
    character(*)                ,intent(out) :: fnNext ! next file name
    integer      ,optional      ,intent(out) :: rc     ! return code

    ! local variables
    integer                :: rCode   ! return code
    integer                :: n      ! loop index
    logical                :: found  ! file name found?
    character(*),parameter :: subName = '(shr_stream_getNextFileName) '
    character(*),parameter :: F00   = "('(shr_stream_getNextFileName) ',8a)"
    !-------------------------------------------------------------------------------

    rCode = 0

    !--- locate input file in the stream's list of files ---
    found = .false.
    do n = 1,strm%nFiles
       if ( trim(fn) == trim(strm%file(n)%name)) then
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       rCode = 1
       write(strm%logunit,F00) "ERROR: input file name is not in stream: ",trim(fn)
       call shr_sys_abort(subName//"ERROR: file name not in stream: "//trim(fn))
    end if

    !--- get next file name ---
    n = n+1  ! next in list
    if (strm%found_lvd .and. strm%found_gvd) then
       if (n > strm%k_gvd)  n = strm%k_lvd ! wrap-around to lvd
    else if (strm%found_lvd ) then
       if (n > strm%nFiles) n = strm%k_lvd ! wrap-around to lvd
    else if (n > strm%nFiles ) then
       n = 1                               ! wrap-around to 1st file
    end if

    fnNext = trim(strm%file(n)%name)
    if ( present(rc) ) rc = rCode

  end subroutine shr_stream_getNextFileName

  !===============================================================================
  subroutine shr_stream_getPrevFileName(strm, fn, fnPrev,rc)

    ! Returns previous file name in sequence

    ! !input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: strm   ! data stream
    character(*)                ,intent(in)  :: fn     ! file name
    character(*)                ,intent(out) :: fnPrev ! preciding file name
    integer      ,optional      ,intent(out) :: rc     ! return code

    !--- local ---
    integer                :: rCode ! return code
    integer                :: n     ! loop index
    logical                :: found ! file name found?
    character(*),parameter :: subName = '(shr_stream_getPrevFileName) '
    character(*),parameter :: F00   = "('(shr_stream_getPrevFileName) ',8a)"

    !-------------------------------------------------------------------------------
    ! Note: will wrap-around data loop if lvd & gvd are known
    ! otherwise may return file name = "unknown"
    !-------------------------------------------------------------------------------

    rCode = 0

    !--- locate input file in the stream's list of files ---
    found = .false.
    do n = 1,strm%nFiles
       if ( trim(fn) == trim(strm%file(n)%name)) then
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       rCode = 1
       write(strm%logunit,F00) "ERROR: input file name is not in stream: ",trim(fn)
       call shr_sys_abort(subName//"ERROR: file name not in stream: "//trim(fn))
    end if

    !--- get previous file name ---
    n = n-1  ! previous in list
    if (strm%found_lvd .and. strm%found_gvd) then
       if ( n < strm%k_lvd) n = strm%k_gvd ! do wrap-around ---
    end if
    if (n>0) then
       fnPrev = trim(strm%file(n)%name)
    else
       fnPrev = "unknown "
    end if
    if ( present(rc) ) rc = rCode

  end subroutine shr_stream_getPrevFileName

  !===============================================================================
  subroutine shr_stream_getNFiles(strm,nfiles)

    ! Returns number of input files in stream

    ! input/output parameters:
    type(shr_stream_streamType),intent(in)  :: strm      ! data stream
    integer                    ,intent(out) :: nfiles    ! number of input files in stream
    !-------------------------------------------------------------------------------

    nfiles = strm%nfiles

  end subroutine shr_stream_getNFiles

  !===============================================================================
  subroutine shr_stream_restWrite(strm, fileName, caseName, nstrms)

    ! Write stream data to a restart file.

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: strm(:)    ! data streams
    character(*)                ,intent(in)  :: fileName   ! name of restart file
    character(*)                ,intent(in)  :: caseName   ! case name
    integer, optional           ,intent(in)  :: nstrms     ! number of streams

    ! local variables
    integer       :: nStreams    ! number of streams
    integer       :: k,n         ! generic loop index
    character( 8) :: dStr        ! F90 wall clock date str yyyymmdd
    character(10) :: tStr        ! F90 wall clock time str hhmmss.sss
    character(CS) :: str         ! generic text string
    integer       :: nUnit       ! a file unit number
    integer       :: nt          ! number of time samples
    integer       :: logunit
    character(*),parameter :: subName = '(shr_stream_restWrite) '
    character(*),parameter :: F00   = "('(shr_stream_restWrite) ',16a) "
    character(*),parameter :: F01   = "('(shr_stream_restWrite) ',a,i5,a,5a) "
    character(*),parameter :: F02   = "('(shr_stream_restWrite) ',a,i5,a,5i8) "
    character(*),parameter :: F03   = "('(shr_stream_restWrite) ',a,i5,a,5l3) "
    !-------------------------------------------------------------------------------

    ! stdout unit
    logunit = strm(1)%logunit ! all streams have the same logunit

    ! error checks
    if (present(nstrms)) then
       if (size(strm) < nstrms) then
          write(logunit,F02) "ERROR: nstrms too large for strm",size(strm),nstrms
          call shr_sys_abort(subname//": ERROR: nstrms too large for strm")
       endif
       nStreams = nstrms
    else
       nStreams = size(strm)
    endif
    do k = 1,nStreams
       if (.not. strm(k)%init) then  ! has stream been initialized?
          write(logunit,F01) "ERROR: can't write uninitialized stream to a restart file, k = ",k
          call shr_sys_abort(subName//": ERROR: given uninitialized stream")
       end if
    end do

    ! write restart data
    open(nUnit,file=trim(fileName),form="unformatted",action="write")
    write(nUnit) nStreams
    do k = 1,nStreams
       write(nUnit) strm(k)%init                  ! has stream been initialized?
       write(nUnit) strm(k)%nFiles                ! number of data files for this stream
       do n = 1,strm(k)%nFiles                    ! data specific to each stream file...
          write(nUnit) strm(k)%file(n)%name       ! the stream file name
          write(nUnit) strm(k)%file(n)%haveData   ! has stream t-coord data been read in?
          write(nUnit) strm(k)%file(n)%nt         ! size of stream time dimension
          if (strm(k)%file(n)%haveData) then      ! ie. if arrays have been allocated
             write(nUnit) strm(k)%file(n)%date(:) ! t-coord date: yyyymmdd
             write(nUnit) strm(k)%file(n)%secs(:) ! t-coord secs: elapsed on date
          end if
       end do
       write(nUnit) strm(k)%offset                ! time axis offset
       write(nUnit) strm(k)%k_lvd                 ! file        of least valid date
       write(nUnit) strm(k)%n_lvd                 !      sample of least valid date
       write(nUnit) strm(k)%found_lvd             ! T <=> k_lvd,n_lvd have been set
       write(nUnit) strm(k)%k_gvd                 ! file        of greatest valid date
       write(nUnit) strm(k)%n_gvd                 !      sample of greatest valid date
       write(nUnit) strm(k)%found_gvd             ! T <=> k_gvd,n_gvd have been set
    end do
    close(nUnit)

    ! write diagnostic log output
    call date_and_time(dStr,tStr)
    write(logunit,F00) "stream restart file created : ",&
         dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '//tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
    write(logunit,F01) "number of streams ",nStreams
    do k = 1,nStreams
       nt = strm(k)%file(1)%nt
       write(logunit,F02) "* stream ",k," offset               = ",strm(k)%offset
       write(logunit,F03) "* stream ",k," init                 = ",strm(k)%init
       write(logunit,F01) "* stream ",k," first file name      = ",trim(strm(k)%file(1)%name)
       write(logunit,F03) "* stream ",k," first file have data = ",strm(k)%file(1)%haveData
       write(logunit,F02) "* stream ",k," first file ntimes    = ",strm(k)%file(1)%nt
       if (strm(k)%file(1)%haveData .and. debug > 0) then
          write(logunit,F02) "* stream ",k," first file date secs = ", strm(k)%file(1)%date(1), strm(k)%file(1)%secs(1)
          write(logunit,F02) "* stream ",k," last  file date secs = ", strm(k)%file(1)%date(nt), strm(k)%file(1)%secs(nt)
       endif
    end do

  end subroutine  shr_stream_restWrite

  !===============================================================================
  subroutine shr_stream_restRead(strm, fileName, nstrms)

    ! read stream data from a restart file

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(inout) :: strm(:)  ! vector of data streams
    character(*)                ,intent(in)    :: fileName ! name of restart file
    integer                     ,intent(in)    :: nstrms   ! number of streams in strm

    ! local variables
    integer         :: rCode                          ! return code
    integer         :: nStreams                       ! number of streams
    integer         :: k,n                            ! generic loop index
    character(CS)   :: str                            ! generic text string
    integer         :: nUnit                          ! a file unit number
    integer         :: offset_input                   ! input integer
    integer         :: nt                             ! size of time dimension
    character(CS)   :: tInterpAlgo                    ! for backwards compatability
    character(CL)   :: name                           ! local variables
    integer         :: nFiles                         ! local variables
    integer         :: k_lvd, n_lvd, k_gvd, n_gvd     ! local variables
    logical         :: found_lvd, found_gvd, haveData ! local variables
    integer,pointer :: date(:),secs(:)                ! local variables
    logical         :: readok                         ! read of restarts ok
    integer         :: logunit
    character(*),parameter :: subName = '(shr_stream_restRead) '
    character(*),parameter :: F00   = "('(shr_stream_restRead) ',16a) "
    character(*),parameter :: F01   = "('(shr_stream_restRead) ',a,i5,a,5a) "
    character(*),parameter :: F02   = "('(shr_stream_restRead) ',a,i5,a,5i8) "
    character(*),parameter :: F03   = "('(shr_stream_restRead) ',a,i5,a,5l3) "
    character(*),parameter :: F04   = "('(shr_stream_restRead) ',a,4i8) "
    character(*),parameter :: F05   = "('(shr_stream_restRead) ',a,2i8,6a) "
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! read the  data
    !----------------------------------------------------------------------------

    logunit = strm(1)%logunit ! all streams have the same logunit

    open(newunit=nUnit,file=trim(fileName),form="unformatted",status="old",action="read", iostat=rCode)
    write(logunit,F01)'reading stream restart info from '//trim(filename)
    if ( rCode /= 0 )then
       call shr_sys_abort(subName//": ERROR: error opening file: "//trim(fileName) )
    end if

    ! number of streams
    read(nUnit) nStreams
    if (nstrms /= nStreams) then
       write(logunit,F02) "ERROR: nstrms ne nStreams on restart",nstrms,' ',nStreams
       call shr_sys_abort(subname//": ERROR: nstrms ne nStreams on restart")
    endif
    nStreams = nstrms
    write(logunit,F01) "Number of streams on restart ",nStreams

    ! loop over streams
    do k = 1,nStreams
       ! has stream been initialized?
       read(nUnit) strm(k)%init
       if (.not. strm(k)%init) then
          write(logunit,F01) "ERROR: uninitialized stream in restart file, k = ",k
          call shr_sys_abort(subName//": ERROR: reading uninitialized stream")
       end if
       readok = .true.

       ! don't overwrite these from input - make local variables fo rinput
       read(nUnit) nFiles
       write(logunit,F02) "Number of files on stream ",n," is ",nStreams
       do n = 1,nFiles
          read(nUnit) name       ! the file name
          read(nUnit) haveData   ! has t-coord data been read in?
          read(nUnit) nt         ! size of time dimension
          if (haveData) then     ! ie. if arrays have been allocated
             allocate(date(nt))
             allocate(secs(nt))
             read(nUnit) date(:) ! t-coord date: yyyymmdd
             read(nUnit) secs(:) ! t-coord secs: elapsed on date
             if (strm(k)%nFiles >= n) then
                if (trim(name) == trim(strm(k)%file(n)%name)) then
                   allocate(strm(k)%file(n)%date(nt))
                   allocate(strm(k)%file(n)%secs(nt))
                   write(logunit,F00) "reading time data for stream file ",trim(name)
                   strm(k)%file(n)%nt = nt
                   strm(k)%file(n)%haveData = haveData
                   strm(k)%file(n)%date(1:nt) = date(1:nt)
                   strm(k)%file(n)%secs(1:nt) = secs(1:nt)
                else
                   write(logunit,F05) "WARNING, skipping reading in time data for stream file ",trim(name)
                   readok = .false.
                endif  ! filenames consistent
             endif  ! strm nfiles
             deallocate(date)
             deallocate(secs)
          end if
       end do
       read(nUnit) offset_input ! time axis offset
       read(nUnit) k_lvd        ! file of least valid date
       read(nUnit) n_lvd        ! sample of least valid date
       read(nUnit) found_lvd    ! T <=> k_lvd,n_lvd have been set
       read(nUnit) k_gvd        ! file of greatest valid date
       read(nUnit) n_gvd        ! sample of greatest valid date
       read(nUnit) found_gvd    ! T <=> k_gvd,n_gvd have been set

       if (offset_input /= strm(k)%offset) then
          write(logunit,F04) " ERROR: offset disagrees ",k,strm(k)%offset,offset_input
          write(logunit,F00) "ERRORS Detected ABORTING NOW"
          call shr_sys_abort(subName//": ERRORS Detected ABORTING NOW")
       endif
       if (readok) then
          ! only overwrite if restart read is ok
          write(logunit,F05) "setting k n and found lvd gvd on restart ",k,n,' ',trim(name)
          strm(k)%k_lvd     = k_lvd
          strm(k)%n_lvd     = n_lvd
          strm(k)%found_lvd = found_lvd
          strm(k)%k_gvd     = k_gvd
          strm(k)%n_gvd     = n_gvd
          strm(k)%found_gvd = found_gvd
       endif

       if (debug > 0) then
          write(logunit,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
          write(logunit,F03) "* stream ",k," first have data = ",strm(k)%file(1)%haveData
          write(logunit,F02) "* stream ",k," first nt        = ",strm(k)%file(1)%nt
          if (strm(k)%file(1)%haveData) then
             nt = strm(k)%file(1)%nt
             write(logunit,F02) "* stream ",k," first date secs = ", strm(k)%file(1)%date(1),strm(k)%file(1)%secs(1)
             write(logunit,F02) "* stream ",k," last  date secs = ", strm(k)%file(1)%date(nt),strm(k)%file(1)%secs(nt)
          end if
       endif

    end do
    close(nUnit)

  end subroutine shr_stream_restRead

  !===============================================================================
  subroutine shr_stream_dataDump(strm)

    ! Dump all data to stdout for debugging

    ! input/output parameters:
    type(shr_stream_streamType),intent(in) :: strm      ! data stream

    !----- local -----
    integer   :: k ! generic loop index
    integer   :: logunit
    character(*),parameter :: F00   = "('(shr_stream_dataDump) ',8a)"
    character(*),parameter :: F01   = "('(shr_stream_dataDump) ',a,3i5)"
    character(*),parameter :: F02   = "('(shr_stream_dataDump) ',a,365i9.8)"
    character(*),parameter :: F03   = "('(shr_stream_dataDump) ',a,365i6)"
    !-------------------------------------------------------------------------------

    logunit = strm%logunit

    if (debug > 0) then
       write(logunit,F00) "dump internal data for debugging..."
       write(logunit,F01) "nFiles        = ", strm%nFiles
       do k=1,strm%nFiles
          write(logunit,F01) "data for file k = ",k
          write(logunit,F00)    "* file(k)%name    = ", trim(strm%file(k)%name)
          if ( strm%file(k)%haveData ) then
             write(logunit,F01) "* file(k)%nt      = ", strm%file(k)%nt
             write(logunit,F02) "* file(k)%date(:) = ", strm%file(k)%date(:)
             write(logunit,F03) "* file(k)%Secs(:) = ", strm%file(k)%secs(:)
          else
             write(logunit,F00) "* time coord data not read in yet for this file"
          end if
       end do
       write(logunit,F01) "yearF/L/A    = ", strm%yearFirst,strm%yearLast,strm%yearAlign
       write(logunit,F01) "offset       = ", strm%offset
       write(logunit,F00) "taxMode      = ", trim(strm%taxMode)
       write(logunit,F00) "meshfile     = ", trim(strm%meshfile)
    end if

  end subroutine shr_stream_dataDump

  !===============================================================================
  subroutine shr_stream_getData(stream, index, filename)

    ! Returns full pathname of stream data file (nt)

    type(shr_stream_streamType) , intent(in)  :: stream
    integer                     , intent(in)  :: index
    character(len=*)            , intent(out) :: filename
    !-------------------------------------------------------------------------------

    filename = trim(stream%file(index)%name)

  end subroutine shr_stream_getData

end module dshr_stream_mod
