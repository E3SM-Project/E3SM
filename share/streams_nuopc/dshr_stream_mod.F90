module dshr_stream_mod

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
  !
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
  ! -------------------------------------------------------------------------------

  use shr_kind_mod   , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod    , only : shr_sys_abort
  use shr_const_mod  , only : shr_const_cday
  use shr_string_mod , only : shr_string_leftalign_and_convert_tabs
  use shr_string_mod , only : shr_string_parseCFtunit
  use shr_string_mod , only : shr_string_listGetName
  use shr_cal_mod    , only : shr_cal_noleap
  use shr_cal_mod    , only : shr_cal_date2ymd
  use shr_cal_mod    , only : shr_cal_ymd2date
  use shr_cal_mod    , only : shr_cal_calendarName
  use shr_cal_mod    , only : shr_cal_advDate
  use shr_cal_mod    , only : shr_cal_advdateint
  use shr_log_mod    , only : s_logunit  => shr_log_Unit
  use shr_log_mod    , only : OOBMsg => shr_log_OOBMsg
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use pio            , only : file_desc_t
  use netcdf
  use perf_mod

  implicit none
  private ! default private

  ! !PUBLIC TYPES:
  public :: shr_stream_streamType        ! stream data type with private components
  public :: shr_stream_fileType

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: shr_stream_init_from_xml
  public :: shr_stream_init_from_fortran ! initial stream type
  public :: shr_stream_default           ! set default values
  public :: shr_stream_parseInput        ! extract fileName,yearAlign, etc. from a string
  public :: shr_stream_findBounds        ! return lower/upper bounding date info
  public :: shr_stream_getMeshFileName   ! return stream filename
  public :: shr_stream_getModelFieldList ! return model field name list
  public :: shr_stream_getStreamFieldList! return stream file field name list
  public :: shr_stream_getFileFieldName  ! return k-th input-file field name
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
  type shr_stream_fileType
     character(CL)         :: name = shr_stream_file_null ! the file name (full pathname)
     logical               :: haveData = .false.          ! has t-coord data been read in?
     integer               :: nt = 0                      ! size of time dimension
     integer  ,allocatable :: date(:)                     ! t-coord date: yyyymmdd
     integer  ,allocatable :: secs(:)                     ! t-coord secs: elapsed on date
  end type shr_stream_fileType

  type shr_stream_streamType
     !private                                          ! no public access to internal components
     !--- input data file names and data ---
     logical           :: init                         ! has stream been initialized?
     integer  ,pointer :: initarr(:) => null()         ! surrogate for init flag
     integer           :: nFiles                       ! number of data files
     integer           :: nvars

     ! data used for time interpolation  - obtained from shr_strdata namelist
     integer           :: yearFirst                    ! first year to use in t-axis (yyyymmdd)
     integer           :: yearLast                     ! last  year to use in t-axis (yyyymmdd)
     integer           :: yearAlign                    ! align yearFirst with this model year
     character(CS)     :: taxMode                      ! cycling option for time axis
     character(CS)     :: tInterpAlgo                  ! Algorithm to use for time interpolation
     character(CS)     :: readMode
     real(r8)          :: dtlimit
     character(CS)     :: mapalgo
     character(CS)     :: mapmask

     ! data used for time interpolation - obtained from stream txt file
     integer           :: offset                       ! offset in seconds of stream data

     ! data used for time interpolation - obtained by reading first stream data file
     character(CS)     :: calendar                     ! stream calendar

     ! stream data metadata - obtained from stream txt file
     character(CL)     :: meshFileName                 ! filename for mesh for all fields on stream (full pathname)
     type(shr_stream_fileType), allocatable :: file(:) ! filenames of stream data files (full pathname)
     character(CXX)    :: fldListFile                  ! field list: file's  field names
     character(CXX)    :: fldListModel                 ! field list: model's field names

     ! useful for quicker searching ---
     integer           :: k_lvd,n_lvd                  ! file/sample of least valid date
     logical           :: found_lvd                    ! T <=> k_lvd,n_lvd have been set
     integer           :: k_gvd,n_gvd                  ! file/sample of greatest valid date
     logical           :: found_gvd                    ! T <=> k_gvd,n_gvd have been set

     ! for keeping files open
     logical           :: fileopen                     ! is current file open
     character(CL)     :: currfile                     ! current filename
     type(file_desc_t) :: currpioid                    ! current pio file desc
     type(stream_data_variable), allocatable :: varlist(:)
  end type shr_stream_streamType

  type stream_data_variable
     character(CS) :: nameinfile
     character(CS) :: nameinmodel
  end type stream_data_variable


  !----- parameters -----
  real(R8)         , parameter :: spd = shr_const_cday ! seconds per day
  integer          , parameter :: initarr_size = 3     ! size of initarr
  integer          , save      :: debug = 0            ! edit/turn-on for debug write statements
  character(len=*) , parameter :: sourcefile = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  ! Complete the dynamic vector definition.

  subroutine shr_stream_init_from_xml(xmlfilename, streamdat, mastertask, rc)
    use FoX_DOM
    type(shr_stream_streamType), intent(inout), pointer :: streamdat(:)
    character(len=*), intent(in) :: xmlfilename
    logical, intent(in) :: mastertask
    integer, intent(out) :: rc

    type(Node), pointer :: Sdoc, p, streamnode
    type(NodeList), pointer :: streamlist, filelist, varlist
    character(len=CL) :: tmpstr
    integer :: i, n, nstrms
    nstrms = 0

    if(mastertask) then
       Sdoc => parseFile(xmlfilename, iostat=rc)
       if(rc /= 0) then
          call shr_sys_abort("Could not open file "//trim(xmlfilename))
       endif
       streamlist => getElementsByTagname(Sdoc, "stream_info")
       nstrms = getLength(streamlist)
       allocate(streamdat(nstrms))
       do i=1, nstrms
          streamnode => item(streamlist, i-1)
          p => item(getElementsByTagname(streamnode, "taxmode"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%taxmode)
          else
             streamdat(i)%taxmode = "cycle"
          endif
          p => item(getElementsByTagname(streamnode, "mapalgo"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%mapalgo)
          else
             streamdat(i)%mapalgo = "bilinear"
          endif
          p => item(getElementsByTagname(streamnode, "mapmask"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%mapmask)
          else
             streamdat(i)%taxmode = "dstmask"
          endif

          p => item(getElementsByTagname(streamnode, "tInterpAlgo"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%tInterpAlgo)
          else
             streamdat(i)%tInterpAlgo = "unused"
          endif
          p => item(getElementsByTagname(streamnode, "readMode"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%readMode)
          else
             streamdat(i)%readMode = "single"
          endif

          p=> item(getElementsByTagname(streamnode, "yearFirst"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearFirst)
          else
             streamdat(i)%yearFirst = 1
          endif
          p=> item(getElementsByTagname(streamnode, "yearLast"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearLast)
          else
             streamdat(i)%yearLast = 30 ! What should this be by default?
          endif
          call extractDataContent(p, streamdat(i)%yearLast)

          p=> item(getElementsByTagname(streamnode, "yearAlign"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%yearAlign)
          else
             streamdat(i)%yearAlign = 1
          endif

          p=> item(getElementsByTagname(streamnode, "dtlimit"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%dtlimit)
          else
             streamdat(i)%dtlimit = 1.0e30
          endif

          p=> item(getElementsByTagname(streamnode, "stream_offset"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%offset)
          else
             streamdat(i)%offset = 0
          endif

          p=> item(getElementsByTagname(streamnode, "stream_mesh_file"), 0)
          if(associated(p)) then
             call extractDataContent(p, streamdat(i)%meshfilename)
          else
             call shr_sys_abort("mesh file name must be provided")
          endif

          p => item(getElementsByTagname(streamnode, "stream_data_files"), 0)
          if(.not. associated(p)) then
             call shr_sys_abort("stream data files must be provided")
          endif
          filelist => getElementsByTagname(p,"file")
          streamdat(i)%nfiles = getLength(filelist)

          allocate(streamdat(i)%file( streamdat(i)%nfiles))
          do n=1, streamdat(i)%nfiles
             p => item(filelist, n-1)
             call extractDataContent(p, streamdat(i)%file(n)%name)
          enddo

          p => item(getElementsByTagname(streamnode, "stream_data_variables"), 0)
          varlist => getElementsByTagname(p, "var")
          streamdat(i)%nvars = getLength(varlist)
          allocate(streamdat(i)%varlist(streamdat(i)%nvars))
          do n=1, streamdat(i)%nvars
             p => item(varlist, n-1)
             call extractDataContent(p, tmpstr)
             streamdat(i)%varlist(n)%nameinfile = tmpstr(1:index(tmpstr, " "))
             streamdat(i)%varlist(n)%nameinmodel = tmpstr(index(trim(tmpstr), " ", .true.)+1:)
          enddo

          call shr_stream_getCalendar(streamdat(i), 1, streamdat(i)%calendar)


       enddo
       call destroy(Sdoc)
    endif
    call broadcast_streamdata(nstrms, streamdat, mastertask)

    ! initialize flag that stream has been set
    do i=1, nstrms
       call shr_stream_setInit(streamdat(i))
    enddo

  end subroutine shr_stream_init_from_xml

  subroutine broadcast_streamdata(nstrms, streamdat, mastertask)
    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMBroadCast
    integer, intent(inout) :: nstrms
    type(shr_stream_streamType), intent(inout), pointer :: streamdat(:)
    logical, intent(in) :: mastertask
    ! broadcast the contents of streamdat from master to all tasks
    type(ESMF_VM) :: vm
    integer :: tmp(6)
    integer :: i
    integer :: n
    integer :: rc
    real(r8) :: rtmp(1)

    call ESMF_VMGetCurrent(vm, rc=rc)
    tmp(1) = nstrms
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    nstrms = tmp(1)

    if(.not. mastertask) then
       allocate(streamdat(nstrms))
    endif
    do i=1,nstrms
       tmp(1) = streamdat(i)%nfiles
       tmp(2) = streamdat(i)%nvars
       tmp(3) = streamdat(i)%yearFirst
       tmp(4) = streamdat(i)%yearLast
       tmp(5) = streamdat(i)%yearAlign
       tmp(6) = streamdat(i)%offset
       call ESMF_VMBroadCast(vm, tmp, 6, 0, rc=rc)
       streamdat(i)%nfiles = tmp(1)
       streamdat(i)%nvars = tmp(2)
       streamdat(i)%yearFirst = tmp(3)
       streamdat(i)%yearLast = tmp(4)
       streamdat(i)%yearAlign = tmp(5)
       streamdat(i)%offset = tmp(6)
       if(.not. mastertask) then
          allocate(streamdat(i)%file(streamdat(i)%nfiles))
          allocate(streamdat(i)%varlist(streamdat(i)%nvars))
       endif
       do n=1,streamdat(i)%nfiles
          call ESMF_VMBroadCast(vm, streamdat(i)%file(n)%name, CL, 0, rc=rc)
       enddo
       do n=1,streamdat(i)%nvars
          call ESMF_VMBroadCast(vm, streamdat(i)%varlist(n)%nameinfile, CS, 0, rc=rc)
          call ESMF_VMBroadCast(vm, streamdat(i)%varlist(n)%nameinmodel, CS, 0, rc=rc)
       enddo

       call ESMF_VMBroadCast(vm, streamdat(i)%meshfilename, CL, 0, rc=rc)
       call ESMF_VMBroadCast(vm, streamdat(i)%taxmode,      CS, 0, rc=rc)
       call ESMF_VMBroadCast(vm, streamdat(i)%readmode,     CS, 0, rc=rc)
       call ESMF_VMBroadCast(vm, streamdat(i)%tinterpAlgo,  CS, 0, rc=rc)
       call ESMF_VMBroadCast(vm, streamdat(i)%mapalgo,      CS, 0, rc=rc)
       call ESMF_VMBroadCast(vm, streamdat(i)%mapmask,      CS, 0, rc=rc)
       rtmp(1) = streamdat(i)%dtlimit
       call ESMF_VMBroadCast(vm, rtmp, 1, 0, rc=rc)
       streamdat(i)%dtlimit = rtmp(1)

    enddo
  end subroutine broadcast_streamdata


  !===============================================================================

  subroutine shr_stream_init_from_fortran(strm, meshfile, &
       yearFirst, yearLast, yearAlign, offset, taxmode, &
       fldlistFile, fldListModel, fileNames)

    ! --------------------------------------------------------
    ! set values of stream datatype independent of a reading in a stream text file
    ! this is used to initialize a stream directly from fortran interface
    ! --------------------------------------------------------

    ! input/output variables
    type(shr_stream_streamType) ,intent(inout) :: strm         ! data stream
    character(*)                ,intent(in)    :: meshFile     ! full pathname to stream mesh file
    integer                     ,intent(in)    :: yearFirst    ! first year to use
    integer                     ,intent(in)    :: yearLast     ! last  year to use
    integer                     ,intent(in)    :: yearAlign    ! align yearFirst with this model year
    integer                     ,intent(in)    :: offset       ! offset in seconds of stream data
    character(*)                ,intent(in)    :: taxMode      ! time axis mode
    character(*)                ,intent(in)    :: fldListFile  ! file field names, colon delim list
    character(*)                ,intent(in)    :: fldListModel ! model field names, colon delim list
    character(*)                ,intent(in)    :: filenames(:) ! stream data filenames (full pathnamesa)

    ! local variables
    integer                   :: n
    character(CS)             :: calendar ! stream calendar
    type(shr_stream_fileType) :: tempFile ! File being constructed.
!    type(fileVector)          :: fileVec  ! Vector used to construct file array.
    character(*),parameter    :: subName = '(shr_stream_init_from_fortran) '
    ! --------------------------------------------------------

    ! set default values for stream
    call shr_stream_default(strm)

    ! overwrite default values
    strm%yearFirst    = yearFirst
    strm%yearLast     = yearLast
    strm%yearAlign    = yearAlign
    strm%offset       = offset
    strm%taxMode      = trim(taxMode)
    strm%meshFileName = trim(meshFile)
    strm%fldListFile  = trim(fldListFile)
    strm%fldListModel = trim(fldListModel)

    ! create a linked list of data files in stream
    do n = 1,size(filenames)
       ! ignore null file names.
       if (trim(filenames(n)) /= trim(shr_stream_file_null)) then
          tempFile%name = trim(filenames(n))
!          call fileVec%push_back(tempFile)
       endif
    enddo

    ! True size after throwing out null names.
!    strm%nFiles = fileVec%vsize()
!    call fileVec%move_out(strm%file)

    ! get initial calendar value
    call shr_stream_getCalendar(strm,1,calendar)
    strm%calendar = trim(calendar)

    ! initialize flag that stream has been set
    call shr_stream_setInit(strm)

  end subroutine shr_stream_init_from_fortran

  !===============================================================================
  subroutine shr_stream_set( strm, yearFirst, yearLast, yearAlign, offset, taxMode,  &
       meshFileName, fldListFile, fldListModel, filenames, rc)

    !-------------------------------------------------------------------------------
    ! set or override stream settings
    !-------------------------------------------------------------------------------

    ! !input/output parameters:
    type(shr_stream_streamType) ,intent(inout) :: strm         ! data stream
    integer      ,optional      ,intent(in)    :: yearFirst    ! first year to use
    integer      ,optional      ,intent(in)    :: yearLast     ! last  year to use
    integer      ,optional      ,intent(in)    :: yearAlign    ! align yearFirst with this model year
    integer      ,optional      ,intent(in)    :: offset       ! offset in seconds of stream data
    character(*) ,optional      ,intent(in)    :: taxMode      ! time axis mode
    character(*) ,optional      ,intent(in)    :: meshFileName ! stream mesh file
    character(*) ,optional      ,intent(in)    :: fldListFile  ! file field names, colon delim list
    character(*) ,optional      ,intent(in)    :: fldListModel ! model field names, colon delim list
    character(*) ,optional      ,intent(in)    :: filenames(:) ! input filenames
    integer      ,optional      ,intent(out)   :: rc           ! return code

    ! local variables
    integer                   :: n
    character(CL)             :: calendar ! stream calendar
    type(shr_stream_fileType) :: tempFile ! File being constructed.
!    type(fileVector)          :: fileVec  ! Vector used to construct file array.
    character(*),parameter    :: subName = '(shr_stream_set) '
    character(*),parameter    :: F00   = "('(shr_stream_set) ',8a)"
    character(*),parameter    :: F01   = "('(shr_stream_set) ',1a,i6)"
    !-------------------------------------------------------------------------------

    call shr_stream_default(strm)

    if ( present(rc) ) rc = 0

    if (present(yearFirst    )) strm%yearFirst    = yearFirst
    if (present(yearLast     )) strm%yearLast     = yearLast
    if (present(yearAlign    )) strm%yearAlign    = yearAlign
    if (present(offset       )) strm%offset       = offset
    if (present(taxMode      )) strm%taxMode      = trim(taxMode)
    if (present(fldListFile  )) strm%fldListFile  = trim(fldListFile)
    if (present(fldListModel )) strm%fldListModel = trim(fldListModel)
    if (present(meshFileName )) strm%meshFileName = trim(meshFileName)

    if (present(filenames)) then
       do n = 1,size(filenames)
          ! Ignore null file names.
          if (trim(filenames(n)) /= trim(shr_stream_file_null)) then
             tempFile%name = trim(filenames(n))
!             call fileVec%push_back(tempFile)
          endif
       enddo
       ! True size after throwing out null names.
!       strm%nFiles = fileVec%vsize()
!       call fileVec%move_out(strm%file)
    endif

    ! get initial calendar value
    call shr_stream_getCalendar(strm,1,calendar)
    strm%calendar = trim(calendar)

    call shr_stream_setInit(strm)

  end subroutine shr_stream_set

  !===============================================================================
  subroutine shr_stream_default(strm)

    !-----------------------------------------------------------------------------
    ! set default values for everything in stream
    !-----------------------------------------------------------------------------

    ! input/output variables:
    type(shr_stream_streamType) ,intent(inout) :: strm      ! data stream
    !-------------------------------------------------------------------------------

    call shr_stream_clearInit(strm)

    if (allocated(strm%file)) deallocate(strm%file)

    strm%meshFileName = ' '
    strm%nFiles       = 0
    strm%yearFirst    = 0
    strm%yearLast     = 0
    strm%yearAlign    = 0
    strm%offset       = 0
    strm%taxMode      = trim(shr_stream_taxis_cycle)
    strm%k_lvd        = -1
    strm%n_lvd        = -1
    strm%found_lvd    = .false.
    strm%k_gvd        = -1
    strm%n_gvd        = -1
    strm%found_gvd    = .false.
    strm%fileopen     = .false.
    strm%currfile     = ''
    strm%fldListFile  = ' '
    strm%fldListModel = ' '
    strm%calendar     = shr_cal_noleap

  end subroutine shr_stream_default
  !===============================================================================

  subroutine shr_stream_readUpToTag(nUnit,tag,optionalTag,rc)

    !----- input/output -----
    integer,intent(in ) :: nUnit       ! i/o unit to read from
    character(*)        ,intent(in ) :: tag         ! string to search for
    logical, optional   ,intent(in ) :: optionalTag ! this is an optional tag
    integer,intent(out) :: rc          ! return code

    !----- local -----
    character(CL)           :: str   ! temp char string
    logical                          :: localOptionalTag ! local version of optionalTag

    !----- formats -----
    character(*),parameter :: subName = '(shr_stream_readUpToTag) '
    character(*),parameter :: F00   = "('(shr_stream_readUpToTag) ',8a)"

    !-------------------------------------------------------------------------------
    ! Note: does not rewind to start of file
    !-------------------------------------------------------------------------------

    rc = 1
    localOptionalTag = .false.
    if (present(optionalTag)) localOptionalTag = optionalTag
    do while (.true.)
       read(nUnit,'(a)',END=999) str
       str = adjustL(str)
       if (str(1:len_trim(adjustL(tag))) == trim(adjustL(tag))) then
          rc = 0
          exit
       end if
    end do
999 continue

    if (rc /= 0 .and. .not. localOptionalTag ) then
       write(s_logunit,F00) "ERROR: tag not found: ",trim(tag)
       call shr_sys_abort(subName//"ERROR: tag not found")
    end if

  end subroutine shr_stream_readUpToTag

  !===============================================================================
  subroutine shr_stream_parseInput(str,fileName,yearAlign,yearFirst,yearLast,rc)

    !-------------------------------------------------------------------------------
    ! shr_stream_parseInput -- extract fileName,yearAlign, etc. from a string
    !-------------------------------------------------------------------------------

    ! input/output parameters:
    character(*)       ,intent(in)  :: str       ! string to parse
    character(*)       ,intent(out) :: fileName  ! file name
    integer            ,intent(out) :: yearFirst ! first year to use
    integer            ,intent(out) :: yearLast  ! last  year to use
    integer            ,intent(out) :: yearAlign ! align yearFirst with this model year
    integer  ,optional ,intent(out) :: rc        ! return code

    ! local variables
    integer                :: n        ! generic index
    character(CL)          :: str2     ! temp work string
    character(*),parameter :: F00   = "('(shr_stream_parseInput) ',8a)"
    character(*),parameter :: F01   = "('(shr_stream_parseInput) ',a,3i10)"
    !-------------------------------------------------------------------------------

    if (debug>1) write(s_logunit,F00) "str       = ",trim(str)

    str2 = adjustL(str)
    n    = index(str2," ")
    fileName = str2(:n)
    read(str2(n:),*) yearAlign,yearFirst,yearLast

    if (debug>1) then
       write(s_logunit,F00) "fileName  = ",trim(fileName)
       write(s_logunit,F01) "yearAlign = ",yearAlign
       write(s_logunit,F01) "yearFirst = ",yearFirst
       write(s_logunit,F01) "yearLast  = ",yearLast
    end if

    if (present(rc)) rc = 0

  end subroutine shr_stream_parseInput

  !===============================================================================
  subroutine shr_stream_findBounds(strm,mDateIn, secIn, &
       mDateLB, dDateLB, secLB, n_lb, fileLB,  mDateUB, dDateUB, secUB, n_ub, fileUB)

    !    Given a stream and a model date, find time coordinates of the upper and
    !    lower time bounds surrounding the models date.  Returns the model date,
    !    data date, elasped seconds, time index, and file names associated with
    !    these upper and lower time bounds.

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

    !----- local -----
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

    if (debug>0) write(s_logunit,F02) "DEBUG: ---------- enter ------------------"

    if ( .not. shr_stream_isInit(strm)) then
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
       write(s_logunit,*) trim(subName),' ERROR: illegal taxMode = ',trim(strm%taxMode)
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
       write(s_logunit,*) trim(subName),' ERROR: dyear lt zero = ',dYear
       call shr_sys_abort(trim(subName)//' ERROR: dyear lt one')
    endif

    dDateIn = dYear*10000 + modulo(mDateIn,10000) ! mDateIn mapped to range of data years
    rDateIn = dDateIn + secIn/spd                 ! dDateIn + fraction of a day
    if(debug>0) then
       write(s_logunit,*) 'tcx fbd1 ',mYear,dYear,dDateIn,rDateIn
       write(s_logunit,*) 'tcx fbd2 ',yrFirst,yrLast,yrAlign,nYears
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
          write(s_logunit,F00)  "ERROR: LVD not found, all data is before yearFirst"
          call shr_sys_abort(trim(subName)//" ERROR: LVD not found, all data is before yearFirst")
       else
          !--- LVD is in or beyond yearFirst, verify it is not beyond yearLast ---
          if ( dDateL <= strm%file(strm%k_lvd)%date(strm%n_lvd) ) then
             write(s_logunit,F00)  "ERROR: LVD not found, all data is after yearLast"
             call shr_sys_abort(trim(subName)//" ERROR: LVD not found, all data is after yearLast")
          end if
       end if
       if (debug>1 ) then
          if (strm%found_lvd) write(s_logunit,F01) "DEBUG: found LVD = ",strm%file(k)%date(n)
       end if
    end if

    if (strm%found_lvd) then
       k = strm%k_lvd
       n = strm%n_lvd
       rDatelvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! LVD date + frac day
    else
       write(s_logunit,F00)  "ERROR: LVD not found yet"
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
       write(s_logunit,*) 'tcx fbd3 ',rDateIn,rDatelvd,rDategvd
    endif

    !-----------------------------------------------------------
    ! dateIn < rDatelvd
    !   limit -> abort
    !   extend -> use lvd value, set LB to 00000101
    !   cycle -> lvd is UB, gvd is LB, shift mDateLB by -nYears
    !-----------------------------------------------------------

    if (rDateIn < rDatelvd) then
       if (limit) then
          write(s_logunit,*)  trim(subName)," ERROR: limit on and rDateIn lt rDatelvd",rDateIn,rDatelvd
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
          !         write(s_logunit,*)'tcx fb1 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
          !         call shr_sys_flush(s_logunit)
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
                      if (debug>1 ) write(s_logunit,F01) "DEBUG: found GVD ",strm%file(k)%date(n)
                      exit B
                   end if
                end do
             end do B
          end if

          if (.not. strm%found_gvd) then
             write(s_logunit,F00)  "ERROR: GVD not found1"
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
          !         write(s_logunit,*)'tcx fb2 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
          !         call shr_sys_flush(s_logunit)
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
          write(s_logunit,*) trim(subName)," ERROR: limit on and rDateIn gt rDategvd",rDateIn,rDategvd
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
          !         write(s_logunit,*)'tcx fb3 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
          !         call shr_sys_flush(s_logunit)
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
          !         write(s_logunit,*)'tcx fb4 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
          !         call shr_sys_flush(s_logunit)
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
                write(s_logunit,*) trim(subName)," ERROR: limit on and rDateIn gt rDategvd",rDateIn,rDategvd
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
                !               write(s_logunit,*)'tcx fb5 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
                !               call shr_sys_flush(s_logunit)
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
                !               write(s_logunit,*)'tcx fb6 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
                !               call shr_sys_flush(s_logunit)
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
                   !                  write(s_logunit,*)'tcx fb7 ',n_lb,mDateLB,secLB,n_ub,mDateUB,secUB
                   !                  call shr_sys_flush(s_logunit)
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
          write(s_logunit,F01) "Incompatable date and secs sizes",size(strm%file(k)%date),size(strm%file(k)%secs)
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
          ! write(s_logunit,*) 'debug ',n,strm%offset,din,sin,dout,sout
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
      if (debug>1 ) write(s_logunit,F01) "checking t-coordinate data   for file k =",k

      if ( .not. strm%file(k)%haveData) then
         rc = 1
         write(s_logunit,F01) "Don't have data for file ",k
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
                  if (debug>1 ) write(s_logunit,F01) "comparing with previous file for file k =",k
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
                  if (debug>1 ) write(s_logunit,F01) "comparing with next     file for file k =",k
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
               write(s_logunit,F01) "ERROR: calendar dates must be increasing"
               write(s_logunit,F02) "date(n), date(n+1) = ",date1,date2
               call shr_sys_abort(subName//"ERROR: calendar dates must be increasing")
            else if ( date1 == date2 ) then
               if ( secs1 >= secs2 ) then
                  rc = 1
                  write(s_logunit,F01) "ERROR: elapsed seconds on a date must be strickly increasing"
                  write(s_logunit,F02) "secs(n), secs(n+1) = ",secs1,secs2
                  call shr_sys_abort(subName//"ERROR: elapsed seconds must be increasing")
               end if
            end if
            if ( secs1 < 0 .or. spd < secs1 ) then
               rc = 1
               write(s_logunit,F01) "ERROR: elapsed seconds out of valid range [0,spd]"
               write(s_logunit,F02) "secs(n) = ",secs1
               call shr_sys_abort(subName//"ERROR: elapsed seconds out of range")
            end if
         end if
      end do

      if (debug>0) write(s_logunit,F01) "data is OK (non-decreasing)  for file k =",k
    end subroutine verifyTCoord

  end subroutine shr_stream_readTCoord

  !===============================================================================
  subroutine shr_stream_getMeshFileName(stream, filename)

    ! Return stream mesh filename

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    character(len=*)            ,intent(out) :: filename
    !-------------------------------------------------------------------------------

    filename = stream%meshFileName

  end subroutine shr_stream_getMeshFileName

  !===============================================================================
  subroutine shr_stream_getModelFieldList(stream, list)

    ! Get list of file fields

    !input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    character(*)                ,intent(out) :: list(:)    ! field list
    !-------------------------------------------------------------------------------
    integer :: i
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
  subroutine shr_stream_getFileFieldName(stream, k, name)

    ! Get name of k-th field in list

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: stream  ! stream in question
    integer                     ,intent(in)  :: k       ! index of field
    character(*)                ,intent(out) :: name    ! k-th name in list
    !-------------------------------------------------------------------------------

    call shr_string_listGetName(stream%fldListFile, k, name)

  end subroutine shr_stream_getFileFieldName

  !===============================================================================
  subroutine shr_stream_getCalendar(strm,k,calendar)

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

    n = len_trim(lcal)
    if (ichar(lcal(n:n)) == 0 ) lcal(n:n) = ' '
    call shr_string_leftalign_and_convert_tabs(lcal)
    calendar = trim(shr_cal_calendarName(trim(lcal)))
    rCode = nf90_close(fid)
    if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_close')

  end subroutine shr_stream_getCalendar

  !===============================================================================
  subroutine shr_stream_getCurrFile(strm,fileopen,currfile,currpioid)

    ! returns current file information

    ! input/output parameters:
    type(shr_stream_streamType),intent(in)  :: strm     ! data stream
    logical           ,optional,intent(out) :: fileopen ! file open flag
    character(*)      ,optional,intent(out) :: currfile ! current filename
    type(file_desc_t) ,optional,intent(out) :: currpioid ! current pioid
    !-------------------------------------------------------------------------------

    if (present(fileopen  )) fileopen = strm%fileopen
    if (present(currfile  )) currfile = strm%currfile
    if (present(currpioid )) currpioid = strm%currpioid

  end subroutine shr_stream_getCurrFile

  !===============================================================================
  subroutine shr_stream_setCurrFile(strm,fileopen,currfile,currpioid)

    ! set current file information

    ! input/output parameters:
    type(shr_stream_streamType),intent(inout) :: strm     ! data stream
    logical           ,optional,intent(in) :: fileopen ! file open flag
    character(*)      ,optional,intent(in) :: currfile ! current filename
    type(file_desc_t) ,optional,intent(in) :: currpioid ! current pioid
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
       write(s_logunit,F00) "ERROR: input file name is not in stream: ",trim(fn)
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
       write(s_logunit,F00) "ERROR: input file name is not in stream: ",trim(fn)
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
  subroutine shr_stream_restWrite(strm, fileName, caseName, caseDesc, nstrms, rc)

    ! Write stream data to a restart file.

    ! input/output parameters:
    type(shr_stream_streamType) ,intent(in)  :: strm(:)    ! data streams
    character(*)                ,intent(in)  :: fileName   ! name of restart file
    character(*)                ,intent(in)  :: caseName   ! case name
    character(*)                ,intent(in)  :: caseDesc   ! case description
    integer, optional           ,intent(in)  :: nstrms     ! number of streams
    integer,optional            ,intent(out) :: rc         ! return code

    ! local variables
    integer       :: rCode       ! return code
    integer       :: nStreams    ! number of streams
    integer       :: k,n         ! generic loop index
    character( 8) :: dStr        ! F90 wall clock date str yyyymmdd
    character(10) :: tStr        ! F90 wall clock time str hhmmss.sss
    character(CS) :: str         ! generic text string
    integer       :: nUnit       ! a file unit number
    integer       :: nt          ! number of time samples
    character(CS) :: tInterpAlgo ! for backwards compatability
    character(*),parameter :: subName = '(shr_stream_restWrite) '
    character(*),parameter :: F00   = "('(shr_stream_restWrite) ',16a) "
    character(*),parameter :: F01   = "('(shr_stream_restWrite) ',a,i5,a,5a) "
    character(*),parameter :: F02   = "('(shr_stream_restWrite) ',a,i5,a,5i8) "
    character(*),parameter :: F03   = "('(shr_stream_restWrite) ',a,i5,a,5l3) "
    !-------------------------------------------------------------------------------

    rCode = 0
    tInterpAlgo = 'unused'

    if (present(nstrms)) then
       if (size(strm) < nstrms) then
          write(s_logunit,F02) "ERROR: nstrms too large for strm",size(strm),nstrms
          call shr_sys_abort(subname//": ERROR: nstrms too large for strm")
       endif
       nStreams = nstrms
    else
       nStreams = size(strm)
    endif
    call date_and_time(dStr,tStr)

    ! log info to stdout
    if (debug > 0) then
       write(s_logunit,F00) "case name        : ",trim(caseName)
       write(s_logunit,F00) "case description : ",trim(caseDesc)
       write(s_logunit,F00) "File created     : ",&
            dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '//tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
       write(s_logunit,F01) "Number of streams ",nStreams
    end if

    !----------------------------------------------------------------------------
    ! write the  data
    !----------------------------------------------------------------------------

    open(newunit=nUnit,file=trim(fileName),form="unformatted",action="write")

    read(nUnit) nStreams
    if (present(nstrms)) then
       if (nstrms /= nStreams) then
          write(s_logunit,F02) "ERROR: nstrms ne nStreams on restart",nstrms,' ',nStreams
          call shr_sys_abort(subname//": ERROR: nstrms ne nStreams on restart")
       endif
       nStreams = nstrms
    endif

    str = "case name        : "//caseName
    write(nUnit) str
    str = "case description : "//caseDesc
    write(nUnit) str
    str = 'File created     : '&
         //dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '//tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
    write(nUnit) str

    write(nUnit) nStreams
    do k = 1,nStreams
       if (.not. shr_stream_isInit(strm(k))) then  ! has stream been initialized?
          rCode = 1
          write(s_logunit,F01) "ERROR: can't write uninitialized stream to a restart file, k = ",k
          call shr_sys_abort(subName//": ERROR: given uninitialized stream")
       end if

       write(nUnit) strm(k)%init         ! has stream been initialized?
       write(nUnit) strm(k)%nFiles       ! number of data files

       if (debug > 0) then
          write(s_logunit,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
          write(s_logunit,F03) "* stream ",k," first have data = ",strm(k)%file(1)%haveData
          write(s_logunit,F02) "* stream ",k," first nt        = ",strm(k)%file(1)%nt
       end if

       nt = strm(k)%file(1)%nt
       if (strm(k)%file(1)%haveData .and. debug > 0) then
          write(s_logunit,F02) "* stream ",k," first date secs = ", strm(k)%file(1)%date(1), strm(k)%file(1)%secs(1)
          write(s_logunit,F02) "* stream ",k," last  date secs = ", strm(k)%file(1)%date(nt), strm(k)%file(1)%secs(nt)
       endif
       do n=1,strm(k)%nFiles                      ! data specific to each file...
          write(nUnit) strm(k)%file(n)%name       ! the file name
          write(nUnit) strm(k)%file(n)%haveData   ! has t-coord data been read in?
          write(nUnit) strm(k)%file(n)%nt         ! size of time dimension
          if (strm(k)%file(n)%haveData) then      ! ie. if arrays have been allocated
             write(nUnit) strm(k)%file(n)%date(:) ! t-coord date: yyyymmdd
             write(nUnit) strm(k)%file(n)%secs(:) ! t-coord secs: elapsed on date
          end if
       end do

       write(nUnit) strm(k)%yearFirst      ! first year to use in t-axis (yyyymmdd)
       write(nUnit) strm(k)%yearLast       ! last  year to use in t-axis (yyyymmdd)
       write(nUnit) strm(k)%yearAlign      ! align yearFirst with this model year
       write(nUnit) strm(k)%offset         ! time axis offset
      !write(nUnit) strm(k)%taxMode        ! time axis cycling mode

       write(nUnit) strm(k)%k_lvd          ! file of least valid date
       write(nUnit) strm(k)%n_lvd          ! sample of least valid date
       write(nUnit) strm(k)%found_lvd      ! T <=> k_lvd,n_lvd have been set
       write(nUnit) strm(k)%k_gvd          ! file of greatest valid date
       write(nUnit) strm(k)%n_gvd          ! sample of greatest valid date
       write(nUnit) strm(k)%found_gvd      ! T <=> k_gvd,n_gvd have been set

       write(nUnit) strm(k)%fldListFile    ! field list: file's  field names
       write(nUnit) strm(k)%fldListModel   ! field list: model's field names
       write(nUnit) tInterpAlgo            ! unused
       write(nUnit) strm(k)%meshFileName   ! mesh filename

    end do

    close(nUnit)
    if ( present(rc) ) rc = rCode

  end subroutine  shr_stream_restWrite

  !===============================================================================
  subroutine shr_stream_restRead(strm,fileName,nstrms,rc)

    ! read stream data from a restart file
    ! Either shr_stream_init xor shr_stream_restRead must be called
    ! Do not call both routines.

    ! input/output parameters:
    type(shr_stream_streamType)  ,intent(inout) :: strm(:)  ! vector of data streams
    character(*)                 ,intent(in)    :: fileName ! name of restart file
    integer,optional,intent(in)    :: nstrms   ! number of streams in strm
    integer,optional,intent(out)   :: rc       ! return code

    !--- local ---
    integer         :: rCode                          ! return code
    integer         :: nStreams                       ! number of streams
    integer         :: k,n                            ! generic loop index
    character(CS)   :: str                            ! generic text string
    integer         :: nUnit                          ! a file unit number
    integer         :: inpi                           ! input integer
    character(CXX)  :: inpcx                          ! input char
    character(CL)   :: inpcl                          ! input char
    character(CS)   :: inpcs                          ! input char
    integer         :: nt                             ! size of time dimension
    character(CS)   :: tInterpAlgo                    ! for backwards compatability
    character(CL)   :: name                           ! local variables
    integer         :: nFiles                         ! local variables
    integer         :: k_lvd, n_lvd, k_gvd, n_gvd     ! local variables
    logical         :: found_lvd, found_gvd, haveData ! local variables
    integer,pointer :: date(:),secs(:)                ! local variables
    logical         :: abort                          ! abort the restart read
    logical         :: readok                         ! read of restarts ok

    !--- formats ---
    character(*),parameter :: subName = '(shr_stream_restRead) '
    character(*),parameter :: F00   = "('(shr_stream_restRead) ',16a) "
    character(*),parameter :: F01   = "('(shr_stream_restRead) ',a,i5,a,5a) "
    character(*),parameter :: F02   = "('(shr_stream_restRead) ',a,i5,a,5i8) "
    character(*),parameter :: F03   = "('(shr_stream_restRead) ',a,i5,a,5l3) "
    character(*),parameter :: F04   = "('(shr_stream_restRead) ',a,4i8) "
    character(*),parameter :: F05   = "('(shr_stream_restRead) ',a,2i8,6a) "
    !-------------------------------------------------------------------------------

    rCode = 0
    tInterpAlgo = 'unused'
    abort = .false.
    inpcl = ' '

    !----------------------------------------------------------------------------
    ! read the  data
    !----------------------------------------------------------------------------

    open(newunit=nUnit,file=trim(fileName),form="unformatted",status="old",action="read", iostat=rCode)
    if ( rCode /= 0 )then
       call shr_sys_abort(subName//": ERROR: error opening file: "//trim(fileName) )
    end if

    read(nUnit) str         ! case name
    if (debug > 0) write(s_logunit,F00) trim(str)
    read(nUnit) str         ! case description
    if (debug > 0) write(s_logunit,F00) trim(str)
    read(nUnit) str         ! file creation date
    if (debug > 0) write(s_logunit,F00) trim(str)

    read(nUnit) nStreams
    if (present(nstrms)) then
       if (nstrms /= nStreams) then
          write(s_logunit,F02) "ERROR: nstrms ne nStreams on restart",nstrms,' ',nStreams
          call shr_sys_abort(subname//": ERROR: nstrms ne nStreams on restart")
       endif
       nStreams = nstrms
    endif
    if (debug > 0) write(s_logunit,F01) "Number of streams ",nStreams

    do k = 1,nStreams
       read(nUnit) strm(k)%init         ! has stream been initialized?
       if (.not. strm(k)%init) then
          rCode = 1
          write(s_logunit,F01) "ERROR: uninitialized stream in restart file, k = ",k
          call shr_sys_abort(subName//": ERROR: reading uninitialized stream")
       end if
       call shr_stream_setInit(strm(k))

       readok = .true.

       ! tcraig, don't overwrite these from input
       read(nUnit) nFiles       ! number of data files

       do n=1,nFiles                     ! data specific to each file...
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
                   write(s_logunit,F05) "reading time axis for stream restart filename ",k,n, &
                        ' ',trim(name),' ',trim(strm(k)%file(n)%name)
                   strm(k)%file(n)%nt = nt
                   strm(k)%file(n)%haveData = haveData
                   allocate(strm(k)%file(n)%date(nt))
                   allocate(strm(k)%file(n)%secs(nt))
                   strm(k)%file(n)%date(1:nt) = date(1:nt)
                   strm(k)%file(n)%secs(1:nt) = secs(1:nt)
                else
                   write(s_logunit,F05) "WARNING, skip time axis for stream restart filename ",k,n,&
                        ' ',trim(name),' ',trim(strm(k)%file(n)%name)
                   readok = .false.
                endif  ! filenames consistent
             endif  ! strm nfiles
             deallocate(date)
             deallocate(secs)
          end if
       end do

       if (debug > 0) then
          write(s_logunit,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
          write(s_logunit,F03) "* stream ",k," first have data = ",strm(k)%file(1)%haveData
          write(s_logunit,F02) "* stream ",k," first nt        = ",strm(k)%file(1)%nt
       end if

       if (strm(k)%file(1)%haveData) then
          nt = strm(k)%file(1)%nt
          if (debug > 0) then
             write(s_logunit,F02) "* stream ",k," first date secs = ", strm(k)%file(1)%date(1),strm(k)%file(1)%secs(1)
             write(s_logunit,F02) "* stream ",k," last  date secs = ", strm(k)%file(1)%date(nt),strm(k)%file(1)%secs(nt)
          end if
       endif

       ! offset is the only field that should not change here for time axis
       read(nUnit) inpi   ! first year to use in t-axis (yyyymmdd)
       read(nUnit) inpi   ! last year to use in t-axis (yyyymmdd)
       read(nUnit) inpi   ! align year to use in t-axis (yyyymmdd)
       read(nUnit) inpi   ! time axis offset
       if (inpi /= strm(k)%offset) then
          write(s_logunit,F04) " ERROR: offset disagrees ",k,strm(k)%offset,inpi
          write(s_logunit,F00) "ERRORS Detected ABORTING NOW"
          call shr_sys_abort(subName//": ERRORS Detected ABORTING NOW")
       endif

       ! read(nUnit) strm(k)%taxMode      ! time axis cycling mode
       read(nUnit) k_lvd        ! file of least valid date
       read(nUnit) n_lvd        ! sample of least valid date
       read(nUnit) found_lvd    ! T <=> k_lvd,n_lvd have been set
       read(nUnit) k_gvd        ! file of greatest valid date
       read(nUnit) n_gvd        ! sample of greatest valid date
       read(nUnit) found_gvd    ! T <=> k_gvd,n_gvd have been set

       ! only overwrite if restart read is ok
       if (readok) then
          write(s_logunit,F05) "setting k n and found lvd gvd on restart ",k,n,' ',trim(name)
          strm(k)%k_lvd     = k_lvd
          strm(k)%n_lvd     = n_lvd
          strm(k)%found_lvd = found_lvd
          strm(k)%k_gvd     = k_gvd
          strm(k)%n_gvd     = n_gvd
          strm(k)%found_gvd = found_gvd
       endif

       ! don't overwrite these from input
       read(nUnit) inpcx !  fldListFile   - field list: file's  field names
       read(nUnit) inpcx !  fldListModel  - field list: model's field names
       read(nUnit) inpcs !  tInterpAlgo   -  unused
       read(nUnit) inpcl !  meshFileName  -  mesh filename

    end do

    close(nUnit)
    if ( present(rc) ) rc = rCode

  end subroutine shr_stream_restRead

  !===============================================================================
  subroutine shr_stream_dataDump(strm)

    ! Dump all data to stdout for debugging

    ! input/output parameters:
    type(shr_stream_streamType),intent(in) :: strm      ! data stream

    !----- local -----
    integer   :: k ! generic loop index
    character(*),parameter :: F00   = "('(shr_stream_dataDump) ',8a)"
    character(*),parameter :: F01   = "('(shr_stream_dataDump) ',a,3i5)"
    character(*),parameter :: F02   = "('(shr_stream_dataDump) ',a,365i9.8)"
    character(*),parameter :: F03   = "('(shr_stream_dataDump) ',a,365i6)"
    !-------------------------------------------------------------------------------

    if (debug > 0) then
       write(s_logunit,F00) "dump internal data for debugging..."
       write(s_logunit,F01) "nFiles        = ", strm%nFiles
       do k=1,strm%nFiles
          write(s_logunit,F01) "data for file k = ",k
          write(s_logunit,F00)    "* file(k)%name    = ", trim(strm%file(k)%name)
          if ( strm%file(k)%haveData ) then
             write(s_logunit,F01) "* file(k)%nt      = ", strm%file(k)%nt
             write(s_logunit,F02) "* file(k)%date(:) = ", strm%file(k)%date(:)
             write(s_logunit,F03) "* file(k)%Secs(:) = ", strm%file(k)%secs(:)
          else
             write(s_logunit,F00) "* time coord data not read in yet for this file"
          end if
       end do
       write(s_logunit,F01) "yearF/L/A    = ", strm%yearFirst,strm%yearLast,strm%yearAlign
       write(s_logunit,F01) "offset       = ", strm%offset
       write(s_logunit,F00) "taxMode      = ", trim(strm%taxMode)

       write(s_logunit,F00) "fldListFile  = ", trim(strm%fldListFile)
       write(s_logunit,F00) "fldListModel = ", trim(strm%fldListModel)
       write(s_logunit,F00) "meshFileName = ", trim(strm%meshFileName)
    end if

  end subroutine shr_stream_dataDump

  !===============================================================================
  logical function shr_stream_isInit(strm,rc)

    ! Checks if stream is initialized

    ! input/output parameters:
    type(shr_stream_streamType),  intent(in)    :: strm
    integer,optional,intent(out)   :: rc
    !-------------------------------------------------------------------------------

    shr_stream_isInit = .false.
    if (size(strm%initarr) == initarr_size) shr_stream_isInit = .true.
    if (present(rc)) rc = 0

  end function shr_stream_isInit

  !===============================================================================
  subroutine shr_stream_setInit(strm,rc)

    ! Sets stream init flag to TRUE

    ! input/output parameters:
    type(shr_stream_streamType),  intent(inout) :: strm
    integer,optional,intent(out)   :: rc

    !--- local ---
    integer :: ier
    !-------------------------------------------------------------------------------

    strm%init = .true.
    deallocate(strm%initarr,stat=ier)
    allocate(strm%initarr(initarr_size))
    if (present(rc)) rc = 0

  end subroutine shr_stream_setInit

  !===============================================================================
  subroutine shr_stream_clearInit(strm,rc)

    !  Checks if stream is initialized

    ! !input/output parameters:
    type(shr_stream_streamType),  intent(inout) :: strm
    integer,optional,intent(out)   :: rc

    !--- local ---
    integer :: ier
    !-------------------------------------------------------------------------------

    strm%init = .true.
    deallocate(strm%initarr,stat=ier)
    allocate(strm%initarr(initarr_size + 5))
    if (present(rc)) rc = 0

  end subroutine shr_stream_clearInit

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
