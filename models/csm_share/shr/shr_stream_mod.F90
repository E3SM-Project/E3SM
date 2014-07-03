!===============================================================================
! SVN $Id: shr_stream_mod.F90 43642 2013-01-31 19:22:00Z sacks $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_130528/shr/shr_stream_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_stream_mod -- Data type and methods to manage input data streams.
!
! !DESCRIPTION:
!    A "data stream" is a sequence of input files where each file contains the 
!    same set of data fields and all the data fields are on the same grid.
!    The sequence of input data files provides an uninterupted time series of 
!    data.
!
!    A stream data type stores information about one data stream, including the
!    range of data date years to use and how data dates align with model dates.  
!
!    Given a model date, this module can return data dates that are upper and 
!    lower time bounds around the given model date and the names of the files
!    containing those dates.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Apr-13 - B. Kauffman - moved code from dshr to shr 
!     2005-Apr-01 - B. Kauffman - first functional version of findBounds
!     2004-Dec-xx - B. Kauffman - initial module
!
! !INTERFACE: ------------------------------------------------------------------

module shr_stream_mod

   use shr_sys_mod    ! shared system calls
   use shr_kind_mod   ! kinds for strong typing
   use shr_const_mod  ! shared constants (including seconds per day)
   use shr_string_mod ! string & list methods
   use shr_mpi_mod    ! shared mpi
   use shr_file_mod   ! file methods
   use shr_cal_mod    ! calendar methods

   use shr_log_mod, only : s_loglev   => shr_log_Level
   use shr_log_mod, only : s_logunit  => shr_log_Unit
   use perf_mod

   implicit none

   private ! default private

! !PUBLIC TYPES:

   public :: shr_stream_streamType        ! stream data type with private components
   public :: shr_stream_fileType 

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_stream_init              ! initialize a stream
   public :: shr_stream_set               ! set stream values
   public :: shr_stream_default           ! set default values
   public :: shr_stream_parseInput        ! extract fileName,yearAlign, etc. from a string
   public :: shr_stream_findBounds        ! return lower/upper bounding date info
   public :: shr_stream_getFileFieldList  ! return input-file field name list
   public :: shr_stream_getModelFieldList ! return model      field name list
   public :: shr_stream_getFileFieldName  ! return k-th input-file field name
   public :: shr_stream_getModelFieldName ! return k-th model      field name list
   public :: shr_stream_getFirstFileName  ! return the 1st file name in stream
   public :: shr_stream_getNextFileName   ! return next file in sequence
   public :: shr_stream_getPrevFileName   ! return previous file in sequence
   public :: shr_stream_getFilePath       ! return file path
   public :: shr_stream_getDataSource     ! return the stream's meta data
   public :: shr_stream_getDomainInfo     ! return the stream's domain info data
   public :: shr_stream_getFile           ! acquire file, return name of file to open
   public :: shr_stream_getNFiles         ! get the number of files in a stream
   public :: shr_stream_getCalendar       ! get the stream calendar
   public :: shr_stream_dataDump          ! internal stream data for debugging
   public :: shr_stream_restWrite         ! write a streams restart file
   public :: shr_stream_restRead          ! read  a streams restart file
   public :: shr_stream_setDebug          ! set internal shr_stream debug level
   public :: shr_stream_setAbort          ! set internal shr_stream abort flag
   public :: shr_stream_getDebug          ! get internal shr_stream debug level
   public :: shr_stream_isInit            ! check if stream is initialized
!   public :: shr_stream_bcast             ! broadcast a stream (untested)

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   character(SHR_KIND_CS),parameter,public :: shr_stream_taxis_cycle  = 'cycle'
   character(SHR_KIND_CS),parameter,public :: shr_stream_taxis_extend = 'extend'
   character(SHR_KIND_CS),parameter,public :: shr_stream_taxis_limit  = 'limit'
   character(SHR_KIND_CS),parameter,public :: shr_stream_file_null  = 'not_set'

   !--- a useful derived type to use inside shr_stream_streamType ---
   type shr_stream_fileType
      character(SHR_KIND_CL)         :: name        ! the file name
      logical                        :: haveData    ! has t-coord data been read in?
      integer  (SHR_KIND_IN)         :: nt          ! size of time dimension
      integer  (SHR_KIND_IN),pointer :: date(:) => null()     ! t-coord date: yyyymmdd
      integer  (SHR_KIND_IN),pointer :: secs(:) => null()     ! t-coord secs: elapsed on date
   end type shr_stream_fileType

   !--- hard-coded array dims ~ could allocate these at run time ---
   integer(SHR_KIND_IN),parameter :: nFileMax = 1000  ! max number of files

   type shr_stream_streamType
      !private                                    ! no public access to internal components
      !--- input data file names and data ---
      logical                   :: init           ! has stream been initialized?
      integer  (SHR_KIND_IN),pointer :: initarr(:) => null()! surrogate for init flag
      integer  (SHR_KIND_IN)    :: nFiles         ! number of data files
      character(SHR_KIND_CS)    :: dataSource     ! meta data identifying data source
      character(SHR_KIND_CL)    :: filePath       ! remote location of data files
      type(shr_stream_fileType) :: file(nFileMax) ! data specific to each file

      !--- specifies how model dates align with data dates ---
      integer(SHR_KIND_IN)      :: yearFirst      ! first year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearLast       ! last  year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearAlign      ! align yearFirst with this model year
      integer(SHR_KIND_IN)      :: offset         ! offset in seconds of stream data
      character(SHR_KIND_CS)    :: taxMode        ! cycling option for time axis

      !--- useful for quicker searching ---
      integer(SHR_KIND_IN) :: k_lvd,n_lvd         ! file/sample of least valid date
      logical              :: found_lvd           ! T <=> k_lvd,n_lvd have been set
      integer(SHR_KIND_IN) :: k_gvd,n_gvd         ! file/sample of greatest valid date
      logical              :: found_gvd           ! T <=> k_gvd,n_gvd have been set

      !--- stream data not used by stream module itself ---
      character(SHR_KIND_CXX):: fldListFile       ! field list: file's  field names
      character(SHR_KIND_CXX):: fldListModel      ! field list: model's field names
      character(SHR_KIND_CL) :: domFilePath       ! domain file: file path of domain file
      character(SHR_KIND_CL) :: domFileName       ! domain file: name
      character(SHR_KIND_CS) :: domTvarName       ! domain file: time-dim var name
      character(SHR_KIND_CS) :: domXvarName       ! domain file: x-dim var name
      character(SHR_KIND_CS) :: domYvarName       ! domain file: y-dim var ame
      character(SHR_KIND_CS) :: domAreaName       ! domain file: area  var name
      character(SHR_KIND_CS) :: domMaskName       ! domain file: mask  var name

      character(SHR_KIND_CS) :: tInterpAlgo       ! Algorithm to use for time interpolation
      character(SHR_KIND_CL) :: calendar          ! stream calendar
   end type shr_stream_streamType

   !----- parameters -----
   real(SHR_KIND_R8)   ,parameter :: spd = SHR_CONST_CDAY ! seconds per day
   integer(SHR_KIND_IN),parameter :: initarr_size = 3     ! size of initarr
   integer(SHR_KIND_IN),save :: debug = 0        ! edit/turn-on for debug write statements
   logical             ,save :: doabort = .true. ! flag if abort on error

!===============================================================================
contains
!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_stream_init -- initialize stream datatype, read description text file
!
! !DESCRIPTION:
!
! !REMARKS:
!    should input be via standard Fortran namelist?
!
! !REVISION HISTORY:
!     2007-Sep-17 - B. Kauffman - reworked wrt new streams.txt format
!     2005-Feb-03 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_init(strm,infoFile,yearFirst,yearLast,yearAlign,taxMode,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)    ,intent(out) :: strm      ! data stream
   character(*)                   ,intent(in)  :: infoFile  ! file with stream info, must read
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearFirst ! first year to use 
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearLast  ! last  year to use 
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearAlign ! align yearFirst with this model year
   character(*)          ,optional,intent(in)  :: taxMode   ! time axis cycling option
   integer  (SHR_KIND_IN),optional,intent(out) :: rc        ! return code
 
!EOP

   !----- local -----
   integer  (SHR_KIND_IN) :: n             ! generic index
   character(SHR_KIND_CL) :: str           ! string to parse from input data file
   integer  (SHR_KIND_IN) :: int           ! integer to parse from input data file
   character(SHR_KIND_CL) :: subStr        ! sub-string of interest
   integer  (SHR_KIND_IN) :: nUnit         ! file i/o unit number
   character(SHR_KIND_CS) :: startTag      ! input file start tag
   character(SHR_KIND_CS) ::   endTag      ! input file   end tag
   character(SHR_KIND_CS) :: fldNameFile   ! field name in data file field list
   character(SHR_KIND_CS) :: fldNameModel  ! field name in model     field list
   character(SHR_KIND_CXX):: fldListFile   ! list of data file fields, colon delim list
   character(SHR_KIND_CXX):: fldListModel  ! list of model     fields, colon delim list
   character(SHR_KIND_CL) :: calendar      ! stream calendar
   integer  (SHR_KIND_IN) :: rCode, rCode2 ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_init) ' 
   character(*),parameter :: F00   = "('(shr_stream_init) ',8a)" 

!-------------------------------------------------------------------------------
! notes: 
! * should this use standard namelist input?
! * needs more robust error checking
! o yearFirst,yearLast,yearAlign are provided by calling routine
! o parse infoFile for remaining, except for...
! o fileNT,fileDates, & fileSecs, which are initially set to -1, but but are replaced with
!   valid values as each file is opened for the first time
!-------------------------------------------------------------------------------

   rCode = 0
   write(s_logunit,F00) 'Reading file ',trim(infoFile)

   call shr_stream_default(strm)

   strm%yearFirst        = yearFirst
   strm%yearLast         = yearLast
   strm%yearAlign        = yearAlign
   if (present(taxMode)) then
      strm%taxMode       = trim(taxMode)
   endif

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading data source'
   !-----------------------------------------------------------------------------  

   nUnit = shr_file_getUnit() ! get unused unit number

   !--- find start tag ---
   startTag =  "<dataSource>"
     endTag = "</dataSource>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ',iostat=rCode)
   if (rCode /= 0) goto 999
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode)
   if (rCode /= 0) goto 999

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   strm%dataSource = str
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * format = ', trim(strm%dataSource)

   close(nUnit)
   call shr_file_freeUnit(nUnit)

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading field data variable names'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<fieldInfo>"
     endTag = "</fieldInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode)
   if (rCode /= 0) goto 999
   startTag =  "<variableNames>"
     endTag = "</variableNames>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode)
   if (rCode /= 0) goto 999

   !--- read data ---
   n=0
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      call shr_string_leftAlign(str)
      n=n+1
      if (str(1:len_trim(endTag)) == trim(endTag)) exit
      fldNameFile  = ""
      fldNameModel = ""
      read(str,*,iostat=rCode) fldNameFile,fldNameModel
      if (len_trim(fldNameFile)==0 .or. len_trim(fldNameModel)==0 ) then
         rCode = 1
         write(s_logunit,F00) "ERROR: reading field names"
         write(s_logunit,F00) '* fldNameFile  = ',trim(fldNameFile)
         write(s_logunit,F00) '* fldNameModel = ',trim(fldNameModel)
         call shr_stream_abort(subName//"ERROR: reading field names")
      end if
      if (n==1) then
         strm%fldListFile  = trim(fldNameFile )
         strm%fldListModel = trim(fldNameModel)
      else  
         strm%fldListFile  = trim(strm%fldListFile ) // ":" // trim(fldNameFile )
         strm%fldListModel = trim(strm%fldListModel) // ":" // trim(fldNameModel)
      end if
   end do
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * file  field list = ',trim(strm%fldListFile )
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * model field list = ',trim(strm%fldListModel)
   if (n==0) then
      rCode = 1
      write(s_logunit,F00) "ERROR: no input field names"
      call shr_stream_abort(subName//"ERROR: no input field names")
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading time-interpolation alogrithm '
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<fieldInfo>"
     endTag = "</fieldInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading offset'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<fieldInfo>"
     endTag = "</fieldInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<offset>"
     endTag = "</offset>"
   call shr_stream_readUpToTag(nUnit,startTag,optionalTag=.true.,rc=rCode2)
   if (rCode2 == 0) then
      !--- read data ---
      read(nUnit,*,END=999) int
      strm%offset = int
   else
      strm%offset = 0
   end if
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * offset ',strm%offset

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading data file path'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<fieldInfo>"
     endTag = "</fieldInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<filePath>"
     endTag = "</filePath>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   n = len_trim(str)
   if (n>0 .and. str(n:n) /= '/') str(n+1:n+2) = "/ " ! must have trailing slash
   if (n==0) str = "./ "                              ! null path => ./
   strm%FilePath = str
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * data file path = ', trim(strm%FilePath)

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading field data file names'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<fieldInfo>"
     endTag = "</fieldInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<fileNames>"
     endTag = "</fileNames>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   !--- read data ---
   n=0
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      call shr_string_leftAlign(str)
      if (str(1:len_trim(endTag)) == trim(endTag)) exit
      n=n+1
      if (n > nFileMax) then
         rCode = 1
         write(s_logunit,F00) "ERROR: exceeded max number of files"
         call shr_stream_abort(subName//"ERROR: exceeded max number of files")
         if ( present(rc) ) rc = rCode
         close(nUnit) 
         call shr_file_freeUnit(nUnit)
         return
      end if
      strm%file(n)%name = str
      if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * ',trim(strm%file(n)%name)
   end do
   strm%nFiles = n
   if (n==0) then
      rCode = 1
      write(s_logunit,F00) "ERROR: no input file names"
      call shr_stream_abort(subName//"ERROR: no input file names")
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading domain data variable names'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<domainInfo>"
     endTag = "</domainInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<variableNames>"
     endTag = "</variableNames>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   !--- read data ---
   n=0
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      call shr_string_leftAlign(str)
      n=n+1
      if (str(1:len_trim(endTag)) == trim(endTag)) exit
      fldNameFile  = ""
      fldNameModel = ""
      read(str,*,iostat=rCode2) fldNameFile,fldNameModel
      if (len_trim(fldNameFile)==0 .or. len_trim(fldNameModel)==0 ) then
         rCode = 1
         write(s_logunit,F00) "ERROR: reading field names"
         write(s_logunit,F00) '* fldNameFile  = ',trim(fldNameFile)
         write(s_logunit,F00) '* fldNameModel = ',trim(fldNameModel)
         call shr_stream_abort(subName//"ERROR: reading field names")
      end if
      if (n==1) then
         fldListFile  = trim(fldNameFile )
         fldListModel = trim(fldNameModel)
      else  
         fldListFile  = trim(fldListFile ) // ":" // trim(fldNameFile )
         fldListModel = trim(fldListModel) // ":" // trim(fldNameModel)
      end if
   end do
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * file  field list = ',trim(fldListFile )
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * model field list = ',trim(fldListModel)

   if (n==0) then
      rCode = 1
      write(s_logunit,F00) "ERROR: no input field names"
      call shr_stream_abort(subName//"ERROR: no input field names")
   else
      !--- get time variable name ---
      n = shr_string_listGetIndexF(fldListModel,"time")
      if (n==0) then
         rCode = 1
         write(s_logunit,F00) "ERROR: no input field names"
         call shr_stream_abort(subName//"ERROR: no time variable name")
      else
         call shr_string_listGetName (fldListFile,n,substr,rc)
         strm%domTvarName = subStr
      endif
      !--- get longitude variable name ---
      n = shr_string_listGetIndexF(fldListModel,"lon")
      if (n==0) then
         rCode = 1
         write(s_logunit,F00) "ERROR: no input field names"
         call shr_stream_abort(subName//"ERROR: no lon variable name")
      else
         call shr_string_listGetName (fldListFile,n,substr,rc)
         strm%domXvarName = subStr
      endif
      !--- get latitude variable name ---
      n = shr_string_listGetIndexF(fldListModel,"lat")
      if (n==0) then
         rCode = 1
         write(s_logunit,F00) "ERROR: no input field names"
         call shr_stream_abort(subName//"ERROR: no lat variable name")
      else
         call shr_string_listGetName (fldListFile,n,substr,rc)
         strm%domYvarName = subStr
      endif
      !--- get area variable name ---
      n = shr_string_listGetIndexF(fldListModel,"area")
      if (n==0) then
!         rCode = 1
!         write(s_logunit,F00) "ERROR: no input field names"
!         call shr_stream_abort(subName//"ERROR: no area variable name")
         strm%domAreaName = 'unknownname'
      else
         call shr_string_listGetName (fldListFile,n,substr,rc)
         strm%domAreaName = subStr
      endif
      !--- get mask variable name ---
      n = shr_string_listGetIndexF(fldListModel,"mask")
      if (n==0) then
!         rCode = 1
!         write(s_logunit,F00) "ERROR: no input field names"
!         call shr_stream_abort(subName//"ERROR: no mask variable name")
         strm%domMaskName = 'unknownname'
      else
         call shr_string_listGetName (fldListFile,n,substr,rc)
         strm%domMaskName = subStr
      endif
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading domain data file path'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   startTag =  "<domainInfo>"
     endTag = "</domainInfo>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<filePath>"
     endTag = "</filePath>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   n = len_trim(str)
   if (n>0 .and. str(n:n) /= '/') str(n+1:n+2) = "/ " ! must have trailing slash
   if (n==0) str = "./ "                              ! null path => ./
   strm%domFilePath = str
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * data file path = ', trim(strm%domFilePath)

   close(nUnit) 
   !-----------------------------------------------------------------------------  
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  reading domain data file name'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   startTag =  "<domainInfo>"
     endTag = "</domainInfo>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if
   startTag =  "<fileNames>"
     endTag = "</fileNames>"
   call shr_stream_readUpToTag(nUnit,startTag,rc=rCode2)
   if (rCode2 /= 0)then
      rCode = rCode2
      goto 999
   end if

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   strm%domFileName = str
   if (debug>0 .and. s_loglev>0) write(s_logunit,F00) '  * ',trim(strm%domFileName)

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   ! get initial calendar value
   !-----------------------------------------------------------------------------  
   call shr_stream_getCalendar(strm,1,calendar)
   strm%calendar = trim(calendar)

   !-----------------------------------------------------------------------------  
   ! normal return or end-of-file problem?
   !-----------------------------------------------------------------------------  
   call shr_stream_setInit(strm)
   if ( present(rc) ) rc = rCode
   call shr_file_freeUnit(nUnit)
   return

999 continue
   write(s_logunit,F00) "ERROR: unexpected end-of-file while reading ",trim(startTag)
   write(s_logunit,F00) " error code = ", rCode
   call shr_stream_abort(subName//"ERROR: unexpected end-of-file")
   close(nUnit) 
   if ( present(rc) ) rc = rCode
   call shr_file_freeUnit(nUnit)
   
end subroutine shr_stream_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_set -- set values of stream datatype
!
! !DESCRIPTION:
!
! !REMARKS:
!    set or override stream settings
!
! !REVISION HISTORY:
!     2010-Apr-20 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_set(strm,yearFirst,yearLast,yearAlign,offset,taxMode, &
                          fldListFile,fldListModel,domFilePath,domFileName, &
                          domTvarName,domXvarName,domYvarName,domAreaName,domMaskName, &
                          filePath,filename,dataSource,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)    ,intent(inout) :: strm      ! data stream
   integer  (SHR_KIND_IN),optional,intent(in)    :: yearFirst ! first year to use 
   integer  (SHR_KIND_IN),optional,intent(in)    :: yearLast  ! last  year to use 
   integer  (SHR_KIND_IN),optional,intent(in)    :: yearAlign ! align yearFirst with this model year
   integer  (SHR_KIND_IN),optional,intent(in)    :: offset    ! offset in seconds of stream data
   character(*)          ,optional,intent(in)    :: taxMode   ! time axis mode
   character(*)          ,optional,intent(in)    :: fldListFile  ! file field names, colon delim list
   character(*)          ,optional,intent(in)    :: fldListModel ! model field names, colon delim list
   character(*)          ,optional,intent(in)    :: domFilePath  ! domain file path
   character(*)          ,optional,intent(in)    :: domFileName  ! domain file name
   character(*)          ,optional,intent(in)    :: domTvarName  ! domain time dim name
   character(*)          ,optional,intent(in)    :: domXvarName  ! domain x dim name
   character(*)          ,optional,intent(in)    :: domYvarName  ! domain y dim nam
   character(*)          ,optional,intent(in)    :: domAreaName  ! domain area name
   character(*)          ,optional,intent(in)    :: domMaskName  ! domain mask name
   character(*)          ,optional,intent(in)    :: filePath   ! path for filenames
   character(*)          ,optional,intent(in)    :: filename(:)  ! input filenames
   character(*)          ,optional,intent(in)    :: dataSource ! comment line
   integer  (SHR_KIND_IN),optional,intent(out)   :: rc        ! return code
 
!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: n
   character(SHR_KIND_CL) :: calendar      ! stream calendar

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_set) ' 
   character(*),parameter :: F00   = "('(shr_stream_set) ',8a)" 
   character(*),parameter :: F01   = "('(shr_stream_set) ',1a,i6)" 

!-------------------------------------------------------------------------------

   call shr_stream_default(strm)

   if ( present(rc) ) rc = 0

   if (present(yearFirst)) then
      strm%yearFirst = yearFirst
   endif
   if (present(yearLast)) then
      strm%yearLast = yearLast
   endif
   if (present(yearAlign)) then
      strm%yearAlign = yearAlign
   endif
   if (present(offset)) then
      strm%offset = offset
   endif
   if (present(taxMode)) then
      strm%taxMode = trim(taxMode)
   endif
   if (present(fldListFile)) then
      strm%fldListFile = trim(fldListFile)
   endif
   if (present(fldListModel)) then
      strm%fldListModel = trim(fldListModel)
   endif
   if (present(domFilePath)) then
      strm%domFilePath = trim(domFilePath)
   endif
   if (present(domFileName)) then
      strm%domFileName = trim(domFileName)
   endif
   if (present(domTvarName)) then
      strm%domTvarName = trim(domTvarName)
   endif
   if (present(domXvarName)) then
      strm%domXvarName = trim(domXvarName)
   endif
   if (present(domYvarName)) then
      strm%domYvarName = trim(domYvarName)
   endif
   if (present(domAreaName)) then
      strm%domAreaName = trim(domAreaName)
   endif
   if (present(domMaskName)) then
      strm%domMaskName = trim(domMaskName)
   endif
   if (present(filePath)) then
      strm%filePath = trim(filePath)
   endif
   if (present(filename)) then
      write(s_logunit,F01) "size of filename = ",size(filename)
      write(s_logunit,F00) "filename = ",filename
      do n = 1,size(filename)
         if (trim(filename(n)) == trim(shr_stream_file_null)) then
            ! ignore it
         else
            if (n > nFileMax) then
               write(s_logunit,F00) "ERROR: exceeded max number of files"
               call shr_stream_abort(subName//"ERROR: exceeded max number of files")
               if ( present(rc) ) rc = 1
               return
            endif
            strm%nFiles = n
            strm%file(n)%name = trim(filename(n))
         endif
      enddo
   endif

   !-----------------------------------------------------------------------------  
   ! get initial calendar value
   !-----------------------------------------------------------------------------  
   call shr_stream_getCalendar(strm,1,calendar)
   strm%calendar = trim(calendar)

   call shr_stream_setInit(strm)

end subroutine shr_stream_set

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_default -- set defaults for stream
!
! !DESCRIPTION:
!
! !REMARKS:
!    set basic default values for streams
!
! !REVISION HISTORY:
!     2010-Oct-20 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_default(strm,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)    ,intent(inout) :: strm      ! data stream
   integer  (SHR_KIND_IN),optional,intent(out)   :: rc        ! return code
 
!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: n

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_default) ' 
   character(*),parameter :: F00   = "('(shr_stream_default) ',8a)" 
   character(*),parameter :: F01   = "('(shr_stream_default) ',1a,i6)" 

!-------------------------------------------------------------------------------

   !-----------------------------------------------------------------------------  
   ! set default values for everything in stream
   !-----------------------------------------------------------------------------  
   call shr_stream_clearInit(strm)
   strm%nFiles           = 0  
   strm%dataSource       = 'undefined'
   strm%filePath         = ' '
   do n=1,nFileMax
      strm%file%name     = trim(shr_stream_file_null) 
      strm%file%haveData = .false.
      strm%file%nt       = 0
   !  strm%file%date     = undefined  ! note: unallocated
   !  strm%file%secs     = undefined  ! note: unallocated
   end do

   strm%yearFirst        = 0
   strm%yearLast         = 0
   strm%yearAlign        = 0
   strm%offset           = 0
   strm%taxMode          = trim(shr_stream_taxis_cycle)

   strm%k_lvd            = -1 
   strm%n_lvd            = -1 
   strm%found_lvd        = .false.
   strm%k_gvd            = -1 
   strm%n_gvd            = -1 
   strm%found_gvd        = .false.

   strm%fldListFile      = ' '
   strm%fldListModel     = ' '
   strm%domFilePath      = ' '
   strm%domFileName      = ' '
   strm%domTvarName      = ' '
   strm%domXvarName      = ' '
   strm%domYvarName      = ' '
   strm%domAreaName      = ' '
   strm%domMaskName      = ' '

   strm%calendar         = shr_cal_noleap

   if ( present(rc) ) rc = 0

end subroutine shr_stream_default
!===============================================================================

subroutine shr_stream_readUpToTag(nUnit,tag,optionalTag,rc)

   !----- input/output -----
   integer(SHR_KIND_IN),intent(in ) :: nUnit       ! i/o unit to read from
   character(*)        ,intent(in ) :: tag         ! string to search for
   logical, optional   ,intent(in ) :: optionalTag ! this is an optional tag
   integer(SHR_KIND_IN),intent(out) :: rc          ! return code

   !----- local -----
   character(SHR_KIND_CL)           :: str   ! temp char string
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
      call shr_stream_abort(subName//"ERROR: tag not found")
   end if
   
end subroutine shr_stream_readUpToTag

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_parseInput -- extract fileName,yearAlign, etc. from a string
!
! !DESCRIPTION:
!     shr_stream_parseInput -- extract fileName,yearAlign, etc. from a string
!
! !REMARKS:
!    should input be via standard Fortran namelist?
!
! !REVISION HISTORY:
!     2007-Aug-01 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_parseInput(str,fileName,yearAlign,yearFirst,yearLast,rc)

! !INPUT/OUTPUT PARAMETERS:

   character(*)                   ,intent(in)  :: str       ! string to parse
   character(*)                   ,intent(out) :: fileName  ! file name
   integer  (SHR_KIND_IN)         ,intent(out) :: yearFirst ! first year to use 
   integer  (SHR_KIND_IN)         ,intent(out) :: yearLast  ! last  year to use 
   integer  (SHR_KIND_IN)         ,intent(out) :: yearAlign ! align yearFirst with this model year
   integer  (SHR_KIND_IN),optional,intent(out) :: rc        ! return code
 
!EOP

   !----- local -----
   integer  (SHR_KIND_IN) :: n        ! generic index
   character(SHR_KIND_CL) :: str2     ! temp work string

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_parseInput) ' 
   character(*),parameter :: F00   = "('(shr_stream_parseInput) ',8a)" 
   character(*),parameter :: F01   = "('(shr_stream_parseInput) ',a,3i10)" 

!-------------------------------------------------------------------------------
! notes: 
! - this routine exists largely because of the difficulty of reading file names
!   that include dir paths, ie. containing "/", from char strings 
!   because the "/" is interpreted as an end-of-record.
!-------------------------------------------------------------------------------

   if (debug>1 .and. s_loglev > 0) write(s_logunit,F00) "str       = ",trim(str)

   str2 = adjustL(str)
   n    = index(str2," ")
   fileName = str2(:n)
   read(str2(n:),*) yearAlign,yearFirst,yearLast

   if (debug>1 .and. s_loglev > 0) then
      write(s_logunit,F00) "fileName  = ",trim(fileName)
      write(s_logunit,F01) "yearAlign = ",yearAlign
      write(s_logunit,F01) "yearFirst = ",yearFirst
      write(s_logunit,F01) "yearLast  = ",yearLast
   end if

   if (present(rc)) rc = 0
   
end subroutine shr_stream_parseInput

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_findBounds -- find stream data bounding a model date
!
! !DESCRIPTION:
!    Given a stream and a model date, find time coordinates of the upper and
!    lower time bounds surrounding the models date.  Returns the model date,
!    data date, elasped seconds, time index, and file names associated with
!    these upper and lower time bounds.
!
! !REVISION HISTORY:
!     2009-Sep-01 - T. Craig - modified
!     2005-Apr-01 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_findBounds(strm,mDateIn,        secIn,               &
                                  &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                                  &   mDateUB,dDateUB,secUB,n_ub,fileUB    )

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(inout):: strm    ! data stream to query
   integer(SHR_KIND_IN)       ,intent(in)   :: mDateIn ! model date (yyyymmdd)
   integer(SHR_KIND_IN)       ,intent(in)   ::   secIn ! elapsed sec on model date
   integer(SHR_KIND_IN)       ,intent(out)  :: mDateLB ! model date    of LB
   integer(SHR_KIND_IN)       ,intent(out)  :: dDateLB ! data  date    of LB
   integer(SHR_KIND_IN)       ,intent(out)  ::   secLB ! elap sec      of LB
   integer(SHR_KIND_IN)       ,intent(out)  ::    n_lb ! t-coord index of LB
   character(*)               ,intent(out)  ::  fileLB ! file containing  LB
   integer(SHR_KIND_IN)       ,intent(out)  :: mDateUB ! model date    of UB
   integer(SHR_KIND_IN)       ,intent(out)  :: dDateUB ! data  date    of UB
   integer(SHR_KIND_IN)       ,intent(out)  ::   secUB ! elap sec      of UB
   integer(SHR_KIND_IN)       ,intent(out)  ::    n_ub ! t-coord index of UB
   character(*)               ,intent(out)  ::  fileUB ! file containing  UB

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: fileName      ! string
   integer  (SHR_KIND_IN) :: nt            ! size of a time-coord dimension
   integer  (SHR_KIND_IN) :: dDateIn       ! model date mapped onto a data date
   integer  (SHR_KIND_IN) :: dDateF        ! first date
   integer  (SHR_KIND_IN) :: dDateL        ! last date
   integer  (SHR_KIND_IN) :: n,nf          ! loop index wrt t-coord array within one file
   integer  (SHR_KIND_IN) :: k,kf          ! loop index wrt list of files
   integer  (SHR_KIND_IN) :: k_ub,k_lb     ! file index of U/L bounds
   integer  (SHR_KIND_IN) :: rCode         ! return code

   integer  (SHR_KIND_IN) :: mYear         ! year of model date
   integer  (SHR_KIND_IN) :: yrFirst       ! first year of data loop
   integer  (SHR_KIND_IN) :: yrLast        ! last year of data loop
   integer  (SHR_KIND_IN) :: yrAlign       ! model year that aligns with yearFirst
   integer  (SHR_KIND_IN) :: nYears        ! number of years in data loop
   integer  (SHR_KIND_IN) :: dYear         ! data year corresponding to model year
   integer  (SHR_KIND_IN) :: yy,mm,dd      ! year,month,day
   real     (SHR_KIND_R8) :: rDateIn       ! model dDateIn + secs/(secs per day)
   real     (SHR_KIND_R8) :: rDate1        ! stream dDateIn + secs/(secs per day)
   real     (SHR_KIND_R8) :: rDate2        ! stream dDateIn + secs/(secs per day)
   real     (SHR_KIND_R8) :: rDatelvd      ! lvd dDate + secs/(secs per day)
   real     (SHR_KIND_R8) :: rDategvd      ! gvd dDate + secs/(secs per day)
   logical                :: cycle         ! is cycling on or off
   logical                :: limit         ! is limiting on or off

   !----- formats -----
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

   if (debug>0 .and. s_loglev > 0) write(s_logunit,F02) "DEBUG: ---------- enter ------------------"

   rCode = 0
   if ( .not. shr_stream_isInit(strm)) then
      rCode = 1
      call shr_stream_abort(trim(subName)//" ERROR: trying to find bounds of uninitialized stream")
      return
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
      call shr_stream_abort(trim(subName)//' ERROR: illegal taxMode = '//trim(strm%taxMode))
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
      call shr_stream_abort(trim(subName)//' ERROR: dyear lt one')
   endif

   dDateIn = dYear*10000 + modulo(mDateIn,10000) ! mDateIn mapped to range of data years
   rDateIn = dDateIn + secIn/spd                 ! dDateIn + fraction of a day

!   write(s_logunit,*) 'tcx fbd1 ',mYear,dYear,dDateIn,rDateIn
!   write(s_logunit,*) 'tcx fbd2 ',yrFirst,yrLast,yrAlign,nYears
!   call shr_sys_flush(s_logunit)

   !----------------------------------------------------------------------------
   ! find least valid date (lvd)
   !----------------------------------------------------------------------------

   if (.not. strm%found_lvd) then
A:    do k=1,strm%nFiles
         if (.not. strm%file(k)%haveData) then  
            call shr_stream_readtCoord(strm, k, rCode)
            if ( rCode /= 0 )then
               call shr_stream_abort(trim(subName)//" ERROR: readtCoord1")
               return
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
         rCode = 1
         write(s_logunit,F00)  "ERROR: LVD not found, all data is before yearFirst"
         call shr_stream_abort(trim(subName)//" ERROR: LVD not found, all data is before yearFirst")
      else
         !--- LVD is in or beyond yearFirst, verify it is not beyond yearLast ---
         if ( dDateL <= strm%file(strm%k_lvd)%date(strm%n_lvd) ) then
            rCode = 1
            write(s_logunit,F00)  "ERROR: LVD not found, all data is after yearLast"
            call shr_stream_abort(trim(subName)//" ERROR: LVD not found, all data is after yearLast")
         end if
      end if
      if (debug>1 .and. s_loglev > 0) then
         if (strm%found_lvd) write(s_logunit,F01) "DEBUG: found LVD = ",strm%file(k)%date(n)
      end if 
   end if

   if (strm%found_lvd) then
      k = strm%k_lvd
      n = strm%n_lvd
      rDatelvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! LVD date + frac day
   else
      write(s_logunit,F00)  "ERROR: LVD not found yet"
      call shr_stream_abort(trim(subName)//" ERROR: LVD not found yet")      
   endif

   if (strm%found_gvd) then
      k = strm%k_gvd
      n = strm%n_gvd
      rDategvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! GVD date + frac day
   else
      rDategvd = 99991231.0
   endif

!   write(s_logunit,*) 'tcx fbd3 ',rDateIn,rDatelvd,rDategvd
!   call shr_sys_flush(s_logunit)

   !-----------------------------------------------------------
   ! dateIn < rDatelvd
   !   limit -> abort
   !   extend -> use lvd value, set LB to 00000101
   !   cycle -> lvd is UB, gvd is LB, shift mDateLB by -nYears
   !-----------------------------------------------------------

   if (rDateIn < rDatelvd) then
      if (limit) then
         write(s_logunit,*)  trim(subName)," ERROR: limit on and rDateIn lt rDatelvd",rDateIn,rDatelvd
         call shr_stream_abort(trim(subName)//" ERROR: rDateIn lt rDatelvd limit true")
         return
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
                     call shr_stream_abort(trim(subName)//" ERROR: readtCoord2")
                     return
                  end if
               end if
               !--- start search at greatest date & move toward least date ---
               do n=strm%file(k)%nt,1,-1
                  if ( strm%file(k)%date(n) < dDateL ) then
                     strm%k_gvd = k
                     strm%n_gvd = n
                     strm%found_gvd = .true.
                     rDategvd = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! GVD date + frac day
                     if (debug>1 .and. s_loglev > 0) write(s_logunit,F01) "DEBUG: found GVD ",strm%file(k)%date(n)
                     exit B
                  end if
               end do
            end do B
         end if

         if (.not. strm%found_gvd) then
            write(s_logunit,F00)  "ERROR: GVD not found1"
            call shr_stream_abort(trim(subName)//" ERROR: GVD not found1")
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
         call shr_stream_abort(trim(subName)//" ERROR: rDateIn gt rDategvd limit true")
         return
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
               call shr_stream_abort(trim(subName)//" ERROR: readtCoord3")
               return
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
               call shr_stream_abort(trim(subName)//" ERROR: rDateIn gt rDategvd limit true")
               return
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

   call shr_stream_abort(trim(subName)//' ERROR: findBounds failed')
   return

end subroutine shr_stream_findBounds

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_readTCoord -- read in time coordinates with possible offset
!
! !DESCRIPTION:
!    verify time coordinate data is OK
!
! !REVISION HISTORY:
!     2009-Sep-01 - T. Craig - modified
!     2005-Apr-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_readTCoord(strm,k,rc)

   use netcdf

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(inout) :: strm ! data stream to query
   integer(SHR_KIND_IN)         ,intent(in)    :: k    ! stream index
   integer(SHR_KIND_IN),optional,intent(out)   :: rc   ! return code
 
!EOP

   !----- local -----
   character(SHR_KIND_CL) :: fileName    ! filename to read
   integer(SHR_KIND_IN)   :: nt
   integer(SHR_KIND_IN)   :: num,n
   integer(SHR_KIND_IN)   :: din,dout
   integer(SHR_KIND_IN)   :: sin,sout,offin
   integer(SHR_KIND_IN)   :: lrc
   integer(SHR_KIND_IN)   :: fid,vid,ndims,rcode
   integer(SHR_KIND_IN),allocatable :: dids(:)
   character(SHR_KIND_CS) :: units,calendar
   character(SHR_KIND_CS) :: bunits        ! time units (days,secs,...)
   integer(SHR_KIND_IN)   :: bdate         ! base date: calendar date
   real(SHR_KIND_R8)      :: bsec          ! base date: elapsed secs
   integer(SHR_KIND_IN)   :: ndate         ! calendar date of time value
   real(SHR_KIND_R8)      :: nsec          ! elapsed secs on calendar date
   real(SHR_KIND_R8),allocatable :: tvar(:)
   !----- formats -----
   character(*),parameter :: subname = '(shr_stream_readTCoord) '
   character(*),parameter :: F01   = "('(shr_stream_readTCoord) ',a,2i7)" 

!-------------------------------------------------------------------------------

   lrc = 0

   !--- need to read in this data ---
   call shr_stream_getFile(strm%filePath,strm%file(k)%name,fileName)
   rCode = nf90_open(fileName,nf90_nowrite,fid)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_open file '//trim(filename))
   rCode = nf90_inq_varid(fid,trim(strm%domTvarName),vid)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inq_varid')
   rCode = nf90_inquire_variable(fid,vid,ndims=ndims)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_variable1')
   allocate(dids(ndims))
   rCode = nf90_inquire_variable(fid,vid,dimids=dids)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_variable2')
   rCode = nf90_inquire_dimension(fid,dids(1),len=nt)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_inquire_dimension')
   deallocate(dids)

   allocate(strm%file(k)%date(nt),strm%file(k)%secs(nt))
   strm%file(k)%nt = nt

   units = ' '
   calendar = ' '
   rCode = nf90_get_att(fid, vid, 'units', units)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_att units')
   rCode = nf90_inquire_attribute(fid, vid, 'calendar')
   if (rCode == nf90_noerr) then
      rCode = nf90_get_att(fid, vid, 'calendar', calendar)
      if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_att calendar')
   else
      calendar = trim(shr_cal_noleap)
   endif
   n = len_trim(units)
   if (ichar(units(n:n)) == 0 ) units(n:n) = ' '
   n = len_trim(calendar)
   if (ichar(calendar(n:n)) == 0 ) calendar(n:n) = ' '
   call shr_string_leftalign(units)
   call shr_string_leftalign(calendar)
   call shr_string_parseCFtunit(units,bunits,bdate,bsec)
   strm%calendar = trim(shr_cal_calendarName(trim(calendar)))

   allocate(tvar(nt))
   rcode = nf90_get_var(fid,vid,tvar)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_get_var')
   rCode = nf90_close(fid)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_close')
   do n = 1,nt
      call shr_cal_advDate(tvar(n),bunits,bdate,bsec,ndate,nsec,calendar)
      strm%file(k)%date(n) = ndate
      strm%file(k)%secs(n) = nint(nsec)
   enddo
   deallocate(tvar)

   if (strm%offset /= 0) then
      if (size(strm%file(k)%date) /= size(strm%file(k)%secs)) then
!          rc = 1
          write(s_logunit,F01) "Incompatable date and secs sizes",size(strm%file(k)%date),size(strm%file(k)%secs)
          call shr_sys_abort()
      endif
      num = size(strm%file(k)%date)
      offin = strm%offset
      do n = 1,num
         din = strm%file(k)%date(n)
         sin = strm%file(k)%secs(n)
         call shr_cal_advDateInt(offin,'seconds',din,sin,dout,sout,calendar)
!        write(s_logunit,*) 'tcx debug rtc1 ',n,strm%offset,din,sin,dout,sout
         strm%file(k)%date(n) = dout
         strm%file(k)%secs(n) = sout
      enddo
   endif

   strm%file(k)%haveData = .true.
   call shr_stream_verifyTCoord(strm,k,lrc) ! check new t-coord data

   if (present(rc)) then
      rc = lrc
   endif

end subroutine shr_stream_readTCoord

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_verifyTCoord -- verify time coordinate data is OK
!
! !DESCRIPTION:
!    verify time coordinate data is OK
!
! !REVISION HISTORY:
!     2005-Apr-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_verifyTCoord(strm,k,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in) :: strm  ! data stream
   integer(SHR_KIND_IN)                   :: k     ! index of file to check
   integer(SHR_KIND_IN)                   :: rc    ! return code
 
!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n           ! generic loop index
   integer(SHR_KIND_IN)   :: nt          ! size of t-dimension
   integer(SHR_KIND_IN)   :: date1,secs1 ! date and seconds for a    time coord
   integer(SHR_KIND_IN)   :: date2,secs2 ! date and seconds for next time coord
   logical                :: checkIt     ! have data / do comparison

   !----- formats -----
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
   if (debug>1 .and. s_loglev > 0) write(s_logunit,F01) "checking t-coordinate data   for file k =",k

   if ( .not. strm%file(k)%haveData) then
      rc = 1
      write(s_logunit,F01) "Don't have data for file ",k
      call shr_stream_abort(subName//"ERROR: can't check -- file not read.")
      return
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
               if (debug>1 .and. s_loglev > 0) write(s_logunit,F01) "comparing with previous file for file k =",k
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
               if (debug>1 .and. s_loglev > 0) write(s_logunit,F01) "comparing with next     file for file k =",k
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
            call shr_stream_abort(subName//"ERROR: calendar dates must be increasing")
            return
         else if ( date1 == date2 ) then
            if ( secs1 >= secs2 ) then
               rc = 1
               write(s_logunit,F01) "ERROR: elapsed seconds on a date must be strickly increasing"
               write(s_logunit,F02) "secs(n), secs(n+1) = ",secs1,secs2
               call shr_stream_abort(subName//"ERROR: elapsed seconds must be increasing")
               return
            end if
         end if
         if ( secs1 < 0 .or. spd < secs1 ) then
            rc = 1
            write(s_logunit,F01) "ERROR: elapsed seconds out of valid range [0,spd]"
            write(s_logunit,F02) "secs(n) = ",secs1
            call shr_stream_abort(subName//"ERROR: elapsed seconds out of range")
            return
         end if
      end if
   end do

   if (debug>0 .and. s_loglev > 0) write(s_logunit,F01) "data is OK (non-decreasing)  for file k =",k

end subroutine shr_stream_verifyTCoord

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getFileFieldList -- Get list of file fields
!
! !DESCRIPTION:
!     Get list of file fields
!     \newline
!     call shr\_stream\_getFileFieldList(stream,list,rc)
!
! !REVISION HISTORY:
!     2005-May-10 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getFileFieldList(stream,list,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: stream  ! stream in question
   character(*)                 ,intent(out) :: list    ! field list
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_getFileFieldList) '
   character(*),parameter :: F00   = "('(shr_stream_getFileFieldList) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   list = stream%fldListFile

   if (present(rc)) rc = rCode

end subroutine shr_stream_getFileFieldList

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getModelFieldList -- Get list of file fields
!
! !DESCRIPTION:
!     Get list of file fields
!     \newline
!     call shr\_stream\_getModelFieldList(stream,list,rc)
!
! !REVISION HISTORY:
!     2005-May-10 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getModelFieldList(stream,list,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: stream  ! stream in question
   character(*)                 ,intent(out) :: list    ! field list
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_getModelFieldList) '
   character(*),parameter :: F00   = "('(shr_stream_getModelFieldList) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   list = stream%fldListModel

   if (present(rc)) rc = rCode

end subroutine shr_stream_getModelFieldList

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getFileFieldName -- Get name of k-th field in list
!
! !DESCRIPTION:
!     Get name of k-th field in list
!     \newline
!     call shr\_stream\_getFileFieldName(stream,k,name,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getFileFieldName(stream,k,name,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: stream  ! stream in question
   integer(SHR_KIND_IN)         ,intent(in)  :: k       ! index of field
   character(*)                 ,intent(out) :: name    ! k-th name in list
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_getFileFieldName) '
   character(*),parameter :: F00   = "('(shr_stream_getFileFieldName) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   call shr_string_listGetName(stream%fldListFile,k,name,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_stream_getFileFieldName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getModelFieldName -- Get name of k-th field in list
!
! !DESCRIPTION:
!     Get name of k-th field in list
!     \newline
!     call shr\_stream\_getModelFieldName(stream,k,name,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getModelFieldName(stream,k,name,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: stream  ! stream in question
   integer(SHR_KIND_IN)         ,intent(in)  :: k       ! index of field
   character(*)                 ,intent(out) :: name    ! k-th name in list
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_getModelFieldName) '
   character(*),parameter :: F00   = "('(shr_stream_getModelFieldName) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   call shr_string_listGetName(stream%fldListModel,k,name,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_stream_getModelFieldName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getFilePath -- return file path
!
! !DESCRIPTION:
!    Returns file path.
!
! !REVISION HISTORY:
!     2005-Nov-23 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getFilepath(strm,path)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm      ! data stream
   character(*)               ,intent(out) :: path      ! file path
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   path = strm%filePath

end subroutine shr_stream_getFilePath

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getDataSource -- return data source meta data
!
! !DESCRIPTION:
!    Returns data source meta data.
!
! !REVISION HISTORY:
!     2005-Feb-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getDataSource(strm,str)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm      ! data stream
   character(*)               ,intent(out) :: str       ! meta data
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   str = strm%dataSource

end subroutine shr_stream_getDataSource

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getCalendar -- return calendar name
!
! !DESCRIPTION:
!    Returns calendar name
!
! !REVISION HISTORY:
!     2010-Oct-11 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getCalendar(strm,k,calendar)

   use netcdf

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm     ! data stream
   integer(SHR_KIND_IN)       ,intent(in)  :: k        ! file to query
   character(*)               ,intent(out) :: calendar ! calendar name
 
!EOP

   integer(SHR_KIND_IN)   :: fid, vid, n
   character(SHR_KIND_CL) :: fileName,strmfile,lcal
   integer(SHR_KIND_IN)   :: rCode
   character(*),parameter :: subName = '(shr_stream_getCalendar) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lcal = ' '
   calendar = ' '
   if (k > strm%nfiles) call shr_sys_abort(subname//' ERROR: k gt nfiles')
   strmfile = strm%file(k)%name
   call shr_stream_getFile(strm%filePath,strmfile,fileName)
   rCode = nf90_open(fileName,nf90_nowrite,fid)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_open file '//trim(filename))
   rCode = nf90_inq_varid(fid,trim(strm%domTvarName),vid)
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
   call shr_string_leftalign(lcal)
   calendar = trim(shr_cal_calendarName(trim(lcal)))
   rCode = nf90_close(fid)
   if (rcode /= nf90_noerr) call shr_sys_abort(subname//' ERROR: nf90_close')

   return

end subroutine shr_stream_getCalendar

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getDomainInfo -- return domain information
!
! !DESCRIPTION:
!    Returns domain information data.
!
! !REVISION HISTORY:
!     2005-Mar-13 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getDomainInfo(strm,filePath,fileName,timeName,lonName,latName,maskName,areaName)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm     ! data stream
   character(*)               ,intent(out) :: filePath ! domain file path
   character(*)               ,intent(out) :: fileName ! domain file name
   character(*)               ,intent(out) :: timeName ! domain time var name
   character(*)               ,intent(out) ::  lonName ! domain lon  var name
   character(*)               ,intent(out) ::  latName ! domain lat  var name
   character(*)               ,intent(out) :: maskName ! domain mask var name
   character(*)               ,intent(out) :: areaName ! domain area var name
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   filePath = strm%domFilePath
   fileName = strm%domFileName
   timeName = strm%domTvarName 
    lonName = strm%domXvarName 
    latName = strm%domYvarName 
   maskName = strm%domMaskName 
   areaName = strm%domAreaName 

end subroutine shr_stream_getDomainInfo

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getFile -- Acquire file, return name of file to open
!
! !DESCRIPTION:
!     Acquire file (if necessary) and return name of file to open
!     \newline
!     call shr\_stream\_getFile(path,fileName,localFileName,rc)
!
! !REVISION HISTORY:
!     2007-Aug-24 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getFile(filePath,fileName,localFile,rc)

   use shr_file_mod, only: shr_file_queryPrefix, shr_file_noPrefix

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)    :: filePath  ! file path
   character(*)                 ,intent(inout) :: fileName  ! file name
   character(*)        ,optional,intent(out)   :: localFile ! name of acquired file
   integer(SHR_KIND_IN),optional,intent(out)   :: rc        ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: localFn ! name of acquired file
   integer  (SHR_KIND_IN) :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_getFile) '
   character(*),parameter :: F00   = "('(shr_stream_getFile) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
! - this routine reflects an added stream file handling requirement: 
!   for files on an nfs-mounted file system (available via unix cp), 
!   there are two options...
!   1) read the file without making a local copy: read path/file
!   2) copy path/file to file, and then read file
! - the shr_file_get/put file name format is used to select the option:
!   using shr_file_queryPrefix -- if recognized prefix found -- do shr_file_get
!   otherwise use the file in place.
! - if   optional argument localFile is present 
!   then fileName is unaltered and localFile is the file to be read
!   else fileName is altered and contains the name of the file to be read
! - this routine is somewhat awkward but reduces redundant code
!-------------------------------------------------------------------------------

   rCode = 0

   if ( shr_file_queryPrefix(filePath) /= shr_file_noPrefix ) then
      localFn = fileName
      call shr_file_get(rCode,localFn, trim(filePath)//fileName)
   else                              ! don't copy file, read original file
      localFn = trim(filePath)//fileName
   end if

   if (debug>0 .and. s_loglev > 0) then
      write(s_logunit,F00) "DEBUG: remote file : ",trim(filePath)//trim(fileName) 
      write(s_logunit,F00) "DEBUG: local  file : ",trim(localFn)
   end if

   if (.not. present(localFile)) fileName  = localFn ! clobber input fileName
   if (      present(localFile)) localFile = localFn ! don't clobber fileName

   if (present(rc)) rc = rCode

end subroutine shr_stream_getFile

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getFirstFileName -- returns first file name 
!
! !DESCRIPTION:
!    Returns first file name in stream.
!
! !REVISION HISTORY:
!     2005-Feb-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getFirstFileName(strm,file,path)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm      ! data stream
   character(*)               ,intent(out) :: file      ! file name
   character(*),optional      ,intent(out) :: path      ! file path
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(path)) path = strm%filePath
   file = strm%file(1)%name

end subroutine shr_stream_getFirstFileName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getNextFileName -- returns next file name in sequence
!
! !DESCRIPTION:
!    Returns next file name in sequence
!
! !REVISION HISTORY:
!     2005-Nov-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

   subroutine shr_stream_getNextFileName(strm,fn,fnNext,path,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)    :: strm   ! data stream
   character(*)               ,intent(in)    :: fn     ! file name
   character(*)               ,intent(out)   :: fnNext ! next file name
   character(*),optional      ,intent(out)   :: path   ! file path
   integer(SHR_KIND_IN),optional,intent(out) :: rc     ! return code
 
!EOP

   !--- local ---
   integer  (SHR_KIND_IN) :: rCode   ! return code
   integer(SHR_KIND_IN) :: n      ! loop index
   logical              :: found  ! file name found?

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_getNextFileName) ' 
   character(*),parameter :: F00   = "('(shr_stream_getNextFileName) ',8a)" 

!-------------------------------------------------------------------------------
! Note: will wrap-around data loop if lvd & gvd are known
! otherwise may return file name = "unknown"
!-------------------------------------------------------------------------------

   rCode = 0
   if (present(path)) path = strm%filePath

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
      call shr_stream_abort(subName//"ERROR: file name not in stream: "//trim(fn))
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
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getPrevFileName -- returns previous file name in sequence
!
! !DESCRIPTION:
!    Returns previous file name in sequence
!
! !REVISION HISTORY:
!     2005-Nov-18 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

   subroutine shr_stream_getPrevFileName(strm,fn,fnPrev,path,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: strm   ! data stream
   character(*)                 ,intent(in)  :: fn     ! file name
   character(*)                 ,intent(out) :: fnPrev ! preciding file name
   character(*),optional        ,intent(out) :: path   ! file path
   integer(SHR_KIND_IN),optional,intent(out) :: rc     ! return code
 
!EOP

   !--- local ---
   integer  (SHR_KIND_IN) :: rCode ! return code
   integer(SHR_KIND_IN)   :: n     ! loop index
   logical                :: found ! file name found?

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_getPrevFileName) ' 
   character(*),parameter :: F00   = "('(shr_stream_getPrevFileName) ',8a)" 

!-------------------------------------------------------------------------------
! Note: will wrap-around data loop if lvd & gvd are known
! otherwise may return file name = "unknown"
!-------------------------------------------------------------------------------

   rCode = 0
   if (present(path)) path = strm%filePath

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
      call shr_stream_abort(subName//"ERROR: file name not in stream: "//trim(fn))
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
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getNFiles -- returns number of input files in stream
!
! !DESCRIPTION:
!    Returns number of input files in stream
!
! !REVISION HISTORY:
!     2010-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getNFiles(strm,nfiles)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm      ! data stream
   integer(SHR_KIND_IN)       ,intent(out) :: nfiles    ! number of input files in stream
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nfiles = strm%nfiles

end subroutine shr_stream_getNFiles

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_restWrite -- write stream data to a restart file
!
! !DESCRIPTION:
!    Write stream data to a restart file.
!
! !REVISION HISTORY:
!     2005-Nov-21 -- B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine  shr_stream_restWrite(strm,fileName,caseName,caseDesc,nstrms,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(in)  :: strm(:)    ! vector of data streams
   character(*)                 ,intent(in)  :: fileName   ! name of restart file
   character(*)                 ,intent(in)  :: caseName   ! case name
   character(*)                 ,intent(in)  :: caseDesc   ! case description
   integer(SHR_KIND_IN),optional,intent(in)  :: nstrms     ! number of streams
   integer(SHR_KIND_IN),optional,intent(out) :: rc         ! return code

!EOP

   !--- local ---
   integer  (SHR_KIND_IN) :: rCode       ! return code
   integer(SHR_KIND_IN)   :: nStreams    ! number of streams 
   integer(SHR_KIND_IN)   :: k,n         ! generic loop index
   character( 8)          :: dStr        ! F90 wall clock date str yyyymmdd
   character(10)          :: tStr        ! F90 wall clock time str hhmmss.sss
   character(SHR_KIND_CS) :: str         ! generic text string
   integer(SHR_KIND_IN)   :: nUnit       ! a file unit number
   integer(SHR_KIND_IN)   :: nt          ! number of time samples
   character(SHR_KIND_CS) :: tInterpAlgo ! for backwards compatability

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_restWrite) '
   character(*),parameter :: F00   = "('(shr_stream_restWrite) ',16a) "
   character(*),parameter :: F01   = "('(shr_stream_restWrite) ',a,i5,a,5a) "
   character(*),parameter :: F02   = "('(shr_stream_restWrite) ',a,i5,a,5i8) "
   character(*),parameter :: F03   = "('(shr_stream_restWrite) ',a,i5,a,5l3) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0
   tInterpAlgo = 'unused'

   if (present(nstrms)) then
      if (size(strm) < nstrms) then
         write(s_logunit,F02) "ERROR: nstrms too large for strm",size(strm),nstrms
         call shr_stream_abort(subname//": ERROR: nstrms too large for strm")
      endif
      nStreams = nstrms
   else
      nStreams = size(strm)
   endif
   call date_and_time(dStr,tStr)

   !--- log info to stdout ---
   if (s_loglev > 0) then
      write(s_logunit,F00) "case name        : ",trim(caseName)
      write(s_logunit,F00) "case description : ",trim(caseDesc)
      write(s_logunit,F00) "File created     : ",dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '  &
                              &     //tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
      write(s_logunit,F01) "Number of streams ",nStreams
   endif

   !----------------------------------------------------------------------------
   ! write the  data
   !----------------------------------------------------------------------------

   nUnit = shr_file_getUnit() ! get an unused unit number
   open(nUnit,file=trim(fileName),form="unformatted",action="write")

   str =        "case name        : "//caseName
   write(nUnit) str
   str =        "case description : "//caseDesc
   write(nUnit) str
   str =        'File created     : '//dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '  &
                              &      //tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
   write(nUnit) str

   write(nUnit) nStreams
   do k = 1,nStreams
      if (.not. shr_stream_isInit(strm(k))) then  ! has stream been initialized?
         rCode = 1
         write(s_logunit,F01) "ERROR: can't write uninitialized stream to a restart file, k = ",k
         call shr_stream_abort(subName//": ERROR: given uninitialized stream")
      end if         

      write(nUnit) strm(k)%init         ! has stream been initialized?
      write(nUnit) strm(k)%nFiles       ! number of data files
      write(nUnit) strm(k)%dataSource   ! meta data identifying data source
      write(nUnit) strm(k)%filePath     ! remote location of files

      if (s_loglev > 0) write(s_logunit,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
      if (s_loglev > 0) write(s_logunit,F03) "* stream ",k," first have data = ",strm(k)%file(1)%haveData
      if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," first nt        = ",strm(k)%file(1)%nt
      nt = strm(k)%file(1)%nt
      if (strm(k)%file(1)%haveData) then
         if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," first date secs = ", &
                                                strm(k)%file(1)%date(1),strm(k)%file(1)%secs(1)
         if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," last  date secs = ", &
                                                strm(k)%file(1)%date(nt),strm(k)%file(1)%secs(nt)
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

      write(nUnit) strm(k)%yearFirst    ! first year to use in t-axis (yyyymmdd)
      write(nUnit) strm(k)%yearLast     ! last  year to use in t-axis (yyyymmdd)
      write(nUnit) strm(k)%yearAlign    ! align yearFirst with this model year
      write(nUnit) strm(k)%offset       ! time axis offset
!     write(nUnit) strm(k)%taxMode      ! time axis cycling mode

      write(nUnit) strm(k)%k_lvd        ! file        of least valid date
      write(nUnit) strm(k)%n_lvd        !      sample of least valid date
      write(nUnit) strm(k)%found_lvd    ! T <=> k_lvd,n_lvd have been set
      write(nUnit) strm(k)%k_gvd        ! file        of greatest valid date
      write(nUnit) strm(k)%n_gvd        !      sample of greatest valid date
      write(nUnit) strm(k)%found_gvd    ! T <=> k_gvd,n_gvd have been set

      write(nUnit) strm(k)%fldListFile  ! field list: file's  field names
      write(nUnit) strm(k)%fldListModel ! field list: model's field names
      write(nUnit) tInterpAlgo          ! unused
      write(nUnit) strm(k)%domFileName  ! domain file: name
      write(nUnit) strm(k)%domFilePath  ! domain file: path
      write(nUnit) strm(k)%domTvarName  ! domain file: time-dim var name
      write(nUnit) strm(k)%domXvarName  ! domain file: x-dim var name
      write(nUnit) strm(k)%domYvarName  ! domain file: y-dim var ame
      write(nUnit) strm(k)%domAreaName  ! domain file: area  var name
      write(nUnit) strm(k)%domMaskName  ! domain file: mask  var name

   end do

   close(nUnit)
   call shr_file_freeUnit(nUnit)
   if ( present(rc) ) rc = rCode
   
end subroutine  shr_stream_restWrite

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_restRead -- read stream data from a restart file
!
! !DESCRIPTION:
!    Read stream data to a restart file.
!    Either shr_stream_init xor shr_stream_restRead must be called
!    Do not call both routines.
!
! !REVISION HISTORY:
!     2005-Nov-21 -- B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine  shr_stream_restRead(strm,fileName,nstrms,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)  ,intent(inout) :: strm(:)  ! vector of data streams
   character(*)                 ,intent(in)    :: fileName ! name of restart file
   integer(SHR_KIND_IN),optional,intent(in)    :: nstrms   ! number of streams in strm
   integer(SHR_KIND_IN),optional,intent(out)   :: rc       ! return code

!EOP

   !--- local ---
   integer  (SHR_KIND_IN) :: rCode       ! return code
   integer(SHR_KIND_IN)   :: nStreams    ! number of streams 
   integer(SHR_KIND_IN)   :: k,n         ! generic loop index
   character(SHR_KIND_CS) :: str         ! generic text string
   integer(SHR_KIND_IN)   :: nUnit       ! a file unit number
   integer(SHR_KIND_IN)   :: inpi        ! input integer
   real(SHR_KIND_R8)      :: inpr        ! input real
   character(SHR_KIND_CXX):: inpcx       ! input char
   character(SHR_KIND_CL) :: inpcl       ! input char
   character(SHR_KIND_CS) :: inpcs       ! input char
   integer(SHR_KIND_IN)   :: nt          ! size of time dimension
   character(SHR_KIND_CS) :: tInterpAlgo ! for backwards compatability
   character(SHR_KIND_CL) :: name                            ! local variables
   integer(SHR_KIND_IN)   :: nFiles                          ! local variables
   integer(SHR_KIND_IN)   :: k_lvd, n_lvd, k_gvd, n_gvd      ! local variables
   logical                :: found_lvd, found_gvd, haveData  ! local variables
   integer(SHR_KIND_IN),pointer :: date(:),secs(:)           ! local variables
   logical                :: abort       ! abort the restart read
   logical                :: readok      ! read of restarts ok

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_restRead) '
   character(*),parameter :: F00   = "('(shr_stream_restRead) ',16a) "
   character(*),parameter :: F01   = "('(shr_stream_restRead) ',a,i5,a,5a) "
   character(*),parameter :: F02   = "('(shr_stream_restRead) ',a,i5,a,5i8) "
   character(*),parameter :: F03   = "('(shr_stream_restRead) ',a,i5,a,5l3) "
   character(*),parameter :: F04   = "('(shr_stream_restRead) ',a,4i8) "
   character(*),parameter :: F05   = "('(shr_stream_restRead) ',a,2i8,6a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0
   tInterpAlgo = 'unused'
   abort = .false.
   inpcl = ' '

   !----------------------------------------------------------------------------
   ! read the  data
   !----------------------------------------------------------------------------

   nUnit = shr_file_getUnit() ! get an unused unit number 
   open(nUnit,file=trim(fileName),form="unformatted",status="old",action="read", iostat=rCode)
   if ( rCode /= 0 )then
      call shr_file_freeUnit(nUnit)
      call shr_stream_abort(subName//": ERROR: error opening file: "//trim(fileName) )
      if ( present(rc) ) rc = rCode
      return
   end if

   read(nUnit) str         ! case name
   if (s_loglev > 0) write(s_logunit,F00) trim(str)  
   read(nUnit) str         ! case description
   if (s_loglev > 0) write(s_logunit,F00) trim(str)  
   read(nUnit) str         ! file creation date
   if (s_loglev > 0) write(s_logunit,F00) trim(str)  

   read(nUnit) nStreams
   if (present(nstrms)) then
      if (nstrms /= nStreams) then
         write(s_logunit,F02) "ERROR: nstrms ne nStreams on restart",nstrms,' ',nStreams
         call shr_stream_abort(subname//": ERROR: nstrms ne nStreams on restart")
      endif
      nStreams = nstrms
   endif
   if (s_loglev > 0) write(s_logunit,F01) "Number of streams ",nStreams

   do k = 1,nStreams
      read(nUnit) strm(k)%init         ! has stream been initialized? 
      if (.not. strm(k)%init) then 
         rCode = 1
         write(s_logunit,F01) "ERROR: uninitialized stream in restart file, k = ",k
         call shr_stream_abort(subName//": ERROR: reading uninitialized stream")
      end if         
      call shr_stream_setInit(strm(k))

      readok = .true.

      ! tcraig, don't overwrite these from input
      read(nUnit) nFiles       ! number of data files
      read(nUnit) inpcs  ! dataSource   ! meta data identifying data source
      read(nUnit) inpcl  ! filePath     ! remote location of files

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

      if (s_loglev > 0) write(s_logunit,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
      if (s_loglev > 0) write(s_logunit,F03) "* stream ",k," first have data = ",strm(k)%file(1)%haveData
      if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," first nt        = ",strm(k)%file(1)%nt
      if (strm(k)%file(1)%haveData) then
         nt = strm(k)%file(1)%nt
         if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," first date secs = ", &
                                                strm(k)%file(1)%date(1),strm(k)%file(1)%secs(1)
         if (s_loglev > 0) write(s_logunit,F02) "* stream ",k," last  date secs = ", &
                                                strm(k)%file(1)%date(nt),strm(k)%file(1)%secs(nt)
      endif

! tcraig, apr 2 2012, offset is the only field that should not change here for time axis
!      read(nUnit) strm(k)%yearFirst     ! last  year to use in t-axis (yyyymmdd)
!      read(nUnit) strm(k)%yearLast     ! last  year to use in t-axis (yyyymmdd)
!      read(nUnit) strm(k)%yearAlign    ! align yearFirst with this model year
!      read(nUnit) strm(k)%offset       ! time axis offset
      read(nUnit) inpi   ! first year to use in t-axis (yyyymmdd)
!      if (inpi /= strm(k)%yearFirst) then
!         write(s_logunit,F04) " ERROR: yearFirst disagrees ",k,strm(k)%yearFirst,inpi
!         abort=.true.
!      endif
      read(nUnit) inpi   ! last year to use in t-axis (yyyymmdd)
!      if (inpi /= strm(k)%yearLast) then
!         write(s_logunit,F04) " ERROR: yearLast disagrees ",k,strm(k)%yearLast,inpi
!         abort=.true.
!      endif
      read(nUnit) inpi   ! align year to use in t-axis (yyyymmdd)
!      if (inpi /= strm(k)%yearAlign) then
!         write(s_logunit,F04) " ERROR: yearAlign disagrees ",k,strm(k)%yearAlign,inpi
!         abort=.true.
!      endif
      read(nUnit) inpi   ! time axis offset
      if (inpi /= strm(k)%offset) then
         write(s_logunit,F04) " ERROR: offset disagrees ",k,strm(k)%offset,inpi
         abort=.true.
      endif

!     read(nUnit) strm(k)%taxMode      ! time axis cycling mode

      read(nUnit) k_lvd        ! file        of least valid date
      read(nUnit) n_lvd        !      sample of least valid date
      read(nUnit) found_lvd    ! T <=> k_lvd,n_lvd have been set
      read(nUnit) k_gvd        ! file        of greatest valid date
      read(nUnit) n_gvd        !      sample of greatest valid date
      read(nUnit) found_gvd    ! T <=> k_gvd,n_gvd have been set
      ! tcraig, april 2012, only overwrite if restart read is ok
      if (readok) then
         write(s_logunit,F05) "setting k n and found lvd gvd on restart ",k,n,' ',trim(name)
         strm(k)%k_lvd     = k_lvd
         strm(k)%n_lvd     = n_lvd
         strm(k)%found_lvd = found_lvd
         strm(k)%k_gvd     = k_gvd
         strm(k)%n_gvd     = n_gvd
         strm(k)%found_gvd = found_gvd
      endif

      ! tcraig, april 2012, don't overwrite these from input
      read(nUnit) inpcx !  fldListFile  ! field list: file's  field names
      read(nUnit) inpcx !  fldListModel ! field list: model's field names
      read(nUnit) inpcs !  tInterpAlgo          ! unused
      read(nUnit) inpcl !  domFileName  ! domain file: name
      read(nUnit) inpcl !  domFilePath  ! domain file: path
      read(nUnit) inpcs !  domTvarName  ! domain file: time-dim var name
      read(nUnit) inpcs !  domXvarName  ! domain file: x-dim var name
      read(nUnit) inpcs !  domYvarName  ! domain file: y-dim var ame
      read(nUnit) inpcs !  domAreaName  ! domain file: area  var name
      read(nUnit) inpcs !  domMaskName  ! domain file: mask  var name

   end do

   if (abort) then
      write(s_logunit,F00) "ERRORS Detected ABORTING NOW"
      call shr_stream_abort(subName//": ERRORS Detected ABORTING NOW")
   endif

   close(nUnit)
   call shr_file_freeUnit(nUnit)
   if ( present(rc) ) rc = rCode

end subroutine  shr_stream_restRead

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_dataDump -- dump all data to stdout for debugging
!
! !DESCRIPTION:
!    Dump all data to stdout for debugging
!
! !REVISION HISTORY:
!     2005-Mar-23 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_dataDump(strm)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in) :: strm      ! data stream
 
!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: k ! generic loop index

   !----- formats -----
   character(*),parameter :: subName = '(shr_stream_dataDump) ' 
   character(*),parameter :: F00   = "('(shr_stream_dataDump) ',8a)" 
   character(*),parameter :: F01   = "('(shr_stream_dataDump) ',a,3i5)" 
   character(*),parameter :: F02   = "('(shr_stream_dataDump) ',a,365i9.8)" 
   character(*),parameter :: F03   = "('(shr_stream_dataDump) ',a,365i6)" 

!-------------------------------------------------------------------------------
! notes: 
!-------------------------------------------------------------------------------

   if (s_loglev <= 0) return

   write(s_logunit,F00) "dump internal data for debugging..."

   !-----------------------------------------------------------------------------  
   ! dump internal data
   !-----------------------------------------------------------------------------  
   write(s_logunit,F01) "nFiles        = ", strm%nFiles            
   write(s_logunit,F00) "filePath      = ", trim(strm%filePath)
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
   write(s_logunit,F00) "domFileName  = ", trim(strm%domFileName)
   write(s_logunit,F00) "domFilePath  = ", trim(strm%domFilePath)
   write(s_logunit,F00) "domTvarName  = ", trim(strm%domTvarName)
   write(s_logunit,F00) "domXvarName  = ", trim(strm%domXvarName)
   write(s_logunit,F00) "domYvarName  = ", trim(strm%domYvarName)
   write(s_logunit,F00) "domAreaName  = ", trim(strm%domAreaName)
   write(s_logunit,F00) "domMaskName  = ", trim(strm%domMaskName)

end subroutine shr_stream_dataDump

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_setDebug -- Set local debug level
!
! !DESCRIPTION:
!    Set local/internal debug level, 0 = production
!    \newline
!    General Usage: call shr\_stream\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_setDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(in) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_setDebug) '
   character(*),parameter :: F00   = "('(shr_stream_setDebug) ',a) "
   character(*),parameter :: F01   = "('(shr_stream_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   debug = level
   if (s_loglev > 0) write(s_logunit,F01) "debug level reset to ",level

end subroutine shr_stream_setDebug

!==============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_getDebug -- return local/internal debug level
!
! !DESCRIPTION:
!    Return internal debug level, 0 = production
!    \newline
!    General Usage: call shr\_stream\_getDebug(level)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_getDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(out) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_getDebug) '
   character(*),parameter :: F00   = "('(shr_stream_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   level = debug

end subroutine shr_stream_getDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_setAbort -- Set abort level
!
! !DESCRIPTION:
!    Set local/internal abort level, .true. = production
!    \newline
!    General Usage: call shr\_stream\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2008-May-28  - E. Kluzek - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_setAbort(flag)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,intent(in) :: flag

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_setAbort) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   doabort = flag

end subroutine shr_stream_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_abort -- Call abort and end
!
! !DESCRIPTION:
!    Local interface for shr_stream abort calls
!    General Usage: call shr\_stream\_abort(msg)
!
! !REVISION HISTORY:
!     2008-May-28  - E. Kluzek - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_abort( msg )

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   character(len=*), optional, intent(IN) :: msg  ! Message to describe error

!EOP

   character(SHR_KIND_CL) :: lmsg

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_abort) '
   character(*),parameter :: F00   = "('(shr_stream_abort) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lmsg = ' '
  if (present(msg)) lmsg= msg

  if (doabort) then
    call shr_sys_abort(lmsg)
  else
    write(s_logunit,F00) trim(lmsg)
  endif

end subroutine shr_stream_abort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_isInit - checks if stream is initialized
!
! !DESCRIPTION:
!    Checks if stream is initialized
!
! !REVISION HISTORY:
!     2010-Oct-22 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

logical function shr_stream_isInit(strm,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),  intent(in)    :: strm
   integer(SHR_KIND_IN),optional,intent(out)   :: rc

!EOP

   !--- local ---
   character(*),parameter :: subName = "(shr_stream_isInit)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_stream_isInit = .false.
   if (size(strm%initarr) == initarr_size) then
      shr_stream_isInit = .true.
   endif

   if (present(rc)) rc = 0

end function shr_stream_isInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_setInit - Sets stream init flag to TRUE
!
! !DESCRIPTION:
!    Checks if stream is initialized
!
! !REVISION HISTORY:
!     2010-Oct-22 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_stream_setInit(strm,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),  intent(inout) :: strm
   integer(SHR_KIND_IN),optional,intent(out)   :: rc

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: ier
   character(*),parameter :: subName = "(shr_stream_setInit)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   strm%init = .true.
   deallocate(strm%initarr,stat=ier)
   allocate(strm%initarr(initarr_size))

   if (present(rc)) rc = 0

end subroutine shr_stream_setInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_clearInit - Sets stream init flag to TRUE
!
! !DESCRIPTION:
!    Checks if stream is initialized
!
! !REVISION HISTORY:
!     2010-Oct-22 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_stream_clearInit(strm,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),  intent(inout) :: strm
   integer(SHR_KIND_IN),optional,intent(out)   :: rc

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: ier
   character(*),parameter :: subName = "(shr_stream_clearInit)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   strm%init = .true.
   deallocate(strm%initarr,stat=ier)
   allocate(strm%initarr(initarr_size + 5))

   if (present(rc)) rc = 0

end subroutine shr_stream_clearInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_stream_bcast -- bcast stream
!
! !DESCRIPTION:
!    Return internal debug level, 0 = production
!    \newline
!    General Usage: call shr\_stream\_bcast(level)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_stream_bcast(stream,comm,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),  intent(inout) :: stream
   integer(SHR_KIND_IN),         intent(in)    :: comm
   integer(SHR_KIND_IN),optional,intent(out)   :: rc

!EOP

   !--- locals ---
   integer :: n,nt
   integer :: pid

   !--- formats ---
   character(*),parameter :: subName = '(shr_stream_bcast) '
   character(*),parameter :: F00   = "('(shr_stream_bcast) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if ( present(rc) ) rc = 0
   call shr_mpi_commRank(comm,pid,subName)

   call shr_mpi_bcast(stream%init        ,comm,subName)
   call shr_mpi_bcast(stream%nFiles      ,comm,subName)
   call shr_mpi_bcast(stream%dataSource  ,comm,subName)
   call shr_mpi_bcast(stream%filePath    ,comm,subName)
   call shr_mpi_bcast(stream%yearFirst   ,comm,subName)
   call shr_mpi_bcast(stream%yearLast    ,comm,subName)
   call shr_mpi_bcast(stream%yearAlign   ,comm,subName)
   call shr_mpi_bcast(stream%offset      ,comm,subName)
   call shr_mpi_bcast(stream%taxMode     ,comm,subName)
   call shr_mpi_bcast(stream%k_lvd       ,comm,subName)
   call shr_mpi_bcast(stream%n_lvd       ,comm,subName)
   call shr_mpi_bcast(stream%found_lvd   ,comm,subName)
   call shr_mpi_bcast(stream%k_gvd       ,comm,subName)
   call shr_mpi_bcast(stream%n_gvd       ,comm,subName)
   call shr_mpi_bcast(stream%found_gvd   ,comm,subName)
   call shr_mpi_bcast(stream%fldListFile ,comm,subName)
   call shr_mpi_bcast(stream%fldListModel,comm,subName)
   call shr_mpi_bcast(stream%domFileName ,comm,subName)
   call shr_mpi_bcast(stream%domFilePath ,comm,subName)
   call shr_mpi_bcast(stream%domTvarName ,comm,subName)
   call shr_mpi_bcast(stream%domXvarName ,comm,subName)
   call shr_mpi_bcast(stream%domYvarName ,comm,subName)
   call shr_mpi_bcast(stream%domMaskName ,comm,subName)
   call shr_mpi_bcast(stream%calendar    ,comm,subName)

   do n = 1,stream%nFiles
      call shr_mpi_bcast(stream%file(n)%name    ,comm,subName)
      call shr_mpi_bcast(stream%file(n)%haveData,comm,subName)
      call shr_mpi_bcast(stream%file(n)%nt      ,comm,subName)
      nt = stream%file(n)%nt
      if (pid /= 0) allocate(stream%file(n)%date(nt),stream%file(n)%secs(nt))
      call shr_mpi_bcast(stream%file(n)%date    ,comm,subName)
      call shr_mpi_bcast(stream%file(n)%secs    ,comm,subName)
   enddo

end subroutine shr_stream_bcast

!===============================================================================
end module shr_stream_mod
!===============================================================================

