!===============================================================================
! SVN $Id: shr_stream_mod.F90 1112 2006-06-02 20:59:21Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_stream_mod.F90 $
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
   use shr_ncread_mod ! shared netCDF file reading module

   implicit none

   private ! default private

! !PUBLIC TYPES:

   public :: shr_stream_streamType        ! stream data type with private components

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_stream_init              ! initialize a stream
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
   public :: shr_stream_dataDump          ! internal stream data for debugging
   public :: shr_stream_restWrite         ! write a streams restart file
   public :: shr_stream_restRead          ! read  a streams restart file
   public :: shr_stream_setDebug          ! set internal shr_stream debug level
   public :: shr_stream_getDebug          ! get internal shr_stream debug level

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !--- a useful derived type to use inside shr_stream_streamType ---
   type shr_stream_fileType
      character(SHR_KIND_CS)         :: name               ! the file name
      logical                        :: haveData = .false. ! has t-coord data been read in?
      integer  (SHR_KIND_IN)         :: nt                 ! size of time dimension
      integer  (SHR_KIND_IN),pointer :: date(:)            ! t-coord date: yyyymmdd
      integer  (SHR_KIND_IN),pointer :: secs(:)            ! t-coord secs: elapsed on date
   end type shr_stream_fileType

   !--- hard-coded array dims ~ could allocate these at run time ---
   integer(SHR_KIND_IN),parameter :: nFileMax = 1000  ! max number of files

   type shr_stream_streamType
      private                                     ! no public access to internal components
      !--- input data file names and data ---
      logical                   :: init = .false. ! has stream been initialized?
      integer  (SHR_KIND_IN)    :: nFiles = 0     ! number of data files
      character(SHR_KIND_CS)    :: dataSource     ! meta data identifying data source
      character(SHR_KIND_CL)    :: filePath       ! remote location of files
      type(shr_stream_fileType) :: file(nFileMax) ! data specific to each file

      !--- specifies how model dates align with data dates ---
      integer(SHR_KIND_IN)      :: yearFirst      ! first year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearLast       ! last  year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearAlign      ! align yearFirst with this model year

      !--- useful for quicker searching ---
      integer(SHR_KIND_IN) :: k_lvd,n_lvd         ! file/sample of least valid date
      logical              :: found_lvd = .false. ! T <=> k_lvd,n_lvd have been set
      integer(SHR_KIND_IN) :: k_gvd,n_gvd         ! file/sample of greatest valid date
      logical              :: found_gvd = .false. ! T <=> k_gvd,n_gvd have been set

      !--- stream data not used by stream module itself ---
      character(SHR_KIND_CX) :: fldListFile       ! field list: file's  field names
      character(SHR_KIND_CX) :: fldListModel      ! field list: model's field names
      character(SHR_KIND_CL) :: domFileName       ! domain file: name
      character(SHR_KIND_CS) :: domTvarName       ! domain file: time-dim var name
      character(SHR_KIND_CS) :: domXvarName       ! domain file: x-dim var name
      character(SHR_KIND_CS) :: domYvarName       ! domain file: y-dim var ame
      character(SHR_KIND_CS) :: domAreaName       ! domain file: area  var name
      character(SHR_KIND_CS) :: domMaskName       ! domain file: mask  var name
   end type shr_stream_streamType

   !----- parameters -----
   real(SHR_KIND_R8)   ,parameter :: spd = SHR_CONST_CDAY ! seconds per day

   integer(SHR_KIND_IN),save :: debug = 0  ! edit/turn-on for debug write statements

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
!     2005-Feb-03 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_init(strm,infoFile,yearFirst,yearLast,yearAlign,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType)    ,intent(out) :: strm      ! data stream
   character(*)                   ,intent(in)  :: infoFile  ! file with stream info, must read
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearFirst ! first year to use 
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearLast  ! last  year to use 
   integer  (SHR_KIND_IN)         ,intent(in)  :: yearAlign ! align yearFirst with this model year
   integer  (SHR_KIND_IN),optional,intent(out) :: rc        ! return code
 
!EOP

   !----- local -----
   integer  (SHR_KIND_IN) :: n            ! generic index
   character(SHR_KIND_CL) :: str          ! string to parse from input data file
   character(SHR_KIND_CL) :: subStr       ! sub-string of interest
   integer  (SHR_KIND_IN) :: nUnit        ! file i/o unit number
   character(SHR_KIND_CS) :: startTag     ! input file start tag
   character(SHR_KIND_CS) ::   endTag     ! input file   end tag
   character(SHR_KIND_CX) :: fldNameFile  ! field list of data file fields
   character(SHR_KIND_CX) :: fldNameModel ! field list of model     fields
   integer  (SHR_KIND_IN) :: rCode        ! return code

   !----- formats -----
   character(*),parameter :: subName = "('shr_stream_init') " 
   character(*),parameter :: F00     = "('(shr_stream_init) ',8a)" 

!-------------------------------------------------------------------------------
! notes: 
! * should this use standard namelist input?
! * needs more robust error checking
! o yearFirst,yearLast,yearAlign are provided by calling routine
! o parse infoFile for remaining, except for...
! o fileNT,fileDates, & fileSecs, which are initially set to -1, but but are replaced with
!   valid values as each file is opened for the first time
!-------------------------------------------------------------------------------

   write(6,F00) 'Reading file ',trim(infoFile)

   !-----------------------------------------------------------------------------  
   ! set default values for everything
   !-----------------------------------------------------------------------------  
   strm%nFiles           = 0  
   do n=1,nFileMax
      strm%file%name     = 'not_set' 
      strm%file%nt       = 0
      strm%file%haveData = .false.
   !  strm%file%date     = undefined  ! note: unallocated
   !  strm%file%secs     = undefined  ! note: unallocated
   end do
   strm%yearFirst        = yearFirst
   strm%yearLast         = yearLast
   strm%yearAlign        = yearAlign
   strm%fldListFile      = ''
   strm%fldListModel     = ''
   strm%filePath         = ' '
   strm%k_lvd            = -1 
   strm%n_lvd            = -1 
   strm%found_lvd        = .false.
   strm%k_gvd            = -1 
   strm%n_gvd            = -1 
   strm%found_gvd        = .false.

   !-----------------------------------------------------------------------------  
   if (debug>0) write(6,F00) '  reading data source (meta data)'
   !-----------------------------------------------------------------------------  

   nUnit = shr_sys_ioUnit() ! get unused unit number

   !--- find start tag ---
   startTag =  "<dataSource>"
     endTag = "</dataSource>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rCode)

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   strm%dataSource = str
   if (debug>0) write(6,F00) '  * format = ', trim(strm%dataSource)

   close(nUnit)

   !-----------------------------------------------------------------------------  
   if (debug>0) write(6,F00) '  reading field names'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   startTag =  "<fieldNames>"
     endTag = "</fieldNames>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rCode)

   !--- read data ---
   n=0
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      call shr_string_leftAlign(str)
      n=n+1
      if (str(1:len_trim(endTag)) == trim(endTag)) exit
      fldNameFile  = ""
      fldNameModel = ""
      read(str,*) fldNameFile,fldNameModel
      if (len_trim(fldNameFile)==0 .or. len_trim(fldNameModel)==0 ) then
         write(6,F00) "ERROR: reading field names"
         write(6,F00) '* fldNameFile  = ',trim(fldNameFile)
         write(6,F00) '* fldNameModel = ',trim(fldNameModel)
         call shr_sys_abort(subName//"ERROR: reading field names")
      end if
      if (n==1) then
         strm%fldListFile  = trim(fldNameFile )
         strm%fldListModel = trim(fldNameModel)
      else  
         strm%fldListFile  = trim(strm%fldListFile ) // ":" // trim(fldNameFile )
         strm%fldListModel = trim(strm%fldListModel) // ":" // trim(fldNameModel)
      end if
   end do
   if (debug>0) write(6,F00) '  * file  field list = ',trim(strm%fldListFile )
   if (debug>0) write(6,F00) '  * model field list = ',trim(strm%fldListModel)
   if (n==0) then
      write(6,F00) "ERROR: no input field names"
      call shr_sys_abort(subName//"ERROR: no input field names")
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0) write(6,F00) '  reading file names'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   startTag =  "<fileNames>"
     endTag = "</fileNames>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rCode)

   !--- read data ---
   n=0
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      call shr_string_leftAlign(str)
      if (str(1:len_trim(endTag)) == trim(endTag)) exit
      n=n+1
      if (n > nFileMax) then
         write(6,F00) "ERROR: exceeded max number of files"
         call shr_sys_abort(subName//"ERROR: exceeded max number of files")
      end if
      strm%file(n)%name = str
      if (debug>0) write(6,F00) '  * ',trim(strm%file(n)%name)
   end do
   strm%nFiles = n
   if (n==0) then
      write(6,F00) "ERROR: no input file names"
      call shr_sys_abort(subName//"ERROR: no input file names")
   end if

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0) write(6,F00) '  reading file path'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   startTag =  "<filePath>"
     endTag = "</filePath>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rCode)
   if (rCode /= 0) goto 999

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   call shr_string_leftAlign(str)
   n = len_trim(str)
   if (n>0 .and. str(n:n) /= '/') str(n+1:n+2) = "/ " ! must have trailing slash
   strm%filePath = str
   if (debug>0) write(6,F00) '  * path = ', trim(strm%filePath)

   close(nUnit) 

   !-----------------------------------------------------------------------------  
   if (debug>0) write(6,F00) '  reading domain file info'
   !-----------------------------------------------------------------------------  

   !--- find start tag ---
   startTag =  "<domainInfo>"
     endTag = "</domainInfo>"
   open(nUnit,file=infoFile,STATUS='OLD',FORM='FORMATTED',ACTION='READ')
   call shr_stream_readUpToTag(nUnit,startTag,rCode)
   if (rCode /= 0) goto 999

   !--- read data ---
   read(nUnit,'(a)',END=999) str
   startTag =  "<fileName>"
     endTag = "</fileName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domFileName = subStr
   if (debug>0) write(6,F00) '  * domain file  name = ', trim(strm%domFileName)

   read(nUnit,'(a)',END=999) str
   startTag =  "<tVarName>"
     endTag = "</tVarName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domTvarName = subStr
   if (debug>0) write(6,F00) '  * time dim     name = ', trim(strm%domTvarName)

   read(nUnit,'(a)',END=999) str
   startTag =  "<xVarName>"
     endTag = "</xVarName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domXvarName = subStr
   if (debug>0) write(6,F00) '  * domain x-dim name = ', trim(strm%domXvarName)

   read(nUnit,'(a)',END=999) str
   startTag =  "<yVarName>"
     endTag = "</yVarName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domYvarName = subStr
   if (debug>0) write(6,F00) '  * domain y-dim name = ', trim(strm%domYvarName)

   read(nUnit,'(a)',END=999) str
   startTag =  "<areaName>"
     endTag = "</areaName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domAreaName = subStr
   if (debug>0) write(6,F00) '  * domain area  name = ', trim(strm%domAreaName)

   read(nUnit,'(a)',END=999) str
   startTag =  "<maskName>"
     endTag = "</maskName>"
   call shr_string_betweenTags(str,startTag,endTag,subStr)
   call shr_string_leftAlign(subStr)
   strm%domMaskName = subStr
   if (debug>0) write(6,F00) '  * domain mask  name = ', trim(strm%domMaskName)

   close(nUnit)

   !-----------------------------------------------------------------------------  
   ! normal return or end-of-file problem?
   !-----------------------------------------------------------------------------  
   strm%init = .true. ! stream has been initialized
   return

999 continue
   write(6,F00) "ERROR: unexpected end-of-file while reading ",trim(startTag)
   call shr_sys_abort(subName//"ERROR: unexpected end-of-file")
   
end subroutine shr_stream_init

!===============================================================================
!===============================================================================

subroutine shr_stream_readUpToTag(nUnit,tag,rc)

   !----- input/output -----
   integer(SHR_KIND_IN),intent(in ) :: nUnit ! i/o unit to read from
   character(*)        ,intent(in ) :: tag   ! string to search for
   integer(SHR_KIND_IN),intent(out) :: rc    ! return code

   !----- local -----
   character(SHR_KIND_CL)           :: str   ! temp char string

   !----- formats -----
   character(*),parameter :: subName = "('shr_stream_readUpToTag') " 
   character(*),parameter :: F00  = "('(shr_stream_readUpToTag) ',8a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rc = 1
   do while (.true.)
      read(nUnit,'(a)',END=999) str
      if (str(1:len_trim(tag)) == trim(tag)) then
        rc = 0
        exit
      end if
   end do

999 continue

   if (rc /= 0) then
      write(6,F00) "ERROR: tag not found: ",trim(tag)
      call shr_sys_abort(subName//"ERROR: tag not found")
   end if
   
end subroutine shr_stream_readUpToTag

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
!     2005-Apr-01 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_findBounds(strm,mDateIn,        secIn,               &
                                  &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                                  &   mDateUB,dDateUB,secUB,n_ub,fileUB,rc )

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

   integer  (SHR_KIND_IN),intent(out),optional :: rc ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: fileName      ! string
   integer  (SHR_KIND_IN) :: nt            ! size of a time-coord dimension
   integer  (SHR_KIND_IN) :: dDateIn       ! model date mapped onto a data date
   integer  (SHR_KIND_IN) :: dDate,sec     ! temporary yymmdd & sec of some date
   integer  (SHR_KIND_IN) :: n             ! loop index wrt t-coord array within one file
   integer  (SHR_KIND_IN) :: k             ! loop index wrt list of files
   integer  (SHR_KIND_IN) :: k_ub,k_lb     ! file index of U/L bounds
   integer  (SHR_KIND_IN) :: rCode         ! return code

   integer  (SHR_KIND_IN) :: mYear         ! year of model date
   integer  (SHR_KIND_IN) :: year0         ! first year of data loop
   integer  (SHR_KIND_IN) :: nYears        ! number of years in data loop
   integer  (SHR_KIND_IN) :: oYear         ! model year that aligns with year0
   integer  (SHR_KIND_IN) :: dYear         ! data year corresponding to model year
   logical                :: caseA,caseB   ! flags special cases A & B
   real     (SHR_KIND_R8) :: rDateIn       ! model  dDateIn + secs/(secs per day)
   real     (SHR_KIND_R8) :: rDate         ! stream dDate   + secs/(secs per day)

   !----- formats -----
   character(*),parameter :: subName = "('shr_stream_findBoundDates') " 
   character(*),parameter :: F00  = "('(shr_stream_findBounds) ',8a)" 
   character(*),parameter :: F01  = "('(shr_stream_findBounds) ',a,i9.8,a)" 
   character(*),parameter :: F02  = "('(shr_stream_findBounds) ',a,2i9.8,i6,i5,1x,a)" 
   character(*),parameter :: F03  = "('(shr_stream_findBounds) ',a,i4)" 
   character(*),parameter :: F04  = "('(shr_stream_findBounds) ',2a,i4)" 

!-------------------------------------------------------------------------------
! Purpose:
!   1) take the model date, map it into the data date range
!   2) find the upper and lower bounding data dates
!   3) return the bounding data and model dates, file names, & t-coord indicies 
!
!  Note: special cases A & B:
!  ...after converting the model date to a data date...
!  * case A: model date is prior to the least    valid date (LVD)
!  * case B: model date is after to the greatest valid date (GVD)
!  Both cases imply the model date is in the wrap-around region of the data loop, 
!  hence the GVD is the lower bound and the LVD is the upper bound and special 
!  care needs to be taken wrt associating model dates with data dates.
!-------------------------------------------------------------------------------

   if (debug>0) write(6,F02) "DEBUG: ---------- enter ------------------"
   if (present(rc)) rc = 0

   !----------------------------------------------------------------------------
   ! convert/map the model year/date into a data year/date
   ! note: these values will be needed later to convert data year to model year
   !----------------------------------------------------------------------------
   mYear  = mDateIn/10000                        ! assumes/require F90 truncation
   year0  = strm%yearFirst                       ! 1st year in data sequence
   nYears = strm%yearLast - strm%yearFirst + 1   ! number of years in data sequence
   oYear  = strm%yearAlign                       ! model year corresponding to year0
   dYear  = year0 + modulo(mYear-oYear,nYears)   ! current data year 
   dDateIn = dYear*10000 + modulo(mDateIn,10000) ! mDateIn mapped to range of data years
   rDateIn = dDateIn + secIn/spd                 ! dDateIn + fraction of a day
   caseA = .false.                               ! flags special case A
   caseB = .false.                               ! flags special case B

   if (debug>1) write(6,F02) "DEBUG: model date, data date = ", mDateIn,dDateIn

   !----------------------------------------------------------------------------
   ! find least valid date (lvd)
   !----------------------------------------------------------------------------
   if (.not. strm%found_lvd) then
      if (debug>0) write(6,F00) "DEBUG: find least valid date"
      dDate = strm%yearFirst * 10000 + 101 ! 1st date in valid range
A:    do k=1,strm%nFiles
         if (.not. strm%file(k)%haveData) then  
            !--- need to read in this data ---
            fileName = trim(strm%filePath)//trim(strm%file(k)%name)
            call shr_ncread_varDimSizes(fileName,strm%domTvarName,nt)
            allocate(strm%file(k)%date(nt),strm%file(k)%secs(nt))
            strm%file(k)%nt = nt
            call shr_ncread_tCoord(fileName         , &
            &                      strm%domTvarName , &
            &                      strm%file(k)%date, &
            &                      strm%file(k)%secs  )
            strm%file(k)%haveData = .true.
            call shr_stream_verifyTCoord(strm,k) ! check new t-coord data
         end if
         do n=1,strm%file(k)%nt
            if ( dDate <= strm%file(k)%date(n) ) then
               !--- found least valid date ---
               strm%k_lvd = k
               strm%n_lvd = n
               strm%found_lvd = .true.
               exit A
            end if
         end do 
      end do A
      if (debug>1) then
         if (      strm%found_lvd) write(6,F01) "DEBUG: found LVD = ",strm%file(k)%date(n)
         if (.not. strm%found_lvd) write(6,F01) "DEBUG: LVD not found."
      end if 
      !--- verify there is at least one valid date in data stream ---
      dDate = (strm%yearLast + 1) * 10000 + 101 ! 1st date beyond valid range
      if ( dDate <= strm%file(strm%k_lvd)%date(strm%n_lvd) ) then
         write(6,F00) "ERROR: no valid dates in input stream"
         call shr_sys_abort("ERROR: no valid dates in input stream")
      end if
   end if

   !----------------------------------------------------------------------------
   ! special case A: wrap-around point -- GVD is lower-bound, LVD is upper-bound
   !----------------------------------------------------------------------------
   k = strm%k_lvd
   n = strm%n_lvd
   rDate = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! LVD date + frac day
   if ( rDateIn < rDate ) then
      caseA = .true.
      if (debug>1) write(6,F03) "DEBUG: special case A"

      !--- find greatest valid date (GVD) ---
      if (.not. strm%found_gvd) then
         dDate = (strm%yearLast + 1) * 10000 + 101 ! 1st date beyond valid range
         !--- start search at last file & move toward first file ---
B:       do k=strm%nFiles,1,-1
            !--- read data for file number k ---
            if (.not. strm%file(k)%haveData) then  
               !--- need to read in this data ---
               fileName = trim(strm%filePath)//trim(strm%file(k)%name)
               call shr_ncread_varDimSizes(fileName,strm%domTvarName,nt)
               allocate(strm%file(k)%date(nt),strm%file(k)%secs(nt))
               strm%file(k)%nt = nt
               call shr_ncread_tCoord(fileName         , &
               &                      strm%domTvarName , &
               &                      strm%file(k)%date, &
               &                      strm%file(k)%secs  )
               strm%file(k)%haveData = .true.
               call shr_stream_verifyTCoord(strm,k) ! check new t-coord data
            end if
            !--- start search at greatest date & move toward least date ---
            do n=strm%file(k)%nt,1,-1
               if ( strm%file(k)%date(n) < dDate ) then
                  strm%k_gvd = k
                  strm%n_gvd = n
                  strm%found_gvd = .true.
                  if (debug>1) write(6,F01) "DEBUG: found GVD ",strm%file(k)%date(n)
                  exit B
               end if
            end do
         end do B
      end if
            
      !--- GVD is lower-bound, LVD is upper-bound ---
      k_lb = strm%k_gvd
      n_lb = strm%n_gvd
      k_ub = strm%k_lvd
      n_ub = strm%n_lvd

   end if

   if (caseA) goto 999  ! case A => already know UB/LB, skip the search
                        ! could eliminate goto by creating a big if-then block

   !----------------------------------------------------------------------------
   ! search: assume lvd is *not* the UB (*not* special case A)
   !----------------------------------------------------------------------------
C: do k=strm%k_lvd,strm%nFiles
      !--- debug info ---
      if (debug>1) write(6,F03) "DEBUG: searching, k=",k
      if (debug>1) call shr_sys_flush(6)
      !--- read data for file number k ---
      if (.not. strm%file(k)%haveData) then  
         !--- need to read in this data ---
         fileName = trim(strm%filePath)//trim(strm%file(k)%name)
         call shr_ncread_varDimSizes(fileName,strm%domTvarName,nt)
         allocate(strm%file(k)%date(nt),strm%file(k)%secs(nt))
         strm%file(k)%nt = nt
         call shr_ncread_tCoord(fileName         , &
         &                      strm%domTvarName , &
         &                      strm%file(k)%date, &
         &                      strm%file(k)%secs  )
         strm%file(k)%haveData = .true.
         call shr_stream_verifyTCoord(strm,k) ! check new t-coord data
      end if
      !--- examine t-coords for file k ---
      n     = strm%file(k)%nt                                 ! last t-index in file
      rDate = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! last date + frac day
      if ( rDate <= rDateIn ) then
         !--- *all* dates in file k are lower-bounds ---
         k_lb = k
         n_lb = strm%file(k)%nt
         !--- ? special case B: upper-bound wrap-around ---
         if ( k == strm%nFiles) then
            k_ub = strm%k_lvd
            n_ub = strm%n_lvd
            caseB = .true.
            if (debug>1) write(6,F03) "DEBUG: special case B-1"
         end if
      else
         !--- the greatest lower-bound is in file k, find it ---
         do n=1,strm%file(k)%nt
            rDate = strm%file(k)%date(n) + strm%file(k)%secs(n)/spd ! date + frac day
            if ( rDate <= rDateIn ) then
               !--- found another/greater lower-bound ---
               k_lb = k
               n_lb = n
            else
               !--- found the least upper-bound ---
               if (debug>1) write(6,F01) "DEBUG: found LUB ",strm%file(k)%date(n)
               k_ub = k
               n_ub = n
               !--- is the lower-bound the greatest valid date? ---
               if (.not. strm%found_gvd) then
                  dDate = (strm%yearLast + 1)*10000 + 0101
                  if (dDate <= strm%file(k_ub)%date(n_ub)) then
                     strm%found_gvd = .true.
                     strm%k_gvd     = k_lb
                     strm%n_gvd     = n_lb
                     if (debug>1) write(6,F01) "DEBUG: found GVD ",strm%file(k_lb)%date(n_lb)
                  end if
               end if
               !--- ? special case B: upper-bound wrap-around ---
               if (strm%found_gvd) then
                  if (k_lb == strm%k_gvd .and. n_lb == strm%n_gvd ) then
                     k_ub = strm%k_lvd
                     n_ub = strm%n_lvd
                     caseB = .true.
                     if (debug>1) write(6,F03) "DEBUG: special case B-2"
                  end if
               end if
               exit C
            end if
         end do
      end if
   end do C

999 continue

   !----------------------------------------------------------------------------
   ! set values for output args, convert from data dates to model dates
   !----------------------------------------------------------------------------
   dDateLB = strm%file(k_lb)%date(n_lb)
   mDateLB = dDateLB + (mYear-dYear)*10000
     secLB = strm%file(k_lb)%secs(n_lb)
    fileLB = strm%file(k_lb)%name

   dDateUB = strm%file(k_ub)%date(n_ub)
   mDateUB = dDateUB + (mYear-dYear)*10000
     secUB = strm%file(k_ub)%secs(n_ub)
    fileUB = strm%file(k_ub)%name

   if (caseA) mDateLB = mDateLB - nYears*10000 ! adjust LB wrt data loop wrap-around
   if (caseB) mDateUB = mDateUB + nYears*10000 ! adjust UB wrt data loop wrap-around
   if (caseA .and. mDateLB<101) then
      write(6,F00) "WARNING: correct lower-bound date has negative year (unsupported)", &
      &            " reset to yy/mm/dd = 00/01/01   sec = 0"
      mDateLB = 0101  
        secLB = 0
   endif

   !----------------------------------------------------------------------------
   ! print debug info?
   !----------------------------------------------------------------------------
   if (debug>0) then
      write(6,F02) "DEBUG: LB mDate,dDate,sec,index,file =",mDateLB,dDateLB,secLB,n_lb,trim(fileLB)
      write(6,F02) "DEBUG: mDateIN,dDateIN,secIN         =",mDateIn,dDateIn,secIn
      write(6,F02) "DEBUG: UB mDate,dDate,sec,index,file =",mDateUB,dDateUB,secUB,n_ub,trim(fileUB)
      if (caseA) write(6,F03) "DEBUG: note ~ special case A"
      if (caseB) write(6,F03) "DEBUG: note ~ special case B"
      caseA = .true.
      write(6,F02) "DEBUG: ---------- exit -------------------"
   end if

   if (present(rc)) rc = rCode

end subroutine shr_stream_findBounds

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

subroutine shr_stream_verifyTCoord(strm,k)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in) :: strm  ! data stream
   integer(SHR_KIND_IN)                   :: k     ! index of file to check
 
!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n           ! generic loop index
   integer(SHR_KIND_IN)   :: nt          ! size of t-dimension
   integer(SHR_KIND_IN)   :: date1,secs1 ! date and seconds for a    time coord
   integer(SHR_KIND_IN)   :: date2,secs2 ! date and seconds for next time coord
   logical                :: checkIt     ! have data / do comparison

   !----- formats -----
   character(*),parameter :: subName = "('shr_stream_verifyTCoord') " 
   character(*),parameter :: F00     = "('(shr_stream_verifyTCoord) ',8a)" 
   character(*),parameter :: F01     = "('(shr_stream_verifyTCoord) ',a,2i7)" 
   character(*),parameter :: F02     = "('(shr_stream_verifyTCoord) ',a,2i9.8)" 

!-------------------------------------------------------------------------------
! Notes:
!   o checks that dates are increasing (must not decrease)
!   o does not check for valid dates (eg. day=0 & month = 13 are "OK")
!   o checks that secs are strictly increasing within any one day
!   o checks that 0 <= secs <= spd (seconds per day)
!   o checks all dates from one file plus last date of previous file and 
!     first date of next file
!-------------------------------------------------------------------------------

   if (debug>1) write(6,F01) "checking t-coordinate data   for file k =",k

   if ( .not. strm%file(k)%haveData) then
      write(6,F01) "Don't have data for file ",k
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
               if (debug>1) write(6,F01) "comparing with previous file for file k =",k
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
               if (debug>1) write(6,F01) "comparing with next     file for file k =",k
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
            write(6,F01) "ERROR: calendar dates must be increasing"
            write(6,F02) "date(n), date(n+1) = ",date1,date2
            call shr_sys_abort(subName//"ERROR: calendar dates must be increasing")
         else if ( date1 == date2 ) then
            if ( secs1 >= secs2 ) then
               write(6,F01) "ERROR: elapsed seconds on a date must be strickly increasing"
               write(6,F02) "secs(n), secs(n+1) = ",secs1,secs2
               call shr_sys_abort(subName//"ERROR: elapsed seconds must be increasing")
            end if
         end if
         if ( secs1 < 0 .or. spd < secs1 ) then
            write(6,F01) "ERROR: elapsed seconds out of valid range [0,spd]"
            write(6,F02) "secs(n) = ",secs1
            call shr_sys_abort(subName//"ERROR: elapsed seconds out of range")
         end if
      end if
   end do

   if (debug>0) write(6,F01) "data is OK (non-decreasing)  for file k =",k

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
   character(*),parameter :: subName =   "(shr_stream_getFileFieldList)"
   character(*),parameter :: F00     = "('(shr_stream_getFileFieldList) ',4a)"

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
   character(*),parameter :: subName =   "(shr_stream_getModelFieldList)"
   character(*),parameter :: F00     = "('(shr_stream_getModelFieldList) ',4a)"

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
   character(*),parameter :: subName =   "(shr_stream_getFileFieldName)"
   character(*),parameter :: F00     = "('(shr_stream_getFileFieldName) ',4a)"

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
   character(*),parameter :: subName =   "(shr_stream_getModelFieldName)"
   character(*),parameter :: F00     = "('(shr_stream_getModelFieldName) ',4a)"

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
! !IROUTINE: shr_stream_getDomainInfo -- return domain information
!
! !DESCRIPTION:
!    Returns domain information data.
!
! !REVISION HISTORY:
!     2005-Mar-13 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------  

subroutine shr_stream_getDomainInfo(strm,fileName,timeName,lonName,latName,maskName,areaName)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm     ! data stream
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

   fileName = trim(strm%filePath)//trim(strm%domFileName)
   timeName = strm%domTvarName 
    lonName = strm%domXvarName 
    latName = strm%domYvarName 
   maskName = strm%domMaskName 
   areaName = strm%domAreaName 

end subroutine shr_stream_getDomainInfo

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

subroutine shr_stream_getFirstFileName(strm,str)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm      ! data stream
   character(*)               ,intent(out) :: str       ! file name
 
!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   str = trim(strm%filePath) // strm%file(1)%name

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

   subroutine shr_stream_getNextFileName(strm,fn,fnNext,path)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm   ! data stream
   character(*)               ,intent(in)  :: fn     ! file name
   character(*)               ,intent(out) :: fnNext ! next file name
   character(*),optional      ,intent(out) :: path   ! file path
 
!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: n      ! loop index
   logical              :: found  ! file name found?

   !--- formats ---
   character(*),parameter :: subName =  "('shr_stream_getNextFileName') " 
   character(*),parameter :: F00     = "('(shr_stream_getNextFileName) ',8a)" 

!-------------------------------------------------------------------------------
! Note: will wrap-around data loop if lvd & gvd are known
! otherwise may return file name = "unknown"
!-------------------------------------------------------------------------------

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
      write(6,F00) "ERROR: input file name is not in stream: ",trim(fn)
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

   subroutine shr_stream_getPrevFileName(strm,fn,fnPrev,path)

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),intent(in)  :: strm   ! data stream
   character(*)               ,intent(in)  :: fn     ! file name
   character(*)               ,intent(out) :: fnPrev ! preciding file name
   character(*),optional      ,intent(out) :: path   ! file path
 
!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: n     ! loop index
   logical              :: found ! file name found?

   !--- formats ---
   character(*),parameter :: subName =  "('shr_stream_getPrevFileName') " 
   character(*),parameter :: F00     = "('(shr_stream_getPrevFileName) ',8a)" 

!-------------------------------------------------------------------------------
! Note: will wrap-around data loop if lvd & gvd are known
! otherwise may return file name = "unknown"
!-------------------------------------------------------------------------------

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
      write(6,F00) "ERROR: input file name is not in stream: ",trim(fn)
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

end subroutine shr_stream_getPrevFileName

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

subroutine  shr_stream_restWrite(strm,fileName,caseName,caseDesc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),pointer     :: strm(:)    ! vector of data streams
   character(SHR_KIND_CL)     ,intent(in)  :: fileName   ! name of restart file
   character(*)               ,intent(in)  :: caseName   ! case name
   character(*)               ,intent(in)  :: caseDesc   ! case description

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: nStreams    ! number of streams 
   integer(SHR_KIND_IN)   :: k,n         ! generic loop index
   character( 8)          :: dStr        ! F90 wall clock date str yyyymmdd
   character(10)          :: tStr        ! F90 wall clock time str hhmmss.sss
   character(SHR_KIND_CS) :: str         ! generic text string
   integer(SHR_KIND_IN)   :: nUnit       ! a file unit number

   !--- formats ---
   character(*),parameter :: subName =  "('shr_stream_restWrite') "
   character(*),parameter :: F00     = "('(shr_stream_restWrite) ',16a) "
   character(*),parameter :: F01     = "('(shr_stream_restWrite) ',a,i5,2a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nStreams = size(strm)     
   call date_and_time(dStr,tStr)

   !--- log info to stdout ---
   write(6,F00) "case name        : ",trim(caseName)
   write(6,F00) "case description : ",trim(caseDesc)
   write(6,F00) "File created     : ",dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '  &
                              &     //tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
   write(6,F01) "Number of streams ",nStreams


   !----------------------------------------------------------------------------
   ! write the  data
   !----------------------------------------------------------------------------

   nUnit = shr_sys_ioUnit() ! get an unused unit number
   open(nUnit,file=fileName,form="unformatted",action="write")

   str =        "case name        : "//caseName
   write(nUnit) str
   str =        "case description : "//caseDesc
   write(nUnit) str
   str =        'File created     : '//dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '  &
                              &      //tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)
   write(nUnit) str

   write(nUnit) nStreams
   do k = 1,nStreams
      if (.not. strm(k)%init) then  ! has stream been initialized?
         write(6,F01) "ERROR: can't write uninitialized stream to a restart file, k = ",k
         call shr_sys_abort(subName//": ERROR: given uninitialized stream")
      end if         

      write(nUnit) strm(k)%init         ! has stream been initialized?
      write(nUnit) strm(k)%nFiles       ! number of data files
      write(nUnit) strm(k)%dataSource   ! meta data identifying data source
      write(nUnit) strm(k)%filePath     ! remote location of files

      write(6,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)
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

      write(nUnit) strm(k)%k_lvd        ! file        of least valid date
      write(nUnit) strm(k)%n_lvd        !      sample of least valid date
      write(nUnit) strm(k)%found_lvd    ! T <=> k_lvd,n_lvd have been set
      write(nUnit) strm(k)%k_gvd        ! file        of greatest valid date
      write(nUnit) strm(k)%n_gvd        !      sample of greatest valid date
      write(nUnit) strm(k)%found_gvd    ! T <=> k_gvd,n_gvd have been set

      write(nUnit) strm(k)%fldListFile  ! field list: file's  field names
      write(nUnit) strm(k)%fldListModel ! field list: model's field names
      write(nUnit) strm(k)%domFileName  ! domain file: name
      write(nUnit) strm(k)%domTvarName  ! domain file: time-dim var name
      write(nUnit) strm(k)%domXvarName  ! domain file: x-dim var name
      write(nUnit) strm(k)%domYvarName  ! domain file: y-dim var ame
      write(nUnit) strm(k)%domAreaName  ! domain file: area  var name
      write(nUnit) strm(k)%domMaskName  ! domain file: mask  var name

   end do

   close(nUnit)
   
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

subroutine  shr_stream_restRead(strm,fileName)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType),pointer     :: strm(:)  ! vector of data streams
   character(*)               ,intent(in)  :: fileName ! name of restart file

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: nStreams    ! number of streams 
   integer(SHR_KIND_IN)   :: k,n         ! generic loop index
   integer(SHR_KIND_IN)   :: nt          ! size of time dimension
   character(SHR_KIND_CS) :: str         ! generic text string
   integer(SHR_KIND_IN)   :: nUnit       ! a file unit number

   !--- formats ---
   character(*),parameter :: subName =  "('shr_stream_restRead') "
   character(*),parameter :: F00     = "('(shr_stream_restRead) ',16a) "
   character(*),parameter :: F01     = "('(shr_stream_restRead) ',a,i5,2a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! write the  data
   !----------------------------------------------------------------------------

   nUnit = shr_sys_ioUnit() ! get an unused unit number 
   open(nUnit,file=fileName,form="unformatted",status="old",action="read")

   read(nUnit) str         ! case name
   write(6,F00) trim(str)  
   read(nUnit) str         ! case description
   write(6,F00) trim(str)  
   read(nUnit) str         ! file creation date
   write(6,F00) trim(str)  

   read(nUnit) nStreams
   write(6,F01) "Number of streams ",nStreams

   allocate(strm(nStreams)) ! NOTE allocation 

   do k = 1,nStreams
      read(nUnit) strm(k)%init         ! has stream been initialized? 
      if (.not. strm(k)%init) then 
         write(6,F01) "ERROR: uninitialized stream in restart file, k = ",k
         call shr_sys_abort(subName//": ERROR: reading uninitialized stream")
      end if         

      read(nUnit) strm(k)%nFiles       ! number of data files
      read(nUnit) strm(k)%dataSource   ! meta data identifying data source
      read(nUnit) strm(k)%filePath     ! remote location of files

      do n=1,strm(k)%nFiles                     ! data specific to each file...
         read(nUnit) strm(k)%file(n)%name       ! the file name
         read(nUnit) strm(k)%file(n)%haveData   ! has t-coord data been read in?
         read(nUnit) strm(k)%file(n)%nt         ! size of time dimension
         nt = strm(k)%file(n)%nt
         if (strm(k)%file(n)%haveData) then     ! ie. if arrays have been allocated
            allocate(strm(k)%file(n)%date(nt))
            allocate(strm(k)%file(n)%secs(nt))
            read(nUnit) strm(k)%file(n)%date(:) ! t-coord date: yyyymmdd
            read(nUnit) strm(k)%file(n)%secs(:) ! t-coord secs: elapsed on date
         end if
      end do

      write(6,F01) "* stream ",k," first file name = ",trim(strm(k)%file(1)%name)

      read(nUnit) strm(k)%yearFirst    ! first year to use in t-axis (yyyymmdd)
      read(nUnit) strm(k)%yearLast     ! last  year to use in t-axis (yyyymmdd)
      read(nUnit) strm(k)%yearAlign    ! align yearFirst with this model year

      read(nUnit) strm(k)%k_lvd        ! file        of least valid date
      read(nUnit) strm(k)%n_lvd        !      sample of least valid date
      read(nUnit) strm(k)%found_lvd    ! T <=> k_lvd,n_lvd have been set
      read(nUnit) strm(k)%k_gvd        ! file        of greatest valid date
      read(nUnit) strm(k)%n_gvd        !      sample of greatest valid date
      read(nUnit) strm(k)%found_gvd    ! T <=> k_gvd,n_gvd have been set

      read(nUnit) strm(k)%fldListFile  ! field list: file's  field names
      read(nUnit) strm(k)%fldListModel ! field list: model's field names
      read(nUnit) strm(k)%domFileName  ! domain file: name
      read(nUnit) strm(k)%domTvarName  ! domain file: time-dim var name
      read(nUnit) strm(k)%domXvarName  ! domain file: x-dim var name
      read(nUnit) strm(k)%domYvarName  ! domain file: y-dim var ame
      read(nUnit) strm(k)%domAreaName  ! domain file: area  var name
      read(nUnit) strm(k)%domMaskName  ! domain file: mask  var name

   end do

   close(nUnit)
   
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
   character(*),parameter :: subName = "('shr_stream_dataDump') " 
   character(*),parameter :: F00     = "('(shr_stream_dataDump) ',8a)" 
   character(*),parameter :: F01     = "('(shr_stream_dataDump) ',a,3i5)" 
   character(*),parameter :: F02     = "('(shr_stream_dataDump) ',a,365i9.8)" 
   character(*),parameter :: F03     = "('(shr_stream_dataDump) ',a,365i6)" 

!-------------------------------------------------------------------------------
! notes: 
!-------------------------------------------------------------------------------

   write(6,F00) "dump internal data for debugging..."

   !-----------------------------------------------------------------------------  
   ! dump internal data
   !-----------------------------------------------------------------------------  
   write(6,F01) "nFiles        = ", strm%nFiles            
   write(6,F00) "filePath      = ", trim(strm%filePath)
   do k=1,strm%nFiles
      write(6,F01) "data for file k = ",k
      write(6,F00)    "* file(k)%name    = ", trim(strm%file(k)%name)
      if ( strm%file(k)%haveData ) then
         write(6,F01) "* file(k)%nt      = ", strm%file(k)%nt
         write(6,F02) "* file(k)%date(:) = ", strm%file(k)%date(:)    
         write(6,F03) "* file(k)%Secs(:) = ", strm%file(k)%secs(:)    
      else
         write(6,F00) "* time coord data not read in yet for this file"
      end if
   end do
   write(6,F01) "yearF/L/A    = ", strm%yearFirst,strm%yearLast,strm%yearAlign         
   write(6,F00) "fldListFile  = ", trim(strm%fldListFile)
   write(6,F00) "fldListModel = ", trim(strm%fldListModel)
   write(6,F00) "domFileName  = ", trim(strm%domFileName)
   write(6,F00) "domTvarName  = ", trim(strm%domTvarName)
   write(6,F00) "domXvarName  = ", trim(strm%domXvarName)
   write(6,F00) "domYvarName  = ", trim(strm%domYvarName)
   write(6,F00) "domAreaName  = ", trim(strm%domAreaName)
   write(6,F00) "domMaskName  = ", trim(strm%domMaskName)

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
   character(*),parameter :: subName =  "('shr_stream_setDebug') "
   character(*),parameter :: F00     = "('(shr_stream_setDebug) ',a) "
   character(*),parameter :: F01     = "('(shr_stream_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   debug = level
   write(6,F01) "debug level reset to ",level

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
   character(*),parameter :: subName =  "('shr_stream_getDebug') "
   character(*),parameter :: F00     = "('(shr_stream_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   level = debug

end subroutine shr_stream_getDebug

!===============================================================================
end module shr_stream_mod
!===============================================================================

