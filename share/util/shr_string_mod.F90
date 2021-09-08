! !MODULE: shr_string_mod -- string and list methods
!
! !DESCRIPTION:
!    General string and specific list method.  A list is a single string
!    that is delimited by a character forming multiple fields, ie,
!    character(len=*) :: mylist = "t:s:u1:v1:u2:v2:taux:tauy"
!    The delimiter is called listDel in this module, is default ":",
!    but can be set by a call to shr_string_listSetDel.
!
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

module shr_string_mod

  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   ! F90 kinds
  use shr_sys_mod    ! shared system calls
  use shr_timer_mod, only : shr_timer_get, shr_timer_start, shr_timer_stop
  use shr_log_mod,   only : errMsg    => shr_log_errMsg
  use shr_log_mod,   only : s_loglev  => shr_log_Level
  use shr_log_mod,   only : s_logunit => shr_log_Unit

  implicit none
  private

  ! !PUBLIC TYPES:

  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_string_countChar       ! Count number of char in string, fn
  public :: shr_string_toUpper         ! Convert string to upper-case
  public :: shr_string_toLower         ! Convert string to lower-case
  public :: shr_string_getParentDir    ! For a pathname get the parent directory name
  public :: shr_string_lastIndex       ! Index of last substr in str
  public :: shr_string_endIndex        ! Index of end of substr in str
  public :: shr_string_leftalign_and_convert_tabs ! remove leading white space and convert all tabs to spaces
  public :: shr_string_convert_tabs    ! Convert all tabs to spaces
  public :: shr_string_alphanum        ! remove all non alpha-numeric characters
  public :: shr_string_betweenTags     ! get the substring between the two tags
  public :: shr_string_parseCFtunit    ! parse CF time units
  public :: shr_string_clean           ! Set string to all white space

  public :: shr_string_listIsValid     ! test for a valid "list"
  public :: shr_string_listGetNum      ! Get number of fields in list, fn
  public :: shr_string_listGetIndex    ! Get index of field
  public :: shr_string_listGetIndexF   ! function version of listGetIndex
  public :: shr_string_listGetName     ! get k-th field name
  public :: shr_string_listIntersect   ! get intersection of two field lists
  public :: shr_string_listUnion       ! get union of two field lists
  public :: shr_string_listDiff        ! get set difference of two field lists
  public :: shr_string_listMerge       ! merge two lists to form third
  public :: shr_string_listAppend      ! append list at end of another
  public :: shr_string_listPrepend     ! prepend list in front of another
  public :: shr_string_listSetDel      ! Set field delimiter in lists
  public :: shr_string_listGetDel      ! Get field delimiter in lists
  public :: shr_string_listFromSuffixes! return colon delimited field list
  ! given array of suffixes and a base string
  public :: shr_string_listCreateField ! return colon delimited field list
  ! given number of fields N and a base string
  public :: shr_string_listAddSuffix   ! add a suffix to every field in a field list
  public :: shr_string_setAbort        ! set local abort flag
  public :: shr_string_setDebug        ! set local debug flag

  ! !PUBLIC DATA MEMBERS:

  ! no public data members

  !EOP

  character(len=1)    ,save :: listDel  = ":"    ! note single exec implications
  character(len=2)    ,save :: listDel2 = "::"   ! note single exec implications
  logical             ,save :: doabort  = .true.
  integer(SHR_KIND_IN),save :: debug    = 0

  !===============================================================================
contains
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_countChar -- Count number of occurances of a character
  !
  ! !DESCRIPTION:
  !  count number of occurances of a single character in a string
  !     \newline
  !     n = shr\_string\_countChar(string,character)
  !
  ! !REVISION HISTORY:
  !     2005-Feb-28 - First version from dshr_bundle
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_countChar(str,char,rc)


    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: str   ! string to search
    character(1)        ,intent(in)           :: char  ! char to search for
    integer(SHR_KIND_IN),intent(out),optional :: rc    ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN) :: count    ! counts occurances of char
    integer(SHR_KIND_IN) :: n        ! generic index
    integer(SHR_KIND_IN) :: t01 = 0  ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_countChar) "
    character(*),parameter :: F00     = "('(shr_string_countChar) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    count = 0
    do n = 1, len_trim(str)
       if (str(n:n) == char) count = count + 1
    end do
    shr_string_countChar = count

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_countChar

  !===============================================================================
  !BOP ===========================================================================
  ! !IROUTINE: shr_string_toUpper -- Convert string to upper case
  !
  ! !DESCRIPTION:
  !     Convert the input string to upper-case.
  !     Use achar and iachar intrinsics to ensure use of ascii collating sequence.
  !
  ! !REVISION HISTORY:
  !     2005-Dec-20 - Move CAM version over to shared code.
  !
  ! !INTERFACE: ------------------------------------------------------------------

  function shr_string_toUpper(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: str      ! String to convert to upper case
    character(len=len(str))      :: shr_string_toUpper

    !----- local -----
    integer(SHR_KIND_IN) :: i             ! Index
    integer(SHR_KIND_IN) :: aseq          ! ascii collating sequence
    integer(SHR_KIND_IN) :: LowerToUpper  ! integer to convert case
    character(len=1)     :: ctmp          ! Character temporary
    integer(SHR_KIND_IN) :: t01 = 0       ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_toUpper) "
    character(*),parameter :: F00     = "('(shr_string_toUpper) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    LowerToUpper = iachar("A") - iachar("a")

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
            ctmp = achar(aseq + LowertoUpper)
       shr_string_toUpper(i:i) = ctmp
    end do

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_toUpper

  !===============================================================================
  !BOP ===========================================================================
  ! !IROUTINE: shr_string_toLower -- Convert string to lower case
  !
  ! !DESCRIPTION:
  !     Convert the input string to lower-case.
  !     Use achar and iachar intrinsics to ensure use of ascii collating sequence.
  !
  ! !REVISION HISTORY:
  !     2006-Apr-20 - Creation
  !
  ! !INTERFACE: ------------------------------------------------------------------
  function shr_string_toLower(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: str      ! String to convert to lower case
    character(len=len(str))      :: shr_string_toLower

    !----- local -----
    integer(SHR_KIND_IN) :: i            ! Index
    integer(SHR_KIND_IN) :: aseq         ! ascii collating sequence
    integer(SHR_KIND_IN) :: UpperToLower ! integer to convert case
    character(len=1)     :: ctmp         ! Character temporary
    integer(SHR_KIND_IN) :: t01 = 0      ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_toLower) "
    character(*),parameter :: F00     = "('(shr_string_toLower) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    UpperToLower = iachar("a") - iachar("A")

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
            ctmp = achar(aseq + UpperToLower)
       shr_string_toLower(i:i) = ctmp
    end do

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_toLower

  !===============================================================================
  !BOP ===========================================================================
  ! !IROUTINE: shr_string_getParentDir -- For pathname get the parent directory name
  !
  ! !DESCRIPTION:
  !     Get the parent directory name for a pathname.
  !
  ! !REVISION HISTORY:
  !     2006-May-09 - Creation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  function shr_string_getParentDir(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: str      ! String to convert to lower case
    character(len=len(str))      :: shr_string_getParentDir

    !----- local -----
    integer(SHR_KIND_IN) :: i       ! Index
    integer(SHR_KIND_IN) :: nlen    ! Length of string
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_getParentDir) "
    character(*),parameter :: F00     = "('(shr_string_getParentDir) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    nlen = len_trim(str)
    if ( str(nlen:nlen) == "/" ) nlen = nlen - 1
    i = index( str(1:nlen), "/", back=.true. )
    if ( i == 0 )then
       shr_string_getParentDir = str
    else
       shr_string_getParentDir = str(1:i-1)
    end if

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_getParentDir

  !===============================================================================
  !BOP ===========================================================================
  !
  !
  ! !IROUTINE: shr_string_lastIndex -- Get index of last substr within string
  !
  ! !DESCRIPTION:
  !  Get index of last substr within string
  !     \newline
  !     n = shr\_string\_lastIndex(string,substring)
  !
  ! !REVISION HISTORY:
  !     2005-Feb-28 - First version from dshr_domain
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_lastIndex(string,substr,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: string ! string to search
    character(*)        ,intent(in)           :: substr ! sub-string to search for
    integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_lastIndex) "
    character(*),parameter :: F00     = "('(shr_string_lastIndex) ',4a)"

    !-------------------------------------------------------------------------------
    ! Note:
    ! - "new" F90 back option to index function makes this home-grown solution obsolete
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    shr_string_lastIndex = index(string,substr,.true.)

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_lastIndex

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_endIndex -- Get the ending index of substr within string
  !
  ! !DESCRIPTION:
  !  Get the ending index of substr within string
  !     \newline
  !     n = shr\_string\_endIndex(string,substring)
  !
  ! !REVISION HISTORY:
  !     2005-May-10 - B. Kauffman, first version.
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_endIndex(string,substr,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: string ! string to search
    character(*)        ,intent(in)           :: substr ! sub-string to search for
    integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

    !EOP

    !--- local ---
    integer(SHR_KIND_IN)   :: i       ! generic index
    integer(SHR_KIND_IN)       :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_endIndex) "
    character(*),parameter :: F00     = "('(shr_string_endIndex) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! * returns zero if substring not found, uses len_trim() intrinsic
    ! * very similar to: i = index(str,substr,back=.true.)
    ! * do we need this function?
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    i = index(trim(string),trim(substr))
    if ( i == 0 ) then
       shr_string_endIndex = 0  ! substr is not in string
    else
       shr_string_endIndex = i + len_trim(substr) - 1
    end if

    !  -------------------------------------------------------------------
    !  i = index(trim(string),trim(substr),back=.true.)
    !  if (i == len(string)+1) i = 0
    !  shr_string_endIndex = i
    !  -------------------------------------------------------------------

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_endIndex

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_leftalign_and_convert_tabs -- remove leading white space and
  ! convert tabs to spaces
  !
  ! !DESCRIPTION:
  !    Remove leading white space (spaces and tabs) and convert tabs to spaces
  !    This even converts tabs in the middle or at the end of the string to spaces
  !     \newline
  !     call shr\_string\_leftalign_and_convert_tabs(string)
  !
  ! !REVISION HISTORY:
  !     2005-Apr-28 - B. Kauffman - First version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_leftalign_and_convert_tabs(str,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(inout)          :: str
    integer(SHR_KIND_IN),intent(out)  ,optional :: rc   ! return code

    !EOP

    !----- local ----
    integer(SHR_KIND_IN) :: t01 = 0 ! timer
    character, parameter :: tab_char = char(9)

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_leftalign_and_convert_tabs) "
    character(*),parameter :: F00     = "('(shr_string_leftalign_and_convert_tabs) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    ! First convert tabs to white space in the string
    str = shr_string_convert_tabs(str, rc)

    ! Now remove the leading white space
    str = adjustL(str)

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_leftalign_and_convert_tabs

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_convert_tabs -- convert all tabs to spaces
  !
  ! !DESCRIPTION:
  !    Convert all tabs to spaces in the given string
  !
  ! !REVISION HISTORY:
  !     2017-May- - M. Vertenstein
  !
  ! !INTERFACE: ------------------------------------------------------------------

  function shr_string_convert_tabs(str_input,rc) result(str_output)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(len=*)    ,intent(in)             :: str_input
    integer(SHR_KIND_IN),intent(out)  ,optional :: rc   ! return code
    character(len=len(str_input))               :: str_output
    !EOP

    !----- local ----
    integer(SHR_KIND_IN) :: inlength, i ! temporaries

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_remove_tabs) "
    character(*),parameter :: F00     = "('(shr_string_remove_tabs) ',4a)"

    ! note that tab is achar(9)
    inlength = len(str_input)
    str_output = ''
    do i = 1, inlength
       if (str_input(i:i) == achar(9)) then
          str_output(i:i) = ' '
       else
          str_output(i:i) = str_input(i:i)
       end if
    end do

    if (present(rc)) rc = 0

  end function shr_string_convert_tabs

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_alphanum -- remove non alpha numeric characters
  !
  ! !DESCRIPTION:
  !    Remove all non alpha numeric characters from string
  !     \newline
  !     call shr\_string\_alphanum(string)
  !
  ! !REVISION HISTORY:
  !     2005-Aug-01 - T. Craig - First version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_alphanum(str,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(inout)          :: str
    integer(SHR_KIND_IN),intent(out)  ,optional :: rc   ! return code

    !EOP

    !----- local ----
    integer(SHR_KIND_IN) :: n,icnt ! counters
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_alphaNum) "
    character(*),parameter :: F00     = "('(shr_string_alphaNum) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    icnt = 0
    do n=1,len_trim(str)
       if ((str(n:n) >= 'a' .and. str(n:n) <= 'z') .or.  &
            (str(n:n) >= 'A' .and. str(n:n) <= 'Z') .or.  &
            (str(n:n) >= '0' .and. str(n:n) <= '9')) then
          icnt = icnt + 1
          str(icnt:icnt) = str(n:n)
       endif
    enddo
    do n=icnt+1,len(str)
       str(n:n) = ' '
    enddo

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_alphanum

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_betweenTags -- Get the substring between the two tags.
  !
  ! !DESCRIPTION:
  !    Get the substring found between the start and end tags.
  !    \newline
  !    call shr\_string\_betweenTags(string,startTag,endTag,substring,rc)
  !
  ! !REVISION HISTORY:
  !     2005-May-11 - B. Kauffman, first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_betweenTags(string,startTag,endTag,substr,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)  :: string      ! string to search
    character(*)        ,intent(in)  :: startTag    ! start tag
    character(*)        ,intent(in)  :: endTag      ! end tag
    character(*)        ,intent(out) :: substr      ! sub-string between tags
    integer(SHR_KIND_IN),intent(out),optional :: rc ! retrun code

    !EOP

    !--- local ---
    integer(SHR_KIND_IN)   :: iStart  ! substring start index
    integer(SHR_KIND_IN)   :: iEnd    ! substring end   index
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN)       :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_betweenTags) "
    character(*),parameter :: F00     = "('(shr_string_betweenTags) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! * assumes the leading/trailing white space is not part of start & end tags
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    iStart = shr_string_endIndex(string,trim(adjustL(startTag))) ! end of start tag
    iEnd   =               index(string,trim(adjustL(endTag  ))) ! start of end tag

    rCode = 0
    substr = ""

    if (iStart < 1) then
       if (s_loglev > 0) then
          write(s_logunit,F00) "ERROR: can't find start tag in string"
          write(s_logunit,F00) "ERROR: start tag = ",trim(startTag)
          write(s_logunit,F00) "ERROR: string    = ",trim(string)
       endif
       rCode = 1
    else if (iEnd < 1) then
       if (s_loglev > 0) then
          write(s_logunit,F00) "ERROR: can't find end tag in string"
          write(s_logunit,F00) "ERROR: end   tag = ",trim(  endTag)
          write(s_logunit,F00) "ERROR: string    = ",trim(string)
       endif
       rCode = 2
    else if ( iEnd <= iStart) then
       if (s_loglev > 0) then
          write(s_logunit,F00) "ERROR: start tag not before end tag"
          write(s_logunit,F00) "ERROR: start tag = ",trim(startTag)
          write(s_logunit,F00) "ERROR: end   tag = ",trim(  endTag)
          write(s_logunit,F00) "ERROR: string    = ",trim(string)
       endif
       rCode = 3
    else if ( iStart+1 == iEnd ) then
       substr = ""
       if (s_loglev > 0) write(s_logunit,F00) "WARNING: zero-length substring found in ",trim(string)
    else
       substr = string(iStart+1:iEnd-1)
       if (len_trim(substr) == 0 .and. s_loglev > 0) &
            & write(s_logunit,F00) "WARNING: white-space substring found in ",trim(string)
    end if

    if (present(rc)) rc = rCode

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_betweenTags

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_parseCFtunit -- Parse CF time unit
  !
  ! !DESCRIPTION:
  !  Parse CF time unit into a delta string name and a base time in yyyymmdd
  !  and seconds (nearest integer actually).
  !     \newline
  !     call shr\_string\_parseCFtunit(string,substring)
  !     \newline
  !  Input string is like "days since 0001-06-15 15:20:45.5 -6:00"
  !    - recognizes "days", "hours", "minutes", "seconds"
  !    - must have at least yyyy-mm-dd, hh:mm:ss.s is optional
  !    - expects a "since" in the string
  !    - ignores time zone part
  !
  ! !REVISION HISTORY:
  !     2005-May-15 - T. Craig - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_parseCFtunit(string,unit,bdate,bsec,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: string ! string to search
    character(*)        ,intent(out)          :: unit   ! delta time unit
    integer(SHR_KIND_IN),intent(out)          :: bdate  ! base date yyyymmdd
    real(SHR_KIND_R8)   ,intent(out)          :: bsec   ! base seconds
    integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

    !EOP

    !--- local ---
    integer(SHR_KIND_IN)   :: i,i1,i2          ! generic index
    character(SHR_KIND_CL) :: tbase            ! baseline time
    character(SHR_KIND_CL) :: lstr             ! local string
    integer(SHR_KIND_IN)   :: yr,mo,da,hr,min  ! time stuff
    real(SHR_KIND_R8)      :: sec              ! time stuff
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_parseCFtunit) "
    character(*),parameter :: F00     = "('(shr_string_parseCFtunit) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! o assume length of CF-1.0 time attribute char string  < SHR_KIND_CL
    !   This is a reasonable assumption.
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    unit = 'none'
    bdate = 0
    bsec = 0.0_SHR_KIND_R8

    i = shr_string_lastIndex(string,'days ')
    if (i > 0) unit = 'days'
    i = shr_string_lastIndex(string,'hours ')
    if (i > 0) unit = 'hours'
    i = shr_string_lastIndex(string,'minutes ')
    if (i > 0) unit = 'minutes'
    i = shr_string_lastIndex(string,'seconds ')
    if (i > 0) unit = 'seconds'

    if (trim(unit) == 'none') then
       write(s_logunit,F00) ' ERROR time unit unknown'
       call shr_string_abort(subName//' time unit unknown')
    endif

    i = shr_string_lastIndex(string,' since ')
    if (i < 1) then
       write(s_logunit,F00) ' ERROR since does not appear in unit attribute for time '
       call shr_string_abort(subName//' no since in attr name')
    endif
    tbase = trim(string(i+6:))
    call shr_string_leftalign_and_convert_tabs(tbase)

    if (debug > 0 .and. s_logunit > 0) then
       write(s_logunit,*) trim(subName)//' '//'unit '//trim(unit)
       write(s_logunit,*) trim(subName)//' '//'tbase '//trim(tbase)
    endif

    yr=0; mo=0; da=0; hr=0; min=0; sec=0
    i1 = 1

    i2 = index(tbase,'-') - 1
    if(i2<0) goto 200
    lstr = tbase(i1:i2)

    read(lstr,*,ERR=200,END=200) yr
    tbase = tbase(i2+2:)
    call shr_string_leftalign_and_convert_tabs(tbase)

    i2 = index(tbase,'-') - 1
    if(i2<0) goto 200
    lstr = tbase(i1:i2)
    read(lstr,*,ERR=200,END=200) mo
    tbase = tbase(i2+2:)
    call shr_string_leftalign_and_convert_tabs(tbase)

    i2 = index(tbase,' ') - 1
    if(i2<0) i2= len_trim(tbase)
    lstr = tbase(i1:i2)
    read(lstr,*,ERR=200,END=200) da
    tbase = tbase(i2+2:)
    call shr_string_leftalign_and_convert_tabs(tbase)

    i2 = index(tbase,':') - 1
    if(i2<0) i2=len_trim(tbase)
    lstr = tbase(i1:i2)
    read(lstr,*,ERR=200,END=100) hr
    tbase = tbase(i2+2:)
    call shr_string_leftalign_and_convert_tabs(tbase)

    i2 = index(tbase,':') - 1
    if(i2<0) i2=len_trim(tbase)
    lstr = tbase(i1:i2)
    read(lstr,*,ERR=200,END=100) min
    tbase = tbase(i2+2:)
    call shr_string_leftalign_and_convert_tabs(tbase)

    i2 = index(tbase,' ') - 1
    if(i2<0) i2=len_trim(tbase)
    lstr = tbase(i1:i2)
    read(lstr,*,ERR=200,END=100) sec

100 continue
    if (debug > 0 .and. s_loglev > 0) write(s_logunit,*) trim(subName),'ymdhms:',yr,mo,da,hr,min,sec

    bdate = abs(yr)*10000 + mo*100 + da
    if (yr < 0) bdate = -bdate
    bsec = real(hr*3600 + min*60,SHR_KIND_R8) + sec

    if (present(rc)) rc = 0

    if (debug>1) call shr_timer_stop (t01)
    return

200 continue
    write(s_logunit,F00) 'ERROR 200 on char num read '
    call shr_string_abort(subName//' ERROR on char num read')
    if (debug>1) call shr_timer_stop (t01)
    return

  end subroutine shr_string_parseCFtunit

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_clean -- Clean a string, set it to "blank"
  !
  ! !DESCRIPTION:
  !     Clean a string, set it to blank
  !     \newline
  !     call shr\_string\_clean(string,rc)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_clean(string,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(inout) :: string  ! list/string
    integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_clean) "
    character(*),parameter :: F00     = "('(shr_string_clean) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    rCode = 0
    string = '       '
    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_clean

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listIsValid -- determine whether string is a valid list
  !
  ! !DESCRIPTION:
  !     Determine whether string is a valid list
  !     \newline
  !     logical_var = shr\_string\_listIsValid(list,rc)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - B. Kauffman
  !
  ! !INTERFACE: ------------------------------------------------------------------

  logical function shr_string_listIsValid(list,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list    ! list/string
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    integer  (SHR_KIND_IN) :: nChar   ! lenth of list
    integer  (SHR_KIND_IN) :: rCode   ! return code
    integer  (SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listIsValid) "
    character(*),parameter :: F00     = "('(shr_string_listIsValid) ',4a)"

    !-------------------------------------------------------------------------------
    ! check that the list conforms to the list format
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    rCode = 0
    shr_string_listIsValid = .true.

    nChar = len_trim(list)
    if (nChar < 1) then                           ! list is an empty string
       rCode = 1
    else if (    list(1:1)     == listDel  ) then ! first char is delimiter
       rCode = 2
    else if (list(nChar:nChar) == listDel  ) then ! last  char is delimiter
       rCode = 3
    else if (index(trim(list)," " )     > 0) then ! white-space in a field name
       rCode = 4
    else if (index(trim(list),listDel2) > 0) then ! found zero length field
       rCode = 5
    end if

    if (rCode /= 0) then
       shr_string_listIsValid = .false.
       if (s_loglev > 0) write(s_logunit,F00) "WARNING: invalid list = ",trim(list)
    endif

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_listIsValid

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listGetName -- Get name of k-th field in list
  !
  ! !DESCRIPTION:
  !     Get name of k-th field in list
  !     \newline
  !     call shr\_string\_listGetName(list,k,name,rc)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - B. Kauffman
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listGetName(list,k,name,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list    ! list/string
    integer(SHR_KIND_IN)         ,intent(in)  :: k       ! index of field
    character(*)                 ,intent(out) :: name    ! k-th name in list
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: i,n   ! generic indecies
    integer(SHR_KIND_IN)   :: kFlds   ! number of fields in list
    integer(SHR_KIND_IN)   :: i0,i1   ! name = list(i0:i1)
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN)   :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listGetName) "
    character(*),parameter :: F00     = "('(shr_string_listGetName) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    rCode = 0

    !--- check that this is a valid list ---
    if (.not. shr_string_listIsValid(list,rCode) ) then
       write(s_logunit,F00) "ERROR: invalid list = ",trim(list)
       call shr_string_abort(subName//" ERROR: invalid list = "//trim(list))
    end if

    !--- check that this is a valid index ---
    kFlds = shr_string_listGetNum(list)
    if (k<1 .or. kFlds<k) then
       write(s_logunit,*) subName,"ERROR: invalid index = ",k
       write(s_logunit,*) subName,"ERROR:          list = ",trim(list)
       call shr_string_abort(subName//" ERROR: invalid index")
    end if

    !--- start with whole list, then remove fields before and after desired field ---
    i0 = 1
    i1 = len_trim(list)

    !--- remove field names before desired field ---
    do n=2,k
       i = index(list(i0:i1),listDel)
       i0 = i0 + i
    end do

    !--- remove field names after desired field ---
    if ( k < kFlds ) then
       i = index(list(i0:i1),listDel)
       i1 = i0 + i - 2
    end if

    !--- copy result into output variable ---
    name = list(i0:i1)//"   "

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listGetName

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listIntersect -- Get intersection of two field lists
  !
  ! !DESCRIPTION:
  !     Get intersection of two fields lists, write into third list
  !     \newline
  !     call shr\_string\_listIntersect(list1,list2,listout)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listIntersect(list1,list2,listout,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list1   ! list/string
    character(*)                 ,intent(in)  :: list2   ! list/string
    character(*)                 ,intent(out) :: listout ! list/string
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: nf,n1,n2 ! counters
    character(SHR_KIND_CS) :: name     ! field name
    integer(SHR_KIND_IN)   :: rCode    ! return code
    integer(SHR_KIND_IN)   :: t01 = 0  ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listIntersect) "
    character(*),parameter :: F00     = "('(shr_string_listIntersect) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    rCode = 0

    nf = shr_string_listGetNum(list1)
    call shr_string_clean(listout)
    do n1 = 1,nf
       call shr_string_listGetName(list1,n1,name,rCode)
       n2 = shr_string_listGetIndexF(list2,name)
       if (n2 > 0) then
          call shr_string_listAppend(listout,name)
       endif
    enddo

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listIntersect

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listUnion -- Get union of two field lists
  !
  ! !DESCRIPTION:
  !     Get union of two fields lists, write into third list
  !     \newline
  !     call shr\_string\_listUnion(list1,list2,listout)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listUnion(list1,list2,listout,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list1   ! list/string
    character(*)                 ,intent(in)  :: list2   ! list/string
    character(*)                 ,intent(out) :: listout ! list/string
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: nf,n1,n2 ! counters
    character(SHR_KIND_CS) :: name     ! field name
    integer(SHR_KIND_IN)   :: rCode    ! return code
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listUnion) "
    character(*),parameter :: F00     = "('(shr_string_listUnion) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)
    rCode = 0

    call shr_string_clean(listout)

    nf = shr_string_listGetNum(list1)
    do n1 = 1,nf
       call shr_string_listGetName(list1,n1,name,rCode)
       n2 = shr_string_listGetIndexF(listout,name)
       if (n2 < 1) then
          call shr_string_listAppend(listout,name)
       endif
    enddo

    nf = shr_string_listGetNum(list2)
    do n1 = 1,nf
       call shr_string_listGetName(list2,n1,name,rCode)
       n2 = shr_string_listGetIndexF(listout,name)
       if (n2 < 1) then
          call shr_string_listAppend(listout,name)
       endif
    enddo

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listUnion

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listDiff -- Get set difference of two field lists
  !
  ! !DESCRIPTION:
  !     Get set difference of two fields lists, write into third list
  !     \newline
  !     call shr\_string\_listDiff(list1,list2,listout)
  !     \newline
  !     listout will contain all elements in list1 but not in list2
  !
  ! !REVISION HISTORY:
  !     2015-April-24 - W. Sacks
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listDiff(list1,list2,listout,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list1   ! list/string
    character(*)                 ,intent(in)  :: list2   ! list/string
    character(*)                 ,intent(out) :: listout ! list/string
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: num_fields, index1, index2
    character(SHR_KIND_CS) :: name     ! field name
    integer(SHR_KIND_IN)   :: rCode    ! return code
    integer(SHR_KIND_IN)   :: t01 = 0  ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listDiff) "

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    rCode = 0

    num_fields = shr_string_listGetNum(list1)
    call shr_string_clean(listout)
    do index1 = 1,num_fields
       call shr_string_listGetName(list1,index1,name,rCode)
       index2 = shr_string_listGetIndexF(list2,name)
       if (index2 <= 0) then
          call shr_string_listAppend(listout,name)
       endif
    enddo

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listDiff


  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listMerge -- Merge lists two list to third
  !
  ! !DESCRIPTION:
  !     Merge two list to third
  !     \newline
  !     call shr\_string\_listMerge(list1,list2,listout)
  !     call shr\_string\_listMerge(list1,list2,list1)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listMerge(list1,list2,listout,rc)

    implicit none
    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)  :: list1   ! list/string
    character(*)                 ,intent(in)  :: list2   ! list/string
    character(*)                 ,intent(out) :: listout ! list/string
    integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

    !EOP

    !----- local -----
    character(len=len(list1)) :: l1
    character(len=len(list2)) :: l2
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN)   :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listMerge) "
    character(*),parameter :: F00     = "('(shr_string_listMerge) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! - no input or output string should be longer than SHR_KIND_CX
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)
    rCode = 0

    call shr_string_clean(l1)
    call shr_string_clean(l2)
    call shr_string_clean(listout)
    l1 = trim(list1)
    l2 = trim(list2)
    call shr_string_leftalign_and_convert_tabs(l1,rCode)
    call shr_string_leftalign_and_convert_tabs(l2,rCode)
    if (len_trim(l1)+len_trim(l2)+1 > len(listout)) &
         call shr_string_abort(subName//'ERROR: output list string not large enough')
    if (len_trim(l1) == 0) then
       listout = trim(l2)
    else
       listout = trim(l1)//":"//trim(l2)
    endif

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listMerge

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listAppend -- Append one list to another
  !
  ! !DESCRIPTION:
  !     Append one list to another
  !     \newline
  !     call shr\_string\_listAppend(list,listadd)
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listAppend(list,listadd,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(inout) :: list   ! list/string
    character(*)                 ,intent(in)    :: listadd ! list/string
    integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

    !EOP

    !----- local -----
    character(SHR_KIND_CX) :: l1      ! local string
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN)   :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listAppend) "
    character(*),parameter :: F00     = "('(shr_string_listAppend) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! - no input or output string should be longer than SHR_KIND_CX
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)
    rCode = 0

    !--- make sure temp string is large enough ---
    if (len(l1) < len_trim(listAdd)) then
       call shr_string_abort(subName//'ERROR: temp string not large enough')
    end if

    call shr_string_clean(l1)
    l1 = trim(listadd)
    call shr_string_leftalign_and_convert_tabs(l1,rCode)
    if (len_trim(list)+len_trim(l1)+1 > len(list)) &
         call shr_string_abort(subName//'ERROR: output list string not large enough')
    if (len_trim(list) == 0) then
       list = trim(l1)
    else
       list = trim(list)//":"//trim(l1)
    endif

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listAppend

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listPrepend -- Prepend one list to another
  !
  ! !DESCRIPTION:
  !     Prepend one list to another
  !     \newline
  !     call shr\_string\_listPrepend(listadd,list)
  !     \newline
  !     results in listadd:list
  !
  ! !REVISION HISTORY:
  !     2005-May-05 - T. Craig
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listPrepend(listadd,list,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)                 ,intent(in)    :: listadd ! list/string
    character(*)                 ,intent(inout) :: list   ! list/string
    integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

    !EOP

    !----- local -----
    character(SHR_KIND_CX) :: l1      ! local string
    integer(SHR_KIND_IN)   :: rCode   ! return code
    integer(SHR_KIND_IN)   :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listPrepend) "
    character(*),parameter :: F00     = "('(shr_string_listPrepend) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! - no input or output string should be longer than SHR_KIND_CX
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)
    rCode = 0

    !--- make sure temp string is large enough ---
    if (len(l1) < len_trim(listAdd)) then
       call shr_string_abort(subName//'ERROR: temp string not large enough')
    end if

    call shr_string_clean(l1)
    l1 = trim(listadd)
    call shr_string_leftalign_and_convert_tabs(l1,rCode)
    call shr_string_leftalign_and_convert_tabs(list,rCode)
    if (len_trim(list)+len_trim(l1)+1 > len(list)) &
         call shr_string_abort(subName//'ERROR: output list string not large enough')
    if (len_trim(l1) == 0) then
       list = trim(list)
    else
       list = trim(l1)//":"//trim(list)
    endif

    if (present(rc)) rc = rCode
    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listPrepend

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listGetIndexF -- Get index of field in string
  !
  ! !DESCRIPTION:
  !     Get index of field in string
  !     \newline
  !     k = shr\_string\_listGetIndex(str,"taux")
  !
  ! !REVISION HISTORY:
  !     2005-Feb-28 - B. Kauffman and J. Schramm - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_listGetIndexF(string,fldStr)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*),intent(in) :: string   ! string
    character(*),intent(in) :: fldStr   ! name of field

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)    :: k        ! local index variable
    integer(SHR_KIND_IN)    :: rc       ! error code
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listGetIndexF) "
    character(*),parameter :: F00     = "('(shr_string_listGetIndexF) ',4a)"

    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    call shr_string_listGetIndex(string,fldStr,k,print=.false.,rc=rc)
    shr_string_listGetIndexF = k

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_listGetIndexF

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listGetIndex -- Get index of field in string
  !
  ! !DESCRIPTION:
  !     Get index of field in string
  !     \newline
  !     call shr\_string\_listGetIndex(str,"taux",k,rc)
  !
  ! !REVISION HISTORY:
  !     2005-Feb-28 - B. Kauffman and J. Schramm - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listGetIndex(string,fldStr,kFld,print,rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: string  ! string
    character(*)        ,intent(in)           :: fldStr  ! name of field
    integer(SHR_KIND_IN),intent(out)          :: kFld    ! index of field
    logical             ,intent(in) ,optional :: print   ! print switch
    integer(SHR_KIND_IN),intent(out),optional :: rc      ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: n                ! index for colon position
    integer(SHR_KIND_IN)   :: k                ! index for field name position
    integer(SHR_KIND_IN)   :: nFields          ! number of fields in a string
    integer(SHR_KIND_IN)   :: i0,i1            ! fldStr == string(i0,i1) ??
    integer(SHR_KIND_IN)   :: j0,j1            ! fldStr == string(j0,j1) ??
    logical                :: found            ! T => field found in fieldNames
    logical                :: lprint           ! local print flag
    integer(SHR_KIND_IN)   :: t01 = 0          ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listGetIndex) "
    character(*),parameter :: F00     = "('(shr_string_listGetIndex) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! - searching from both ends of the list at the same time seems to be 20% faster
    !   but I'm not sure why (B. Kauffman, Feb 2007)
    ! - I commented out sanity check to a little gain speed (B. Kauffman, Mar 2007)
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)
    if (present(rc)) rc = 0

    lprint = .false.
    if (present(print)) lprint = print

    !--- confirm proper size of input data ---
    if (len_trim(fldStr) < 1) then
       if (lprint) write(s_logunit,F00) "ERROR: input field name has 0 length"
       call shr_string_abort(subName//"invalid field name")
    end if

    !--- search for field name in string's list of fields ---
    found   = .false.
    kFld    = 0
    i0      =  1  ! ?? fldStr == string(i0:i1) ??
    i1      = -1
    j0      = -1  ! ?? fldStr == string(j0:j1) ??
    j1      =  len_trim(string)
    nFields = shr_string_listGetNum(string)
    do k = 1,nFields
       !--------------------------------------------------------
       ! search from end of list to end of list
       !--------------------------------------------------------
       !--- get end index of of field number k ---
       n = index(string(i0:len_trim(string)),listDel)
       if (n > 0) then
          i1 = i0 + n - 2       ! *not*  the  last field name in fieldNames
       else
          i1 = len_trim(string) ! this is the last field name in fieldNames
       endif
       !--- sanity check ---
       !  if ((k <nFields .and. n<1) .or. (k==nFields .and. n>0)) then
       !     call shr_string_abort(subName//"ERROR: wrong string%nf ?")
       !  end if
       !--- is it a match? ---
       if (trim(fldStr) == string(i0:i1)) then
          found = .true.
          kFld = k
          exit
       endif
       i0 = i1 + 2 ! start index for next iteration
       !--------------------------------------------------------
       ! search from end of list to start of list
       !--------------------------------------------------------
       !--- get start index of field number (nFields + 1 - k ) ---
       n = index(string(1:j1),listDel,back=.true.)
       j0 = n + 1 ! n==0 => the first field name in fieldNames
       !--- sanity check ---
       !  if ((k <nFields .and. n<1) .or. (k==nFields .and. n>0)) then
       !     call shr_string_abort(subName//"ERROR: wrong string%nf ?")
       !  end if
       !--- is it a match? ---
       if (trim(fldStr) == string(j0:j1)) then
          found = .true.
          kFld = nFields + 1 - k
          exit
       endif
       j1 = j0 - 2 ! end index for next iteration
       !--------------------------------------------------------
       ! exit if all field names have been checked
       !--------------------------------------------------------
       if (2*k >= nFields) exit
    end do

    !--- not finding a field is not a fatal error ---
    if (.not. found) then
       kFld = 0
       if (lprint .and. s_loglev > 0) write(s_logunit,F00) "FYI: field ",trim(fldStr)," not found in list ",trim(string)
       if (present(rc)) rc = 1
    end if

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listGetIndex

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listGetNum -- get number of fields in a string list
  !
  ! !DESCRIPTION:
  !  return number of fields in string list
  !
  ! !REVISION HISTORY:
  !     2005-Apr-28 - T. Craig - First version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_listGetNum(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*),intent(in) :: str   ! string to search

    !EOP

    !----- local -----
    integer(SHR_KIND_IN) :: count    ! counts occurances of char
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_listGetNum) "
    character(*),parameter :: F00     = "('(shr_string_listGetNum) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    shr_string_listGetNum = 0

    if (len_trim(str) > 0) then
       count = shr_string_countChar(str,listDel)
       shr_string_listGetNum = count + 1
    endif

    if (debug>1) call shr_timer_stop (t01)

  end function shr_string_listGetNum

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listSetDel -- Set list delimiter character
  !
  ! !DESCRIPTION:
  !     Set field delimiter character in lists
  !     \newline
  !     call shr\_string\_listSetDel(":")
  !
  ! !REVISION HISTORY:
  !     2005-Apr-30  - T. Craig - first prototype
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listSetDel(cflag)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(len=1),intent(in) :: cflag

    !EOP

    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !--- formats ---
    character(*),parameter :: subName =   "(shr_string_listSetDel) "
    character(*),parameter :: F00     = "('(shr_string_listSetDel) ',a) "

    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) 'changing listDel from '//trim(listDel)//' to '//trim(cflag)
    listDel = trim(cflag)
    listDel2 = listDel//listDel

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listSetDel

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_listGetDel -- Get list delimiter character
  !
  ! !DESCRIPTION:
  !     Get field delimiter character in lists
  !     \newline
  !     call shr\_string\_listGetDel(del)
  !
  ! !REVISION HISTORY:
  !     2005-May-15  - T. Craig - first prototype
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_listGetDel(del)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*),intent(out) :: del

    !EOP

    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !--- formats ---
    character(*),parameter :: subName =   "(shr_string_listGetDel) "
    character(*),parameter :: F00     = "('(shr_string_listGetDel) ',a) "

    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    del = trim(listDel)

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_listGetDel

  !===============================================================================
  !
  ! shr_string_listFromSuffixes
  !
  !   Returns a string of colon delimited fields given an array of suffixes and a base string
  !
  !        given suffixes = ['_s1', '_s2', '_s3'] and strBase = 'foo', returns:
  !            'foo_s1:foo_s2:foo_s3'
  !
  !===============================================================================
  function shr_string_listFromSuffixes( suffixes, strBase ) result ( retString )

    character(len=*), intent(in)  :: suffixes(:)
    character(len=*), intent(in)  :: strBase
    character(len=:), allocatable :: retString

    integer :: nfields
    integer :: i
    integer(SHR_KIND_IN) :: t01 = 0     ! timer

    character(len=*), parameter :: subName = "(shr_string_listFromSuffixes) "

    !-------------------------------------------------------------------------------

    if ( debug > 1 .and. t01 < 1 ) call shr_timer_get( t01,subName )
    if ( debug > 1 ) call shr_timer_start( t01 )

    nfields = size(suffixes)
    retString = trim(strBase) // suffixes(1)
    do i = 2, nfields
       retString = trim(retString) // ':' // trim(strBase) // suffixes(i)
    end do

    if ( debug > 1 ) call shr_timer_stop ( t01 )

  end function shr_string_listFromSuffixes

  !===============================================================================
  !
  ! shr_string_listCreateField
  !
  !   Returns a string of colon delimited fields for use in shr_strdata_create
  !   arguments, fldListFile and fldListModel.
  !   Use to create actual args for shr_strdata_create (fldListFile and
  !   flidListModel).
  !
  !   This works for numFields up to 999.  Modify the string write if you want
  !   more range.
  !
  !   retString = shr_string_listCreateField(numFields, strBase)
  !        given numFields = 5 and strBase = LAI, returns:
  !            LAI_1:LAI_2:LAI_3:LAI_4:LAI_5
  !
  !===============================================================================
  function shr_string_listCreateField( numFields, strBase ) result ( retString )

    implicit none

    integer(SHR_KIND_IN), intent(in) :: numFields   ! number of fields
    character(len=*)    , intent(in) :: strBase     ! input string base
    character(SHR_KIND_CXX)          :: retString   ! colon delimited field list

    integer                          :: idx         ! index for looping over numFields
    integer(SHR_KIND_IN)             :: t01 = 0     ! timer
    character(SHR_KIND_CX)           :: tmpString   ! temporary
    character(SHR_KIND_CX)           :: intAsChar   ! temporary
    character(1), parameter          :: colonStr = ':'
    character(1), parameter          :: underStr = '_'

    !--- formats ---
    character(*),parameter :: subName = "(shr_string_listCreateField) "
    character(*),parameter :: F00     = "('(shr_string_listCreateField) ',a) "

    !-------------------------------------------------------------------------------

    if ( debug > 1 .and. t01 < 1 ) call shr_timer_get( t01,subName )
    if ( debug > 1 ) call shr_timer_start( t01 )

    !
    ! this assert isn't that accurate since it counts all integers as being one
    ! digit, but it should catch most errors and under rather than overestimates
    !
    SHR_ASSERT_FL( ( ( ( len(strBase) + 3 ) * numFields ) <= 1024 ) , __FILE__, __LINE__)

    retString = ''
    do idx = 1,numFields

       ! reset temps per numField
       intAsChar = ''
       tmpString = ''

       ! string conversion based on 1,2,3 digits
       if ( idx < 10 ) then
          write(intAsChar, "(I1)") idx
       else if ( idx >= 10 .and. idx < 100 ) then
          write(intAsChar, "(I2)") idx
       else
          write(intAsChar, "(I3)") idx
       end if

       tmpString = trim(StrBase)//trim(underStr)//trim(intAsChar)

       if ( idx > 1 ) then
          tmpString = trim(colonStr)//trim(tmpString)
       end if

       retString = trim(retString)//trim(tmpString)

    end do

    if ( debug > 1 ) call shr_timer_stop ( t01 )

  end function shr_string_listCreateField

  !===============================================================================
  !
  ! shr_string_listAddSuffix
  !
  !   Given an existing list and a suffix, returns a new list with that suffix added to the
  !   end of every field in the list.
  !
  !   call shr_string_listAddSuffix('a:b:c', '00', new_list)
  !     gives new_list = 'a00:b00:c00'
  !
  !===============================================================================
  subroutine shr_string_listAddSuffix(list, suffix, new_list)

    implicit none

    character(len=*), intent(in)  :: list
    character(len=*), intent(in)  :: suffix
    character(len=*), intent(out) :: new_list

    integer :: num_fields
    integer :: field_num
    character(SHR_KIND_CS) :: this_field
    character(len(this_field) + len(suffix)) :: this_field_with_suffix
    character(len(new_list)) :: temp_list

    num_fields = shr_string_listGetNum(list)
    new_list = ' '

    do field_num = 1, num_fields
       call shr_string_listGetName(list, field_num, this_field)
       this_field_with_suffix = trim(this_field) // suffix
       temp_list = new_list
       call shr_string_listMerge(temp_list, this_field_with_suffix, new_list)
    end do
  end subroutine shr_string_listAddSuffix

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_setAbort -- Set local shr_string abort flag
  !
  ! !DESCRIPTION:
  !     Set local shr_string abort flag, true = abort, false = print and continue
  !     \newline
  !     call shr\_string\_setAbort(.false.)
  !
  ! !REVISION HISTORY:
  !     2005-Apr-30  - T. Craig - first prototype
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_setAbort(flag)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    logical,intent(in) :: flag

    !EOP

    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !--- formats ---
    character(*),parameter :: subName =   "(shr_string_setAbort) "
    character(*),parameter :: F00     = "('(shr_string_setAbort) ',a) "

    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    if (debug > 0 .and. s_loglev > 0) then
       if (flag) then
          write(s_logunit,F00) 'setting abort to true'
       else
          write(s_logunit,F00) 'setting abort to false'
       endif
    endif

    doabort = flag

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_setAbort

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_setDebug -- Set local shr_string debug level
  !
  ! !DESCRIPTION:
  !     Set local shr_string debug level, 0 = production
  !     \newline
  !     call shr\_string\_setDebug(2)
  !
  ! !REVISION HISTORY:
  !     2005-Apr-30  - T. Craig - first prototype
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_string_setDebug(iFlag)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in) :: iFlag ! requested debug level

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: t01 = 0 ! timer

    !--- formats ---
    character(*),parameter :: subName =   "(shr_string_setDebug) "
    character(*),parameter :: F00     = "('(shr_string_setDebug) ',a) "
    character(*),parameter :: F01     = "('(shr_string_setDebug) ',a,i3,a,i3) "

    !-------------------------------------------------------------------------------
    ! NTOE: write statement can be expensive if called many times.
    !-------------------------------------------------------------------------------

    if (iFlag>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (iFlag>1) call shr_timer_start(t01)

    !   if (s_loglev > 0) write(s_logunit,F01) 'changing debug level from ',debug,' to ',iflag
    debug = iFlag

    if (iFlag>1) call shr_timer_stop (t01)

  end subroutine shr_string_setDebug

  !===============================================================================
  !===============================================================================

  subroutine shr_string_abort(string)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*),optional,intent(in) :: string

    !EOP

    integer(SHR_KIND_IN)  :: t01 = 0 ! timer

    !--- local ---
    character(SHR_KIND_CX) :: lstring
    character(*),parameter :: subName =   "(shr_string_abort)"
    character(*),parameter :: F00     = "('(shr_string_abort) ',a)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    ! - no input or output string should be longer than SHR_KIND_CX
    !-------------------------------------------------------------------------------

    if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    if (debug>1) call shr_timer_start(t01)

    lstring = ''
    if (present(string)) lstring = string

    if (doabort) then
       call shr_sys_abort(trim(lstring))
    else
       write(s_logunit,F00) ' no abort:'//trim(lstring)
    endif

    if (debug>1) call shr_timer_stop (t01)

  end subroutine shr_string_abort

  !===============================================================================
  !===============================================================================

end module shr_string_mod
