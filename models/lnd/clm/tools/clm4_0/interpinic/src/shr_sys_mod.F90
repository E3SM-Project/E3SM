!===============================================================================
! SVN $Id: shr_sys_mod.F90 28978 2011-06-27 20:37:05Z jedwards $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_110803/shr/shr_sys_mod.F90 $
!===============================================================================

MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none

! PUBLIC: Public interfaces

   private

   public :: shr_sys_system  ! make a system call
   public :: shr_sys_chdir   ! change current working dir
   public :: shr_sys_getenv  ! get an environment variable
   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_irtc    ! returns real-time clock tick
   public :: shr_sys_sleep   ! have program sleep for a while
   public :: shr_sys_flush   ! flush an i/o buffer

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_system(str,rcode)

   IMPLICIT none

   !----- arguments ---
   character(*)        ,intent(in)  :: str    ! system/shell command string
   integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code

   !----- functions -----
#if (defined CRAY) || (defined UNICOSMP)
   integer(SHR_KIND_IN),external    :: ishell ! function to envoke shell command
#endif
#if (defined OSF1 || defined SUNOS || (defined LINUX && !defined __GFORTRAN__ && !defined CATAMOUNT))
   integer(SHR_KIND_IN),external    :: system ! function to envoke shell command
#endif

   !----- local -----
#if (defined CATAMOUNT)
   character(2*SHR_KIND_CL) :: file1            ! one or two filenames
   character(  SHR_KIND_CL) :: file2            ! 2nd file name
   integer(SHR_KIND_IN)     :: iloc             ! index/location within a string
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_system) '
   character(*),parameter :: F00     = "('(shr_sys_system) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
! NOTE: 
! - for Catamount (Cray, pheonix at ORNL) there is no system call -- workarounds 
!   exist only for simple "rm" and "cp" commands
!-------------------------------------------------------------------------------


#if (defined CRAY) || (defined UNICOSMP)

   rcode=ishell(str)

#elif (defined IRIX64 || defined NEC_SX)

   rcode = 0
   call system(str)

#elif (defined AIX)

   call system(str,rcode)

#elif (defined OSF1 || defined SUNOS || defined __GFORTRAN__ || (defined LINUX && !defined CATAMOUNT))

   rcode = system(str)

#elif (defined CATAMOUNT)
   if (str(1:3) == 'rm ') then
      call unlink(str(4:))
      if (s_loglev > 0) write(s_logunit,F00) 'CATAMOUNT unlink ',trim(str(4:))
      rcode = 0
   elseif (str(1:3) == 'mv ') then
      file1 = str(4:)
      iloc = index(file1,' ') + 3
      if (iloc < 6) then
         if (s_loglev > 0) write(s_logunit,*) 'CATAMOUNT mv error ',trim(str),iloc
         rcode = -1
      else
         file1 = str(4:iloc)
         file2 = str(iloc+1:)
         call rename(trim(file1),trim(file2))
         if (s_loglev > 0) write(s_logunit,F00) 'CATAMOUNT rename ',trim(file1)," ",trim(file2)
         rcode = 0
      endif
   else
      rcode = -1
   endif

#else

   write(s_logunit,F00) 'ERROR: no implementation of system call for this architecture'
   call shr_sys_abort(subName//'no implementation of system call for this architecture')

#endif

END SUBROUTINE shr_sys_system

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_chdir(path, rcode)

   IMPLICIT none

   !----- arguments -----
   character(*)        ,intent(in)  :: path    ! chdir to this dir
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenpath ! length of path
#if (defined AIX || defined OSF1 || defined SUNOS || (defined LINUX && !defined __GFORTRAN__) || defined NEC_SX)
   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_chdir) '
   character(*),parameter :: F00     = "('(shr_sys_chdir) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)

   call pxfchdir(path, lenpath, rcode)

#elif (defined AIX)

   rcode = chdir(%ref(path(1:lenpath)//'\0'))

#elif (defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)

   rcode=chdir(path(1:lenpath))

#else

   write(s_logunit,F00) 'ERROR: no implementation of chdir for this architecture'
   call shr_sys_abort(subname//'no implementation of chdir for this machine')

#endif

END SUBROUTINE shr_sys_chdir

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(*)        ,intent(in)  :: name    ! env var name
   character(*)        ,intent(out) :: val     ! env var value
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenname ! length of env var name
   integer(SHR_KIND_IN)             :: lenval  ! length of env var value
   character(SHR_KIND_CL)           :: tmpval  ! temporary env var value

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_getenv) '
   character(*),parameter :: F00     = "('(shr_sys_getenv) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenname=len_trim(name)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)

   call pxfgetenv(name, lenname, val, lenval, rcode)

#elif (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)

   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0         ) rcode = 1
   if (len_trim(val) >  SHR_KIND_CL) rcode = 2

#else

   write(s_logunit,F00) 'ERROR: no implementation of getenv for this architecture'
   call shr_sys_abort(subname//'no implementation of getenv for this machine')

#endif

END SUBROUTINE shr_sys_getenv

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_abort(string,rc)

   IMPLICIT none

   character(*)        ,optional :: string  ! error message string
   integer(SHR_KIND_IN),optional :: rc      ! error code

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_abort) '
   character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   call shr_sys_flush(s_logunit)
   if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
   write(s_logunit,F00) 'WARNING: stopping'
   call shr_sys_flush(s_logunit)
   call abort()
   stop

END SUBROUTINE shr_sys_abort

!===============================================================================
!===============================================================================

integer(SHR_KIND_I8) FUNCTION shr_sys_irtc( rate )

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_I8), optional :: rate

   !----- local -----
   integer(SHR_KIND_IN)      :: count
   integer(SHR_KIND_IN)      :: count_rate
   integer(SHR_KIND_IN)      :: count_max
   integer(SHR_KIND_IN),save :: last_count   = -1
   integer(SHR_KIND_I8),save :: count_offset =  0

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_irtc) '
   character(*),parameter :: F00     = "('(shr_sys_irtc) ',4a)"

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!-------------------------------------------------------------------------------

   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   shr_sys_irtc = count

   !--- adjust for clock wrap-around ---
   if ( last_count /= -1 ) then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   shr_sys_irtc = shr_sys_irtc + count_offset
   last_count   = count

END FUNCTION shr_sys_irtc

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_sleep(sec)

   IMPLICIT none

   !----- arguments -----
   real   (SHR_KIND_R8),intent(in) :: sec  ! number of seconds to sleep

   !----- local -----
   integer(SHR_KIND_IN) :: isec   ! integer number of seconds
   integer(SHR_KIND_IN) :: rcode  ! return code
   character(90)        :: str    ! system call string

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_sleep) '
   character(*),parameter :: F00     = "('(shr_sys_sleep) ',4a)"
   character(*),parameter :: F10     = "('sleep ',i8 )"

!-------------------------------------------------------------------------------
! PURPOSE: Sleep for approximately sec seconds
!-------------------------------------------------------------------------------

   isec = nint(sec)

   if (isec < 0) then
      if (s_loglev > 0) write(s_logunit,F00) 'ERROR: seconds must be > 0, sec=',sec
   else if (isec == 0) then
      ! Don't consider this an error and don't call system sleep
   else
#if defined(CATAMOUNT)
      call sleep(isec)
#else
      write(str,FMT=F10) isec
      call shr_sys_system( str, rcode )
#endif
   endif

END SUBROUTINE shr_sys_sleep

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_flush) '
   character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

#if (defined IRIX64 || defined CRAY || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX || defined UNICOSMP)

   call flush(unit)

#elif (defined AIX)

   call flush_(unit)

#else

   if (s_loglev > 0) write(s_logunit,F00) 'WARNING: no implementation of flush for this architecture'

#endif

END SUBROUTINE shr_sys_flush

!===============================================================================
!===============================================================================

END MODULE shr_sys_mod
