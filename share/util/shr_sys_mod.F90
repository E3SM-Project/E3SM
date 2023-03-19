!===============================================================================
! SVN $Id: shr_sys_mod.F90 66411 2014-12-19 22:40:08Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_150116/shr/shr_sys_mod.F90 $
!===============================================================================

! Currently supported by all compilers
#define HAVE_GET_ENVIRONMENT
#define HAVE_SLEEP

! Except this combination?
#if defined CPRPGI && defined CNL
#undef HAVE_GET_ENVIRONMENT
#endif

#if defined CPRNAG
#define HAVE_EXECUTE
#endif

MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use shr_abort_mod, only: shr_sys_abort => shr_abort_abort
   use shr_abort_mod, only: shr_sys_backtrace => shr_abort_backtrace

#ifdef CPRNAG
   ! NAG does not provide these as intrinsics, but it does provide modules
   ! that implement commonly used POSIX routines.
   use f90_unix_dir, only: chdir
   use f90_unix_proc, only: abort, sleep
#endif

   implicit none

! PUBLIC: Public interfaces

   private

   public :: shr_sys_system  ! make a system call
   public :: shr_sys_chdir   ! change current working dir
   public :: shr_sys_getenv  ! get an environment variable
   public :: shr_sys_irtc    ! returns real-time clock tick
   public :: shr_sys_sleep   ! have program sleep for a while
   public :: shr_sys_flush   ! flush an i/o buffer

   ! Imported from shr_abort_mod and republished with renames. Other code that wishes to
   ! use these routines should use these shr_sys names rather than directly using the
   ! routines from shr_abort_abort. (This is for consistency with older code, from when
   ! these routines were defined in shr_sys_mod.)
   public :: shr_sys_abort     ! abort a program
   public :: shr_sys_backtrace ! print a backtrace, if possible

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
#if (defined LINUX && !defined CPRGNU)
   integer(SHR_KIND_IN),external    :: system ! function to envoke shell command
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_system) '
   character(*),parameter :: F00     = "('(shr_sys_system) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------
   rcode = 0
#ifdef HAVE_EXECUTE
   call execute_command_line(str,exitstat=rcode)   ! Intrinsic as of F2008
#else
#if (defined AIX)

   call system(str,rcode)

#elif (defined CPRGNU || defined LINUX)

   rcode = system(str)

#else

   write(s_logunit,F00) 'ERROR: no implementation of system call for this architecture'
   call shr_sys_abort(subName//'no implementation of system call for this architecture')
#endif
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
#if (defined AIX || (defined LINUX && !defined CPRGNU && !defined CPRNAG) || defined CPRINTEL)
   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_chdir) '
   character(*),parameter :: F00     = "('(shr_sys_chdir) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)

#if (defined AIX)

   rcode = chdir(%ref(path(1:lenpath)//'\0'))

#elif (defined Darwin || (defined LINUX && !defined CPRNAG))

   rcode=chdir(path(1:lenpath))

#elif (defined CPRNAG)

   call chdir(path(1:lenpath), errno=rcode)

#else
   rcode=-999
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
#ifndef HAVE_GET_ENVIRONMENT
   integer(SHR_KIND_IN)             :: lenname ! length of env var name
   integer(SHR_KIND_IN)             :: lenval  ! length of env var value
   character(SHR_KIND_CL)           :: tmpval  ! temporary env var value
#endif
   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_getenv) '
   character(*),parameter :: F00     = "('(shr_sys_getenv) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------

!$OMP master


#ifdef HAVE_GET_ENVIRONMENT
   call get_environment_variable(name=name,value=val,status=rcode)  ! Intrinsic in F2003
#else
   lenname=len_trim(name)
#if (defined AIX || defined LINUX)

   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0         ) rcode = 1
   if (len_trim(val) >  SHR_KIND_CL) rcode = 2

#else

   write(s_logunit,F00) 'ERROR: no implementation of getenv for this architecture'
   call shr_sys_abort(subname//'no implementation of getenv for this machine')

#endif
#endif
!$OMP end master

END SUBROUTINE shr_sys_getenv

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
!$OMP THREADPRIVATE (last_count, count_offset)

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_irtc) '
   character(*),parameter :: F00     = "('(shr_sys_irtc) ',4a)"

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!
! This function is not intended to measure elapsed time between
! multi-threaded regions with different numbers of threads. However,
! use of the threadprivate declaration does guarantee accurate
! measurement per thread within a single multi-threaded region as
! long as the number of threads is not changed dynamically during
! execution within the multi-threaded region.
!
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
#ifndef HAVE_SLEEP
   integer(SHR_KIND_IN) :: rcode  ! return code
   character(90)        :: str    ! system call string
#endif
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
#ifdef HAVE_SLEEP
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

   !----- local -----
   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_flush) '
   character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!
! This is probably no longer needed; the "flush" statement is supported by
! all compilers that CESM supports for years now.
!
!-------------------------------------------------------------------------------
!$OMP SINGLE
   flush(unit)
!$OMP END SINGLE
!
! The following code was originally present, but there's an obvious issue.
! Since shr_sys_flush is usually used to flush output to a log, when it
! returns an error, does it do any good to print that error to the log?
!
!   if (ierr > 0) then
!      write(s_logunit,*) subname,' Flush reports error: ',ierr
!   endif
!

END SUBROUTINE shr_sys_flush

!===============================================================================
!===============================================================================

END MODULE shr_sys_mod
