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

   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

   use shr_kind_mod  ! defines real & integer kinds
   use shr_mpi_mod   ! wraps MPI layer
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

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
   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_irtc    ! returns real-time clock tick
   public :: shr_sys_sleep   ! have program sleep for a while
   public :: shr_sys_flush   ! flush an i/o buffer
   public :: shr_sys_backtrace   ! print a backtrace, if possible

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

SUBROUTINE shr_sys_abort(string,rc)

   IMPLICIT none

   character(*)        ,optional :: string  ! error message string
   integer(SHR_KIND_IN),optional :: rc      ! error code

   !----- local -----
   logical              :: flag

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_abort) '
   character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

   ! Local version of the string.
   ! (Gets a default value if string is not present.)
   character(len=shr_kind_cx) :: local_string

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   if (present(string)) then
      local_string = trim(string)
   else
      local_string = "Unknown error submitted to shr_sys_abort."
   end if

   call print_error_to_logs("ERROR", local_string)

   call shr_sys_backtrace()

   call shr_mpi_initialized(flag)

   if (flag) then
      if (present(rc)) then
         call shr_mpi_abort(trim(local_string),rc)
      else
         call shr_mpi_abort(trim(local_string))
      endif
   endif

  ! A compiler's abort method may print a backtrace or do other nice
  ! things, but in fact we can rarely leverage this, because MPI_Abort
  ! usually sends SIGTERM to the process, and we don't catch that signal.
   call abort()

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
   integer(SHR_KIND_IN) :: ierr  ! error code

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
   flush(unit)
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

subroutine shr_sys_backtrace()

  ! This routine uses compiler-specific facilities to print a backtrace to
  ! error_unit (standard error, usually unit 0).

#if defined(CPRIBM)

  ! This theoretically should be in xlfutility, but using it from that
  ! module doesn't seem to always work.
  interface
     subroutine xl_trbk()
     end subroutine xl_trbk
  end interface

  call xl__trbk()

#elif defined(CPRGNU) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8 ))

  ! gfortran 4.8 and later implement this intrinsic. We explicitly call it
  ! out as such to make sure that it really is available, just in case the
  ! CPP logic above screws up.
  intrinsic :: backtrace

  call backtrace()

#elif defined(CPRINTEL)

  ! tracebackqq uses optional arguments, so *must* have an explicit
  ! interface.
  use ifcore, only: tracebackqq

  ! An exit code of -1 is a special value that prevents this subroutine
  ! from aborting the run.
  call tracebackqq(user_exit_code=-1)

#else

  ! Currently we have no means to request a backtrace from the NAG runtime,
  ! even though it is capable of emitting backtraces itself, if you use the
  ! "-gline" option.

  ! Similarly, PGI has a -traceback option, but no user interface for
  ! requesting a backtrace to be printed.

#endif

  flush(error_unit)

end subroutine shr_sys_backtrace

!===============================================================================
!===============================================================================

!
! This routine prints error messages to s_logunit (which is standard output
! for most tasks in CESM) and also to standard error if s_logunit is a
! file.
!
! It also flushes these output units.
!
subroutine print_error_to_logs(error_type, message)
  character(len=*), intent(in) :: error_type, message

  integer, allocatable :: log_units(:)

  integer :: i

  if (s_logunit == output_unit .or. s_logunit == error_unit) then
     ! If the log unit number is standard output or standard error, just
     ! print to that.
     allocate(log_units(1), source=[s_logunit])
  else
     ! Otherwise print the same message to both the log unit and standard
     ! error.
     allocate(log_units(2), source=[error_unit, s_logunit])
  end if

  do i = 1, size(log_units)
     write(log_units(i),*) trim(error_type), ": ", trim(message)
     flush(log_units(i))
  end do

end subroutine print_error_to_logs

!===============================================================================
!===============================================================================

END MODULE shr_sys_mod
