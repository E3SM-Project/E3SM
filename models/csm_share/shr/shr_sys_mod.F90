!===============================================================================
! SVN $Id: shr_sys_mod.F90 46802 2013-05-07 18:05:27Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_130528/shr/shr_sys_mod.F90 $
!===============================================================================

#if defined CPRIBM || defined CPRPGI || defined CPRINTEL || defined __GFORTRAN__ || defined CPRCRAY
#define HAVE_GET_ENVIRONMENT
#define HAVE_SLEEP
#define HAVE_FLUSH
#endif
#if defined CPRPATHSCALE
#define HAVE_GET_ENVIRONMENT
#define HAVE_SLEEP
#endif
#if defined CPRNAG
#define HAVE_GET_ENVIRONMENT
#define HAVE_FLUSH
#define HAVE_EXECUTE
#endif


#if defined CPRPGI && defined CNL
#undef HAVE_GET_ENVIRONMENT
#endif 


MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_mpi_mod   ! wraps MPI layer
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
! PURPOSE: an architecture independent system call
! NOTE: 
! - for Catamount (Cray, pheonix at ORNL) there is no system call -- workarounds 
!   exist only for simple "rm" and "cp" commands
!-------------------------------------------------------------------------------
   rcode = 0
#ifdef HAVE_EXECUTE
   call execute_command_line(str,exitstat=rcode)   ! Intrinsic as of F2008
#else
#if (defined CRAY) || (defined UNICOSMP)

   rcode=ishell(str)

#elif (defined IRIX64 || defined NEC_SX)


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
#if (defined AIX || defined OSF1 || defined SUNOS || (defined LINUX && !defined __GFORTRAN__) || defined NEC_SX || defined CPRINTEL)
   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_chdir) '
   character(*),parameter :: F00     = "('(shr_sys_chdir) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)

   call pxfchdir(path, lenpath, rcode)

#elif (defined AIX)

   rcode = chdir(%ref(path(1:lenpath)//'\0'))

#elif (defined OSF1 || defined SUNOS || defined NEC_SX || defined Darwin || (defined LINUX && !defined CPRNAG))

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
! PURPOSE: an architecture independent system call
!-------------------------------------------------------------------------------

!$OMP master


#ifdef HAVE_GET_ENVIRONMENT
   call get_environment_variable(name=name,value=val,status=rcode)  ! Intrinsic in F2003
#else
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
   integer(SHR_KIND_IN) :: ierr
   logical              :: flag

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_abort) '
   character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   call shr_sys_flush(s_logunit)
   if (present(string)) then
      if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
   end if
   write(s_logunit,F00) 'WARNING: calling shr_mpi_abort() and stopping'
   call shr_sys_flush(s_logunit)

   ! Make sure we always print to stdout as well; this way the abort
   ! message can always be found in the CESM log.
   if ( s_logunit/= 6 ) then
      if (present(string)) then
         if (len_trim(string) > 0)  write(6,F00) 'ERROR: '//trim(string)
      end if
      write(6,F00) 'WARNING: calling shr_mpi_abort() and stopping'
      call shr_sys_flush(6)
   end if

   call shr_mpi_initialized(flag)

#if defined(CPRIBM)
   close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
#if !defined(BGL) && !defined(BGP)
   call xl__trbk()
#endif
#endif

   if (flag) then
     if (present(string).and.present(rc)) then
       call shr_mpi_abort(trim(string),rc)
     elseif (present(string)) then
       call shr_mpi_abort(trim(string))
     elseif (present(rc)) then
       call shr_mpi_abort(rcode=rc)
     else
       call shr_mpi_abort()
     endif
   endif
   call shr_sys_flush(s_logunit)
#ifndef CPRNAG
   call abort()
#endif
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
!-------------------------------------------------------------------------------
#ifdef HAVE_FLUSH
   flush(unit,iostat=ierr)   !  F2003
   if (ierr > 0) then
!      write(s_logunit,*) subname,' Flush reports error: ',ierr
   endif
#else
   call flush(unit)   !
#endif

END SUBROUTINE shr_sys_flush

!===============================================================================
!===============================================================================

END MODULE shr_sys_mod
