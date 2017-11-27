MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
!  use shr_mpi_mod   ! wraps MPI layer

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
   public :: shr_sys_ioUnit  ! returns an unused fortran i/o unit number

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_system(str,rcode)

   IMPLICIT none

   !----- arguments ---
   character(len=*)    ,intent(in)  :: str    ! system/shell command string
   integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code

#if (defined CRAY) || (defined UNICOSMP)
   !----- functions -----
   integer(SHR_KIND_IN),external    :: ishell ! function to envoke shell command
#endif
#if (defined OSF1 || defined SUNOS || (defined LINUX && !defined G95))
   !----- functions -----
   integer(SHR_KIND_IN),external    :: system ! function to envoke shell command
#endif

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

#if (defined CRAY) || (defined UNICOSMP)
   rcode=ishell(str)
#endif

#if (defined IRIX64 || defined NEC_SX)
   rcode = 0
   call system(str)
#endif

#if (defined AIX)
   call system(str,rcode)
#endif

#if (defined OSF1 || defined SUNOS || defined LINUX)
   rcode = system(str)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(*,*) '(shr_sys_system) ERROR: no implementation for this architecture'
   call shr_sys_abort('no system routine on this machine')
#endif

END SUBROUTINE shr_sys_system

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_chdir(path, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: path    ! chdir to this dir
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenpath ! length of path
#if (defined AIX || defined OSF1 || defined SUNOS || (defined LINUX && !defined G95) || defined NEC_SX)
   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)
   call pxfchdir(path, lenpath, rcode)
#endif

#if (defined AIX)
   rcode=chdir(%ref(path(1:lenpath)//'\0'))
#endif

#if (defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
   rcode=chdir(path(1:lenpath))
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(*,*) '(shr_sys_chdir) ERROR: no implementation for this architecture'
   call shr_sys_abort('no implementation of chdir for this machine')
#endif

END SUBROUTINE shr_sys_chdir

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: name    ! env var name
   character(len=*)    ,intent(out) :: val     ! env var value
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenname ! length of env var name
   integer(SHR_KIND_IN)             :: lenval  ! length of env var value
   character(len=512)               :: tmpval  ! temporary env var value

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenname=len_trim(name)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)
   call pxfgetenv(name, lenname, val, lenval, rcode)
#endif

#if (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0 ) rcode = 1
   if (len_trim(val) > 512 ) rcode = 2
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(*,*) '(shr_sys_getenv) ERROR: no implementation for this architecture'
   call shr_sys_abort('no implementation of getenv for this machine')
#endif

END SUBROUTINE shr_sys_getenv

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_abort(string,rcode)

   IMPLICIT none

   character*(*),optional :: string    ! error message string
   integer(SHR_KIND_IN),optional :: rcode   ! error code

   !----- local -----
   integer(SHR_KIND_IN) :: ierr
   logical              :: flag

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_sys_abort) ',a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   call shr_sys_flush(6)
   if (len_trim(string) > 0) write(6,F00) 'ERROR: '//trim(string)
!  write(6,F00) 'WARNING: calling shr_mpi_abort() and stopping'
   write(6,F00) 'WARNING:                 abort() and stopping'
   call shr_sys_flush(6)
!  call shr_mpi_initialized(flag)
!  if (flag) then
!    if (present(string).and.present(rcode)) then
!      call shr_mpi_abort(trim(string),rcode)
!    elseif (present(string)) then
!      call shr_mpi_abort(trim(string))
!    elseif (present(rcode)) then
!      call shr_mpi_abort(rcode=rcode)
!    else
!      call shr_mpi_abort()
!    endif
!  endif
   call shr_sys_flush(6)
   call abort()
   stop

END SUBROUTINE shr_sys_abort

!===============================================================================
!===============================================================================

integer(SHR_KIND_I8) FUNCTION shr_sys_irtc( rate )

   IMPLICIT none
   !----- optional output argument -----
   integer(SHR_KIND_I8), optional :: rate

   !----- local -----
   integer(SHR_KIND_IN)          :: count
   integer(SHR_KIND_IN)          :: count_rate
   integer(SHR_KIND_IN)          :: count_max

   integer(SHR_KIND_IN),save :: last_count = -1
   integer(SHR_KIND_I8),save :: count_offset = 0

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!-------------------------------------------------------------------------------
   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   shr_sys_irtc = count
!
! System clock is a 24-hour clock, so must check each time pass over midnight
!
   if ( last_count /= -1 )then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   shr_sys_irtc = shr_sys_irtc + count_offset
   last_count = count

END FUNCTION shr_sys_irtc

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_sleep(sec)

   IMPLICIT none

   !----- input -----
   real   (shr_kind_r8),intent(in) :: sec  ! number of seconds to sleep

   !----- local -----
   integer(shr_kind_in) :: isec   ! integer number of seconds
   integer(shr_kind_in) :: rcode  ! return code
   character(len=90) :: sleep_var

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('sleep ',i8 )"

   save

!-------------------------------------------------------------------------------
! PURPOSE: Sleep for approximately sec seconds
!-------------------------------------------------------------------------------

   isec = nint(sec)

   if (isec < 0) then
      write(6,*) 'ERROR: seconds must be > 0, sec=',sec
   else if (isec == 0) then
!     Don't consider this an error and don't call system sleep
   else
      write(sleep_var,FMT=F00) isec
      call shr_sys_system( sleep_var, rcode )
   endif

END SUBROUTINE shr_sys_sleep

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

#if (defined IRIX64 || defined CRAY || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX || defined UNICOSMP)
   call flush(unit)
#endif
#if (defined AIX)
   call flush_(unit)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(*,*) '(shr_sys_flush) WARNING: no implementation for this architecture'
#endif

END SUBROUTINE shr_sys_flush

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_sys_iounit -- returns an unused fortran i/o unit number
!
! !DESCRIPTION:
!    Returns an unused fortran i/o unit number.
!
! !REVISION HISTORY:
!     2005-Nov-21 -- B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

integer(SHR_KIND_IN) function shr_sys_ioUnit()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !-- no input, output via function value

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: nUnit       ! a file unit number
   logical                :: open        ! true if file unit is open

   !--- formats ---
   character(*),parameter :: subName =  "('shr_sys_ioUnit') "
   character(*),parameter :: F00     = "('(shr_sys_ioUnit) ',16a) "

!-------------------------------------------------------------------------------
! find an unused unit number
!-------------------------------------------------------------------------------

   do nUnit=10,99
      inquire(nUnit,opened=open)
      if (.not. open) then
         shr_sys_ioUnit = nUnit
         exit
      end if
   end do

   if (open) then
      write(6,F00) "ERROR: couldn't find an unused fortran unit number"
      call shr_sys_abort(subName//": all units open?")
   end if

end function shr_sys_ioUnit

!===============================================================================

END MODULE shr_sys_mod
