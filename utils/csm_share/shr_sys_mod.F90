!===============================================================================
! SVN $Id: shr_sys_mod.F90 3880 2007-04-10 20:59:37Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_sys_mod.F90 $
!===============================================================================

MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_mpi_mod   ! wraps MPI layer

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
   character(*)        ,intent(in)  :: str    ! system/shell command string
   integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code

   !----- functions -----
#if (defined CRAY) || (defined UNICOSMP)
   integer(SHR_KIND_IN),external    :: ishell ! function to envoke shell command
#endif
#if (defined OSF1 || defined SUNOS || (defined LINUX && !defined G95) || (defined LINUX && !defined CATAMOUNT))
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
#endif

#if (defined IRIX64 || defined NEC_SX)
   rcode = 0
   call system(str)
#endif

#if (defined AIX)
   call system(str,rcode)
#endif

#if (defined OSF1 || defined SUNOS || defined LINUX && !defined CATAMOUNT)
   rcode = system(str)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(6,F00) 'ERROR: no implementation for this architecture'
   call shr_sys_abort(subName//'no implementation for this architecture')
#endif

#if (defined CATAMOUNT)
   if (str(1:3) == 'rm ') then
      call unlink(str(4:))
      write(6,F00) 'CATAMOUNT unlink ',trim(str(4:))
      rcode = 0
   elseif (str(1:3) == 'mv ') then
      file1 = str(4:)
      iloc = index(file1,' ') + 3
      if (iloc < 6) then
         write(6,*) 'CATAMOUNT mv error ',trim(str),iloc
         rcode = -1
      else
         file1 = str(4:iloc)
         file2 = str(iloc+1:)
         call rename(trim(file1),trim(file2))
         write(6,F00) 'CATAMOUNT rename ',trim(file1)," ",trim(file2)
         rcode = 0
      endif
   else
      rcode = -1
   endif
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
#if (defined AIX || defined OSF1 || defined SUNOS || (defined LINUX && !defined G95) || defined NEC_SX)
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
#endif

#if (defined AIX)
   rcode=chdir(%ref(path(1:lenpath)//'\0'))
#endif

#if (defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
   rcode=chdir(path(1:lenpath))
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(6,F00) 'ERROR: no implementation for this architecture'
   call shr_sys_abort('no implementation of chdir for this machine')
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
#endif

#if (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0         ) rcode = 1
   if (len_trim(val) >  SHR_KIND_CL) rcode = 2
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(6,F00) 'ERROR: no implementation for this architecture'
   call shr_sys_abort('no implementation of getenv for this machine')
#endif

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

   call shr_sys_flush(6)
   if (len_trim(string) > 0) write(6,F00) 'ERROR: '//trim(string)
   write(6,F00) 'WARNING: calling shr_mpi_abort() and stopping'
   call shr_sys_flush(6)
   call shr_mpi_initialized(flag)
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
   call shr_sys_flush(6)
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
      write(6,F00) 'ERROR: seconds must be > 0, sec=',sec
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
#endif
#if (defined AIX)
   call flush_(unit)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
   write(6,F00) 'WARNING: no implementation for this architecture'
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
!===============================================================================

END MODULE shr_sys_mod
