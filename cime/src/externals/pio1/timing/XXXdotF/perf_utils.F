      module perf_utils

C-----------------------------------------------------------------------
C
C Purpose: This module supplies the csm_share and CAM utilities
C          needed by perf_mod.F90 (when the csm_share and CAM utilities
C          are not available).
C
C Author:  P. Worley, October 2007
C
C $Id$
C
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C- module boilerplate --------------------------------------------------
C-----------------------------------------------------------------------
      implicit none
C     Make the default access private
      private
      save

C-----------------------------------------------------------------------
C Public interfaces ----------------------------------------------------
C-----------------------------------------------------------------------
      public perfutils_setunit
      public shr_sys_abort
      public shr_mpi_barrier
      public shr_file_getUnit
      public shr_file_freeUnit
      public find_group_name
      public to_lower
      public shr_mpi_bcast

      interface shr_mpi_bcast ; module procedure
     &  shr_mpi_bcastl0,
     &  shr_mpi_bcasti0
      end interface

C-----------------------------------------------------------------------
C Private interfaces ---------------------------------------------------
C-----------------------------------------------------------------------
      private shr_sys_flush
      private shr_mpi_chkerr
      private shr_mpi_abort

C-----------------------------------------------------------------------
C- include statements --------------------------------------------------
C-----------------------------------------------------------------------
#include <mpif.h>
#include "gptl.inc"

C-----------------------------------------------------------------------
C Public data ---------------------------------------------------------
C-----------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! precision/kind constants (from csm_share/shr/shr_kind_mod.F90)
   !----------------------------------------------------------------------------
C  8 byte real
      integer,parameter,public :: SHR_KIND_R8 = selected_real_kind(12)
C  native integer
      integer,parameter,public :: SHR_KIND_IN = kind(1)
C  long char
      integer,parameter,public :: SHR_KIND_CL = 256
C extra-long char
      integer,parameter,public :: SHR_KIND_CX = 512

C-----------------------------------------------------------------------
C Private data ---------------------------------------------------------
C-----------------------------------------------------------------------

C  default
      integer, parameter :: def_pu_logunit = 6
C  unit number for log output
      integer, private   :: pu_logunit = def_pu_logunit

C=======================================================================
      contains
C=======================================================================

C
C========================================================================
C
      subroutine perfutils_setunit(LogUnit)
C-----------------------------------------------------------------------
C Purpose:  Set log unit number.
C Author:   P. Worley
C-----------------------------------------------------------------------
C---------------------------Input arguments-----------------------------
C
C  Unit number for log output
      integer(SHR_KIND_IN), intent(IN) :: LogUnit
C-----------------------------------------------------------------------
      pu_logunit = LogUnit
C
      return
C
      end subroutine perfutils_setunit

C============== Routines from csm_share/shr/shr_sys_mod.F90 ============
C=======================================================================

      SUBROUTINE shr_sys_abort(string)

      IMPLICIT none

C  error message string
      character(*)        ,optional :: string

C  !----- local -----
      integer(SHR_KIND_IN) :: ierr
      logical              :: flag

C  !----- formats -----
      character(*),parameter :: subName =   '(shr_sys_abort) '
      character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

C-------------------------------------------------------------------------------
C PURPOSE: consistent stopping mechanism
C (dumbed down from original shr_sys_mod.F90 version for use in perf_mod)
C-------------------------------------------------------------------------------

      call shr_sys_flush(pu_logunit)

      if ( present(string) ) then
         if (len_trim(string) > 0) then
            write(pu_logunit,*) trim(subName),' ERROR: ',trim(string)
         else
            write(pu_logunit,*) trim(subName),' ERROR '
         endif
      else
         write(pu_logunit,*) trim(subName),' ERROR '
      endif

      write(pu_logunit,F00) 'WARNING: calling mpi_abort() and stopping'
      call shr_sys_flush(pu_logunit)
      call mpi_abort(MPI_COMM_WORLD,0,ierr)
      call shr_sys_flush(pu_logunit)
      call abort()

      stop

      END SUBROUTINE shr_sys_abort

C===============================================================================
C===============================================================================

      SUBROUTINE shr_sys_flush(unit)

      IMPLICIT none

C  !----- arguments -----
C  flush output buffer for this unit
      integer(SHR_KIND_IN) :: unit

C  !----- formats -----
      character(*),parameter :: subName =   '(shr_sys_flush) '
      character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"

C-------------------------------------------------------------------------------
C PURPOSE: an architecture independant system call
C-------------------------------------------------------------------------------

#if (defined IRIX64 || defined CRAY || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX || defined UNICOSMP)
      call flush(unit)
#endif
#if (defined AIX)
      call flush_(unit)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX && !defined NEC_SX && !defined UNICOSMP)
Cpw if (s_loglev > 0) write(pu_logunit,F00) 'WARNING: no implementation of flush for this architecture'
#endif

      END SUBROUTINE shr_sys_flush

C===============================================================================

C================== Routines from csm_share/shr/shr_mpi_mod.F90 ===============
C===============================================================================

      SUBROUTINE shr_mpi_chkerr(rcode,string)

      IMPLICIT none

C  !----- arguments ---
C  input MPI error code
      integer(SHR_KIND_IN), intent(in) :: rcode
C  message
      character(*),         intent(in) :: string

C  !----- local ---
      character(*),parameter           :: subName = '(shr_mpi_chkerr) '
      character(MPI_MAX_ERROR_STRING)  :: lstring
      integer(SHR_KIND_IN)             :: len
      integer(SHR_KIND_IN)             :: ierr

C-------------------------------------------------------------------------------
C PURPOSE: layer on MPI error checking
C-------------------------------------------------------------------------------

      if (rcode /= MPI_SUCCESS) then
        call MPI_ERROR_STRING(rcode,lstring,len,ierr)
        write(pu_logunit,*) trim(subName),":",lstring(1:len)
        call shr_mpi_abort(string,rcode)
      endif

      END SUBROUTINE shr_mpi_chkerr

C===============================================================================
C===============================================================================

      SUBROUTINE shr_mpi_abort(string,rcode)

      IMPLICIT none

C  !----- arguments ---
C  message
      character(*),optional,intent(in)   :: string
C  optional code
      integer,optional,intent(in)        :: rcode

C  !----- local ---
      character(*),parameter             :: subName = '(shr_mpi_abort) '
      integer(SHR_KIND_IN)               :: ierr

C-------------------------------------------------------------------------------
C PURPOSE: MPI abort
C-------------------------------------------------------------------------------

      if ( present(string) .and. present(rcode) ) then
         write(pu_logunit,*) trim(subName),":",trim(string),rcode
      endif
      call MPI_ABORT(MPI_COMM_WORLD,rcode,ierr)

      END SUBROUTINE shr_mpi_abort

C===============================================================================
C===============================================================================

      SUBROUTINE shr_mpi_barrier(comm,string)

      IMPLICIT none

C  !----- arguments ---
      integer,intent(in)               :: comm
C  message
      character(*),optional,intent(in) :: string

C  !----- local ---
      character(*),parameter           :: subName = '(shr_mpi_barrier) '
      integer(SHR_KIND_IN)             :: ierr

C-------------------------------------------------------------------------------
C PURPOSE: MPI barrier
C-------------------------------------------------------------------------------

      call MPI_BARRIER(comm,ierr)
      if (present(string)) then
        call shr_mpi_chkerr(ierr,subName//trim(string))
      else
        call shr_mpi_chkerr(ierr,subName)
      endif

      END SUBROUTINE shr_mpi_barrier

C===============================================================================
C===============================================================================

      SUBROUTINE shr_mpi_bcasti0(vec,comm,string)

      IMPLICIT none

C  !----- arguments ---
C  vector of 1
      integer(SHR_KIND_IN), intent(inout):: vec
C  mpi communicator
      integer(SHR_KIND_IN), intent(in)   :: comm
C  message
      character(*),optional,intent(in)   :: string

C  !----- local ---
      character(*),parameter           :: subName = '(shr_mpi_bcasti0) '
      integer(SHR_KIND_IN)             :: ierr
      integer(SHR_KIND_IN)             :: lsize

C-------------------------------------------------------------------------------
C PURPOSE: Broadcast an integer
C-------------------------------------------------------------------------------

      lsize = 1

      call MPI_BCAST(vec,lsize,MPI_INTEGER,0,comm,ierr)
      if (present(string)) then
        call shr_mpi_chkerr(ierr,subName//trim(string))
      else
        call shr_mpi_chkerr(ierr,subName)
      endif

      END SUBROUTINE shr_mpi_bcasti0

C===============================================================================
C===============================================================================

      SUBROUTINE shr_mpi_bcastl0(vec,comm,string)

      IMPLICIT none

C  !----- arguments ---
C  vector of 1
      logical, intent(inout):: vec
C  mpi communicator
      integer(SHR_KIND_IN), intent(in)   :: comm
C  message
      character(*),optional,intent(in)   :: string

C  !----- local ---
      character(*),parameter           :: subName = '(shr_mpi_bcastl0) '
      integer(SHR_KIND_IN)             :: ierr
      integer(SHR_KIND_IN)             :: lsize

C-------------------------------------------------------------------------------
C PURPOSE: Broadcast a logical
C-------------------------------------------------------------------------------

      lsize = 1

      call MPI_BCAST(vec,lsize,MPI_LOGICAL,0,comm,ierr)
      if (present(string)) then
        call shr_mpi_chkerr(ierr,subName//trim(string))
      else
        call shr_mpi_chkerr(ierr,subName)
      endif

      END SUBROUTINE shr_mpi_bcastl0

C===============================================================================

C================== Routines from csm_share/shr/shr_file_mod.F90 ===============
C===============================================================================
CBOP ===========================================================================
C
C !IROUTINE: shr_file_getUnit -- Get a free FORTRAN unit number
C
C !DESCRIPTION: Get the next free FORTRAN unit number.
C
C !REVISION HISTORY:
C     2005-Dec-14 - E. Kluzek - creation
C     2007-Oct-21 - P. Worley - dumbed down for use in perf_mod
C
C !INTERFACE: ------------------------------------------------------------------

      INTEGER FUNCTION shr_file_getUnit ()

      implicit none

CEOP

C  !----- local parameters -----
C  Min unit number to give
      integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10
C  Max unit number to give
      integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99

C  !----- local variables -----
C  loop index
      integer(SHR_KIND_IN)   :: n
C  If unit opened or not
      logical                :: opened

   !----- formats -----
      character(*),parameter :: subName = '(shr_file_getUnit) '
      character(*),parameter :: F00   = "('(shr_file_getUnit) ',A,I4,A)"

C-------------------------------------------------------------------------------
C Notes:
C-------------------------------------------------------------------------------

C  ! --- Choose first available unit other than 0, 5, or 6  ------
      do n=shr_file_minUnit, shr_file_maxUnit
         inquire( n, opened=opened )
         if (n == 5 .or. n == 6 .or. opened) then
            cycle
         end if
         shr_file_getUnit = n
         return
      end do

      call shr_sys_abort( subName//': Error: no available units found' )

      END FUNCTION shr_file_getUnit
C===============================================================================

C===============================================================================
CBOP ===========================================================================
C
C !IROUTINE: shr_file_freeUnit -- Free up a FORTRAN unit number
C
C !DESCRIPTION: Free up the given unit number
C
C !REVISION HISTORY:
C     2005-Dec-14 - E. Kluzek - creation
C     2007-Oct-21 - P. Worley - dumbed down for use in perf_mod
C
C !INTERFACE: ------------------------------------------------------------------

      SUBROUTINE shr_file_freeUnit ( unit)

      implicit none

C !INPUT/OUTPUT PARAMETERS:

C  unit number to be freed
      integer(SHR_KIND_IN),intent(in) :: unit

CEOP

C  !----- local parameters -----
C  Min unit number to give
      integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10
C  Max unit number to give
      integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99

C  !----- formats -----
      character(*), parameter :: subName = '(shr_file_freeUnit) '
      character(*), parameter :: F00 = "('(shr_file_freeUnit) ',A,I4,A)"

C-------------------------------------------------------------------------------
C Notes:
C-------------------------------------------------------------------------------

      if (unit < 0 .or. unit > shr_file_maxUnit) then
!pw   if (s_loglev > 0) write(pu_logunit,F00) 'invalid unit number request:', unit
      else if (unit == 0 .or. unit == 5 .or. unit == 6) then
         call shr_sys_abort(
     &      subName//': Error: units 0, 5, and 6 must not be freed' )
      end if

      return

      END SUBROUTINE shr_file_freeUnit
C===============================================================================

C============= Routines from atm/cam/src/utils/namelist_utils.F90 ==============
C===============================================================================

      subroutine find_group_name(unit, group, status)

C---------------------------------------------------------------------------------------
C Purpose:
C Search a file that contains namelist input for the specified namelist group name.
C Leave the file positioned so that the current record is the first record of the
C input for the specified group.
C
C Method:
C Read the file line by line.  Each line is searched for an '&' which may only
C be preceded by blanks, immediately followed by the group name which is case
C insensitive.  If found then backspace the file so the current record is the
C one containing the group name and return success.  Otherwise return -1.
C
C Author:  B. Eaton, August 2007
C---------------------------------------------------------------------------------------

C  fortran unit attached to file
      integer,          intent(in)  :: unit
C  namelist group name
      character(len=*), intent(in)  :: group
C  0 for success, -1 if group name not found
      integer,          intent(out) :: status

C  ! Local variables

      integer           :: len_grp
C  io status
      integer           :: ios
C  first 80 characters of input record
      character(len=80) :: inrec
C  left adjusted input record
      character(len=80) :: inrec2
      character(len=len(group)) :: lc_group

C  !---------------------------------------------------------------------------

      len_grp = len_trim(group)
      lc_group = to_lower(group)

      ios = 0
      do while (ios <= 0)

         read(unit, '(a)', iostat=ios, end=102) inrec

         if (ios <= 0) then
C           ios < 0  indicates an end of record condition

C         ! look for group name in this record

C         ! remove leading blanks
            inrec2 = adjustl(inrec)

C         ! check for leading '&'
            if (inrec2(1:1) == '&') then

C           ! check for case insensitive group name
               if (trim(lc_group) == to_lower(inrec2(2:len_grp+1))) then

C              ! found group name.  backspace to leave file position at this record
                  backspace(unit)
                  status = 0
                  return

               end if
            end if
         end if

      end do

102   continue  ! end of file processing
      status = -1

      end subroutine find_group_name
C===============================================================================

C================ Routines from atm/cam/src/utils/string_utils.F90 =============
C===============================================================================

      function to_lower(str)

C-----------------------------------------------------------------------
C Purpose:
C Convert character string to lower case.
C
C Method:
C Use achar and iachar intrinsics to ensure use of ascii collating sequence.
C
C Author:  B. Eaton, July 2001
C
C $Id$
C-----------------------------------------------------------------------
      implicit none

C  String to convert to lower case
      character(len=*), intent(in) :: str
      character(len=len(str))      :: to_lower

C Local variables

C  Index
      integer :: i
C  ascii collating sequence
      integer :: aseq
C  integer to convert case
      integer :: upper_to_lower
C  Character temporary
      character(len=1) :: ctmp
C-----------------------------------------------------------------------
      upper_to_lower = iachar("a") - iachar("A")

      do i = 1, len(str)
         ctmp = str(i:i)
         aseq = iachar(ctmp)
         if ( aseq >= iachar("A") .and. aseq <= iachar("Z") )
     &        ctmp = achar(aseq + upper_to_lower)
         to_lower(i:i) = ctmp
      end do

      end function to_lower
C===============================================================================

      end module perf_utils
