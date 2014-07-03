!===============================================================================
! SVN $Id: shr_mpi_mod.F90 59033 2014-04-11 01:55:15Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140509/shr/shr_mpi_mod.F90 $
!===============================================================================

Module shr_mpi_mod

!-------------------------------------------------------------------------------
! PURPOSE: general layer on MPI functions
!-------------------------------------------------------------------------------

   use shr_kind_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none
   private

! PUBLIC: Public interfaces

   public :: shr_mpi_chkerr
   public :: shr_mpi_send
   public :: shr_mpi_recv
   public :: shr_mpi_bcast
   public :: shr_mpi_gathScatVInit
   public :: shr_mpi_gatherV
   public :: shr_mpi_scatterV
   public :: shr_mpi_sum
   public :: shr_mpi_min
   public :: shr_mpi_max
   public :: shr_mpi_commsize
   public :: shr_mpi_commrank
   public :: shr_mpi_initialized
   public :: shr_mpi_abort
   public :: shr_mpi_barrier
   public :: shr_mpi_init
   public :: shr_mpi_finalize

   interface shr_mpi_send ; module procedure &
     shr_mpi_sendi0, &
     shr_mpi_sendi1, &
     shr_mpi_sendr0, &
     shr_mpi_sendr1, &
     shr_mpi_sendr3
   end interface
   interface shr_mpi_recv ; module procedure &
     shr_mpi_recvi0, &
     shr_mpi_recvi1, &
     shr_mpi_recvr0, &
     shr_mpi_recvr1, &
     shr_mpi_recvr3
   end interface
   interface shr_mpi_bcast ; module procedure &
     shr_mpi_bcastc0, &
     shr_mpi_bcastc1, &
     shr_mpi_bcastl0, &
     shr_mpi_bcastl1, &
     shr_mpi_bcasti0, &
     shr_mpi_bcasti1, &
     shr_mpi_bcasti2, &
     shr_mpi_bcastr0, &
     shr_mpi_bcastr1, &
     shr_mpi_bcastr2, &
     shr_mpi_bcastr3
   end interface
   interface shr_mpi_gathScatVInit ; module procedure &
     shr_mpi_gathScatVInitr1
   end interface
   interface shr_mpi_gatherv ; module procedure &
     shr_mpi_gatherVr1
   end interface
   interface shr_mpi_scatterv ; module procedure &
     shr_mpi_scatterVr1
   end interface
   interface shr_mpi_sum ; module procedure &
     shr_mpi_sumi0, &
     shr_mpi_sumi1, &
     shr_mpi_sumb0, &
     shr_mpi_sumb1, &
     shr_mpi_sumr0, &
     shr_mpi_sumr1, &
     shr_mpi_sumr2, &
     shr_mpi_sumr3
   end interface
   interface shr_mpi_min ; module procedure &
     shr_mpi_mini0, &
     shr_mpi_mini1, &
     shr_mpi_minr0, &
     shr_mpi_minr1
   end interface
   interface shr_mpi_max ; module procedure &
     shr_mpi_maxi0, &
     shr_mpi_maxi1, &
     shr_mpi_maxr0, &
     shr_mpi_maxr1
   end interface

#include <mpif.h>         ! mpi library include file

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE shr_mpi_chkerr(rcode,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: rcode  ! input MPI error code
   character(*),         intent(in) :: string ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_chkerr) '
   character(MPI_MAX_ERROR_STRING)  :: lstring
   integer(SHR_KIND_IN)             :: len
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: layer on MPI error checking
!-------------------------------------------------------------------------------

   if (rcode /= MPI_SUCCESS) then
     call MPI_ERROR_STRING(rcode,lstring,len,ierr)
     write(s_logunit,*) trim(subName),":",lstring(1:len)
     call shr_mpi_abort(string,rcode)
   endif

END SUBROUTINE shr_mpi_chkerr

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! send value
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a single integer
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a real scalar
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sendr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(in) :: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to send to
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! mpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sendr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call MPI_SEND(array,lsize,MPI_REAL8,pid,tag,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_sendr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvi0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvi0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvi1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvi1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr0(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr0) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = 1

   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr1(lvec,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(out):: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid      ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag      ! tag
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr1) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_recvr3(array,pid,tag,comm,string)

   IMPLICIT none

   !----- arguments ---
   real   (SHR_KIND_R8), intent(out):: array(:,:,:)  ! in/out local values
   integer(SHR_KIND_IN), intent(in) :: pid           ! pid to recv from
   integer(SHR_KIND_IN), intent(in) :: tag           ! tag
   integer(SHR_KIND_IN), intent(in) :: comm          ! mpi communicator
   character(*),optional,intent(in) :: string        ! message

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_recvr3) '
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(array)

   call MPI_RECV(array,lsize,MPI_REAL8,pid,tag,comm,status,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_recvr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcasti0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastl0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastl0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastl0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastc0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   character(len=*), intent(inout)    :: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastc0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a character string
!-------------------------------------------------------------------------------

   lsize = len(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_CHARACTER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastc0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastc1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   character(len=*), intent(inout)    :: vec(:)   ! 1D vector
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastc1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a character string
!-------------------------------------------------------------------------------

   lsize = size(vec)*len(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_CHARACTER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastc1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr0(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastr0) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a real
!-------------------------------------------------------------------------------

   lsize = 1
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcasti1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastl1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   logical, intent(inout):: vec(:)      ! vector of 1
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastl1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a logical
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastl1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr1(vec,comm,string,pebcast)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(inout):: vec(:)   ! vector 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_bcastr1) '
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
!-------------------------------------------------------------------------------

   lsize = size(vec)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(vec,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr2(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcastr2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcasti2(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   integer,              intent(inout):: arr(:,:) ! array, 2d 
   integer(SHR_KIND_IN), intent(in)   :: comm     ! mpi communicator
   character(*),optional,intent(in)   :: string   ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcasti2) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 2d array of integers
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_INTEGER,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcasti2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_bcastr3(arr,comm,string,pebcast)

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(inout):: arr(:,:,:) ! array, 3d 
   integer(SHR_KIND_IN), intent(in)   :: comm       ! mpi communicator
   character(*),optional,intent(in)   :: string     ! message
   integer(SHR_KIND_IN), optional, intent(in)   :: pebcast  ! bcast pe (otherwise zero)

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr
   integer(SHR_KIND_IN)               :: lsize
   integer(SHR_KIND_IN)               :: lpebcast

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_bcastr3) '

!-------------------------------------------------------------------------------
! PURPOSE: Broadcast a 3d array of reals
!-------------------------------------------------------------------------------

   lsize = size(arr)
   lpebcast = 0
   if (present(pebcast)) lpebcast = pebcast

   call MPI_BCAST(arr,lsize,MPI_REAL8,lpebcast,comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_bcastr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_gathScatvInitr1(comm, rootid, locArr, glob1DArr, globSize, &
                                   displs, string )

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to gather/scatter on
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array of distributed data
   real(SHR_KIND_R8),    pointer      :: glob1DArr(:)  ! Global 1D array of gathered data
   integer(SHR_KIND_IN), pointer      :: globSize(:)   ! Size of each distributed piece
   integer(SHR_KIND_IN), pointer      :: displs(:)     ! Displacements for receive
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: npes          ! Number of MPI tasks
   integer(SHR_KIND_IN)               :: locSize       ! Size of local distributed data
   integer(SHR_KIND_IN), pointer      :: sendSize(:)   ! Size to send for initial gather
   integer(SHR_KIND_IN)               :: i             ! Index
   integer(SHR_KIND_IN)               :: rank          ! Rank of this MPI task
   integer(SHR_KIND_IN)               :: nSize         ! Maximum size to send
   integer(SHR_KIND_IN)               :: ierr          ! Error code
   integer(SHR_KIND_IN)               :: nSiz1D        ! Size of 1D global array
   integer(SHR_KIND_IN)               :: maxSize       ! Maximum size

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_gathScatvInitr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Setup arrays for a gatherv/scatterv operation
!-------------------------------------------------------------------------------

   locSize = size(locarr)
   call shr_mpi_commsize( comm, npes )
   call shr_mpi_commrank( comm, rank )
   allocate( globSize(npes) )
   !
   ! --- Gather the send global sizes from each MPI task -----------------------
   !
   allocate( sendSize(npes) )
   sendSize(:) = 1
   globSize(:) = 1
   call MPI_GATHER( locSize, 1, MPI_INTEGER, globSize, sendSize, &
                    MPI_INTEGER, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif
   deallocate( sendSize )
   !
   ! --- Prepare the displacement and allocate arrays -------------------------
   !
   allocate( displs(npes) )
   displs(1) = 0
   if ( rootid /= rank )then
      maxSize = 1
      globSize = 1
   else
      maxSize = maxval(globSize)
   end if
   nsiz1D  = min(maxSize,globSize(1))
   do i = 2, npes
      nSize = min(maxSize,globSize(i-1))
      displs(i) = displs(i-1) + nSize
      nsiz1D = nsiz1D + min(maxSize,globSize(i))
   end do
   allocate( glob1DArr(nsiz1D) )
   !----- Do some error checking for the root task arrays computed ----
   if ( rootid == rank )then
      if ( nsiz1D /= sum(globSize) ) &
         call shr_mpi_abort( subName//" : Error, size of global array not right" )
      if ( any(displs < 0) .or. any(displs >= nsiz1D) ) &
         call shr_mpi_abort( subName//" : Error, displacement array not right" )
      if ( (displs(npes)+globSize(npes)) /= nsiz1D ) &
         call shr_mpi_abort( subName//" : Error, displacement array values too big" )
   end if

END SUBROUTINE shr_mpi_gathScatvInitr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_gathervr1(locarr, locSize, glob1DArr, globSize, displs, rootid, &
                             comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(in)   :: locArr(:)     ! Local array
   real(SHR_KIND_R8),    intent(inout):: glob1DArr(:)  ! Global 1D array to receive in on
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to send this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to receive each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for receive
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to gather on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_gathervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Gather a 1D array of reals
!-------------------------------------------------------------------------------

   call MPI_GATHERV( locarr, locSize, MPI_REAL8, glob1Darr, globSize, displs, &
                     MPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_gathervr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_scattervr1(locarr, locSize, glob1Darr, globSize, displs, rootid, &
                              comm, string )

   IMPLICIT none

   !----- arguments -----
   real(SHR_KIND_R8),    intent(out)  :: locarr(:)     ! Local array
   real(SHR_KIND_R8),    intent(in)   :: glob1Darr(:)  ! Global 1D array to send from
   integer(SHR_KIND_IN), intent(in)   :: locSize       ! Number to receive this PE
   integer(SHR_KIND_IN), intent(in)   :: globSize(:)   ! Number to send to each PE
   integer(SHR_KIND_IN), intent(in)   :: displs(:)     ! Displacements for send
   integer(SHR_KIND_IN), intent(in)   :: rootid        ! MPI task to scatter on
   integer(SHR_KIND_IN), intent(in)   :: comm          ! mpi communicator
   character(*),optional,intent(in)   :: string        ! message

   !----- local -----
   integer(SHR_KIND_IN)               :: ierr          ! Error code

   !----- formats -----
   character(*),parameter             :: subName = '(shr_mpi_scattervr1) '

!-------------------------------------------------------------------------------
! PURPOSE: Scatter a 1D array of reals
!-------------------------------------------------------------------------------


   call MPI_SCATTERV( glob1Darr, globSize, displs, MPI_REAL8, locarr, locSize, &
                      MPI_REAL8, rootid, comm, ierr )
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_scattervr1


!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumb0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_I8), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_I8), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumb0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumb0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumb1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_I8), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_I8), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumb1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumb1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr2(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:)! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:)! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr2) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr2

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_sumr3(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:,:,:) ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:,:,:) ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_sumr3) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_sumr3

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_mini0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_mini0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_mini0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_mini1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_mini1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_mini1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_minr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_minr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_minr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_minr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_minr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_minr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxi0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec     ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxi0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxi0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxi1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   integer(SHR_KIND_IN), intent(in) :: lvec(:)  ! in/out local values
   integer(SHR_KIND_IN), intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxi1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxi1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxr0(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec     ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec     ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxr0) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxr0

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_maxr1(lvec,gvec,comm,string,all)

   IMPLICIT none

   !----- arguments ---
   real(SHR_KIND_R8),    intent(in) :: lvec(:)  ! in/out local values
   real(SHR_KIND_R8),    intent(out):: gvec(:)  ! in/out global values
   integer(SHR_KIND_IN), intent(in) :: comm     ! mpi communicator
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   !----- local ---
   character(*),parameter           :: subName = '(shr_mpi_maxr1) '
   logical                          :: lall
   character(SHR_KIND_CL)           :: lstring
   integer(SHR_KIND_IN)             :: reduce_type  ! mpi reduction type
   integer(SHR_KIND_IN)             :: lsize
   integer(SHR_KIND_IN)             :: gsize
   integer(SHR_KIND_IN)             :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
!-------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif
   if (present(string)) then
     lstring = trim(subName)//":"//trim(string)
   else
     lstring = trim(subName)
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call shr_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_ALLREDUCE")
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,0,comm,ierr)
     call shr_mpi_chkerr(ierr,trim(lstring)//" MPI_REDUCE")
   endif

END SUBROUTINE shr_mpi_maxr1

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_commsize(comm,size,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: size
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_commsize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI commsize
!-------------------------------------------------------------------------------

   call MPI_COMM_SIZE(comm,size,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_commsize

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_commrank(comm,rank,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   integer,intent(out)                :: rank
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_commrank) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI commrank
!-------------------------------------------------------------------------------

   call MPI_COMM_RANK(comm,rank,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_commrank

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_initialized(flag,string)

   IMPLICIT none

   !----- arguments ---
   logical,intent(out)                :: flag
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_initialized) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI initialized
!-------------------------------------------------------------------------------

   call MPI_INITIALIZED(flag,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_initialized

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_abort(string,rcode)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message
   integer,optional,intent(in)        :: rcode    ! optional code

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_abort) '
   integer(SHR_KIND_IN)               :: ierr
   integer                            :: rc       ! return code

!-------------------------------------------------------------------------------
! PURPOSE: MPI abort
!-------------------------------------------------------------------------------

   if ( present(string) .and. present(rcode) ) then
      write(s_logunit,*) trim(subName),":",trim(string),rcode
   endif
   if ( present(rcode) )then
      rc = rcode
   else
      rc = 1001
   end if
   call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)

END SUBROUTINE shr_mpi_abort

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_barrier(comm,string)

   IMPLICIT none

   !----- arguments ---
   integer,intent(in)                 :: comm
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_barrier) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI barrier
!-------------------------------------------------------------------------------

   call MPI_BARRIER(comm,ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_barrier

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_init(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_init) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI init
!-------------------------------------------------------------------------------

   call MPI_INIT(ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_init

!===============================================================================
!===============================================================================

SUBROUTINE shr_mpi_finalize(string)

   IMPLICIT none

   !----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   !----- local ---
   character(*),parameter             :: subName = '(shr_mpi_finalize) '
   integer(SHR_KIND_IN)               :: ierr

!-------------------------------------------------------------------------------
! PURPOSE: MPI finalize
!-------------------------------------------------------------------------------

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
   if (present(string)) then
     call shr_mpi_chkerr(ierr,subName//trim(string))
   else
     call shr_mpi_chkerr(ierr,subName)
   endif

END SUBROUTINE shr_mpi_finalize

!===============================================================================
!===============================================================================

END MODULE shr_mpi_mod
