!===============================================================================
! SVN $Id: shr_mpi_mod.F90 65839 2014-11-30 14:40:15Z jedwards $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_150116/shr/shr_mpi_mod.F90 $
!===============================================================================

Module bshr_mpi_mod

!-------------------------------------------------------------------------------
! PURPOSE: general layer on MPI functions
!-------------------------------------------------------------------------------

   use bshr_kind_mod, only : SHR_KIND_IN
   use bshr_kind_mod, only : SHR_KIND_R8
   use bshr_kind_mod, only : SHR_KIND_I8
   use bshr_kind_mod, only : SHR_KIND_CL
   use bshr_log_mod, only  : s_loglev  => shr_log_Level
   use bshr_log_mod, only  : s_logunit => shr_log_Unit

   implicit none
   private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
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
     shr_mpi_bcasti80, &
     shr_mpi_bcasti81, &
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

   ! mpi library include file

#include "bmpif.h"

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
   if (rcode > 0) continue
   if (len(string) > 0) continue



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
   if (lvec > 0) continue
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   if (lvec(1) > 0) continue
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   if (lvec > 0) continue
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   if (lvec(1) > 0) continue
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   if (array(1, 1, 1) > 0) continue
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   lvec = 0
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue


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
   lvec(1) = 0
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue


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
   lvec = 0
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue


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
   lvec(1) = 0
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue


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
   array(1, 1, 1) = 0
   if (pid > 0) continue
   if (tag > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue



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
   vec = vec
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



END SUBROUTINE shr_mpi_bcasti0

SUBROUTINE shr_mpi_bcasti80(vec,comm,string,pebcast)

  IMPLICIT none

  !----- arguments ---
  integer(SHR_KIND_I8), intent(inout):: vec      ! vector of 1
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
   vec = vec
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



END SUBROUTINE shr_mpi_bcasti80

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
   vec = vec
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue


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
   vec = vec
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   vec(1) = vec(1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue


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
   vec = vec
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   vec(1) = vec(1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue


END SUBROUTINE shr_mpi_bcasti1

SUBROUTINE shr_mpi_bcasti81(vec,comm,string,pebcast)

  IMPLICIT none

  !----- arguments ---
  integer(SHR_KIND_I8), intent(inout):: vec(:)   ! vector
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
   vec(1) = vec(1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



END SUBROUTINE shr_mpi_bcasti81

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
   vec(1) = vec(1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   vec(1) = vec(1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue

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
   arr(1, 1) = arr(1, 1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   arr(1, 1) = arr(1, 1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   arr(1, 1, 1) = arr(1, 1, 1)
   if (comm > 0) continue
   if (len(string) > 0) continue
   if (pebcast > 0) continue



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
   if (comm > 0) continue
   if (rootid > 0) continue
   glob1DArr(1) = locArr(1)
   if (locSize > 0) continue
   if (globSize(1) > 0) continue
   if (displs(1) > 0) continue
   if (len(string) > 0) continue


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
   glob1DArr(1) = locarr(1)
   if (locSize > 0) continue
   if (globSize(1) > 0) continue
   if (displs(1) > 0) continue
   if (rootid > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue


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
   locarr(1) = glob1DArr(1)
   if (locSize > 0) continue
   if (globSize(1) > 0) continue
   if (displs(1) > 0) continue
   if (rootid > 0) continue
   if (comm > 0) continue
   if (len(string) > 0) continue




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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1,1) = lvec(1,1)
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1,1,1) = lvec(1,1,1)
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue

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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue



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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec = lvec
   if (all) continue

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
   if (len(string) > 0) continue
   if (comm > 0) continue
   gvec(1) = lvec(1)
   if (all) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue
   size = 0

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
   if (len(string) > 0) continue
   if (comm > 0) continue
   rank = 0

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

   if (present(string)) then
      if (len(string) > 0) continue
   endif
   flag = .true.


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
   if (len(string) > 0) continue
   if (rcode > 0) continue


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
   if (len(string) > 0) continue
   if (comm > 0) continue


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
   if (len(string) > 0) continue



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
   if (len(string) > 0) continue


END SUBROUTINE shr_mpi_finalize

!===============================================================================
!===============================================================================

END MODULE bshr_mpi_mod
