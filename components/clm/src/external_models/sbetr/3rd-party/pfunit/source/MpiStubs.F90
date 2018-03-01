!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MpiStubs
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module MpiStubs_mod
   implicit none
   private

   public :: MPI_COMM_WORLD
   public :: MPI_COMM_NULL
   public :: MPI_COMM_SUCCESS
   public :: MPI_Comm_rank
   public :: MPI_Comm_size
   public :: MPI_Comm_dup
   public :: MPI_Comm_group
   public :: MPI_Group_range_incl
   public :: MPI_Comm_create

   integer, parameter :: MPI_COMM_WORLD = -1
   integer, parameter :: MPI_COMM_NULL = -1
   integer, parameter :: MPI_COMM_SUCCESS = 0

   integer :: nextCommunicator = MPI_COMM_WORLD

contains

   subroutine MPI_Comm_rank(comm, rank, ier)
      integer, intent(in) :: comm
      integer, intent(out) :: rank
      integer, intent(out) :: ier

      rank = 0
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Comm_rank

   subroutine MPI_Comm_size(comm, size, ier)
      integer, intent(in) :: comm
      integer, intent(out) :: size
      integer, intent(out) :: ier

      size = 1
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Comm_size

   subroutine MPI_Comm_dup(comm, newComm, ier)
      integer, intent(in) :: comm
      integer, intent(out) :: newComm
      integer, intent(out) :: ier


      newComm = newCommunicator()
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Comm_dup

   subroutine MPI_Comm_group(comm, group, ier)
      integer, intent(in) :: comm
      integer, intent(out) :: group
      integer, intent(out) :: ier

      group = 0
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Comm_group

   subroutine MPI_Group_range_incl(group, n, ranges, newGroups, ier)
      integer, intent(in) :: group
      integer, intent(in) :: n
      integer, intent(in) :: ranges(3,n)
      integer, intent(out) :: newGroups(n)
      integer, intent(out) :: ier

      newGroups = 0
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Group_range_incl

   subroutine MPI_Comm_create(comm, group, newComm, ier)
      integer, intent(in) :: comm
      integer, intent(in) :: group
      integer, intent(out) :: newComm
      integer, intent(out) :: ier

      newComm = newCommunicator()
      ier = MPI_COMM_SUCCESS
      
   end subroutine MPI_Comm_create

   integer function newCommunicator()
      nextCommunicator = nextCommunicator + 1
      newCommunicator = nextCommunicator
   end function newCommunicator
   
end module MpiStubs_mod
