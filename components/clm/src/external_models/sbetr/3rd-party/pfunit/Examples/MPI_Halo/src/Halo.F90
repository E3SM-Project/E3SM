module Halo_mod
   implicit none

   include 'mpif.h'

contains

   ! Fills first and last column of array with interior data from
   ! neigbboring process.  Assumes periodic domain.
   subroutine haloFill(array, communicator)
      real, intent(in) :: array(:,0:)
      integer, intent(in) :: communicator
      
      integer :: jm, im
      integer :: PE, PE_left, PE_right, NPES
      integer :: status(MPI_STATUS_SIZE)
      integer :: ier

      integer :: tagA = 10
      integer :: tagB = 11

      im = size(array,1)
      jm = size(array,2) - 2 ! exclude boundaries

      call Mpi_Comm_rank(communicator, PE, ier)
      call Mpi_Comm_size(communicator, NPES, ier)

      PE_left = mod(PE - 1 + NPES, NPES)
      PE_right = mod(PE + 1, NPES)

      ! Send right, receive left
      call mpi_sendrecv( &
           & array(:,jm), im, MPI_REAL, PE_right, tagA, &
           & array(:,0), im, MPI_REAL, PE_left, tagA, &
           & communicator, status, ier)

      ! Send left, receive right
      call mpi_sendrecv( &
           & array(:,1), im, MPI_REAL, PE_left, tagB, &
           & array(:,jm+1), im, MPI_REAL, PE_right, tagB, &
           & communicator, status, ier)

   end subroutine haloFill

end module Halo_mod
