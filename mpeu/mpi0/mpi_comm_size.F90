!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_comm_size -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_comm_size(comm,size,ier)
    use m_mpi0,only : mpi0_initialized
    implicit none
    integer,intent(in) :: comm
    integer,intent(out) :: size
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_comm_size'

  if(.not.mpi0_initialized) call mpi_init(ier)

    size=1
    ier=0
end subroutine mpi_comm_size
