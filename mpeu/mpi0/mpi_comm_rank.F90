!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_comm_rank -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_comm_rank(comm,rank,ier)
    use m_mpi0,only : mpi0_initialized
    implicit none
    integer,intent(in) :: comm
    integer,intent(out) :: rank
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_comm_rank'

  if(.not.mpi0_initialized) call mpi_init(ier)

    rank=0
    ier=0
end subroutine mpi_comm_rank
