!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_comm_free -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_comm_free(comm,ier)
    use m_mpi0,only : mpi0_initialized
    implicit none
    integer,intent(inout) :: comm
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_comm_free'

  if(.not.mpi0_initialized) call mpi_init(ier)

    comm=-1
    ier=0

end subroutine mpi_comm_free
