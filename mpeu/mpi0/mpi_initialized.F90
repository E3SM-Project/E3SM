!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_initialized -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_initialized(flag,ier)
    use m_mpi0,only : mpi0_initialized
    implicit none
    logical,intent(out) :: flag
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_initialized'

    flag=mpi0_initialized
    ier=0
end subroutine mpi_initialized
