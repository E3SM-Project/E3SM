!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_abort -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_abort(comm,errorcode,ier)
    implicit none
    integer,intent(in) :: comm
    integer,intent(in) :: errorcode
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_abort'

    call exit(errorcode)
    ier=0
end subroutine mpi_abort

