!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_bcast -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mpi_bcast(buffer,count,datatype,root,comm,ierror)
      use m_mpi0,only : mpi0_initialized
      implicit none
      integer,dimension(*)  :: buffer	! intent is (in or out)
      integer,intent(in)    :: count
      integer,intent(in)    :: datatype
      integer,intent(in)    :: root
      integer,intent(in)    :: comm
      integer,intent(out)   :: ierror

! !REVISION HISTORY:
! 	29Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_bcast'

  if(.not.mpi0_initialized) call mpi_init(ierror)

  ierror=0
end subroutine mpi_bcast
