!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_wtick -
!
! !DESCRIPTION:
!
! !INTERFACE:

  function MPI_wtick()
    implicit none
    double precision :: MPI_wtick

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_wtick'
  real*8, dimension(1:1) :: tic

    call get_ztick(tic)
    MPI_wtick=tic(1)
end function mpi_wtick
