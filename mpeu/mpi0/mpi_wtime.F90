!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_wtime -
!
! !DESCRIPTION:
!
! !INTERFACE:

  function MPI_wtime()
    implicit none
    double precision :: MPI_wtime

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_wtime'
  real*8, dimension(1:5) :: zts

    call get_zeits(zts)
    MPI_wtime=zts(1)
end function mpi_wtime
