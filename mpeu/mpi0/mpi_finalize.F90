!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_finalize -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_finalize(ier)
    use m_mpi0
    implicit none
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_finalize'

#include "Const_MPItype.H"

  if(.not.mpi0_initialized) call mpi_init(ier)

  mpi0_initialized	=.false.

#ifndef Const_MPItype_
  MPI_INTEGER		=-1
  MPI_REAL		=-1
  MPI_DOUBLE_PRECISION	=-1
  MPI_COMPLEX		=-1
  MPI_DOUBLE_COMPLEX	=-1
  MPI_LOGICAL		=-1
  MPI_CHARACTER		=-1
  MPI_BYTE		=-1
  MPI_2INTEGER		=-1
  MPI_2REAL		=-1
  MPI_2DOUBLE_PRECISION	=-1
  ! MPI_2COMPLEX	=-1	! not supported on IRIX64
  ! MPI_2DOUBLE_COMPLEX	=-1	! not supported on IRIX64
  MPI_INTEGER1		=-1
  MPI_INTEGER2		=-1
  MPI_INTEGER4		=-1
  ! MPI_REAL2		=-1	! not supported on IRIX64
  MPI_REAL4		=-1
  MPI_REAL8		=-1
#endif

  ier=0

end subroutine mpi_finalize
