!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mpi0 - mpi0 flag
!
! !DESCRIPTION:
!
! !INTERFACE:

	module m_mpi0
	implicit none
	private	! except

	public :: mpi0_initialized	! The class data structure
	public :: MPI_INTEGER
	public :: MPI_REAL
	public :: MPI_DOUBLE_PRECISION
	public :: MPI_COMPLEX
	public :: MPI_DOUBLE_COMPLEX
	public :: MPI_LOGICAL
	public :: MPI_CHARACTER
	public :: MPI_BYTE
	public :: MPI_2INTEGER
	public :: MPI_2REAL
	public :: MPI_2DOUBLE_PRECISION
	! public :: MPI_2COMPLEX	! not supported on IRIX64
	! public :: MPI_2DOUBLE_COMPLEX	! not supported on IRIX64
	public :: MPI_INTEGER1
	public :: MPI_INTEGER2
	public :: MPI_INTEGER4
	! public :: MPI_REAL2		! not supported on IRIX64
	public :: MPI_REAL4
	public :: MPI_REAL8

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

	character(len=*),parameter :: myname='m_mpi0'
	include "mpif.h"

	logical,save :: mpi0_initialized=.false.

	end module m_mpi0
