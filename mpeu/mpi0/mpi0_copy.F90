!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi0_copy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mpi0_copy(sbuf,stype,scount,rbuf,rcount,rtype,ier)
      use m_mpi0
      implicit none
      integer,dimension(*),intent(in) :: sbuf
      integer,		   intent(in) :: scount
      integer,		   intent(in) :: stype

      integer,dimension(*),intent(out):: rbuf
      integer,		   intent(in) :: rcount
      integer,		   intent(in) :: rtype

      integer,		   intent(out):: ier

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi0_copy'

  if(stype/=rtype .or. scount/=rcount) then
    ier=-1
    return
  endif

  ier=0
  if(	  stype==MPI_INTEGER		) then
    call copy_INTEGER(sbuf,rbuf,scount)
  elseif( stype==MPI_REAL		) then
    call copy_REAL(sbuf,rbuf,scount)
  elseif( stype==MPI_DOUBLE_PRECISION	) then
    call copy_DOUBLE_PRECISION(sbuf,rbuf,scount)
  elseif( stype==MPI_LOGICAL		) then
    call copy_LOGICAL(sbuf,rbuf,scount)
  elseif( stype==MPI_CHARACTER		) then
    call copy_CHARACTER(sbuf,rbuf,scount)
  elseif( stype==MPI_INTEGER4		) then
    call copy_INTEGER4(sbuf,rbuf,scount)
  elseif( stype==MPI_REAL4		) then
    call copy_REAL4(sbuf,rbuf,scount)
  elseif( stype==MPI_REAL8		) then
    call copy_REAL8(sbuf,rbuf,scount)
  else
    ier=-1
  endif

end subroutine mpi0_copy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_INTEGER - copy INTEGERs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_INTEGER(sbuf,rbuf,n)
      implicit none
      integer,dimension(*),intent(in)  :: sbuf
      integer,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_INTEGER'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_INTEGER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_REAL - copy REALs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_REAL(sbuf,rbuf,n)
      implicit none
      real   ,dimension(*),intent(in)  :: sbuf
      real   ,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_REAL'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_REAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_DOUBLE_PRECISION - copy DOUBLE_PRECISIONs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_DOUBLE_PRECISION(sbuf,rbuf,n)
      implicit none
      double precision,dimension(*),intent(in)  :: sbuf
      double precision,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_DOUBLE_PRECISION'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_DOUBLE_PRECISION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_LOGICAL - copy LOGICALs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_LOGICAL(sbuf,rbuf,n)
      implicit none
      logical,dimension(*),intent(in)  :: sbuf
      logical,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_LOGICAL'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_LOGICAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_CHARACTER - copy CHARACTERs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_CHARACTER(sbuf,rbuf,n)
      implicit none
      character,dimension(*),intent(in)  :: sbuf
      character,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_CHARACTER'
  integer :: i

  do i=1,n
    rbuf(i)=sbuf(i)
  end do

end subroutine copy_CHARACTER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_INTEGER4 - copy INTEGER*4
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_INTEGER4(sbuf,rbuf,n)
      implicit none
      integer*4,dimension(*),intent(in)  :: sbuf
      integer*4,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_INTEGER4'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_INTEGER4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_REAL4 - copy REAL*4
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_REAL4(sbuf,rbuf,n)
      implicit none
      real*4,dimension(*),intent(in)  :: sbuf
      real*4,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_REAL4'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_REAL4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: copy_REAL8 - copy REAL*8
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_REAL8(sbuf,rbuf,n)
      implicit none
      real*8,dimension(*),intent(in)  :: sbuf
      real*8,dimension(*),intent(out) :: rbuf
      integer,intent(in) :: n

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='copy_REAL8'

  rbuf(1:n)=sbuf(1:n)

end subroutine copy_REAL8
