!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: mpi_error_string -
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine MPI_error_string(ierror,cerror,ln,ier)
    implicit none
    integer,intent(in) :: ierror
    character(len=*),intent(out) :: cerror
    integer,intent(out) :: ln
    integer,intent(out) :: ier

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='mpi_error_string'
  character(len=16) :: cint

    write(cint,'(i16)') ierror
    cerror='MPI0: ierror = '//trim(adjustl(cint))
    ln=len_trim(cerror)
    ier=0
end subroutine mpi_error_string
