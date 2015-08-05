!>
!! @file pio_support.F90
!! @brief internal code for compiler workarounds, aborts and debug functions
!!
!! $Revision$
!! $LastChangedDate$
!<
module pio_support
  use pio_kinds
  use iso_c_binding
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
  private
#ifdef NO_MPIMOD
  include 'mpif.h'    ! _EXTERNAL
#endif
  public :: piodie
  public :: CheckMPIreturn
  public :: pio_readdof
  public :: pio_writedof
  public :: replace_c_null

  logical, public :: Debug=.FALSE.
  logical, public :: DebugIO=.FALSE.
  logical, public :: DebugAsync=.FALSE.
  integer,private,parameter :: versno = 1001

  character(len=*), parameter :: modName='pio_support'

contains
  subroutine replace_c_null(istr)
    use iso_c_binding, only : C_NULL_CHAR
    character(len=*),intent(inout) :: istr
    integer :: i
    do i=1,len(istr)
       if(istr(i:i) == C_NULL_CHAR) exit
    end do
    istr(i:len(istr))=''
  end subroutine replace_c_null

  subroutine piodie (file,line, msg, ival1, msg2, ival2, msg3, ival3, mpirank)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Abort the model for abnormal termination
    !
    ! Author: Jim Edwards
    !
    ! Change History
    ! 20070608 R. Loy  added optional args
    !-----------------------------------------------------------------------
    ! $Id$
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: file
    integer,intent(in) :: line
    character(len=*), intent(in), optional :: msg,msg2,msg3
    integer,intent(in),optional :: ival1,ival2,ival3, mpirank

    character(len=*), parameter :: subName=modName//'::pio_die'
    integer :: ierr, myrank=-1
    
    if(present(mpirank)) myrank=mpirank

    if (present(ival3)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
            msg,ival1,msg2,ival2,msg3,ival3
    else if (present(msg3)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
            msg,ival1,msg2,ival2,msg3
    else if (present(ival2)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2,ival2
    else if (present(msg2)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2
    else if (present(ival1)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1
    else if (present(msg)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg
    else
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': (no message)'
    endif


#if defined(CPRXLF) && !defined(BGQ)
  close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
  call xl__trbk()
#endif

    ! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
    ! error code of 257
    call mpi_abort (MPI_COMM_WORLD, 1, ierr)  

#ifdef CPRNAG
    stop
#else
    call abort
#endif


  end subroutine piodie

!=============================================
!  CheckMPIreturn:
!
!      Check and prints an error message
!  if an error occured in a MPI subroutine.
!=============================================
  subroutine CheckMPIreturn(locmesg, errcode, file, line)

     character(len=*), intent(in) :: locmesg
     integer(i4), intent(in) :: errcode
     character(len=*),optional :: file
     integer, intent(in),optional :: line
     character(len=MPI_MAX_ERROR_STRING) :: errorstring

     integer(i4) :: errorlen

     integer(i4) :: ierr
     if (errcode .ne. MPI_SUCCESS) then
        call MPI_Error_String(errcode,errorstring,errorlen,ierr)
        write(*,*) TRIM(ADJUSTL(locmesg))//errorstring(1:errorlen)
        if(present(file).and.present(line)) then
           call piodie(file,line)
        endif
     end if
  end subroutine CheckMPIreturn

!>
!! @defgroup pio_writedof MAPtools 
!!  Fortran interface to write a mapping file
!!
!<

  subroutine pio_writedof (file, gdims, DOF, comm, punit)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Write a DOF to standard format
    !
    ! Author: T Craig
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*),intent(in) :: file
    integer, intent(in) :: gdims(:)
    integer(PIO_OFFSET_KIND)  ,intent(in) :: dof(:)
    integer         ,intent(in) :: comm
    integer,optional,intent(in) :: punit
    integer :: err
    integer :: ndims
    

    interface
       integer(c_int) function PIOc_writemap_from_f90(file, ndims, gdims, maplen, map, f90_comm) &
            bind(C,name="PIOc_writemap_from_f90")
         use iso_c_binding
         character(C_CHAR), intent(in) :: file
         integer(C_INT), value, intent(in) :: ndims
         integer(C_INT), intent(in) :: gdims(*)
         integer(C_SIZE_T), value, intent(in) :: maplen 
         integer(C_SIZE_T), intent(in) :: map(*)
         integer(C_INT), value, intent(in) :: f90_comm
       end function PIOc_writemap_from_f90
    end interface
    ndims = size(gdims)
    err = PIOc_writemap_from_f90(trim(file)//C_NULL_CHAR, ndims, gdims, int(size(dof),C_SIZE_T), dof, comm)

  end subroutine pio_writedof
!>
!! @addtogroup pio_readdof MAPtools
!!  Fortran interface to read a mapping file
!!
!<

  subroutine pio_readdof (file, ndims, gdims, DOF, comm, punit)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Read a DOF to standard format
    !
    ! Author: T Craig
    !
    ! Change History
    ! 
    !-----------------------------------------------------------------------
    ! $Id$
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*),intent(in) :: file
    integer(PIO_OFFSET_KIND),pointer:: dof(:)
    integer         ,intent(in) :: comm
    integer,optional,intent(in) :: punit
    integer, intent(out) :: ndims
    integer, pointer :: gdims(:)
    integer(PIO_OFFSET_KIND) :: maplen
    integer :: ierr
    type(C_PTR) :: tgdims, tmap
    interface
       integer(C_INT) function PIOc_readmap_from_f90(file, ndims, gdims, maplen, map, f90_comm) &
            bind(C,name="PIOc_readmap_from_f90") 
         use iso_c_binding
         character(C_CHAR), intent(in) :: file
         integer(C_INT), intent(out) :: ndims
         type(C_PTR), intent(out) :: gdims
         integer(C_SIZE_T), intent(out) :: maplen
         type(C_PTR) :: map
         integer(C_INT), value, intent(in) :: f90_comm
       end function PIOc_readmap_from_f90
    end interface

    ierr = PIOc_readmap_from_f90(trim(file)//C_NULL_CHAR, ndims, tgdims, maplen, tmap, comm);

    call c_f_pointer(tgdims, gdims, (/ndims/))
    call c_f_pointer(tmap, DOF, (/maplen/))
!    DOF = DOF+1
  end subroutine pio_readdof



end module pio_support
