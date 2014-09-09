module pio_utils
  use pio_types, only : file_desc_t, var_desc_t
  use pio_types, only : pio_int, pio_real, pio_double, pio_char
  use pio_types, only : iotype_netcdf, iotype_pnetcdf, PIO_internal_error
  use pio_types, only : PIO_iotype_netcdf4p, pio_iotype_netcdf4c
  use pio_types, only : PIO_bcast_error 
  use pio_kinds, only : i4, r4, r8
  use pio_support, only : checkmpireturn, piodie, Debug

#ifdef _NETCDF
  use netcdf            ! _EXTERNAL
#endif
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
#ifdef USE_PNETCDF_MOD
  use pnetcdf
#endif
  implicit none
  private
#ifdef NO_MPIMOD
  include 'mpif.h'      ! _EXTERNAL
#endif
#ifdef _PNETCDF
#ifndef USE_PNETCDF_MOD
#include <pnetcdf.inc>   /* _EXTERNAL */
#endif
#endif

  public :: check_netcdf 
  public :: bad_iotype 

  

contains
  subroutine replace_c_null(istr)
    use iso_c_binding, only : C_NULL_CHAR
    character(len=*),intent(inout) :: istr
    integer :: i
    do i=1,len(istr)
       if(istr(i:i) == C_NULL_CHAR) istr(i:i)=''
    end do
  end subroutine replace_c_null



  subroutine check_netcdf(File, status, filestr, line)
    type(file_desc_t), intent(in) :: file
    integer, intent(inout) :: status
    character(len=*), intent(in) :: filestr
    integer, intent(in) :: line

    integer :: mpierr, iotype

!  Three choices for error handling:
!  1: abort on error from any task           PIO_INTERNAL_ERROR
!  2: broadcast an error from io_rank 0      PIO_BCAST_ERROR
!  3: do nothing - allow the user to handle it PIO_RETURN_ERROR
!
    iotype = file%iotype
    
    if(Debug) call mpi_barrier(file%iosystem%union_comm, mpierr)

    select case(iotype)
    case(iotype_pnetcdf)
#ifdef _PNETCDF
       if(file%iosystem%error_handling==PIO_INTERNAL_ERROR) then
          if(status /= nf_noerr) then
             call piodie(filestr,line,trim(nfmpi_strerror(status)))
          end if
       else if(file%iosystem%error_handling==PIO_BCAST_ERROR) then
          call MPI_BCAST(status,1,MPI_INTEGER,file%iosystem%iomaster,File%iosystem%my_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
       end if

#endif
    case(iotype_netcdf,pio_iotype_netcdf4p,pio_iotype_netcdf4c)
#ifdef _NETCDF
       if(status /= nf90_noerr) then
          print *,trim(nf90_strerror(status))
       endif
       if(File%iosystem%error_handling==PIO_INTERNAL_ERROR) then
          if(status /= nf90_noerr) then
             call piodie(filestr,line,trim(nf90_strerror(status)))
          end if
       else if(file%iosystem%error_handling==PIO_BCAST_ERROR) then
          call MPI_BCAST(status,1,MPI_INTEGER,file%iosystem%iomaster,File%iosystem%my_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
       end if
#endif
    end select

  end subroutine check_netcdf



!>
!! @private
!<
  subroutine bad_iotype(iotype,file,line)
    integer iotype
    character(len=*) file
    integer line

#ifndef _PNETCDF
    if (iotype==iotype_pnetcdf) then
       call piodie(file,line,'PNETCDF not enabled in the build')
    endif
#endif
#ifndef _NETCDF
    if (iotype==iotype_netcdf) then
       call piodie(file,line,'NETCDF not enabled in the build')
    endif
#endif
#ifndef _NETCDF4
    if (iotype==PIO_iotype_netcdf4p .or. iotype==pio_iotype_netcdf4c) then
       call piodie(file,line,'NETCDF4 not enabled in the build')
    endif
#endif
    print *,'Invalid iotype, value=',iotype
    call piodie(file,line,'Quitting')

  end subroutine bad_iotype

  

end module pio_utils
