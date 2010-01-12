module pio_utils
  use pio_types, only : file_desc_t, var_desc_t
  use pio_types, only : pio_int, pio_real, pio_double, pio_char
  use pio_types, only : iotype_netcdf, iotype_pnetcdf, PIO_internal_error
  use pio_types, only : PIO_iotype_netcdf4p, pio_iotype_netcdf4c
  use pio_types, only : PIO_bcast_error 
  use pio_kinds, only : i4, r4, r8
  use pio_support, only : checkmpireturn, piodie

#ifdef _NETCDF
  use netcdf            ! _EXTERNAL
#endif

  implicit none
  private

  include 'mpif.h'      ! _EXTERNAL
#ifdef _PNETCDF
#include <pnetcdf.inc>   /* _EXTERNAL */
#endif

  public :: check_netcdf 
  public :: bad_iotype 
#ifdef _PNETCDF
  public :: pnetcdf_version_check
#endif
  

contains

  subroutine check_netcdf(File, status, filestr, line)
    type(file_desc_t), intent(in) :: file
    integer, intent(inout) :: status
    character(len=*), intent(in) :: filestr
    integer, intent(in) :: line

    integer :: mpierr

!  Three choices for error handling:
!  1: abort on error from any task           PIO_INTERNAL_ERROR
!  2: broadcast an error from io_rank 0      PIO_BCAST_ERROR
!  3: do nothing - allow the user to handle it PIO_RETURN_ERROR
!


    if(file%iotype==iotype_pnetcdf) then
#ifdef _PNETCDF
       if(file%iosystem%error_handling==PIO_INTERNAL_ERROR) then
          if(status /= nf_noerr) then
             call piodie(filestr,line,trim(nfmpi_strerror(status)))
          end if
       else if(file%iosystem%error_handling==PIO_BCAST_ERROR) then
          call MPI_BCAST(status,1,MPI_INTEGER,file%iosystem%iomaster,File%iosystem%comp_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
       end if

#endif
    else if(file%iotype==iotype_netcdf) then
#ifdef _NETCDF
       if(File%iosystem%error_handling==PIO_INTERNAL_ERROR) then
          if(status /= nf90_noerr) then
             call piodie(filestr,line,trim(nf90_strerror(status)))
          end if
       else if(file%iosystem%error_handling==PIO_BCAST_ERROR) then
          call MPI_BCAST(status,1,MPI_INTEGER,file%iosystem%iomaster,File%iosystem%comp_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
       end if
#endif
    end if

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

#ifdef _PNETCDF
  subroutine pnetcdf_version_check()
    character(len=80) :: version
    integer :: dot1, dot2, s
    integer :: v1, v2, v3
    integer, parameter :: rv1=1, rv2=1, rv3=0
    logical,save :: pnetcdf_version_okay=.false.
    character(len=80) :: error

    if(pnetcdf_version_okay) return


    write(error,*) 'Pnetcdf version appears to be older than minimum required ',rv1,'.',rv2,'.',rv3

    version = nfmpi_inq_libvers()

    dot1 = index(version,'.')
    dot2 = index(version,'.',.true.)

    s = index(version(1:dot1),' ',.true.)
    read(version(s:dot1-1),'(i)') v1
    read(version(dot1+1:dot2-1),'(i)') v2
    s = dot2+index(version(dot2:),' ')
!   print *,__FILE__,__LINE__,version,dot2,s
    read(version(dot2+1:s-1),'(i)') v3

    if(v1<rv1) then
       call piodie(version,__LINE__,error)
    else if(v1==rv1) then
       if(v2<rv2) then
          call piodie(version,__LINE__,error)
       else if(v2==rv2) then
          if(v3<rv3) then
             call piodie(version,__LINE__,error)
          end if
       end if
    end if
    pnetcdf_version_okay=.true.
  end subroutine pnetcdf_version_check
#endif


end module pio_utils
