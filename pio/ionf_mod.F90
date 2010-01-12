#define _FILE_ "ionf_mod.F90"
module ionf_mod
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf      ! _EXTERNAL
#endif
  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types

  use pio_utils, only: bad_iotype, check_netcdf
  use pio_support, only : Debug, DebugIO, piodie   
#ifdef _NETCDF
  use netcdf            ! _EXTERNAL
#endif
  use pio_support, only : CheckMPIReturn

  implicit none
  private

  include 'mpif.h'      ! _EXTERNAL
#ifdef _PNETCDF
#include <pnetcdf.inc>   /* _EXTERNAL */
#endif
 

   public :: create_nf
   public :: open_nf 
   public :: close_nf 
   public :: sync_nf 

contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! create_nf
  !

  integer function create_nf(File,fname, amode) result(ierr)

    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in)      :: fname
    integer(i4),  intent(in) :: amode
    integer(i4) :: iotype, mpierr
    integer :: nmode

    
    if(iotype==iotype_netcdf4c) then
       nmode=ior(amode,NF90_NETCDF4)
    else
       nmode=amode
    end if
    if(Debug) print *, 'Create file ', fname, 'with amode ',nmode
    ierr=PIO_noerr
    if(File%iosystem%ioproc) then
       iotype = File%iotype 

       select case (iotype) 

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr  = nfmpi_create(File%iosystem%IO_comm,fname,nmode ,File%iosystem%info,File%fh)
! Set default to NOFILL for performance.  
!   pnetcdf is nofill by default and doesn't support a fill mode
!	  ierr = nfmpi_set_fill(File%fh, NF_NOFILL, nmode)
#endif
#ifdef _NETCDF
       case(iotype_netcdf, iotype_netcdf4c)
          ! Only io proc 0 will do writing
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_create(fname, nmode , File%fh)
! Set default to NOFILL for performance.  
             ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
          call MPI_BCAST(File%fh,1,MPI_INTEGER,0, File%iosystem%IO_comm  , mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
       case(iotype_netcdf4p)
          ierr = nf90_create_par(fname, ior(NF90_NETCDF4,nmode), &
               File%iosystem%IO_comm,File%iosystem%info,File%fh) 
          if(ierr==NF90_NOERR) &
               ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
#endif
       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    end if
    call check_netcdf(File, ierr,_FILE_,__LINE__)

    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(File%fh,1,MPI_INTEGER,File%iosystem%IOMaster, File%iosystem%Comp_comm  , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function create_nf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open_nf
  !

  integer function open_nf(File,fname, mode) result(ierr)

    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in)      :: fname
    integer(i4), optional, intent(in) :: mode
    integer(i4) :: iotype, amode , mpierr, ier2

    ierr=PIO_noerr
    if(file%iosystem%ioproc) then
       iotype = File%iotype 

       select case (iotype) 

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          if(present(mode)) then
             amode = mode
          else
             amode = NF_NOWRITE
          end if
          ierr  = nfmpi_open(File%iosystem%IO_comm,fname,amode,File%iosystem%info,File%fh)
          if(Debug) print *, _FILE_,__LINE__,'CFILE open ',file%fh
#endif

#ifdef _NETCDF
       case(iotype_netcdf, iotype_netcdf4c)

          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             if(present(mode)) then
                if(mode == 1) then
                   amode = NF90_WRITE
                else
                   amode = mode
                end if
             else
                amode = NF90_NOWRITE
             end if
             ierr = nf90_open(fname,amode,File%fh)
             ! Set default to NOFILL for performance.  
             if(ierr .eq. NF90_NOERR .and. iand(amode, NF90_WRITE) > 0) then
                ierr = nf90_set_fill(File%fh, NF90_NOFILL, ier2)
             end if
          endif
          if(File%iosystem%num_iotasks>1) then
             call MPI_BCAST(File%fh,1,MPI_INTEGER,0, File%iosystem%IO_comm  , mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
       case(iotype_netcdf4p)
          


#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    end if

    call check_netcdf(File, ierr,_FILE_,__LINE__)


    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(File%fh,1,MPI_INTEGER,File%iosystem%IOMaster, File%iosystem%Comp_comm  , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function open_nf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! close_nf
  !


  integer function close_nf(File) result(ierr)

    type (File_desc_t), intent(inout) :: File

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       if(Debug) print *,_FILE_,__LINE__,'CFILE closing : ',file%fh
       select case (File%iotype) 
#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_close(file%fh)
#endif
#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr= nf90_sync(File%fh)
             ierr= nf90_close(File%fh)
          endif
#endif
       case default
          call bad_iotype(File%iotype,_FILE_,__LINE__)
       end select
    end if
    call check_netcdf(File, ierr,_FILE_,__LINE__)
  end function close_nf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! sync_nf
  !


  integer function sync_nf(File) result(ierr)

    type (File_desc_t), intent(inout) :: File

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       if(Debug) print *,_FILE_,__LINE__,'CFILE syncing : ',file%fh
       select case (File%iotype) 
#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_sync(file%fh)
#endif
#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr= nf90_sync(File%fh)
          endif
#endif
       case default
          call bad_iotype(File%iotype,_FILE_,__LINE__)
       end select
    end if
    call check_netcdf(File, ierr,_FILE_,__LINE__)
  end function sync_nf

end module ionf_mod
