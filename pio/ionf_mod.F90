#define _FILE_ "ionf_mod.F90"
module ionf_mod
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf      ! _EXTERNAL
#endif
  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types
#ifndef _PNETCDF
  use pio_utils, only: bad_iotype, check_netcdf
#else
  use pio_utils, only: bad_iotype, check_netcdf, pnetcdf_version_check
#endif
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

    
    nmode=amode
    if(Debug) print *, 'Create file ', fname, 'with amode ',nmode
    ierr=PIO_noerr
    File%fh=-1

    if(File%iosystem%ioproc) then
       iotype = File%iotype 
       select case (iotype) 
#ifdef _PNETCDF
       case(PIO_iotype_pnetcdf)
          call pnetcdf_version_check()
          ierr  = nfmpi_create(File%iosystem%IO_comm,fname,nmode ,File%iosystem%info,File%fh)
! Set default to NOFILL for performance.  
!   pnetcdf is nofill by default and doesn't support a fill mode
!	  ierr = nfmpi_set_fill(File%fh, NF_NOFILL, nmode)
#endif
#ifdef _NETCDF
#ifdef _NETCDF4
       case(PIO_iotype_netcdf4p)
          if(iand(PIO_64BIT_OFFSET,amode)==PIO_64BIT_OFFSET) then
             nmode = ieor(amode,PIO_64bit_OFFSET)
          else
             nmode=amode
          end if
          nmode = ior(nmode,NF90_NETCDF4)
#ifdef _MPISERIAL
          ierr = nf90_create(fname, nmode , File%fh)
#else
          ierr = nf90_create_par(fname, nmode, File%iosystem%io_comm, File%iosystem%info, File%fh)
#endif
! Set default to NOFILL for performance.  
          if(ierr==NF90_NOERR) &
               ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          call check_netcdf(File, ierr,_FILE_,__LINE__)
       case(PIO_iotype_netcdf4c)
          if(iand(PIO_64BIT_OFFSET,amode)==PIO_64BIT_OFFSET) then
             nmode = ieor(amode,PIO_64bit_OFFSET)
          else
             nmode=amode
          end if
          nmode = ior(nmode,NF90_NETCDF4)

          ! Only io proc 0 will do writing
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_create(fname, nmode , File%fh)
! Set default to NOFILL for performance.  
             if(ierr==NF90_NOERR) &
                  ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
#endif          
       case(PIO_iotype_netcdf)
          ! Only io proc 0 will do writing
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_create(fname, nmode , File%fh)
! Set default to NOFILL for performance.  
             if(ierr==NF90_NOERR) &
                  ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
#endif
       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    end if
    call check_netcdf(File, ierr,_FILE_,__LINE__)

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
    File%fh=-1
    if(file%iosystem%ioproc) then
       call check_file_type(File, fname)
       iotype = File%iotype 
#ifdef _NETCDF
       if(present(mode)) then
          if(mode == 1) then
             amode = NF90_WRITE
          else
             amode = mode
          end if
       else
          amode = NF90_NOWRITE
       end if
#endif
#ifdef _PNETCDF
       if(iotype==PIO_iotype_pnetcdf) then
          call pnetcdf_version_check()
          if(present(mode)) then
             amode = mode
          else
             amode = NF_NOWRITE
          end if
          ierr  = nfmpi_open(File%iosystem%IO_comm,fname,amode,File%iosystem%info,File%fh)

#ifdef _NETCDF
#ifdef _NETCDF4
          if(ierr /= PIO_NOERR) then    ! try hdf5 format
             File%iotype = pio_iotype_netcdf4c
             iotype = pio_iotype_netcdf4c
          end if
#endif
#endif
          if(Debug) print *, _FILE_,__LINE__,'CFILE open ',file%fh
       end if
#endif

#ifdef _NETCDF
#ifdef _NETCDF4
        if(iotype==PIO_iotype_netcdf4p .or. iotype ==pio_iotype_netcdf4c) then
#ifdef _MPISERIAL
           ierr = nf90_open(fname,amode,File%fh)           
#else
           ierr = nf90_open_par(fname, amode, File%iosystem%io_comm, File%iosystem%info, File%fh)
#endif
        end if
#endif

       if(iotype==PIO_iotype_netcdf) then
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_open(fname,amode,File%fh)
             ! Set default to NOFILL for performance.  
             if(ierr .eq. NF90_NOERR .and. iand(amode, NF90_WRITE) > 0) then
                ierr = nf90_set_fill(File%fh, NF90_NOFILL, ier2)
             end if
          endif
#endif
       end if
    end if

    call check_netcdf(File, ierr,_FILE_,__LINE__)

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
       case(PIO_iotype_pnetcdf)
          ierr=nfmpi_close(file%fh)
#endif
#ifdef _NETCDF
       case(PIO_iotype_netcdf, pio_iotype_netcdf4c, pio_iotype_netcdf4p)
          if (File%fh>0) then
             ierr= nf90_sync(File%fh)
             ierr= nf90_close(File%fh)
          endif
#endif
       case default
          call bad_iotype(File%iotype,_FILE_,__LINE__)
       end select
    end if
    file%fh=-1
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
       case(PIO_iotype_pnetcdf)
          ierr=nfmpi_sync(file%fh)
#endif
#ifdef _NETCDF
       case(PIO_iotype_netcdf)
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

  subroutine check_file_type(File, filename) 
    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in) :: filename
    character(len=4) :: magic
    integer :: fh, mpierr, reclength=4, i, eof 
    logical :: UNITOK, UNITOP

!   Check format of existing files opened to read.

    inquire(file=filename, exist=UNITOK) 
    if(.not. UNITOK) return

    magic='fail'
   
    if(File%iosystem%ioproc) then      
       if(File%iosystem%io_rank==0) then
!  Find a unique unit number to open the file
          do fh=12,99
             inquire (unit=fh,exist=UNITOK,opened=UNITOP)
             if (UNITOK .and. .not. UNITOP) then
                open (unit = fh,File=filename,access='direct',recl=reclength,&
                     FORM='UNFORMATTED',STATUS='OLD',err=100)
! Read the first 4 bytes and look for the CDF or HDF stamp
                read (fh,rec=1,err=101) magic

                close(fh)
                exit
             endif
          end do
          if(magic(1:3) .eq. 'CDF') then
             ! No need to do anything here
          else if(magic(2:4).eq.'HDF') then
#ifdef _NETCDF4
             if(File%iotype /= PIO_IOTYPE_NETCDF4C .and. &
                  File%iotype /= PIO_IOTYPE_NETCDF4P) then
                if(debug) print *,'Changing file type to netcdf4p'
                File%iotype=pio_iotype_netcdf4c
             end if
#else
             call piodie(__FILE__,__LINE__,'You must link with the netcdf4 '\\&
                  'library built with hdf5 support to read this file',0,filename)
#endif       
          else 
             ! The HDF identifier could be offset further into the file.
             
             open (unit = fh,file=filename,access='direct',recl=reclength,&
                  form='UNFORMATTED',STATUS='OLD',err=100)

             i=128
             eof=0
             do while(eof>=0)
                read (fh,rec=i, iostat=eof, err=101) magic

                if(magic(2:4).eq.'HDF') then
                   if(debug) print *,'Changing file type to netcdf4p'
                   File%iotype=pio_iotype_netcdf4c
                   exit
                end if
                i=i*2
             end do
             close(fh)
             if(eof<0) call piodie(__FILE__,__LINE__,'Unrecognized file format ',0,filename)             
          end if

       end if
       
       call mpi_bcast(file%iotype,1,mpi_integer, 0, file%iosystem%io_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
    return
100 call piodie(__FILE__,__LINE__,'File open error ',0,filename)
101 call piodie(__FILE__,__LINE__,'File read error ',0,filename)



  end subroutine check_file_type




end module ionf_mod
