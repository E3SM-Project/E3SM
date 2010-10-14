#define _FILE_ "ionf_mod.F90"
module ionf_mod
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf      ! _EXTERNAL
#endif
  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types
  use pio_utils, only: bad_iotype, check_netcdf

  use pio_support, only : Debug, DebugIO, piodie, DebugAsync   
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
    integer :: nmode, tmpfh

    nmode=amode
    
    ierr=PIO_noerr
    File%fh=-1

    if(File%iosystem%ioproc) then
       iotype = File%iotype 
       select case (iotype) 
#ifdef _PNETCDF
       case(PIO_iotype_pnetcdf)
          ierr  = nfmpi_create(File%iosystem%IO_comm,fname,nmode ,File%iosystem%info,File%fh)
! Set default to NOFILL for performance.  
!   pnetcdf is nofill by default and doesn't support a fill mode
!	  ierr = nfmpi_set_fill(File%fh, NF_NOFILL, nmode)
#endif
#ifdef _NETCDF
#ifdef _NETCDF4
       case(PIO_iotype_netcdf4p)
!         The 64 bit option is not compatable with hdf5 format files

          if(iand(PIO_64BIT_OFFSET,amode)==PIO_64BIT_OFFSET) then
             nmode = ieor(amode,PIO_64bit_OFFSET)
          else
             nmode=amode
          end if

          nmode = ior(nmode,NF90_NETCDF4)
#ifdef _MPISERIAL
          ierr = nf90_create(fname, nmode , File%fh)
#else
          nmode = ior(nmode,NF90_MPIIO)
          ierr = nf90_create(fname, nmode, File%fh, &
               comm=File%iosystem%io_comm, info=File%iosystem%info)
#endif
! Set default to NOFILL for performance.  
          if(ierr==PIO_NOERR) ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
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
             ierr = nf90_create(fname, amode, File%fh, &
                  info=File%iosystem%info )
! Set default to NOFILL for performance.  
             if(ierr==PIO_NOERR) &
                  ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
#endif          
       case(PIO_iotype_netcdf)
          ! Only io proc 0 will do writing
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_create(fname, nmode , File%fh)
             if(Debug .or. Debugasync) print *,__FILE__,__LINE__,file%fh, ierr
! Set default to NOFILL for performance.  
             if(ierr==NF90_NOERR) &
                  ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
          call mpi_bcast(file%fh,1,mpi_integer, 0, file%iosystem%io_comm, mpierr)
          if(File%iosystem%io_rank > 0) File%fh=-File%fh
#endif
       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
       if(Debug) print *,__FILE__,__LINE__,file%fh,ierr
    end if
    tmpfh = file%fh
    if(Debug.or.DebugAsync) print *,__FILE__,__LINE__,file%fh,ierr
    
    call mpi_bcast(tmpfh,1,mpi_integer, file%iosystem%iomaster, file%iosystem%my_comm, mpierr)

    if(.not. file%iosystem%ioproc) file%fh=-tmpfh
    
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
    integer :: tmpfh


    ierr=PIO_noerr
    File%fh=-1
    if(file%iosystem%ioproc) then
!       This subroutine seems to break pgi compiler for large files.
!       call check_file_type(File, fname)
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
          if(present(mode)) then
             amode = mode
          else
             amode = NF_NOWRITE
          end if
          ierr  = nfmpi_open(File%iosystem%IO_comm,fname,amode,File%iosystem%info,File%fh)

#ifdef _NETCDF
#ifdef _NETCDF4
          if(ierr /= PIO_NOERR) then    ! try hdf5 format
             if(Debug) print *, 'try netcdf4 format'
             File%iotype = pio_iotype_netcdf4p
             iotype = pio_iotype_netcdf4p
          end if
#endif
#endif
       end if
#endif

#ifdef _NETCDF
#ifdef _NETCDF4
        if(iotype==PIO_iotype_netcdf4p .or. iotype ==pio_iotype_netcdf4c) then
#ifdef _MPISERIAL
           ierr = nf90_open(fname,amode,File%fh)           
#else
           ierr = nf90_open(fname,  ior(amode,ior(NF90_NETCDF4,NF90_MPIIO)), File%fh, &
                comm=File%iosystem%io_comm, info=File%iosystem%info)
           if(ierr==nf90_enotnc4 .or. ierr==nf90_einval) then
              ierr = nf90_open(fname, amode, File%fh,info=File%iosystem%info)
              print *,__FILE__,__LINE__,ierr
           end if
#endif
        end if
#endif

       if(iotype==PIO_iotype_netcdf) then
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
             ierr = nf90_open(fname,amode,File%fh)
             if(Debug .or. Debugasync) print *,__FILE__,__LINE__,file%fh, ierr
             ! Set default to NOFILL for performance.  
             if(ierr .eq. NF90_NOERR .and. iand(amode, NF90_WRITE) > 0) then
                ierr = nf90_set_fill(File%fh, NF90_NOFILL, ier2)
             end if
          endif
          call mpi_bcast(file%fh,1,mpi_integer, 0, file%iosystem%io_comm, mpierr)
          if(File%iosystem%io_rank > 0) File%fh=-File%fh
       end if
#endif
    end if

    tmpfh = file%fh
    call mpi_bcast(tmpfh,1,mpi_integer, file%iosystem%iomaster, file%iosystem%my_comm, mpierr)

    if(.not. file%iosystem%ioproc) file%fh=-tmpfh
      
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
             if(ierr==PIO_NOERR) then
                ierr= nf90_close(File%fh)
             else
                if(Debug) print *,__FILE__,__LINE__,ierr
             end if
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
!          print *,__FILE__,__LINE__, magic
          if(magic(1:3) .eq. 'CDF') then
             ! No need to do anything here
          else if(magic(2:4).eq.'HDF') then
#ifdef _NETCDF4
             if(File%iotype /= PIO_IOTYPE_NETCDF4C .and. &
                  File%iotype /= PIO_IOTYPE_NETCDF4P) then
                print *,'Changing file type to netcdf4p'
                File%iotype=pio_iotype_netcdf4c
             end if
#else
             call piodie(__FILE__,__LINE__,'You must link with the netcdf4 ',0,&
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
                   File%iotype=pio_iotype_netcdf4p
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
