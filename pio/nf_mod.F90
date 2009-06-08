#define _FILE_ "nf_mod.F90"
module nf_mod
#ifdef TIMING
  use perf_mod, only : t_startf, t_stopf      ! _EXTERNAL
#endif
  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types

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


  ! Removed definition of PIO_error (not used)
  ! and PIO_noerror (defined in pio_types)
  ! 

  public :: enddef,    &
       def_dim,   &
       def_var   
  !
  !  Attribute functions
  !
  public :: inquire,        &
       inq_attname,    & 
       inq_att,        &
       inq_attlen,     &
       inq_varid,      &
       inq_varname,    &
       inq_vartype,    &
       inq_varndims,   &
       inq_vardimid,   &
       inq_varnatts,   &
       inq_dimid,      &
       inq_dimname,    &
       inq_dimlen,     &
       redef,          &             
       copy_att

  interface def_var
     module procedure &
          def_var0d, &
          def_varmd
  end interface
          

  interface inq_varid
     module procedure inq_varid_vid, &
          inq_varid_vardesc
  end interface

  public :: check_netcdf, bad_iotype


  public :: create_nf,open_nf,close_nf, sync_nf

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
       case(iotype_netcdf)
          ! Only io proc 0 will do writing
          if (File%iosystem%io_rank == 0) then
             ! Stores the ncid in File%fh
!             ierr = nf__create(fname, amode, 2147483647, 46006272,File%fh )
             ierr = nf90_create(fname, nmode , File%fh)
! Set default to NOFILL for performance.  
             ierr = nf90_set_fill(File%fh, NF90_NOFILL, nmode)
          endif
          call MPI_BCAST(File%fh,1,MPI_INTEGER,0, File%iosystem%IO_comm  , mpierr)
          call CheckMPIReturn('nf_mod',mpierr)

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
       case(iotype_netcdf)

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
             if(iand(amode, NF90_WRITE) > 0) then
                ierr = nf90_set_fill(File%fh, NF90_NOFILL, ier2)
             end if
          endif
          if(File%iosystem%num_iotasks>1) then
             call MPI_BCAST(File%fh,1,MPI_INTEGER,0, File%iosystem%IO_comm  , mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
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

    print *,'Invalid iotype, value=',iotype
    call piodie(file,line,'Quitting')

  end subroutine bad_iotype




  !============================================
  !
  !  Inquire function {netcdf,pnetcdf} only
  !============================================

  integer function inquire(File,nDimensions,nVariables,nAttributes,unlimitedDimID) result(ierr)
    type (File_desc_t), intent(in) :: File
    integer, optional, intent(out) :: &
         nDimensions,  &! number of dimensions
         nVariables,   &! number of variables
         nAttributes,  & ! number of global attributes
         unlimitedDimID ! ID of unlimited dimension
    integer :: vals(4)
    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr
    vals(:) = -1
    call mpi_barrier(file%iosystem%comp_comm, mpierr)
    call CheckMPIReturn('nf_mod',mpierr)

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq( File%fh,vals(1),vals(2), &
               vals(3),vals(4))

#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire( File%fh,vals(1),vals(2), &
                  vals(3),vals(4))
          endif

          if(File%iosystem%num_iotasks>1) then
             call MPI_BCAST(vals,4,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif

    call check_netcdf(File, ierr, _FILE_,__LINE__)

    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(vals,4,MPI_INTEGER,File%iosystem%IOMaster, File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

    if(present(nDimensions)) then	
       ndimensions = vals(1)
    endif
    if(present(nVariables)) then	
       nVariables = vals(2)
    endif
    if(present(nAttributes)) then	
       nAttributes = vals(3)
    endif
    if(present(unlimitedDimID)) then	
       unlimitedDimID = vals(4)
    endif

  end function inquire



  !============================================
  ! inq_att:
  !
  !  Inquire function {netcdf,pnetcdf} only
  !============================================

  integer function inq_att(File,varid,name,xtype,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)           :: varid
    character(len=*), intent(in)      :: name
    integer, intent(out)              :: xtype
    integer, intent(out)              :: len !Attribute length



    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    integer(kind=PIO_Offset) :: clen

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_att(File%fh,varid,name,xtype,clen)

          len = INT(clen,kind=i4)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_attribute( File%fh,varid,name, &
                  xtype=xtype,len=len)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(xtype,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
             call MPI_BCAST(len,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(xtype,1,MPI_INTEGER,File%iosystem%IOMaster, File%iosystem%Comp_comm , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
       call MPI_BCAST(len,1,MPI_INTEGER,File%iosystem%IOMaster, File%iosystem%Comp_comm  , mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_att



  !============================================
  ! inq_attlen:
  !
  !  Inquire function {netcdf,pnetcdf} only
  !============================================

  integer function inq_attlen(File,varid,name,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)            :: varid
    character(len=*), intent(in)      :: name
    integer, intent(out),optional     :: len !Attribute length


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    integer(kind=PIO_Offset) :: clen

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          if(present(len)) then 
             ierr=nfmpi_inq_attlen(File%fh,varid,name,clen)
             len = INT(clen,kind=i4)
          endif

#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_attribute( File%fh,varid,name, &
                  len=len)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(len,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(len,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_attlen



  !============================================
  ! inq_attname:
  !
  !  Inquire function {netcdf,pnetcdf} only
  !============================================

  integer function inq_attname(File,varid,attnum,name) result(ierr)

    type (File_desc_t), intent(inout) :: File
    integer(i4), intent(in)           :: varid
    integer, intent(in)              :: attnum !Attribute number
    character(len=*), intent(out)     :: name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_attname(File%fh,varid,attnum,name)

#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inq_attname(File%fh,varid,attnum,name)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(name,len(name),MPI_CHARACTER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(name,len(name),MPI_CHARACTER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_attname



  !============================================
  ! Inq_varid:
  !
  !  Returns the ID of a netCDF variable give its name
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_varid_vid(File,name,varid) result(ierr)

    type (File_desc_t), intent(in)   :: File
    character(len=*), intent(in)     :: name
    integer(i4), intent(out)       :: varid
    integer :: ierr2

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_varid(File%fh,name,varid)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inq_varid(File%fh,name,varid)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(varid,1,MPI_INTEGER,0,File%iosystem%IO_comm,ierr2)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)

    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(varid,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm,ierr2)
    end if
    call MPI_BCAST(ierr,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm,ierr2)
  end function inq_varid_vid

  integer function inq_varid_vardesc(File,name,varDesc) result(ierr)

    type (File_desc_t), intent(in)   :: File
    character(len=*), intent(in)     :: name
    type (Var_desc_t), intent(inout) :: varDesc


    ierr = inq_varid_vid(File, name, vardesc%varid)
    vardesc%rec=-1
  end function inq_varid_vardesc



  !============================================
  ! Inq_varname:
  !
  !  Returns the name of the variable
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_varname(File,varDesc,name) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(inout) :: varDesc
    character(len=*), intent(out)    :: name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_varname(File%fh,varDesc%varid,name)

#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varDesc%varid,name=name)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(name,len(name),MPI_CHARACTER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(name,len(name),MPI_CHARACTER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_varname



  !============================================
  ! Inq_vartype:
  !
  !  Returns the type of the variable
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_vartype(File,varDesc,xtype) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(inout) :: varDesc
    integer(i4), intent(out)    :: xtype


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_vartype(File%fh,varDesc%varid,xtype)

#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varDesc%varid,xtype=xtype)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(xtype,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr, _FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(xtype,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_vartype



  !============================================
  ! Inq_varndims:
  !
  !  Returns the number of dimensions of the variable
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_varndims(File,varDesc,ndims) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(inout) :: varDesc
    integer(i4), intent(out)    :: ndims


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          ierr=nfmpi_inq_varndims(File%fh,varDesc%varid,ndims)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varDesc%varid,ndims=ndims)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(ndims,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File,ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(ndims,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_varndims



  !============================================
  ! Inq_vardimid:
  !
  !  Returns the dimids of the variable as an integer array
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_vardimid(File,varDesc,dimids) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(inout) :: varDesc
    integer(i4), intent(out)    :: dimids(:)


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_vardimid(File%fh,varDesc%varid,dimids)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varDesc%varid,dimids=dimids)
             if(ierr<0) print *,_FILE_,__LINE__,vardesc%varid
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(dimids,size(dimids),MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File,ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(dimids,size(dimids),MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_vardimid



  !============================================
  ! Inq_varnatts:
  !
  !  Returns the number of attributes associated with a variable
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_varnatts(File,varDesc,natts) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(inout) :: varDesc
    integer(i4), intent(out)         :: natts


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_varnatts(File%fh,varDesc%varid,natts)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varDesc%varid,nAtts=natts)
          endif

          call MPI_BCAST(natts,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
          call CheckMPIReturn('nf_mod',mpierr)
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(natts,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_varnatts



  !============================================
  ! inq_dimid:
  !
  !  Returns the ID of a netCDF variabl give its name
  !   {netcdf,pnetcdf} only
  !  We do not want internal error checking on this function.
  !============================================

  integer function inq_dimid(File,name,dimid) result(ierr)

    type (File_desc_t), intent(in) :: File
    character(len=*), intent(in)   :: name
    integer, intent(out)           :: dimid        !dimension ID


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_dimid(File%fh,name,dimid)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inq_dimid(File%fh,name,dimid)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(dimid,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif

    call check_netcdf(File, ierr,_FILE_,__LINE__)

    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(dimid,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_dimid



  !============================================
  ! inq_dimname:
  !
  !  Returns the dimname of a netCDF dimension given its dimid
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_dimname(File,dimid,dimname) result(ierr)

    type (File_desc_t), intent(in) :: File
    integer         , intent(in)   :: dimid
    character(len=*), intent(out)  :: dimname        !dimension name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)
          ierr=nfmpi_inq_dimname(File%fh,dimid,dimname)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_dimension(File%fh,dimid,name=dimname)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(dimname,len(dimname),MPI_CHARACTER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(dimname,len(dimname),MPI_CHARACTER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_dimname



  !============================================
  ! inq_dimlen:
  !
  !  Returns the dimlen of a netCDF dimension given its dimid
  !   {netcdf,pnetcdf} only
  !============================================

  integer function inq_dimlen(File,dimid,dimlen) result(ierr)

    type (File_desc_t), intent(in) :: File
    integer(i4)     , intent(in)   :: dimid
    integer(i4)     , intent(out)  :: dimlen        !dimension name


    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    integer(kind=PIO_OFFSET) :: clen

    iotype = File%iotype
    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          ierr=nfmpi_inq_dimlen(File%fh,dimid,clen)
          dimlen = INT(clen,kind=i4)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)


          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_dimension(File%fh,dimid,len=dimlen)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(dimlen,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(dimlen,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if

  end function inq_dimlen



  !============================================
  ! EndDef:
  !
  !  End definition mode {netcdf,pnetcdf} only
  !============================================

  integer function EndDef(File) result(ierr)
    type (File_desc_t), intent(inout) :: File

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    logical, parameter :: Check = .TRUE.

    iotype = File%iotype

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          ierr=nfmpi_enddef(File%fh)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)



          if (File%iosystem%io_rank==0) then
             ierr=nf90_enddef(File%fh)
          endif

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
  end function Enddef


  !============================================
  ! EndDef:
  !
  !  End definition mode {netcdf,pnetcdf} only
  !============================================

  integer function ReDef(File) result(ierr)
    type (File_desc_t), intent(inout) :: File

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    logical, parameter :: Check = .TRUE.

    iotype = File%iotype

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          ierr=nfmpi_redef(File%fh)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)



          if (File%iosystem%io_rank==0) then
             ierr=nf90_redef(File%fh)
          endif

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
  end function ReDef



  !============================================
  ! def_dim:
  !
  !   Defines the dimension of a variable
  !   {netcdf,pnetcdf} only
  !============================================

  integer function def_dim(File,name,len,dimid) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer(i4), intent(in)         :: len
    integer(i4), intent(out)        :: dimid

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr
    integer(kind=PIO_Offset)  :: clen

    iotype = File%iotype

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          clen = len
          ierr=nfmpi_def_dim(File%fh,name,clen,dimid)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)



          if (File%iosystem%io_rank==0) then
             ierr=nf90_def_dim(ncid=File%fh,name=name,len=len,dimid=dimid)
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(dimid, 1, MPI_INTEGER, 0, File%iosystem%IO_Comm, ierr)
          end if

#endif
       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(dimid, 1, MPI_INTEGER, File%iosystem%IOMaster, File%iosystem%Comp_Comm, ierr)
    end if
  end function def_dim



  !============================================
  ! def_var:
  !
  !   Defines the dimension of a variable
  !   {netcdf,pnetcdf} only
  !============================================
  integer function def_var0d(File,name,type,varDesc) result(ierr)
    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    type (Var_desc_t), intent(inout) :: varDesc
    integer :: len=0, dimids(1)

    ierr = def_varmd(File,name,type,dimids(1:len),vardesc)

  end function def_var0d

  integer function def_varmd(File,name,type,dimids,varDesc) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    integer, intent(in)             :: dimids(:)
    type (Var_desc_t), intent(inout) :: varDesc

    !------------------
    ! Local variables
    !------------------
    integer :: iotype, mpierr,len

    iotype = File%iotype

    ierr=PIO_noerr
    vardesc%rec=-1

    vardesc%type = type

    if(File%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          len = SIZE(dimids)
          ierr=nfmpi_def_var(File%fh,name,type,len,dimids,varDesc%varid)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)



          ! assuming type valid for both pnetcdf and netcdf

          if (File%iosystem%io_rank==0) then
             ierr=nf90_def_var( ncid=File%fh,name=name,xtype=type, &
                  dimids=dimids,varid=varDesc%varid)
             if (Debug) print *, '0: def_var fh=',File%fh, &
                  'name=',name,' id=',varDesc%varid
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(varDesc%varid, 1, MPI_INTEGER, 0, File%iosystem%IO_Comm, ierr)
          end if

#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File, ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(vardesc%varid, 1, MPI_INTEGER, File%iosystem%IOMaster, File%iosystem%Comp_Comm, ierr)
    end if
  end function def_varmd

  integer function copy_att(infile, invarid, name, outfile, outvarid) result(ierr)

    type (File_desc_t), intent(in)  :: inFile, outfile
    character(len=*), intent(in)    :: name
    integer, intent(in) :: invarid, outvarid
    integer :: iotype, mpierr

    ierr=PIO_noerr
    iotype = infile%iotype
    if(InFile%iosystem%IOproc) then
       select case(iotype)

#ifdef _PNETCDF
       case(iotype_pnetcdf)

          ierr = nfmpi_copy_att(infile%fh, invarid, name, &
               outfile%fh, outvarid)
#endif
#ifdef _NETCDF
       case(iotype_netcdf)


          if (inFile%iosystem%io_rank==0) then
             ierr = nf90_copy_att(infile%fh,invarid,name,&
                  outfile%fh,outvarid)     
          end if
#endif    
       end select
    end if
    call check_netcdf(outFile, ierr,_FILE_,__LINE__)
  end function copy_att



end module nf_mod
