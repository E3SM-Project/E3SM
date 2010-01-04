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

  !
  !  Attribute functions
  !
  public :: pio_def_var,   &
       pio_inq_attname,    & 
       pio_inq_att,        &
       pio_inq_attlen,     &
       pio_inq_varid,      &
       pio_inq_varname,    &
       pio_inq_vartype,    &
       pio_inq_varndims,   &
       pio_inq_vardimid,   &
       pio_inq_varnatts,   &
#ifdef _PNETCDF
       pnetcdf_version_check, &
#endif
       pio_inquire_variable
!>
!! \defgroup PIO_def_var
!<
  interface pio_def_var
     module procedure &
          def_var_0d, &
          def_var_md
  end interface

!>
!! \defgroup PIO_inq_varid
!<
  interface pio_inq_varid
     module procedure inq_varid_vid, &
          inq_varid_vardesc
  end interface
!>
!! \defgroup PIO_inq_att
!<
  interface pio_inq_att
     module procedure inq_att_vid, &
          inq_att_vardesc
  end interface

!>
!! \defgroup PIO_inq_attlen
!<
  interface pio_inq_attlen
     module procedure inq_attlen_vid, &
          inq_attlen_vardesc
  end interface

!>
!! \defgroup PIO_inq_attname
!<
  interface pio_inq_attname
     module procedure inq_attname_vid, &
          inq_attname_vardesc
  end interface

!>
!! \defgroup PIO_inq_varname
!<
  interface pio_inq_varname
     module procedure inq_varname_vid, inq_varname_vdesc
  end interface

!>
!! \defgroup PIO_inq_varndims
!<
  interface pio_inq_varndims
     module procedure inq_varndims_vid, inq_varndims_vdesc
  end interface

!>
!! \defgroup PIO_inq_varnatts
!<
  interface pio_inq_varnatts
     module procedure inq_varnatts_vid, inq_varnatts_vdesc
  end interface

!>
!! \defgroup PIO_inq_vardimid
!<
  interface pio_inq_vardimid
     module procedure inq_vardimid_vid, inq_vardimid_vdesc
  end interface

!>
!! \defgroup PIO_inq_vartype
!<
  interface pio_inq_vartype
     module procedure inq_vartype_vid, inq_vartype_vdesc
  end interface

!>
!! \defgroup PIO_inquire_variable
!<
  interface pio_inquire_variable
     module procedure inquire_variable_vid, inquire_variable_vdesc
  end interface

!>
!! @defgroup PIO_def_dim
!<
   public :: PIO_def_dim

!>
!! @defgroup PIO_enddef
!<
  public :: PIO_enddef

!>
!! \defgroup PIO_redef
!<
  public :: PIO_redef

!> 
!! \defgroup PIO_inquire
!<
   public :: PIO_inquire

!>
!! \defgroup PIO_inq_dimid
!<
   public :: PIO_inq_dimid

!>
!! \defgroup PIO_inq_dimname
!<
   public :: PIO_inq_dimname

!>
!! \defgroup PIO_inq_dimlen
!<
   public :: PIO_inq_dimlen

!>
!! \defgroup PIO_inquire_dimension
!<
  public :: PIO_inquire_dimension

!> 
!! \defgroup PIO_copy_att
!<
  public :: PIO_copy_att

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
             if(ierr .eq. NF90_NOERR .and. iand(amode, NF90_WRITE) > 0) then
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

!>
!! @public 
!! @ingroup PIO_inquire
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param nDimensions :  Number of dimensions defined for the netcdf file
!! @param nVariables : Number of variables defined for the netcdf file 
!! @param nAttributes : Number of attributes defined for the netcdf file
!! @param unlimitedDimID : the Unlimited dimension ID
!! @retval ierr @copydoc error_return
!>
  integer function pio_inquire(File,nDimensions,nVariables,nAttributes,unlimitedDimID) result(ierr)
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

  end function pio_inquire

!>
!! @public 
!! @ingroup PIO_inq_att
!! @brief Gets information about attributes
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable identifier
!! @param name : Name of the attribute 
!! @param xtype : The type of attribute
!! @param len : The length of the attribute 
!! @retval ierr @copydoc error_return
!>
  integer function inq_att_vid(File,varid,name,xtype,len) result(ierr)


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
  end function inq_att_vid


!>
!! @public 
!! @ingroup PIO_inq_att
!! @brief  Gets information about attributes
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : Name of the attribute 
!! @param xtype : The type of attribute
!! @param len : The length of the attribute 
!! @retval ierr @copydoc error_return
!>
  integer function inq_att_vardesc(File,vardesc,name,xtype,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    type(var_desc_t), intent(in)           :: vardesc
    character(len=*), intent(in)      :: name
    integer, intent(out)              :: xtype
    integer, intent(out)              :: len !Attribute length

    ierr = pio_inq_att(file, vardesc%varid, name, xtype, len)

  end function inq_att_vardesc

!>
!! @public 
!! @ingroup PIO_inq_attlen
!! @brief Gets the attribute length 
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : attribute id
!! @param name : name of attribute
!! @param len : Length of attribute
!! @retval ierr @copydoc error_return
!>
  integer function inq_attlen_vid(File,varid,name,len) result(ierr)

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

  end function inq_attlen_vid

!>
!! @public 
!! @ingroup PIO_inq_attlen
!! @brief  Gets the attribute length 
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : name of attribute
!! @param len : Length of attribute
!! @retval ierr @copydoc error_return
!>
  integer function inq_attlen_vardesc(File,vardesc,name,len) result(ierr)

    type (File_desc_t), intent(inout) :: File
    type (Var_desc_t), intent(in)            :: vardesc
    character(len=*), intent(in)      :: name
    integer, intent(out),optional     :: len !Attribute length

    ierr = pio_inq_attlen(file, vardesc%varid, name, len)

  end function inq_attlen_vardesc

!> 
!! @public 
!! @ingroup PIO_inq_attname
!! @brief Returns the name of a netcdf attribute 
!! @details
!! @param File @copydoc file_desc_t
!! @param varid :  The variable ID 
!! @param attnum : Attribute number returned from function ????
!! @param name   : Name of the returned attribute
!! @retval ierr @copydoc error_return
!<
  integer function inq_attname_vid(File,varid,attnum,name) result(ierr)

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

  end function inq_attname_vid

!> 
!! @public 
!! @ingroup PIO_inq_attname
!! @brief  Returns the name of a netcdf attribute.
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t 
!! @param attnum : Attribute number returned from function ????
!! @param name   : Name of the returned attribute
!! @retval ierr @copydoc error_return
!<
  integer function inq_attname_vardesc(File,vardesc,attnum,name) result(ierr)
    type (File_desc_t), intent(inout) :: File
    type(var_desc_t), intent(in)           :: vardesc
    integer, intent(in)              :: attnum !Attribute number
    character(len=*), intent(out)     :: name

    ierr = pio_inq_attname(file, vardesc%varid, attnum, name)

  end function inq_attname_vardesc

!> 
!! @public 
!! @ingroup PIO_inq_varid
!! @brief  Returns the ID of a netcdf variable given its name 
!! @details
!! @param File @copydoc file_desc_t
!! @param name : Name of the returned attribute
!! @param varid : variable ID
!! @retval ierr @copydoc error_return
!<
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

!> 
!! @public 
!! @ingroup PIO_inq_varid
!! @brief Returns the ID of a netcdf variable given its name 
!! @details
!! @param File @copydoc file_desc_t
!! @param name   : Name of the returned attribute
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function inq_varid_vardesc(File,name,vardesc) result(ierr)

    type (File_desc_t), intent(in)   :: File
    character(len=*), intent(in)     :: name
    type (Var_desc_t), intent(inout) :: vardesc

    ierr = pio_inq_varid(File, name, vardesc%varid)
    vardesc%rec=-1
  end function inq_varid_vardesc

!>
!! @public 
!! @ingroup PIO_inq_varname
!! @brief Get the name associated with a variable
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param name : The name of the netcdf variable.
!! @retval ierr @copydoc error_return
!>
  integer function inq_varname_vdesc(File,vardesc,name) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in)    :: vardesc
    character(len=*), intent(out)    :: name
    
    ierr = pio_inq_varname(file,vardesc%varid,name)

  end function inq_varname_vdesc

!>
!! @public 
!! @ingroup PIO_inq_varname
!! @brief Get the name associated with a variable
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id.
!! @param name : The name of the netcdf variable.
!! @retval ierr @copydoc error_return
!>
  integer function inq_varname_vid(File,varid,name) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
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
          ierr=nfmpi_inq_varname(File%fh,varid,name)

#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,name=name)
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
  end function inq_varname_vid

!>
!! @public 
!! @ingroup PIO_inq_varndims
!! @brief Gets the number of dimension associated with a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The variable identifier
!! @param ndims : The number of dimensions for the variable 
!! @retval ierr @copydoc error_return
!>
  integer function inq_varndims_vid(File,varid,ndims) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
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
          ierr=nfmpi_inq_varndims(File%fh,varid,ndims)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,ndims=ndims)
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
  end function inq_varndims_vid

!>
!! @public 
!! @ingroup PIO_inq_varndims
!! @brief Gets the number of dimension associated with a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param ndims : The number of dimensions for the variable 
!! @retval ierr @copydoc error_return
!>
  integer function inq_varndims_vdesc(File,vardesc,ndims) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: ndims

    ierr = pio_inq_varndims(File, vardesc%varid, ndims)
  end function inq_varndims_vdesc

!>
!! @public 
!! @ingroup PIO_inq_vartype
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id
!! @param type : The type of variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_vartype_vid(File,varid,type) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer, intent(in) :: varid
    integer(i4), intent(out)    :: type


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
          ierr=nfmpi_inq_vartype(File%fh,varid,type)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)

          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,xtype=type)
          endif

          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(type,1,MPI_INTEGER,0,File%iosystem%IO_comm, mpierr)
             call CheckMPIReturn('nf_mod',mpierr)
          end if
#endif

       case default
          call bad_iotype(iotype,_FILE_,__LINE__)

       end select
    endif
    call check_netcdf(File,ierr,_FILE_,__LINE__)
    if(File%iosystem%num_tasks>File%iosystem%num_iotasks) then
       call MPI_BCAST(type,1,MPI_INTEGER,File%iosystem%IOMaster,File%iosystem%Comp_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
  end function inq_vartype_vid

!>
!! @public 
!! @ingroup PIO_inq_vartype
!! @brief Gets metadata information for netcdf file.
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param type : The type of variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_vartype_vdesc(File,vardesc,type) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: type

    ierr = pio_inq_vartype(File, vardesc%varid, type)
  end function inq_vartype_vdesc

!>
!! @public 
!! @ingroup PIO_inq_vardimid
!! @brief returns the dimids of the variable as an interger array
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The variable id
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @retval ierr @copydoc error_return
!>
  integer function inq_vardimid_vid(File,varid,dimids) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer,            intent(in) :: varid
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
          ierr=nfmpi_inq_vardimid(File%fh,varid,dimids)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,dimids=dimids)
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
  end function inq_vardimid_vid

!>
!! @public 
!! @ingroup PIO_inq_vardimid
!! @brief returns the dimids of the variable as an interger array
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @retval ierr @copydoc error_return
!>
  integer function inq_vardimid_vdesc(File,vardesc,dimids) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in) :: vardesc
    integer(i4), intent(out)    :: dimids(:)


    ierr = pio_inq_vardimid(File, vardesc%varid, dimids)
  end function inq_vardimid_vdesc

!>
!! @public 
!! @ingroup PIO_inq_varnatts
!! @brief Returns the number of attributes associated with a varaible
!! @details
!! @param File @copydoc file_desc_t
!! @param varid : The netcdf variable id
!! @param natts : The number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_varnatts_vid(File,varid,natts) result(ierr)

    type (File_desc_t), intent(in)   :: File
    integer           , intent(in) :: varid
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
          ierr=nfmpi_inq_varnatts(File%fh,varid,natts)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)
          if (File%iosystem%io_rank==0) then
             ierr=nf90_inquire_variable(File%fh,varid,nAtts=natts)
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
  end function inq_varnatts_vid

!>
!! @public 
!! @ingroup PIO_inq_varnatts
!! @brief Returns the number of attributes associated with a varaible
!! @details
!! @param File @copydoc file_desc_t
!! @param vardesc @copydoc var_desc_t
!! @param natts : The number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inq_varnatts_vdesc(File,vardesc,natts) result(ierr)

    type (File_desc_t), intent(in)   :: File
    type (Var_desc_t), intent(in)    :: vardesc
    integer(i4), intent(out)         :: natts


    ierr = pio_inq_varnatts(file, vardesc%varid, natts)
  end function inq_varnatts_vdesc

!>
!! @public 
!! @ingroup PIO_inq_dimid
!! @brief Returns the netcdf dimension id for the name.
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the netcdf dimension.
!! @param dimid : The netcdf dimension id.
!! @retval ierr @copydoc error_return
!!
!! Note that we do not want internal error checking for this funtion.
!>
  integer function pio_inq_dimid(File,name,dimid) result(ierr)

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

  end function pio_inq_dimid

!>
!! @public 
!! @ingroup PIO_inq_dimname
!! @brief Gets the name of a dimension given its ID
!! @details
!! @param File @copydoc file_desc_t
!! @param dimid : The netcdf dimension id.
!! @param dimname : The name associated with the netcdf dimension id.
!! @retval ierr @copydoc error_return
!>
  integer function pio_inq_dimname(File,dimid,dimname) result(ierr)

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

  end function pio_inq_dimname

!>
!! @public 
!! @ingroup PIO_inq_dimlen
!! @brief Returns the extent of a netCDF dimension 
!! @details
!! @param File @copydoc file_desc_t
!! @param dimid : The netcdf dimension.
!! @param dimlen : The extent of the netcdf dimension.
!! @retval ierr @copydoc error_return
!>
  integer function pio_inq_dimlen(File,dimid,dimlen) result(ierr)

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

  end function pio_inq_dimlen

!> 
!! @public
!! @ingroup PIO_enddef
!! @brief Exits netcdf define mode.
!! @details
!! @param File @copydoc file_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function PIO_enddef(File) result(ierr)
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
  end function PIO_enddef

!> 
!! @public
!! @ingroup PIO_redef
!! @brief Re-enters netcdf define mode.   
!! @details 
!! @warning Entering and leaving netcdf define mode causes a file sync operation to 
!!          occur, these operations can be very expensive in parallel systems.   We 
!!          recommend structuring your code to minimize calls to this function.
!! @param File @copydoc file_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function PIO_redef(File) result(ierr)
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
  end function PIO_redef

!> 
!! @public
!! @ingroup PIO_def_dim
!! @brief Defines the netcdf dimension
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the dimension to define
!! @param len :  The size of the dimension
!! @param dimid : The dimension identifier
!<
  integer function PIO_def_dim(File,name,len,dimid) result(ierr)

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
  end function PIO_def_dim

!> 
!! @public 
!! @ingroup PIO_def_var
!! @brief Defines a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the variable to define
!! @param type : The type of variable 
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function def_var_0d(File,name,type,vardesc) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    type (Var_desc_t), intent(inout) :: vardesc
    integer :: len=0, dimids(1)

    ierr = def_var_md(File,name,type,dimids(1:len),vardesc)

  end function def_var_0d

!> 
!! @public
!! @ingroup PIO_def_var
!! @brief Defines the a netcdf variable
!! @details
!! @param File @copydoc file_desc_t
!! @param name : The name of the variable to define
!! @param type : The type of variable 
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param vardesc @copydoc var_desc_t
!! @retval ierr @copydoc error_return
!<
  integer function def_var_md(File,name,type,dimids,vardesc) result(ierr)

    type (File_desc_t), intent(in)  :: File
    character(len=*), intent(in)    :: name
    integer, intent(in)             :: type
    integer, intent(in)             :: dimids(:)
    type (Var_desc_t), intent(inout) :: vardesc

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
          ierr=nfmpi_def_var(File%fh,name,type,len,dimids,vardesc%varid)
#endif

#ifdef _NETCDF
       case(iotype_netcdf)



          ! assuming type valid for both pnetcdf and netcdf

          if (File%iosystem%io_rank==0) then
             ierr=nf90_def_var( ncid=File%fh,name=name,xtype=type, &
                  dimids=dimids,varid=vardesc%varid)
             if (Debug) print *, '0: def_var fh=',File%fh, &
                  'name=',name,' id=',vardesc%varid
          endif
          if(File%iosystem%num_tasks==File%iosystem%num_iotasks) then
             call MPI_BCAST(vardesc%varid, 1, MPI_INTEGER, 0, File%iosystem%IO_Comm, ierr)
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
  end function def_var_md

!>
!! @public
!! @ingroup PIO_copy_att
!! @brief No idea what this function does
!! @details 
!! @param infile @copydoc file_desc_t
!! @param invarid :
!! @param name : 
!! @param outfile :
!! @param outvarid :
!! @retval ierr @copydoc error_return
!<
  integer function pio_copy_att(infile, invarid, name, outfile, outvarid) result(ierr)

    type (File_desc_t), intent(in)  :: infile, outfile
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
  end function pio_copy_att


!>
!! @public 
!! @ingroup PIO_inquire_variable
!! @brief Inquires if a NetCDF variable is present and returns its attributes  
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param varid : The netcdf variable ID.
!! @param name : The name of the variable
!! @param xtype : The type of the variable
!! @param ndims : The number of dimensions for the variable.
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param natts : Number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inquire_variable_vid(ncid, varid, name, xtype, ndims, dimids, natts) result(ierr)
    type(file_desc_t), intent(in) :: ncid
    integer,                         intent( in) :: varid
    character (len = *),   optional, intent(out) :: name
    integer,               optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, intent(out) :: dimids
    integer,               optional, intent(out) :: natts

    
    if(present(name)) ierr = pio_inq_varname(ncid, varid, name)
    if(present(ndims)) ierr = pio_inq_varndims(ncid, varid, ndims)
    if(present(dimids)) ierr = pio_inq_vardimid(ncid, varid, dimids)
    if(present(natts)) ierr = pio_inq_varnatts(ncid, varid, natts)
    if(present(xtype)) ierr = pio_inq_vartype(ncid, varid, xtype)



  end function inquire_variable_vid

!>
!! @public 
!! @ingroup PIO_inquire_variable
!! @brief Inquires if a NetCDF variable is present and returns its attributes  
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param vardesc @copydoc var_desc_t
!! @param name : The name of the variable
!! @param xtype : The type of the variable
!! @param ndims : The number of dimensions for the variable.
!! @param dimids : The dimension identifier returned by \ref PIO_def_dim
!! @param natts : Number of attributes associated with the variable
!! @retval ierr @copydoc error_return
!>
  integer function inquire_variable_vdesc(ncid, vardesc, name, xtype, ndims, dimids, natts) result(ierr)
    type(file_desc_t),               intent(in) :: ncid
    type(var_desc_t),                intent( in) :: vardesc
    character (len = *),   optional, intent(out) :: name
    integer,               optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, intent(out) :: dimids
    integer,               optional, intent(out) :: natts

    if(present(name)) ierr = pio_inq_varname(ncid, vardesc, name)
    if(present(ndims)) ierr = pio_inq_varndims(ncid, vardesc, ndims)
    if(present(dimids)) ierr = pio_inq_vardimid(ncid, vardesc, dimids)
    if(present(natts)) ierr = pio_inq_varnatts(ncid, vardesc, natts)
    if(present(xtype)) ierr = pio_inq_vartype(ncid, vardesc, xtype)

  end function inquire_variable_vdesc

!>
!! @public 
!! @ingroup PIO_inquire_dimension
!! @brief  Get information about a particular dimension in netcdf file 
!! @details
!! @param ncid : A netcdf file descriptor returned by \ref PIO_openfile or \ref PIO_createfile.
!! @param dimid : The netcdf dimension ID.
!! @param name : The name of the dimension.
!! @param len : The length of the dimesions name.
!! @retval ierr @copydoc error_return
!>
  integer function PIO_inquire_dimension(ncid, dimid, name, len) result(ierr)
    type(file_desc_T),             intent(in)  :: ncid
    integer,                       intent( in) :: dimid
    character (len = *), optional, intent(out) :: name
    integer,             optional, intent(out) :: len

    if(present(len)) ierr = pio_inq_dimlen(ncid, dimid, len)
    if(present(name)) ierr = pio_inq_dimname(ncid, dimid,name)

  end function PIO_inquire_dimension

#ifdef _PNETCDF  
  subroutine pnetcdf_version_check()
    character(len=80) :: version
    integer :: dot1, dot2, s
    integer :: v1, v2, v3
    integer, parameter :: rv1=1, rv2=1, rv3=0
    character(len=80) :: error

    write(error,*) 'Pnetcdf version appears to be older than minimum required ',rv1,'.',rv2,'.',rv3

    version = nfmpi_inq_libvers()

    dot1 = index(version,'.')
    dot2 = index(version,'.',.true.)
    
    s = index(version(1:dot1),' ',.true.)
    read(version(s:dot1-1),'(i)') v1
    read(version(dot1+1:dot2-1),'(i)') v2
    s = index(version(dot2:),' ')
    read(version(dot2+1:s-1),'(s)') v3

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

  end subroutine pnetcdf_version_check
#endif  

end module nf_mod
