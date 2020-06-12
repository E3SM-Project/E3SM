
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Dumpy McDump Face (DMDF)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A convenience tool for dumping data to NetCDF files and reading it.
!! - The tool assumes you'll dump records of data, and every dmdf_write()
!!   routine except the attribute routines automatically appends an
!!   unlimited dimension to the end of the dimension list.
!! - Currently (and probably permenantly) you have to write one file per
!!   MPI task. It's not ideal, but coordinating the same file among multiple
!!   tasks when you don't know who's going to write data when is does not
!!   appear to be something NetCDF was designed to do.
!!   - This should be ameliorated by the fact that at scale, we'll always
!!     be sampling a small proportion of a random uniform distribution to
!!     determine when to write data. Thus, there shouldn't be a massive
!!     number of writes at the same time at least.
!! - There is "logical" data type functionality, which is not natively
!!   supported by NetCDF. I hack it by storing to and from an integer type.
!!   Obviously, once written, there's no way to know whether it was supposed
!!   to be integer or logical, so the user keeps up with that.
!! - There is the ability to keep a file open for multiple writes at a time
!!   since in most cases, the inputs and outputs for a routine involve more
!!   than a single piece of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PUBLIC INTERFACES:
!!
!!   dmdf_write_attr(val,rank,fprefix,aname)
!!   dmdf_read_attr(val,fprefix,aname)
!!   dmdf_write(dat,rank,fprefix,vname       ,first,last) !For scalar values
!!   dmdf_write(dat,rank,fprefix,vname,dnames,first,last) !For array values
!!   dmdf_read(dat,fprefix,vname,ind1,ind2,first,last)
!!
!! PARAMETERS:
!! - first == .true. means open the file
!! - last  == .true. means close the file
!! - dnames are the dimension names in a FORTRAN (/array/)
!! - vname is the variable name
!! - aname is the attribute name
!! - rank is an arbitrary integer ID (which is usually MPI rank)
!! - val and dat are the data
!! - fprefix is the previx of the filename, which becomes: prefix_rank.nc
!!   - The file names use preprended zeros when printing the rank.
!! - ind1 and ind2 are the bounds in terms of the unlimited dimension (i.e.,
!!   (ind2-ind1+1) is the total number of samples being read in.
!!
!! NOTES:
!! - All arrays are assumed shape arrays, the dimensions of which are
!!   gathered by shape()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Typical workflow:
!! - Dump data to file using write routines (one file per MPI rank)
!! - Combine files with 'ncrcat prefix_*.nc prefix.nc' (NCO tool)
!! - Read from the concatenated file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Known issues:
!! - You cannot define a variable and then use a dimension by the same name
!! - You need to know how many records to read ahead of time for now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define _NCERR(x) if ((x) /= NF90_NOERR) then; write(error_string,fmt='(A)') nf90_strerror(x); success = .false.; return; endif
#define _ERR(x) write(error_string,fmt='(A)') x; success = .false.; return
#define _RET_IF_ERR if (.not. success) return


module dmdf
  use netcdf
  implicit none
  private

  integer :: ncid

  character(len=1024), public :: error_string
  logical            , public :: success
  
  interface dmdf_write_attr
    module procedure dmdf_write_attr_real4
    module procedure dmdf_write_attr_real8
    module procedure dmdf_write_attr_int4
    module procedure dmdf_write_attr_int8
    module procedure dmdf_write_attr_char
    module procedure dmdf_write_attr_log
  end interface
  
  interface dmdf_read_attr
    module procedure dmdf_read_attr_real4
    module procedure dmdf_read_attr_real8
    module procedure dmdf_read_attr_int4
    module procedure dmdf_read_attr_int8
    module procedure dmdf_read_attr_char
    module procedure dmdf_read_attr_log
  end interface
  
  interface dmdf_write
    module procedure dmdf_write_real4_scalar
    module procedure dmdf_write_real8_scalar
    module procedure dmdf_write_int4_scalar
    module procedure dmdf_write_int8_scalar
    module procedure dmdf_write_log_scalar

    module procedure dmdf_write_real4_1d
    module procedure dmdf_write_real8_1d
    module procedure dmdf_write_int4_1d
    module procedure dmdf_write_int8_1d
    module procedure dmdf_write_log_1d

    module procedure dmdf_write_real4_2d
    module procedure dmdf_write_real8_2d
    module procedure dmdf_write_int4_2d
    module procedure dmdf_write_int8_2d
    module procedure dmdf_write_log_2d

    module procedure dmdf_write_real4_3d
    module procedure dmdf_write_real8_3d
    module procedure dmdf_write_int4_3d
    module procedure dmdf_write_int8_3d
    module procedure dmdf_write_log_3d

    module procedure dmdf_write_real4_4d
    module procedure dmdf_write_real8_4d
    module procedure dmdf_write_int4_4d
    module procedure dmdf_write_int8_4d
    module procedure dmdf_write_log_4d
  end interface
  
  interface dmdf_read
    module procedure dmdf_read_real4_scalar
    module procedure dmdf_read_real8_scalar
    module procedure dmdf_read_int4_scalar
    module procedure dmdf_read_int8_scalar
    module procedure dmdf_read_log_scalar

    module procedure dmdf_read_real4_1d
    module procedure dmdf_read_real8_1d
    module procedure dmdf_read_int4_1d
    module procedure dmdf_read_int8_1d
    module procedure dmdf_read_log_1d

    module procedure dmdf_read_real4_2d
    module procedure dmdf_read_real8_2d
    module procedure dmdf_read_int4_2d
    module procedure dmdf_read_int8_2d
    module procedure dmdf_read_log_2d

    module procedure dmdf_read_real4_3d
    module procedure dmdf_read_real8_3d
    module procedure dmdf_read_int4_3d
    module procedure dmdf_read_int8_3d
    module procedure dmdf_read_log_3d

    module procedure dmdf_read_real4_4d
    module procedure dmdf_read_real8_4d
    module procedure dmdf_read_int4_4d
    module procedure dmdf_read_int8_4d
    module procedure dmdf_read_log_4d
  end interface

  !dmdf_write_attr(val,rank,fprefix,aname)
  public :: dmdf_write_attr

  !dmdf_read_attr(val,fprefix,aname)
  public :: dmdf_read_attr

  !dmdf_write(dat,rank,fprefix,vname       ,first,last)   !For scalar values
  !dmdf_write(dat,rank,fprefix,vname,dnames,first,last)   !For array values
  public :: dmdf_write

  !dmdf_read(dat,fprefix,vname,ind1,ind2,first,last) !For multiple indices
  !dmdf_read(dat,fprefix,vname,ind      ,first,last) !For a single index
  public :: dmdf_read

  !dmdf_num_records(fprefix,num_records) !For a single index
  public :: dmdf_num_records


contains


  subroutine dmdf_num_records(fprefix,num_records)
    implicit none
    character(len=*), intent(in   ) :: fprefix
    integer         , intent(  out) :: num_records
    integer :: ierr, dimid
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)                            ; _RET_IF_ERR
    _NCERR( nf90_inq_dimid( ncid , 'unlim', dimid ) )
    _NCERR( nf90_inquire_dimension( ncid , dimid , len=num_records ) )
    call close_file(.true.,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_num_records


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! GLOBAL READ ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dmdf_read_real4_scalar(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 1
    real(4)         , intent(  out) :: dat(:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real4_scalar


  subroutine dmdf_read_real8_scalar(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 1
    real(8)         , intent(  out) :: dat(:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real8_scalar


  subroutine dmdf_read_int4_scalar(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 1
    integer(4)      , intent(  out) :: dat(:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int4_scalar


  subroutine dmdf_read_int8_scalar(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 1
    integer(8)      , intent(  out) :: dat(:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int8_scalar


  subroutine dmdf_read_log_scalar(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 1
    logical         , intent(  out) :: dat(:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr, dat_int(ind2-ind1+1)
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat_int,start,count) )
    dat = merge(.true.,.false.,dat_int == 1)
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_log_scalar


  subroutine dmdf_read_real4_1d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 2
    real(4)         , intent(  out) :: dat(:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real4_1d


  subroutine dmdf_read_real8_1d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 2
    real(8)         , intent(  out) :: dat(:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real8_1d


  subroutine dmdf_read_int4_1d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 2
    integer(4)      , intent(  out) :: dat(:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int4_1d


  subroutine dmdf_read_int8_1d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 2
    integer(8)      , intent(  out) :: dat(:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int8_1d


  subroutine dmdf_read_log_1d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 2
    logical         , intent(  out) :: dat(:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr, shp(ndims)
    integer, allocatable :: dat_int(:,:)
    success = .true.
    shp = shape(dat)
    allocate(dat_int(shp(1),shp(2)))
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat_int,start,count) )
    dat = merge(.true.,.false.,dat_int == 1)
    deallocate(dat_int)
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_log_1d


  subroutine dmdf_read_real4_2d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 3
    real(4)         , intent(  out) :: dat(:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real4_2d


  subroutine dmdf_read_real8_2d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 3
    real(8)         , intent(  out) :: dat(:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real8_2d


  subroutine dmdf_read_int4_2d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 3
    integer(4)      , intent(  out) :: dat(:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int4_2d


  subroutine dmdf_read_int8_2d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 3
    integer(8)      , intent(  out) :: dat(:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int8_2d


  subroutine dmdf_read_log_2d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 3
    logical         , intent(  out) :: dat(:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr, shp(ndims)
    integer, allocatable :: dat_int(:,:,:)
    success = .true.
    shp = shape(dat)
    allocate(dat_int(shp(1),shp(2),shp(3)))
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat_int,start,count) )
    dat = merge(.true.,.false.,dat_int == 1)
    deallocate(dat_int)
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_log_2d


  subroutine dmdf_read_real4_3d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 4
    real(4)         , intent(  out) :: dat(:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real4_3d


  subroutine dmdf_read_real8_3d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 4
    real(8)         , intent(  out) :: dat(:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real8_3d


  subroutine dmdf_read_int4_3d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 4
    integer(4)      , intent(  out) :: dat(:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int4_3d


  subroutine dmdf_read_int8_3d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 4
    integer(8)      , intent(  out) :: dat(:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int8_3d


  subroutine dmdf_read_log_3d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 4
    logical         , intent(  out) :: dat(:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr, shp(ndims)
    integer, allocatable :: dat_int(:,:,:,:)
    success = .true.
    shp = shape(dat)
    allocate(dat_int(shp(1),shp(2),shp(3),shp(4)))
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat_int,start,count) )
    dat = merge(.true.,.false.,dat_int == 1)
    deallocate(dat_int)
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_log_3d


  subroutine dmdf_read_real4_4d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 5
    real(4)         , intent(  out) :: dat(:,:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real4_4d


  subroutine dmdf_read_real8_4d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 5
    real(8)         , intent(  out) :: dat(:,:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_real8_4d


  subroutine dmdf_read_int4_4d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 5
    integer(4)      , intent(  out) :: dat(:,:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int4_4d


  subroutine dmdf_read_int8_4d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 5
    integer(8)      , intent(  out) :: dat(:,:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_int8_4d


  subroutine dmdf_read_log_4d(dat,fprefix,vname,ind1,ind2,first,last)
    implicit none
    integer, parameter :: ndims = 5
    logical         , intent(  out) :: dat(:,:,:,:,:)
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    logical         , intent(in   ) :: first
    logical         , intent(in   ) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr, shp(ndims)
    integer, allocatable :: dat_int(:,:,:,:,:)
    success = .true.
    shp = shape(dat)
    allocate(dat_int(shp(1),shp(2),shp(3),shp(4),shp(5)))
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,-1,fprefix,ncid)                                      ; _RET_IF_ERR
    call read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)                 ; _RET_IF_ERR
    _NCERR( nf90_get_var(ncid,varid,dat_int,start,count) )
    dat = merge(.true.,.false.,dat_int == 1)
    deallocate(dat_int)
    call close_file(last,ncid)                                                      ; _RET_IF_ERR
  end subroutine dmdf_read_log_4d


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! GLOBAL ATTRIBUTE READ ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dmdf_read_attr_real4(val,fprefix,aname)
    implicit none
    real(4)         , intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_real4


  subroutine dmdf_read_attr_real8(val,fprefix,aname)
    implicit none
    real(8)         , intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_real8


  subroutine dmdf_read_attr_int4(val,fprefix,aname)
    implicit none
    integer(4)      , intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_int4


  subroutine dmdf_read_attr_int8(val,fprefix,aname)
    implicit none
    integer(8)      , intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_int8


  subroutine dmdf_read_attr_char(val,fprefix,aname)
    implicit none
    character(len=*), intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_char


  subroutine dmdf_read_attr_log(val,fprefix,aname)
    implicit none
    logical         , intent(  out) :: val
    character(len=*), intent(in   ) :: fprefix
    character(len=*), intent(in   ) :: aname
    integer :: ierr, ival
    success = .true.
    call procure_fileid(.true.,-1,fprefix,ncid)               ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_get_att(ncid,NF90_GLOBAL,trim(aname),ival) )
    val = merge(.true.,.false.,ival==1)
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_read_attr_log


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! GLOBAL ATTRIBUTE WRITE ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dmdf_write_attr_real4(val,rank,fprefix,aname)
    implicit none
    real(4)         , intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_real4


  subroutine dmdf_write_attr_real8(val,rank,fprefix,aname)
    implicit none
    real(8)         , intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_real8


  subroutine dmdf_write_attr_int4(val,rank,fprefix,aname)
    implicit none
    integer(4)      , intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_int4


  subroutine dmdf_write_attr_int8(val,rank,fprefix,aname)
    implicit none
    integer(8)      , intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_int8


  subroutine dmdf_write_attr_char(val,rank,fprefix,aname)
    implicit none
    character(len=*), intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),val) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_char


  subroutine dmdf_write_attr_log(val,rank,fprefix,aname)
    implicit none
    logical         , intent(in) :: val
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: aname
    integer :: ierr
    success = .true.
    call procure_fileid(.true.,rank,fprefix,ncid)             ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    _NCERR( nf90_put_att(ncid,NF90_GLOBAL,trim(aname),merge(1,0,val)) )
    call close_file(.true.,ncid)                              ; _RET_IF_ERR
  endsubroutine dmdf_write_attr_log


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! SCALAR WRITE ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dmdf_write_real4_scalar(dat,rank,fprefix,vname,first,last)
    implicit none
    integer, parameter :: ndims = 1
    real(4)         , intent(in) :: dat
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)                ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_FLOAT,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,(/dat/),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real4_scalar


  subroutine dmdf_write_real8_scalar(dat,rank,fprefix,vname,first,last)
    implicit none
    integer, parameter :: ndims = 1
    real(8)         , intent(in) :: dat
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)                ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_DOUBLE,varid)        ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,(/dat/),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real8_scalar


  subroutine dmdf_write_int4_scalar(dat,rank,fprefix,vname,first,last)
    implicit none
    integer, parameter :: ndims = 1
    integer(4)      , intent(in) :: dat
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)                ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,(/dat/),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int4_scalar


  subroutine dmdf_write_int8_scalar(dat,rank,fprefix,vname,first,last)
    implicit none
    integer, parameter :: ndims = 1
    integer(8)      , intent(in) :: dat
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)                ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT64,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,(/dat/),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int8_scalar


  subroutine dmdf_write_log_scalar(dat,rank,fprefix,vname,first,last)
    implicit none
    integer, parameter :: ndims = 1
    logical         , intent(in) :: dat
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)                ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,(/merge(1,0,dat)/),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_log_scalar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ARRAY WRITE ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dmdf_write_real4_1d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 2
    real(4)         , intent(in) :: dat(:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_FLOAT,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real4_1d


  subroutine dmdf_write_real8_1d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 2
    real(8)         , intent(in) :: dat(:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_DOUBLE,varid)        ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real8_1d


  subroutine dmdf_write_int4_1d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 2
    integer(4)      , intent(in) :: dat(:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int4_1d


  subroutine dmdf_write_int8_1d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 2
    integer(8)      , intent(in) :: dat(:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT64,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int8_1d


  subroutine dmdf_write_log_1d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 2
    logical         , intent(in) :: dat(:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,merge(1,0,dat),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_log_1d


  subroutine dmdf_write_real4_2d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 3
    real(4)         , intent(in) :: dat(:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_FLOAT,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real4_2d


  subroutine dmdf_write_real8_2d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 3
    real(8)         , intent(in) :: dat(:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_DOUBLE,varid)        ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real8_2d


  subroutine dmdf_write_int4_2d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 3
    integer(4)      , intent(in) :: dat(:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int4_2d


  subroutine dmdf_write_int8_2d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 3
    integer(8)      , intent(in) :: dat(:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT64,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int8_2d


  subroutine dmdf_write_log_2d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 3
    logical         , intent(in) :: dat(:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,merge(1,0,dat),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_log_2d


  subroutine dmdf_write_real4_3d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 4
    real(4)         , intent(in) :: dat(:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_FLOAT,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real4_3d


  subroutine dmdf_write_real8_3d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 4
    real(8)         , intent(in) :: dat(:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_DOUBLE,varid)        ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real8_3d


  subroutine dmdf_write_int4_3d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 4
    integer(4)      , intent(in) :: dat(:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int4_3d


  subroutine dmdf_write_int8_3d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 4
    integer(8)      , intent(in) :: dat(:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT64,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int8_3d


  subroutine dmdf_write_log_3d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 4
    logical         , intent(in) :: dat(:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,merge(1,0,dat),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_log_3d


  subroutine dmdf_write_real4_4d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 5
    real(4)         , intent(in) :: dat(:,:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_FLOAT,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real4_4d


  subroutine dmdf_write_real8_4d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 5
    real(8)         , intent(in) :: dat(:,:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_DOUBLE,varid)        ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_real8_4d


  subroutine dmdf_write_int4_4d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 5
    integer(4)      , intent(in) :: dat(:,:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int4_4d


  subroutine dmdf_write_int8_4d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 5
    integer(8)      , intent(in) :: dat(:,:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT64,varid)         ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,dat,start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_int8_4d


  subroutine dmdf_write_log_4d(dat,rank,fprefix,vname,dnames,first,last)
    implicit none
    integer, parameter :: ndims = 5
    logical         , intent(in) :: dat(:,:,:,:)
    integer         , intent(in) :: rank
    character(len=*), intent(in) :: fprefix
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: dnames(ndims-1)
    logical         , intent(in) :: first
    logical         , intent(in) :: last
    integer :: dimids(ndims), start(ndims), count(ndims), dsizes(ndims)
    integer :: varid, unlim_len, ierr
    success = .true.
    dsizes(1:ndims-1) = shape(dat)
    dsizes(ndims) = NF90_UNLIMITED
    call procure_fileid(first,rank,fprefix,ncid)                         ; _RET_IF_ERR
    ierr = nf90_redef(ncid)
    call procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)       ; _RET_IF_ERR
    call procure_varid(ndims,ncid,vname,dimids,NF90_INT,varid)           ; _RET_IF_ERR
    call compute_start_count(first,ndims,dsizes,unlim_len,start,count)   ; _RET_IF_ERR
    ierr = nf90_enddef(ncid)
    _NCERR( nf90_put_var(ncid,varid,merge(1,0,dat),start,count) )
    call close_file(last,ncid)                                           ; _RET_IF_ERR
  end subroutine dmdf_write_log_4d


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! REUSED INTERNAL ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !Get var ID, then dimids, then dsizes, then fix start,count
  subroutine read_process(ndims,ncid,vname,ind1,ind2,varid,start,count)
    implicit none
    integer         , intent(in   ) :: ndims
    integer         , intent(in   ) :: ncid
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: ind1
    integer         , intent(in   ) :: ind2
    integer         , intent(  out) :: varid
    integer         , intent(  out) :: start (ndims)
    integer         , intent(  out) :: count (ndims)
    integer :: i, ierr, len, dimids(ndims)
    _NCERR( nf90_inq_varid( ncid , trim(vname) , varid ) )
    _NCERR( nf90_inquire_variable( ncid , varid, dimids=dimids ) )
    do i = 1 , ndims
      start(i) = 1
      _NCERR( nf90_inquire_dimension( ncid , dimids(i) , len=count(i) ) )
    enddo
    start(ndims) = ind1
    count(ndims) = ind2 - ind1 + 1
  end subroutine read_process


  subroutine close_file(last,ncid)
    implicit none
    logical, intent(in) :: last
    integer, intent(in) :: ncid
    if (last) then
      _NCERR( nf90_close(ncid) )
    endif
  end subroutine close_file


  subroutine compute_start_count(first,ndims,dsizes,unlim_len,start,count)
    logical, intent(in   ) :: first
    integer, intent(in   ) :: ndims
    integer, intent(in   ) :: dsizes(ndims)
    integer, intent(in   ) :: unlim_len
    integer, intent(  out) :: start (ndims)
    integer, intent(  out) :: count (ndims)
    integer :: i
    do i = 1 , ndims
      if (dsizes(i) == NF90_UNLIMITED) then
        if (first) then
          start(i) = unlim_len+1
        else
          start(i) = unlim_len
        endif
        count(i) = 1
      else
        start(i) = 1
        count(i) = dsizes(i)
      endif
    enddo
  end subroutine compute_start_count


  subroutine procure_varid(ndims,ncid,vname,dimids,xtype,varid)
    implicit none
    integer         , intent(in   ) :: ndims
    integer         , intent(in   ) :: ncid
    character(len=*), intent(in   ) :: vname
    integer         , intent(in   ) :: dimids(ndims)
    integer         , intent(in   ) :: xtype
    integer         , intent(  out) :: varid
    integer :: i, ierr, xtype_file
    integer, allocatable :: dimids_file(:)
    allocate(dimids_file(ndims))
    !This section procures the variable ID, whether by creating it or finding it.
    ierr = nf90_inq_varid( ncid , trim(vname) , varid )
    if     (ierr == NF90_ENOTVAR) then
      !If the variable doesn't exist, then define it
      _NCERR( nf90_def_var( ncid , trim(vname) , xtype , dimids , varid ) )
    elseif (ierr == NF90_NOERR  ) then
      !The variable already exists. Make sure it has the same dimension IDs.
      ierr = nf90_inquire_variable( ncid , varid, xtype=xtype_file , dimids=dimids_file )
      _NCERR( ierr )
      do i = 1 , ndims
        if (dimids(i) /= dimids_file(i)) then
          _ERR( 'Specified variable dimensions differ from file variable dimensions' )
        endif
      enddo
      if (xtype_file /= xtype) then
        _ERR( 'Passed variable type differs from file variable type' )
      endif
    else
      !Different error. Handle it
      _NCERR( ierr )
    endif
    deallocate(dimids_file)
  end subroutine procure_varid


  subroutine procure_dimids(ndims,ncid,dnames,dsizes,dimids,unlim_len)
    implicit none
    integer         , intent(in   ) :: ndims
    integer         , intent(in   ) :: ncid
    character(len=*), intent(in   ) :: dnames(ndims)
    integer         , intent(in   ) :: dsizes(ndims)
    integer         , intent(  out) :: dimids(ndims)
    integer         , intent(  out) :: unlim_len
    integer :: i, ierr
    integer                       :: len
    !The goal of this pass is to procure dimension IDs, whether by creating them or finding them.
    do i = 1 , ndims-1
      !If the dimension is defined already, get the dimension id
      ierr = nf90_inq_dimid( ncid , trim(dnames(i)) , dimids(i) )
      if     (ierr == NF90_EBADDIM) then
        !If the dimension is not defined, then define it
        _NCERR( nf90_def_dim( ncid , trim(dnames(i)) , dsizes(i) , dimids(i) ) )
      elseif (ierr == NF90_NOERR  ) then
        !If the dimension is defined, then make sure the sizes are the same
        _NCERR( nf90_inquire_dimension( ncid , dimids(i) , len=len ) )
        if (len /= dsizes(i)) then
          _ERR( 'Specified dimension size does not match existing dimension size' )
        endif
      else
        !If there's a different error, then handle it
        _NCERR( ierr )
      endif
    enddo

    call procure_dimid_unlim(ndims,ncid,dimids,unlim_len)
  end subroutine procure_dimids


  subroutine procure_dimid_unlim(ndims,ncid,dimids,unlim_len)
    implicit none
    integer         , intent(in   ) :: ndims
    integer         , intent(in   ) :: ncid
    integer         , intent(  out) :: dimids(ndims)
    integer         , intent(  out) :: unlim_len
    integer :: i, ierr, len
    !If the dimension is defined already, get the dimension id
    ierr = nf90_inq_dimid( ncid , 'unlim', dimids(ndims) )
    if     (ierr == NF90_EBADDIM) then
      !If the dimension is not defined, then define it, and set unlim_len to zero
      _NCERR( nf90_def_dim( ncid , 'unlim' , NF90_UNLIMITED , dimids(ndims) ) )
       unlim_len = 0
    elseif (ierr == NF90_NOERR  ) then
      !If the dimension is defined, get the length
      _NCERR( nf90_inquire_dimension( ncid , dimids(ndims) , len=len ) )
      unlim_len = len
    else
      !If there's a different error, then handle it
      _NCERR( ierr )
    endif
  end subroutine procure_dimid_unlim


  subroutine procure_fileid(first,rank,fprefix,ncid)
    implicit none
    logical         , intent(in   ) :: first
    integer         , intent(in   ) :: rank
    character(len=*), intent(in   ) :: fprefix
    integer         , intent(  out) :: ncid
    logical :: file_exists
    integer :: groupnum
    character(len=256) :: fname
    if (first) then
      if (rank == -1) then
        write(fname,fmt='(A,A)') trim(fprefix)
      else
        write(fname,fmt='(A,A,I0.6,A)') trim(fprefix) , '_' , rank , '.nc'
      endif

      !If the file exists, open for writing. Otherwise, create
      inquire(file=trim(fname), exist=file_exists) 
      if (file_exists) then
        if (rank == -1) then
          _NCERR( nf90_open  ( trim(fname) , NF90_NOWRITE , ncid ) )
        else
          _NCERR( nf90_open  ( trim(fname) , NF90_WRITE   , ncid ) )
        endif
      else
        _NCERR( nf90_create( trim(fname) , NF90_HDF5  , ncid ) )
      endif
    endif
  end subroutine procure_fileid


end module dmdf

