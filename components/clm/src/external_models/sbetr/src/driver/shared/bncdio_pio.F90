module bncdio_pio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
! stand alone version
! revision history
! created by Jinyun Tang, modified from clm3.5 years ago.
! updated by Jinyun Tang, May/10/2014.
! minimum functioning to enable the standalone runs
! open file, create file, read data/write data
! check dimension, check variable

  use netcdf
  use shr_sys_mod,  only : shr_sys_abort
  use shr_kind_mod, only : r8 => shr_kind_r8, i4 => shr_kind_i4
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use shr_assert_mod
  use clm_varctl  , only : single_column, iulog
  use spmdMod     , only : masterproc
  implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  private
  save

  public :: check_ret          ! checks return status of netcdf calls
  public :: check_var          ! determine if variable is on netcdf file
  public :: check_dim          ! validity check on dimension
  public :: check_att          ! check attributes
  public :: ncd_pio_openfile   ! open a file
  public :: ncd_pio_openfile_for_write
  public :: ncd_pio_createfile ! create a new file
  public :: ncd_pio_closefile  ! close a file
  public :: ncd_enddef         ! end define mode
  public :: ncd_putatt         ! put attribute
  public :: ncd_getatt         ! get attribute
  public :: ncd_defdim         ! define dimension
  public :: ncd_inqdid         ! inquire dimension id
  public :: ncd_inqdname       ! inquire dimension name
  public :: ncd_inqdlen        ! inquire dimension length
  public :: ncd_inqfdims       ! inquire file dimnesions
  public :: ncd_defvar         ! define variables
  public :: get_dim_len        ! get dimension length

  public :: ncd_inqvid         ! inquire variable id
  public :: ncd_inqvname       ! inquire variable name
  public :: ncd_putvar         ! write local data
  public :: ncd_getvar         ! read data
  public :: ncd_io


  interface ncd_putvar
    module procedure ncd_putvar_int
    module procedure ncd_putvar_real_sp
    module procedure ncd_putvar_int_1d
    module procedure ncd_putvar_real_sp_1d

    module procedure ncd_putvar_int_2d
    module procedure ncd_putvar_real_sp_2d
    module procedure ncd_putvar_int_3d
    module procedure ncd_putvar_real_sp_3d
  end interface ncd_putvar

  interface ncd_getvar
    module procedure ncd_getvar_int
    module procedure ncd_getvar_real_sp
    module procedure ncd_getvar_int_1d
    module procedure ncd_getvar_real_sp_1d

    module procedure ncd_getvar_int_2d
    module procedure ncd_getvar_real_sp_2d
    module procedure ncd_getvar_int_3d
    module procedure ncd_getvar_real_sp_3d

    module procedure ncd_getvar_real_sp_all_1d
    module procedure ncd_getvar_real_sp_all_2d
    module procedure ncd_getvar_real_sp_all_3d
    module procedure ncd_getvar_real_sp_all_4d
  end interface ncd_getvar

  interface get_dim_len
    module procedure get_dim_len_fl
    module procedure get_dim_len_idn
  end interface get_dim_len

  interface ncd_putatt
     module procedure ncd_putatt_int
     module procedure ncd_putatt_real
     module procedure ncd_putatt_char
  end interface

  interface ncd_getatt
     module procedure ncd_getatt_char
  end interface ncd_getatt

  interface ncd_io
     module procedure ncd_io_char_var0_start_glob

     !DIMS 0,1
     module procedure ncd_io_0d_log_glob
     !DIMS 0,1
     module procedure ncd_io_1d_log_glob

     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_0d_int_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_1d_int_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_2d_int_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_3d_int_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_0d_double_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_1d_double_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_2d_double_glob
     !TYPE int,double
     !DIMS 0,1,2,3
     module procedure ncd_io_3d_double_glob

     !TYPE text
     !DIMS 0,1,2
     module procedure ncd_io_0d_text_glob
     !TYPE text
     !DIMS 0,1,2
     module procedure ncd_io_1d_text_glob
     !TYPE text
     !DIMS 0,1,2
     module procedure ncd_io_2d_text_glob

     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_1d_int
     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_2d_int
     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_3d_int
     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_1d_double
     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_2d_double
     !TYPE int,double
     !DIMS 1,2,3
     module procedure ncd_io_3d_double

     !TYPE logical
     !DIMS 1
     module procedure ncd_io_1d_logical
  end interface

  type, public :: file_desc_t
    integer(i4) :: fh
  end type file_desc_t

  type, public :: Var_desc_t
    integer(i4) :: varID
    integer(i4) :: rec     ! this is a record number or pointer into the unlim dimension of the netcdf file
    integer(i4) :: type
  end type

  type IO_desc_t
    integer :: sth_dump
  end type IO_desc_t

  type iodesc_plus_type
     character(len=64) :: name
     type(IO_desc_t)   :: iodesc
     integer           :: type
     integer           :: ndims
     integer           :: dims(4)
     integer           :: dimids(4)
  end type iodesc_plus_type
  integer, public, parameter :: ncd_int =  nf90_int
  integer,parameter,public :: ncd_float     = nf90_float
  integer,parameter,public :: ncd_double    = nf90_double
  integer,parameter,public :: ncd_char      = nf90_char
  integer,parameter,public :: ncd_global    = nf90_global
  integer,parameter,public :: ncd_write     = nf90_write
  integer,parameter,public :: ncd_nowrite   = nf90_nowrite
  integer,parameter,public :: ncd_clobber   = nf90_clobber
  integer,parameter,public :: ncd_noclobber = nf90_noclobber
  integer,parameter,public :: ncd_nofill    = nf90_nofill
  integer,parameter,public :: ncd_unlimited = nf90_unlimited
  integer, parameter :: ncd_log = -nf90_int
!

! !PRIVATE METHODS:
!

!-----------------------------------------------------------------------

  contains


!

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION:
! Validity check on dimension
!
! !ARGUMENTS:

  implicit none
  type(file_desc_t), intent(in) :: ncid    ! file id
  character(len=*), intent(in) :: dimname   !dimension name
  integer, intent(in) :: value              ! dimension length
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: dimid, dimlen    ! temporaries
  character(len=300) :: msg
  character(len=50) :: name
!-----------------------------------------------------------------------

  call check_ret(nf90_inq_dimid (ncid%fh, trim(dimname), dimid), 'check_dim')

  call check_ret(nf90_inquire_dimension (ncid%fh, dimid, name, dimlen), 'check_dim')

  if (dimlen /= value) then
    write (iulog,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
      ' with expected value ',value,' for variable ', trim(dimname)

    call shr_sys_abort(errMsg(mod_filename,__LINE__))
  end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, vardesc, readvar, print_err)
!
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! USES
  use netcdf
  use abortutils, only : endrun
! !ARGUMENTS:
  implicit none
  type(file_desc_t),      intent(in) :: ncid
  character(len=*)       , intent(in) :: varname   ! Varible name to check
  class(Var_desc_t)     , intent(out) :: vardesc   ! Output variable descriptor
  logical               , intent(out) :: readvar   ! If variable exists or not
  logical, optional      , intent(in) :: print_err ! If should print about error

!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: ret     ! return value
  logical :: log_err ! if should log error
  character(len=255) :: msg
!-----------------------------------------------------------------------

  if ( present(print_err) )then
    log_err = print_err
  else
    log_err = .true.
  end if

  readvar = .true.

  ret = nf90_inq_varid (ncid%fh, varname, vardesc%varid)

  if (ret/=nf90_noerr) then
    readvar = .false.
    if(log_err)then
      write(iulog,*)'CHECK_VAR: variable ',trim(varname), ' is not on initial dataset'
    endif
  end if

  end subroutine check_var
!-----------------------------------------------------------------------
  subroutine check_att(ncid, varid, attrib, att_found)
    !
    ! !DESCRIPTION:
    ! Check if attribute is on file
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    logical           ,intent(out)   :: att_found ! true if the attribute was found
    !
    ! !LOCAL VARIABLES:
    integer :: att_type  ! attribute type
    integer :: att_len   ! attribute length
    integer :: status

    character(len=*), parameter :: subname = 'check_att'
    !-----------------------------------------------------------------------

    att_found = .true.

    status = nf90_inquire_attribute(ncid%fh, varid, trim(attrib), att_type, att_len)
    if (status /= nf90_noerr) then
       att_found = .false.
    end if

  end subroutine check_att

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
  use netcdf
  use abortutils, only : endrun
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: ret
  character(len=*) :: calling

! local variable
  character(len=255) :: msg
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

  if (ret /= nf90_noerr) then
    write(6,*)'netcdf error from ',trim(calling)
    write(msg,*)'netcdf strerror = ',trim(nf90_strerror(ret))
    call endrun(trim(msg)//' check_ret')
  end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar
!
! !INTERFACE:
  subroutine ncd_defvar(ncid, varname, xtype, &
      dim1name, dim2name, dim3name, dim4name, dim5name, &
      long_name, units, cell_method, missing_value, fill_value, &
      imissing_value, ifill_value)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
  use netcdf
  use abortutils, only : endrun
  use shr_kind_mod, only : r8 => shr_kind_r8
! !ARGUMENTS:
  implicit none
  type(file_desc_t)     , intent(inout)  :: ncid                    ! input unit
  character(len=*), intent(in)  :: varname                 ! variable name
  integer         , intent(in)  :: xtype                   ! external type
  character(len=*), intent(in), optional :: dim1name       ! dimension name
  character(len=*), intent(in), optional :: dim2name       ! dimension name
  character(len=*), intent(in), optional :: dim3name       ! dimension name
  character(len=*), intent(in), optional :: dim4name       ! dimension name
  character(len=*), intent(in), optional :: dim5name       ! dimension name
  character(len=*), intent(in), optional :: long_name      ! attribute
  character(len=*), intent(in), optional :: units          ! attribute
  character(len=*), intent(in), optional :: cell_method    ! attribute
  real(r8)        , intent(in), optional :: missing_value  ! attribute for real
  real(r8)        , intent(in), optional :: fill_value     ! attribute for real
  integer         , intent(in), optional :: imissing_value ! attribute for int
  integer         , intent(in), optional :: ifill_value    ! attribute for int
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n              ! indices
  integer :: ndims          ! dimension counter
  integer :: dimid(5)       ! dimension ids
  integer :: varid          ! variable id
  integer :: itmp           ! temporary
  integer :: ncid_local
  integer :: lxtype
  type(file_desc_t) :: ncid_tmp
  character(len=256) :: str ! temporary
  character(len=32) :: subname='ncd_def_var' ! subroutine name
!-----------------------------------------------------------------------

  if (.not. masterproc) return
  !print*,'define:',trim(varname)
    ! Determine dimension ids for variable
  ncid_local = ncid%fh
  dimid(:) = 0
  if( xtype == ncd_log) then
    lxtype = ncd_int
  else
    lxtype = xtype
  endif
  if (present(dim1name)) then
    call check_ret(nf90_inq_dimid(ncid_local, dim1name, dimid(1)), subname//' dim1: '//dim1name)
  endif

  if (present(dim2name)) then
    call check_ret(nf90_inq_dimid(ncid_local, dim2name, dimid(2)), subname//' dim2: '//dim2name)
  endif

  if (present(dim3name)) then
    call check_ret(nf90_inq_dimid(ncid_local, dim3name, dimid(3)), subname//' dim3: '//dim3name)
  endif

  if (present(dim4name)) then
    call check_ret(nf90_inq_dimid(ncid_local, dim4name, dimid(4)), subname//' dim4: '//dim4name)
  endif

  if (present(dim5name)) then
    call check_ret(nf90_inq_dimid(ncid_local, dim5name, dimid(5)), subname//' dim5: '//dim5name)
  endif

  !print*,'Define variable'

  if (present(dim1name)) then
    ndims = 0
    do n = 1, size(dimid)
      if (dimid(n) /= 0) ndims = ndims + 1
    end do
    call check_ret(nf90_def_var(ncid_local, trim(varname), xtype, dimid(1:ndims), varid), subname)
  else
    call check_ret(nf90_def_var(ncid_local, trim(varname), xtype,  varid), subname)
  end if
  !print*,'define att'
  if (present(long_name)) then
    call check_ret(nf90_put_att(ncid_local, varid, 'long_name',  trim(long_name)), subname)
  end if

  if (present(units)) then
    call check_ret(nf90_put_att(ncid_local, varid, 'units',  trim(units)), subname)
  end if

  if (present(cell_method)) then
    str = 'time: ' // trim(cell_method)
    call check_ret(nf90_put_att(ncid_local, varid, 'cell_method',  trim(str)), subname)
  end if

  if (present(fill_value)) then
    call ncd_putatt(ncid, varid, '_FillValue',  fill_value, lxtype)
  end if

  if (present(missing_value)) then
    call ncd_putatt(ncid, varid, 'missing_value', missing_value,lxtype)
  end if

  if (present(ifill_value)) then
    call ncd_putatt(ncid, varid, '_FillValue', ifill_value, lxtype)
  end if

  if (present(imissing_value)) then
    call ncd_putatt(ncid, varid, 'missing_value', imissing_value, lxtype)
  end if

  end subroutine ncd_defvar
!----------------------------------------------------------------------
  subroutine ncd_putvar_int(ncid, varname, rec, data)
  !
  !DESCRIPTION
  !put an integer scalar to file

  use netcdf
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  integer, intent(in) :: data
  character(len=*),intent(in) :: varname
  integer :: ans
  integer :: xtype, ndims, varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data, start = (/rec/)),'ncd_putvar_int')

  end subroutine ncd_putvar_int
!----------------------------------------------------------------------

  subroutine ncd_putvar_real_sp(ncid, varname, rec, data)
  !
  !DESCRIPTION
  !put a real scalar to file
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  real(r8), intent(in) :: data
  character(len=*),intent(in) :: varname
  integer :: ans
  integer :: xtype, ndims, varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)
  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data,  start = (/rec/)),'ncd_putvar_real_sp')

  end subroutine ncd_putvar_real_sp
!----------------------------------------------------------------------
  subroutine ncd_putvar_int_1d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! put 1d integer array to file
  use netcdf

!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  integer, dimension(:), intent(in) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data, start = (/1,rec/)),'ncd_putvar_int_1d')

  end subroutine ncd_putvar_int_1d
!----------------------------------------------------------------------
  subroutine ncd_putvar_real_sp_1d(ncid, varname, rec, data)
  !
  ! DESCRIPTIONS
  ! put 1d real array to file
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  real(r8),dimension(:), intent(in) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data, &
     start = (/1,rec/)),'ncd_putvar_real_sp_1d')

  end subroutine ncd_putvar_real_sp_1d
!*****************************************************************

  subroutine ncd_putvar_int_2d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! put 2d integer array to file
  use netcdf
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  integer, dimension(:,:), intent(in) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data, &
    start = (/1,1,rec/)),'ncd_putvar_int_2d')

  end subroutine ncd_putvar_int_2d
!----------------------------------------------------------------------
  subroutine ncd_putvar_real_sp_2d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! put 2d real array to file
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8

!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  real(r8),dimension(:,:), intent(in) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data,  &
     start = (/1,1,rec/)),'ncd_putvar_real_sp_2d')

  end subroutine ncd_putvar_real_sp_2d
!----------------------------------------------------------------------

  subroutine ncd_putvar_int_3d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! put 3d integer to file
  !
  use netcdf
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  integer, dimension(:,:,:), intent(in) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data,  &
   start = (/1,1,1,rec/)),'ncd_putvar_int_2d')

  end subroutine ncd_putvar_int_3d
!----------------------------------------------------------------------
  subroutine ncd_putvar_real_sp_3d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! put 3d real array to file
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  real(r8),dimension(:,:,:), intent(in) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_put_var(ncid%fh, vardesc%varid, data,  &
     start = (/1,1,1,rec/)),'ncd_putvar_real_sp_2d')

      end subroutine ncd_putvar_real_sp_3d
!*********************************************************************
  subroutine ncd_getvar_int(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read a integer scalar
  use netcdf
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  character(len=*),intent(in) :: varname
  integer, intent(out) :: data
  integer :: ans
  integer :: xtype, ndims, varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data, &
    start = (/rec/)),'ncd_getvar_int')

  end subroutine ncd_getvar_int
  !----------------------------------------------------------------------

    subroutine ncd_getvar_real_sp_all_1d(ncid, varname, data)
    !
    !DESCRIPTION
    ! read a real scalar
    use netcdf
    use shr_kind_mod, only : r8 => shr_kind_r8
  !**********************************************************************
    implicit none
    type(file_desc_t), intent(in) :: ncid
    character(len=*),intent(in) :: varname
    REAL(r8), dimension(:), intent(out) :: data

    integer :: varid
    logical :: readvar
    type(Var_desc_t)  :: vardesc

    call check_var(ncid, trim(varname), vardesc, readvar)

    call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data),'ncd_getvar_real_sp_all_1d')

    end subroutine ncd_getvar_real_sp_all_1d

    !----------------------------------------------------------------------

      subroutine ncd_getvar_real_sp_all_2d(ncid, varname, data)
      !
      !DESCRIPTION
      ! read a real scalar
      use netcdf
      use shr_kind_mod, only : r8 => shr_kind_r8
    !**********************************************************************
      implicit none
      type(file_desc_t), intent(in) :: ncid
      character(len=*),intent(in) :: varname
      REAL(r8), dimension(:,:), intent(out) :: data

      integer :: varid
      logical :: readvar
      type(Var_desc_t)  :: vardesc

      call check_var(ncid, trim(varname), vardesc, readvar)

      call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data),'ncd_getvar_real_sp_all_2d')

      end subroutine ncd_getvar_real_sp_all_2d
      !----------------------------------------------------------------------

        subroutine ncd_getvar_real_sp_all_3d(ncid, varname, data)
        !
        !DESCRIPTION
        ! read a real scalar
        use netcdf
        use shr_kind_mod, only : r8 => shr_kind_r8
      !**********************************************************************
        implicit none
        type(file_desc_t), intent(in) :: ncid
        character(len=*),intent(in) :: varname
        REAL(r8), dimension(:,:,:), intent(out) :: data

        integer :: varid
        logical :: readvar
        type(Var_desc_t)  :: vardesc

        call check_var(ncid, trim(varname), vardesc, readvar)

        call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data),'ncd_getvar_real_sp_all_2d')

        end subroutine ncd_getvar_real_sp_all_3d
        !----------------------------------------------------------------------

          subroutine ncd_getvar_real_sp_all_4d(ncid, varname, data)
          !
          !DESCRIPTION
          ! read a real scalar
          use netcdf
          use shr_kind_mod, only : r8 => shr_kind_r8
        !**********************************************************************
          implicit none
          type(file_desc_t), intent(in) :: ncid
          character(len=*),intent(in) :: varname
          REAL(r8), dimension(:,:,:,:), intent(out) :: data

          integer :: varid
          logical :: readvar
          type(Var_desc_t)  :: vardesc

          call check_var(ncid, trim(varname), vardesc, readvar)

          call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data),'ncd_getvar_real_sp_all_2d')

          end subroutine ncd_getvar_real_sp_all_4d
!----------------------------------------------------------------------

  subroutine ncd_getvar_real_sp(ncid, varname, rec, data)
  !
  !DESCRIPTION
  ! read a real scalar
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  character(len=*),intent(in) :: varname
  REAL(r8), intent(out) :: data

  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data, &
    start = (/rec/)),'ncd_getvar_real_sp')

  end subroutine ncd_getvar_real_sp

!----------------------------------------------------------------------
  subroutine ncd_getvar_int_1d(ncid, varname, rec, data)
  !
  !DESCRIPTION
  !read 1d integer array
  !
  use netcdf
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  integer, dimension(:), intent(out) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
     start = (/1,rec/)),'ncd_getvar_int_1d')

  end subroutine ncd_getvar_int_1d
!----------------------------------------------------------------------
  subroutine ncd_getvar_real_sp_1d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read 1d real array
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  real(r8),dimension(:), intent(out) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
     start = (/1,rec/)),'ncd_getvar_real_sp_1d')

  end subroutine ncd_getvar_real_sp_1d


!----------------------------------------------------------------------
  subroutine ncd_getvar_int_2d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read 2d integer array
  !
  use netcdf

!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  integer, dimension(:,:), intent(out) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
      start = (/1,1,rec/)),'ncd_getvar_int_2d')

  end subroutine ncd_getvar_int_2d
!----------------------------------------------------------------------
  subroutine ncd_getvar_real_sp_2d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read 2d real array
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  real(r8),dimension(:,:), intent(out) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
    start = (/1,1,rec/)),'ncd_getvar_real_sp_2d')

  end subroutine ncd_getvar_real_sp_2d
!----------------------------------------------------------------------
  subroutine ncd_getvar_int_3d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read 3d integer array
  use netcdf

!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer, intent(in) :: rec
  integer, dimension(:,:,:), intent(out) :: data
  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
     start = (/1,1,1,rec/)),'ncd_getvar_int_3d')

  end subroutine ncd_getvar_int_3d
!----------------------------------------------------------------------
  subroutine ncd_getvar_real_sp_3d(ncid, varname, rec, data)
  !
  ! DESCRIPTION
  ! read 3d real array
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8
!**********************************************************************
  implicit none
  type(file_desc_t), intent(in) :: ncid
  integer,           intent(in) :: rec
  real(r8),dimension(:,:,:), intent(out) :: data

  character(len=*), intent(in) :: varname
  integer :: varid
  logical :: readvar
  type(Var_desc_t)  :: vardesc

  call check_var(ncid, trim(varname), vardesc, readvar)

  call check_ret( nf90_get_var(ncid%fh, vardesc%varid, data,  &
    start = (/1,1,1,rec/)),'ncd_getvar_real_sp_2d')

  end subroutine ncd_getvar_real_sp_3d

!------------------------------------------------------
  function get_dim_len_idn(ncid, dimname ) result(dimlen)
!
! !DESCRIPTION:
! get dimension length
!
! !ARGUMENTS:
   use netcdf
   implicit none
   type(file_desc_t), intent(in) :: ncid
   character(len=*), intent(in) :: dimname
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: dimid, dimlen    ! temporaries
   character(len=300) :: msg
   character(len=50) :: name
!-----------------------------------------------------------------------

   call check_ret(nf90_inq_dimid (ncid%fh, trim(dimname), dimid), &
          'check_dim')
   call check_ret(nf90_inquire_dimension (ncid%fh, dimid, name,  &
        dimlen), 'check_dim')


   end function get_dim_len_idn

!-----------------------------------------------------------------------

   function get_dim_len_fl(fname, dim_name)result(ans)
   !
   ! DESCRIPTION
   ! get dimension length
   use netcdf
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: dim_name
   integer :: ncid_local, ans
   type(file_desc_t) :: ncid

   call check_ret(nf90_open(fname,      &
     NF90_NOWRITE, ncid_local),'open file '//trim(fname))
   ncid%fh=ncid_local
   ans = get_dim_len_idn(ncid,trim(dim_name))

   call check_ret(nf90_close(ncid_local),'close file'//trim(fname))
   end function get_dim_len_fl
!-----------------------------------------------------------------------

  subroutine ncd_pio_createfile(file, fname)
    !
    ! !DESCRIPTION:
    ! Create a new NetCDF file with PIO
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: file    ! PIO file descriptor
    character(len=*) ,  intent(in)    :: fname   ! File name to create
    !
    !-----------------------------------------------------------------------

    call check_ret(nf90_create(fname, nf90_clobber, file%fh), 'create file '//fname)

    !write(iulog,*) 'Opened file ', trim(fname),  ' to write', file%fh

  end subroutine ncd_pio_createfile
!-----------------------------------------------------------------------
  subroutine ncd_pio_closefile(file)
    !
    ! !DESCRIPTION:
    ! Close a NetCDF PIO file
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: file   ! PIO file handle to close
    !-----------------------------------------------------------------------

    call check_ret(nf90_close(file%fh), 'ncd_pio_closefile')

  end subroutine ncd_pio_closefile

!-----------------------------------------------------------------------
  subroutine ncd_pio_openfile(file, fname, mode)
    !
    ! !DESCRIPTION:
    ! Open a NetCDF PIO file
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: file   ! Output PIO file handle
    character(len=*)   , intent(in)    :: fname  ! Input fname to open
    integer            , intent(in)    :: mode   ! file mode
    !
    ! !LOCAL VARIABLES:
    integer :: ierr
    !-----------------------------------------------------------------------


    !write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    call check_ret(nf90_open(fname, mode, file%fh),'open file '//trim(fname))


  end subroutine ncd_pio_openfile

  subroutine ncd_pio_openfile_for_write(file, fname)
    !
    ! !DESCRIPTION:
    ! Open a NetCDF PIO file
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: file   ! Output PIO file handle
    character(len=*)   , intent(in)    :: fname  ! Input fname to open
    !
    ! !LOCAL VARIABLES:
    integer :: ierr
    !-----------------------------------------------------------------------


    call ncd_pio_openfile(file, fname, NF90_WRITE)

  end subroutine ncd_pio_openfile_for_write

!-----------------------------------------------------------------------
  subroutine ncd_enddef(ncid)
    !
    ! !DESCRIPTION:
    ! enddef netcdf file
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    !
    ! !LOCAL VARIABLES:
    integer :: status   ! error status
    !-----------------------------------------------------------------------

    call check_ret(nf90_enddef(ncid%fh), 'ncd_enddef')

  end subroutine ncd_enddef

!-----------------------------------------------------------------------

  subroutine ncd_putatt_int(ncid,varid,attrib,value,xtype)
    !
    ! !DESCRIPTION:
    ! put integer attributes
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    integer           ,intent(in)    :: value     ! netcdf attrib value
    integer,optional  ,intent(in)    :: xtype     ! netcdf data type
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_put_att(ncid%fh,varid,trim(attrib),value)

  end subroutine ncd_putatt_int

!-----------------------------------------------------------------------

  subroutine ncd_putatt_real(ncid,varid,attrib,value,xtype)
    !
    ! !DESCRIPTION:
    ! put real attributes
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    real(r8)          ,intent(in)    :: value     ! netcdf attrib value
    integer           ,intent(in)    :: xtype     ! netcdf data type
    !
    ! !LOCAL VARIABLES:
    integer :: status
    real*4  :: value4
    !-----------------------------------------------------------------------

    value4 = value

    if (xtype == nf90_double) then
       status = nf90_put_att(ncid%fh,varid,trim(attrib),value)
    else
       status = nf90_put_att(ncid%fh,varid,trim(attrib),value4)
    endif

  end subroutine ncd_putatt_real

  subroutine ncd_putatt_char(ncid,varid,attrib,value,xtype)
    !
    ! !DESCRIPTION:
    ! put character attributes
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*)  ,intent(in)    :: value     ! netcdf attrib value
    integer,optional  ,intent(in)    :: xtype     ! netcdf data type
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_put_att(ncid%fh,varid,trim(attrib),value)

  end subroutine ncd_putatt_char


  subroutine ncd_getatt_char(ncid,varid,attrib,value)
    !
    ! !DESCRIPTION:
    ! get a character attribute
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*)  ,intent(out)   :: value     ! netcdf attrib value
    !
    ! !LOCAL VARIABLES:
    integer :: status

    character(len=*), parameter :: subname = 'ncd_getatt_char'
    !-----------------------------------------------------------------------

    status = nf90_get_att(ncid%fh,varid,trim(attrib),value)

  end subroutine ncd_getatt_char

    !-----------------------------------------------------------------------

  subroutine ncd_defdim(ncid,attrib,value,dimid)
    !
    ! !DESCRIPTION:
    ! define dimension
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    character(len=*)  , intent(in) :: attrib    ! netcdf attrib
    integer           , intent(in) :: value     ! netcdf attrib value
    integer           , intent(out):: dimid     ! netcdf dimension id
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*), parameter :: subname = 'ncd_defdim'
    !-----------------------------------------------------------------------

    call check_ret(nf90_def_dim(ncid%fh,attrib,value,dimid), subname)

  end subroutine ncd_defdim
    !-----------------------------------------------------------------------

  subroutine ncd_inqdid(ncid,name,dimid,dimexist)
    !
    ! !DESCRIPTION:
    ! inquire on a dimension id
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! netcdf file id
    character(len=*) , intent(in) :: name      ! dimension name
    integer          , intent(out):: dimid     ! dimension id
    logical,optional , intent(out):: dimexist  ! if this dimension exists or not
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------


    status = nf90_inq_dimid(ncid%fh,name,dimid)
    if ( present(dimexist) )then
       if ( status == nf90_noerr)then
          dimexist = .true.
       else
          dimexist = .false.
       end if
    end if


  end subroutine ncd_inqdid


  subroutine ncd_inqdname(ncid,dimid,dname)
    !
    ! !DESCRIPTION:
    ! inquire dim name
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    integer           , intent(in) :: dimid     ! dimension id
    character(len=*)  , intent(out):: dname     ! dimension name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_inquire_dimension(ncid%fh,dimid,dname)

  end subroutine ncd_inqdname


  subroutine ncd_inqdlen(ncid,dimid,len,name)
    !
    ! !DESCRIPTION:
    ! enddef netcdf file
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid       ! netcdf file id
    integer           , intent(inout) :: dimid      ! dimension id
    integer           , intent(out)   :: len        ! dimension len
    character(len=*), optional, intent(in) :: name  ! dimension name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    if ( present(name) )then
       call ncd_inqdid(ncid,name,dimid)
    end if
    len = -1
    status = nf90_inquire_dimension(ncid%fh,dimid,len=len)

  end subroutine ncd_inqdlen



  subroutine ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout):: ncid
    logical           , intent(out)  :: isgrid2d
    integer           , intent(out)  :: ni
    integer           , intent(out)  :: nj
    integer           , intent(out)  :: ns
    !
    ! !LOCAL VARIABLES:
    integer  :: dimid                                ! netCDF id
    integer  :: ier                                  ! error status
    character(len=32) :: subname = 'ncd_inqfdims' ! subroutine name
    !-----------------------------------------------------------------------

    if (single_column) then
       ni = 1
       nj = 1
       ns = 1
       isgrid2d = .true.
       RETURN
    end if

    ni = 0
    nj = 0

    ier = nf90_inq_dimid (ncid%fh, 'lon', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=ni)

    ier = nf90_inq_dimid (ncid%fh, 'lat', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=nj)

    ier = nf90_inq_dimid (ncid%fh, 'lsmlon', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=ni)


    ier = nf90_inq_dimid (ncid%fh, 'lsmlat', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=nj)


    ier = nf90_inq_dimid (ncid%fh, 'ni', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=ni)

    ier = nf90_inq_dimid (ncid%fh, 'nj', dimid)
    if (ier == nf90_noerr) ier = nf90_inquire_dimension(ncid%fh, dimid, len=nj)

    ier = nf90_inq_dimid (ncid%fh, 'gridcell', dimid)
    if (ier == nf90_noerr) then
      ier = nf90_inquire_dimension(ncid%fh, dimid, len=ni)
      if(ier == nf90_noerr) nj = 1
    endif


    if (ni == 0 .or. nj == 0) then
       write(iulog,*) trim(subname),' ERROR: ni,nj = ',ni,nj,' cannot be zero '
       call shr_sys_abort(errMsg(mod_filename, __LINE__))
    end if

    if (nj == 1) then
       isgrid2d = .false.
    else
       isgrid2d = .true.
    end if

    ns = ni*nj

  end subroutine ncd_inqfdims


  subroutine ncd_inqvname(ncid,varid,vname,vardesc)
    !
    ! !DESCRIPTION:
    ! inquire variable name
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
    integer           , intent(in)   :: varid     ! variable id
    character(len=*)  , intent(out)  :: vname     ! variable vname
    type(Var_desc_t)  , intent(inout):: vardesc   ! variable descriptor
    !
    ! !LOCAL VARIABLES:
    integer :: status
    integer :: xtype
    !-----------------------------------------------------------------------

    vname = ''
    status = nf90_inquire_variable(ncid%fh,varid,name=vname, xtype=xtype)
    if(status==nf90_noerr)then
      vardesc%varID = varid
      vardesc%type  = xtype
    endif
  end subroutine ncd_inqvname


  subroutine ncd_inqvid(ncid,name,varid,vardesc,readvar)
    !
    ! !DESCRIPTION:
    ! Inquire on a variable ID
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*)  , intent(in)    :: name      ! variable name
    integer           , intent(out)   :: varid     ! variable id
    type(Var_desc_t)  , intent(out)   :: vardesc   ! variable descriptor
    logical, optional , intent(out)   :: readvar   ! does variable exist
    !
    ! !LOCAL VARIABLES:
    integer :: ret               ! return code
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name
    !-----------------------------------------------------------------------

    if (present(readvar)) then
       readvar = .false.
       ret = nf90_inq_varid(ncid%fh,name,varid)
       if (ret /= nf90_noerr) then
          if (masterproc) write(iulog,*) subname//': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
    else
       ret = nf90_inq_varid(ncid%fh,name,varid)
    endif
    vardesc%varid = varid
  end subroutine ncd_inqvid



  subroutine ncd_io_char_var0_start_glob(vardesc, data, flag, ncid, start )
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global character array with start indices input
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    type(var_desc_t) , intent(in)    :: vardesc          ! local vardesc pointer
    character(len=*) , intent(inout) :: data             ! raw data for this index
    integer          , intent(in)    :: start(:)         ! output bounds
    !
    ! !LOCAL VARIABLES:
    integer :: status               ! error code
    character(len=*),parameter :: subname='ncd_io_char_var0_start_glob'

  end subroutine ncd_io_char_var0_start_glob

  subroutine ncd_io_0d_log_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global integer variable
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncid      ! netcdf file id
    character(len=*)   , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*)   , intent(in)    :: varname   ! variable name
    logical            , intent(inout) :: data ! raw data
    logical, optional  , intent(out)   :: readvar   ! was var read?
    integer, optional  , intent(in)    :: nt        ! time sample index
    logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: start(2), count(2) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    integer           :: idata
    integer, pointer  :: idata1d(:)         ! Temporary integer data to send to file
    character(len=32) :: vname              ! variable error checking
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_0d_log_glob'
  end subroutine ncd_io_0d_log_glob


  subroutine ncd_io_1d_log_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global integer variable
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncid      ! netcdf file id
    character(len=*)   , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*)   , intent(in)    :: varname   ! variable name
    logical            , intent(inout) :: data(:) ! raw data
    logical, optional  , intent(out)   :: readvar   ! was var read?
    integer, optional  , intent(in)    :: nt        ! time sample index
    logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: start(2), count(2) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    integer           :: idata
    integer, pointer  :: idata1d(:)         ! Temporary integer data to send to file
    character(len=32) :: vname              ! variable error checking
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_1d_log_glob'
  end subroutine ncd_io_1d_log_glob

  subroutine ncd_io_0d_int_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    integer(i4)         ,           intent(inout) :: data ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    integer(i4) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_0d_int_glob'
  end subroutine ncd_io_0d_int_glob

  subroutine ncd_io_1d_int_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    integer(i4)         ,           intent(inout) :: data(:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    integer(i4) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_1d_int_glob'
  end subroutine ncd_io_1d_int_glob


  subroutine ncd_io_2d_int_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    integer(i4)         ,           intent(inout) :: data(:,:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    integer(i4) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_2d_int_glob'
  end subroutine ncd_io_2d_int_glob


  subroutine ncd_io_3d_int_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    integer(i4)         ,           intent(inout) :: data(:,:,:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    integer(i4) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_3d_int_glob'
  end subroutine ncd_io_3d_int_glob

  subroutine ncd_io_0d_double_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    real(r8)         ,           intent(inout) :: data ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    real(r8) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_0d_double_glob'
  end subroutine ncd_io_0d_double_glob

  subroutine ncd_io_1d_double_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    real(r8)         ,           intent(inout) :: data(:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    real(r8) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_1d_double_glob'
  end subroutine ncd_io_1d_double_glob

  subroutine ncd_io_2d_double_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    real(r8)         ,           intent(inout) :: data(:,:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    real(r8) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_2d_double_glob'
  end subroutine ncd_io_2d_double_glob

  subroutine ncd_io_3d_double_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    real(r8)         ,           intent(inout) :: data(:,:,:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    logical           :: found              ! if true, found lat/lon dims on file
    character(len=32) :: vname              ! variable error checking
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    real(r8) :: temp(1)
    character(len=*),parameter :: subname='ncd_io_3d_double_glob'
  end subroutine ncd_io_3d_double_glob

  subroutine ncd_io_0d_text_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*)         ,           intent(inout) :: data ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_0d_text_glob'
  end subroutine ncd_io_0d_text_glob

  subroutine ncd_io_1d_text_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*)         ,           intent(inout) :: data(:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_1d_text_glob'
  end subroutine ncd_io_1d_text_glob

  subroutine ncd_io_2d_text_glob(varname, data, flag, ncid, readvar, nt, posNOTonfile)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    type(file_desc_t),         intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*)         ,           intent(inout) :: data(:,:) ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    logical         , optional, intent(in)    :: posNOTonfile ! position is NOT on this file
    !
    ! !LOCAL VARIABLES:
    integer           :: m
    integer           :: varid              ! netCDF variable id
    integer           :: start(4), count(4) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=1)  :: tmpString(128)     ! temp for manipulating output string
    type(var_desc_t)  :: vardesc            ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_2d_text_glob'
  end subroutine ncd_io_2d_text_glob

  subroutine ncd_io_1d_int(varname, data, dim1name, flag, ncid, nt, readvar, cnvrtnan2fill)
    !
    ! !DESCRIPTION:
    ! netcdf I/O for 1d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout)        :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    integer(i4)          , pointer               :: data(:)       ! local decomposition data
    character(len=*) , intent(in)            :: dim1name      ! dimension name
    integer          , optional, intent(in)  :: nt            ! time sample index
    logical          , optional, intent(out) :: readvar       ! true => variable is on initial dataset (read only)
    logical          , optional, intent(in)  :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    !
    ! Local Variables
    character(len=8)                 :: clmlevel   ! clmlevel
    character(len=32)                :: dimname    ! temporary
    integer                          :: n          ! index
    integer                          :: iodnum     ! iodesc num in list
    integer                          :: varid      ! varid
    integer                          :: ndims      ! ndims for var
    integer                          :: ndims_iod  ! ndims iodesc for var
    integer                          :: dims(4)    ! dim sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: start(3)   ! netcdf start index
    integer                          :: count(3)   ! netcdf count index
    integer                          :: status     ! error code
    logical                          :: varpresent ! if true, variable is on tape
    integer                , pointer :: idata(:)   ! Temporary integer data to send to file
    integer                , pointer :: compDOF(:)
    type(iodesc_plus_type) , pointer :: iodesc_plus
    type(var_desc_t)                 :: vardesc
    character(len=*),parameter       :: subname='ncd_io_1d_int' ! subroutine name
  end subroutine ncd_io_1d_int

  subroutine ncd_io_1d_double(varname, data, dim1name, flag, ncid, nt, readvar, cnvrtnan2fill)
    !
    ! !DESCRIPTION:
    ! netcdf I/O for 1d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout)        :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    real(r8)          , pointer               :: data(:)       ! local decomposition data
    character(len=*) , intent(in)            :: dim1name      ! dimension name
    integer          , optional, intent(in)  :: nt            ! time sample index
    logical          , optional, intent(out) :: readvar       ! true => variable is on initial dataset (read only)
    logical          , optional, intent(in)  :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    !
    ! Local Variables
    character(len=8)                 :: clmlevel   ! clmlevel
    character(len=32)                :: dimname    ! temporary
    integer                          :: n          ! index
    integer                          :: iodnum     ! iodesc num in list
    integer                          :: varid      ! varid
    integer                          :: ndims      ! ndims for var
    integer                          :: ndims_iod  ! ndims iodesc for var
    integer                          :: dims(4)    ! dim sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: start(3)   ! netcdf start index
    integer                          :: count(3)   ! netcdf count index
    integer                          :: status     ! error code
    logical                          :: varpresent ! if true, variable is on tape
    integer                , pointer :: idata(:)   ! Temporary integer data to send to file
    integer                , pointer :: compDOF(:)
    type(iodesc_plus_type) , pointer :: iodesc_plus
    type(var_desc_t)                 :: vardesc
    character(len=*),parameter       :: subname='ncd_io_1d_double' ! subroutine name
  end subroutine ncd_io_1d_double

  subroutine ncd_io_1d_logical(varname, data, dim1name, flag, ncid, nt, readvar, cnvrtnan2fill)
    !
    ! !DESCRIPTION:
    ! netcdf I/O for 1d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout)        :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    logical          , pointer               :: data(:)       ! local decomposition data
    character(len=*) , intent(in)            :: dim1name      ! dimension name
    integer          , optional, intent(in)  :: nt            ! time sample index
    logical          , optional, intent(out) :: readvar       ! true => variable is on initial dataset (read only)
    logical          , optional, intent(in)  :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    !
    ! Local Variables
    character(len=8)                 :: clmlevel   ! clmlevel
    character(len=32)                :: dimname    ! temporary
    integer                          :: n          ! index
    integer                          :: iodnum     ! iodesc num in list
    integer                          :: varid      ! varid
    integer                          :: ndims      ! ndims for var
    integer                          :: ndims_iod  ! ndims iodesc for var
    integer                          :: dims(4)    ! dim sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: start(3)   ! netcdf start index
    integer                          :: count(3)   ! netcdf count index
    integer                          :: status     ! error code
    logical                          :: varpresent ! if true, variable is on tape
    integer                , pointer :: idata(:)   ! Temporary integer data to send to file
    integer                , pointer :: compDOF(:)
    type(iodesc_plus_type) , pointer :: iodesc_plus
    type(var_desc_t)                 :: vardesc
    character(len=*),parameter       :: subname='ncd_io_1d_logical' ! subroutine name
  end subroutine ncd_io_1d_logical

  subroutine ncd_io_2d_int(varname, data, dim1name, lowerb2, upperb2, &
       flag, ncid, nt, readvar, switchdim, cnvrtnan2fill)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 2d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer(i4)          , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)  :: switchdim       ! true=> permute dim1 and dim2 for output
    logical, optional, intent(in)  :: cnvrtnan2fill   ! true => convert any NaN's to _FillValue (spval)
  end subroutine ncd_io_2d_int

  subroutine ncd_io_2d_double(varname, data, dim1name, lowerb2, upperb2, &
       flag, ncid, nt, readvar, switchdim, cnvrtnan2fill)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 2d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)          , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)  :: switchdim       ! true=> permute dim1 and dim2 for output
    logical, optional, intent(in)  :: cnvrtnan2fill   ! true => convert any NaN's to _FillValue (spval)
  end subroutine ncd_io_2d_double

  subroutine ncd_io_3d_int(varname, data, dim1name, flag, ncid, nt, readvar)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 3d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer(i4)          , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    !
    ! !LOCAL VARIABLES:
    integer                          :: ndim1,ndim2
    character(len=8)                 :: clmlevel   ! clmlevel
    character(len=32)                :: dimname    ! temporary
    integer                          :: status     ! error status
    integer                          :: ndims      ! ndims total for var
    integer                          :: ndims_iod  ! ndims iodesc for var
    integer                          :: varid      ! varid
    integer                          :: n          ! index
    integer                          :: dims(4)    ! dim sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: iodnum     ! iodesc num in list
    integer                          :: start(5)   ! netcdf start index
    integer                          :: count(5)   ! netcdf count index
    logical                          :: varpresent ! if true, variable is on tape
    type(iodesc_plus_type) , pointer :: iodesc_plus
    type(var_desc_t)                 :: vardesc
    character(len=*),parameter :: subname='ncd_io_3d_int' ! subroutine name
  end subroutine ncd_io_3d_int

  subroutine ncd_io_3d_double(varname, data, dim1name, flag, ncid, nt, readvar)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 3d
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)          , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    !
    ! !LOCAL VARIABLES:
    integer                          :: ndim1,ndim2
    character(len=8)                 :: clmlevel   ! clmlevel
    character(len=32)                :: dimname    ! temporary
    integer                          :: status     ! error status
    integer                          :: ndims      ! ndims total for var
    integer                          :: ndims_iod  ! ndims iodesc for var
    integer                          :: varid      ! varid
    integer                          :: n          ! index
    integer                          :: dims(4)    ! dim sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: iodnum     ! iodesc num in list
    integer                          :: start(5)   ! netcdf start index
    integer                          :: count(5)   ! netcdf count index
    logical                          :: varpresent ! if true, variable is on tape
    type(iodesc_plus_type) , pointer :: iodesc_plus
    type(var_desc_t)                 :: vardesc
    character(len=*),parameter :: subname='ncd_io_3d_double' ! subroutine name
  end subroutine ncd_io_3d_double


  end module bncdio_pio
