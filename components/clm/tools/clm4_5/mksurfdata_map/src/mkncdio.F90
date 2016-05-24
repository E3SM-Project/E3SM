module mkncdio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkncdio
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files, and other useful netcdf operations
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod    , only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save

  private

  public :: check_ret       ! checks return status of netcdf calls
  public :: ncd_defvar      ! define netCDF input variable
  public :: get_dim_lengths ! get dimension lengths of a netcdf variable
!
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
  logical  :: masterproc = .true. ! always use 1 proc
  real(r8) :: spval = 1.e36       ! special value

  public :: nf_open
  public :: nf_close
  public :: nf_write
  public :: nf_sync
  public :: nf_inq_attlen
  public :: nf_inq_dimlen
  public :: nf_inq_dimname
  public :: nf_inq_varid
  public :: nf_inq_varndims 
  public :: nf_inq_vardimid
  public :: nf_get_att_double
  public :: nf_get_att_text
  public :: nf_get_var_double
  public :: nf_get_vara_double
  public :: nf_get_var_int
  public :: nf_get_vara_int
  public :: nf_put_var_double
  public :: nf_put_vara_double
  public :: nf_put_var_int
  public :: nf_put_vara_int
  public :: nf_inq_dimid
  public :: nf_max_name
  public :: nf_max_var_dims
  public :: nf_noerr
!EOP
!-----------------------------------------------------------------------

contains

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
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling), ' rcode = ', ret, &
                 ' error = ', NF_STRERROR(ret)
       call abort()
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
! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: ncid                    ! input unit
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
!
! !LOCAL VARIABLES:
!EOP
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=32) :: subname='NCD_DEFVAR_REAL' ! subroutine name
!-----------------------------------------------------------------------

    if (.not. masterproc) return

    ! Determine dimension ids for variable

    dimid(:) = 0

    if (present(dim1name)) then
       call check_ret(nf_inq_dimid(ncid, dim1name, dimid(1)), subname)
    end if
    if (present(dim2name)) then
       call check_ret(nf_inq_dimid(ncid, dim2name, dimid(2)), subname)
    end if
    if (present(dim3name)) then
       call check_ret(nf_inq_dimid(ncid, dim3name, dimid(3)), subname)
    end if
    if (present(dim4name)) then
       call check_ret(nf_inq_dimid(ncid, dim4name, dimid(4)), subname)
    end if
    if (present(dim5name)) then
       call check_ret(nf_inq_dimid(ncid, dim5name, dimid(5)), subname)
    end if

    ! Define variable

    if (present(dim1name)) then
       ndims = 0
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
       call check_ret(nf_def_var(ncid, trim(varname), xtype, ndims, dimid(1:ndims), varid), subname)
    else
       call check_ret(nf_def_var(ncid, varname, xtype, 0, 0, varid), subname)
    end if
    if (present(long_name)) then
       call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
    end if
    if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call check_ret(nf_put_att_text(ncid, varid, 'cell_method', len_trim(str), trim(str)), subname)
    end if
    if (present(fill_value)) then
       call check_ret(nf_put_att_double(ncid, varid, '_FillValue', xtype, 1, fill_value), subname)
    end if
    if (present(missing_value)) then
       call check_ret(nf_put_att_double(ncid, varid, 'missing_value', xtype, 1, missing_value), subname)
    end if
    if (present(ifill_value)) then
       call check_ret(nf_put_att_int(ncid, varid, '_FillValue', xtype, 1, ifill_value), subname)
    end if
    if (present(imissing_value)) then
       call check_ret(nf_put_att_int(ncid, varid, 'missing_value', xtype, 1, imissing_value), subname)
    end if

  end subroutine ncd_defvar

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_dim_lengths
!
! !INTERFACE:
subroutine get_dim_lengths(ncid, varname, ndims, dim_lengths)
!
! !DESCRIPTION:
! Returns the number of dimensions and an array containing the dimension lengths of a
! variable in an open netcdf file.
!
! Entries 1:ndims in the returned dim_lengths array contain the dimension lengths; the
! remaining entries in that vector are meaningless. The dim_lengths array must be large
! enough to hold all ndims values; if not, the code aborts (this can be ensured by passing
! in an array of length nf_max_var_dims).
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   integer         , intent(in) :: ncid           ! netcdf id of an open netcdf file
   character(len=*), intent(in) :: varname        ! name of variable of interest
   integer         , intent(out):: ndims          ! number of dimensions of variable
   integer         , intent(out):: dim_lengths(:) ! lengths of dimensions of variable
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
   integer :: varid
   integer :: dimids(size(dim_lengths))
   integer :: i
   character(len=*), parameter :: subname = 'get_dim_lengths'
!EOP
!------------------------------------------------------------------------------
   call check_ret(nf_inq_varid(ncid, varname, varid), subname)
   call check_ret(nf_inq_varndims(ncid, varid, ndims), subname)

   if (ndims > size(dim_lengths)) then
      write(6,*) trim(subname), ' ERROR: dim_lengths too small'
      call abort()
   end if

   call check_ret(nf_inq_vardimid(ncid, varid, dimids), subname)

   dim_lengths(:) = 0  ! pre-fill with 0 so we won't have garbage in elements past ndims
   do i = 1, ndims
      call check_ret(nf_inq_dimlen(ncid, dimids(i), dim_lengths(i)), subname)
   end do
 end subroutine get_dim_lengths

end module mkncdio
