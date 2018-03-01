module dimensions
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: info
  public :: maxdims
  public :: is_ewdim, is_nsdim, is_zdim, add_dim
  public :: get_dimlen, get_shape
  public :: ewdim, nsdim, zdim
  public :: ewdimi, nsdimi, zdimi

  include 'netcdf.inc'

  integer, parameter :: maxdims = 4

  type info
    character*(nf_max_name) :: dimname
    integer                 :: dimlen
  end type info

  type(info) :: ewdim(maxdims), nsdim(maxdims), zdim(maxdims)
  type(info) :: ewdimi(maxdims), nsdimi(maxdims), zdimi(maxdims)

  logical init_done
  data init_done /.false./

  character*(nf_max_name) ewnames (maxdims)
  character*(nf_max_name) nsnames (maxdims)
  character*(nf_max_name) znames (maxdims)

  integer :: newnames, nnsnames, nznames

  contains
!-------------------------------------------------------------------------------
  subroutine init_dims

    integer n

    ewnames(:) = ' '
    nsnames(:) = ' '
    znames(:)  = ' '

    ewnames(1) = 'lon'
    ewnames(2) = 'slon'
    nsnames(1) = 'lat'
    nsnames(2) = 'slat'
    znames(1)  = 'lev'
    znames(2)  = 'ilev'

    do n=1,size (ewnames)
      if (ewnames(n) == ' ') exit
    end do
    newnames = n - 1

    do n=1,size (nsnames)
      if (nsnames(n) == ' ') exit
    end do
    nnsnames = n - 1

    do n=1,size (znames)
      if (znames(n) == ' ') exit
    end do
    nznames = n - 1

    ewdim(:)%dimname = ' '
    nsdim(:)%dimname = ' '
    zdim(:)%dimname = ' '

    ewdimi(:)%dimname = ' '
    nsdimi(:)%dimname = ' '
    zdimi(:)%dimname = ' '

    init_done = .true.
    
    return
  end subroutine init_dims
!-------------------------------------------------------------------------------
  logical function is_ewdim (name)
    character*(*) name

    integer n

    if (.not. init_done) call init_dims
    is_ewdim = .false.

    do n=1,newnames
      if (name == ewnames(n)) then
        is_ewdim = .true.
      end if
    end do

    return
  end function is_ewdim
!-------------------------------------------------------------------------------
  logical function is_nsdim (name)
    character*(*) name

    integer n
    
    if (.not. init_done) call init_dims
    is_nsdim = .false.

    do n=1,nnsnames
      if (name == nsnames(n)) then
        is_nsdim = .true.
      end if
    end do

    return
  end function is_nsdim
!-------------------------------------------------------------------------------
  logical function is_zdim (name)
    character*(*) name

    integer n
    
    if (.not. init_done) call init_dims
    is_zdim = .false.

    do n=1,nznames
      if (name == znames(n)) then
        is_zdim = .true.
      end if
    end do

    return
  end function is_zdim
!-------------------------------------------------------------------------------
  subroutine add_dim (arr, dimname, dimlen)
    implicit none

    type(info) :: arr(:)
    character*(nf_max_name) dimname

    integer dimlen
    
    integer n

    n = size (arr,1)
    if (arr(n)%dimname /= ' ') then
      call err_exit ('add_dim: not enough space allocated for dimension array')
    end if

    do n=1,size (arr,1)
      if (arr(n)%dimname == ' ') then
        exit
      else if (arr(n)%dimname == dimname) then
        write(6,*)'add_dim: ', trim(dimname), ' already exists'
        stop 999
      end if
    end do

    arr(n)%dimname = dimname
    arr(n)%dimlen = dimlen

    return
  end subroutine add_dim
!-------------------------------------------------------------------------------
  integer function get_dimlen (arr, dimname)
    type(info) :: arr(:)
    character*(nf_max_name) dimname

    integer n

    do n=1,size (arr,1)
      if (arr(n)%dimname == dimname) then
        get_dimlen = arr(n)%dimlen
        return
      end if
    end do

    write(6,*) 'get_dimlen: dimname ',trim(dimname), ' not found'
    stop 999
    
  end function get_dimlen
!-------------------------------------------------------------------------------
  character*8 function get_shape (ncid, vardids, nvdims, dimnames)

    implicit none

    include 'netcdf.inc'
!
! Input arguments
!
    integer ncid
    integer vardids(nf_max_var_dims) ! variable dimension id's
    integer nvdims
    character*(nf_max_name) dimnames(3)

    dimnames(:) = ' '

    if (nvdims > 0) call wrap_inq_dimname (ncid, vardids(1), dimnames(1))
    if (nvdims > 1) call wrap_inq_dimname (ncid, vardids(2), dimnames(2))
    if (nvdims > 2) call wrap_inq_dimname (ncid, vardids(3), dimnames(3))

    get_shape = 'unknown'

    if (     is_ewdim (dimnames(1)) .and. is_nsdim (dimnames(2)) .and. &
            .not. is_zdim (dimnames(3))) then

      get_shape = 'xy'

    else if (is_ewdim (dimnames(1)) .and. is_nsdim (dimnames(2)) .and. &
             is_zdim (dimnames(3))) then

      get_shape = 'xyz'

    else if (is_ewdim (dimnames(1)) .and. is_zdim  (dimnames(2)) .and. &
             is_nsdim (dimnames(3))) then

      get_shape = 'xzy'

    end if

    return
  end function get_shape
end module dimensions
