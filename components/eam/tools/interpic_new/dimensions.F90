module dimensions

   use shr_kind_mod, only: r8 => shr_kind_r8
   use control,      only: verbose

   implicit none

   private

   public :: info
   public :: maxdims
   public :: is_ewdim, is_nsdim, is_zdim, add_dim, is_ncoldim
   public :: get_dimlen, get_shape
   public :: ewdim, nsdim, zdim, ncoldim
   public :: ewdimi, nsdimi, zdimi
   public :: ncoldimid

   include 'netcdf.inc'

   integer, parameter :: maxdims = 4
   integer            :: ncoldimid = -1

   type info
      character*(nf_max_name) :: dimname
      integer                 :: dimlen
   end type info

   type(info) :: ewdim(maxdims), nsdim(maxdims), zdim(maxdims), ncoldim(maxdims)
   type(info) :: ewdimi(maxdims), nsdimi(maxdims), zdimi(maxdims)

   logical :: init_done = .false.

   character*(nf_max_name) :: ewnames(maxdims)
   character*(nf_max_name) :: nsnames(maxdims)
   character*(nf_max_name) :: znames(maxdims)
   character*(nf_max_name) :: ncolnames(maxdims)

   integer :: newnames, nnsnames, nznames, nncolnames

contains
!-------------------------------------------------------------------------------

subroutine init_dims

   integer :: n

   ewnames(:) = ' '
   nsnames(:) = ' '
   znames(:)  = ' '
   ncolnames(:)  = ' '

   ewnames(1) = 'lon'
   ewnames(2) = 'slon'
   nsnames(1) = 'lat'
   nsnames(2) = 'slat'
   znames(1)  = 'lev'
   znames(2)  = 'ilev'
   ncolnames(1) = 'ncol'

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

   do n=1,size (ncolnames)
      if (ncolnames(n) == ' ') exit
   end do
   nncolnames = n - 1

   ewdim(:)%dimname = ' '
   nsdim(:)%dimname = ' '
   zdim(:)%dimname = ' '
   ncoldim(:)%dimname = ' '

   ewdimi(:)%dimname = ' '
   nsdimi(:)%dimname = ' '
   zdimi(:)%dimname = ' '

   init_done = .true.
    
end subroutine init_dims

!-------------------------------------------------------------------------------

logical function is_ewdim (name)

   character*(*) :: name

   integer :: n

   if (.not. init_done) call init_dims
   is_ewdim = .false.

   do n=1,newnames
      if (name == ewnames(n)) then
         is_ewdim = .true.
      end if
   end do

end function is_ewdim

!-------------------------------------------------------------------------------

logical function is_ncoldim (name)

   character*(*) :: name

   integer :: n

   if (.not. init_done) call init_dims
   is_ncoldim = .false.

   do n=1,nncolnames
      if (name == ncolnames(n)) then
         is_ncoldim = .true.
      end if
   end do

end function is_ncoldim

!-------------------------------------------------------------------------------

logical function is_nsdim (name)

   character*(*) :: name

   integer :: n
    
   if (.not. init_done) call init_dims
   is_nsdim = .false.

   do n=1,nnsnames
      if (name == nsnames(n)) then
         is_nsdim = .true.
      end if
   end do

end function is_nsdim

!-------------------------------------------------------------------------------

logical function is_zdim (name)

   character*(*) :: name

   integer :: n
    
   if (.not. init_done) call init_dims
   is_zdim = .false.

   do n=1,nznames
      if (name == znames(n)) then
         is_zdim = .true.
      end if
   end do

end function is_zdim

!-------------------------------------------------------------------------------

subroutine add_dim (arr, dimname, dimlen)
   
   type(info)              :: arr(:)
   character*(nf_max_name) :: dimname
   integer                 :: dimlen
    
   integer :: n

   n = size (arr,1)
   if (arr(n)%dimname(1:1) /= ' ') then
      call err_exit ('add_dim: not enough space allocated for dimension array')
   end if

   ! Postition array arr at first empty slot.
   ! If dimension name is already in arr then stop with error.
   do n = 1, size(arr, 1)
      if (arr(n)%dimname(1:1) == ' ') then
         exit
      else if (trim(arr(n)%dimname) == trim(dimname)) then
         write(6,*)'add_dim: ', trim(dimname), ' already exists'
         stop 999
      end if
   end do

   ! Add dimension name/size to arr
   arr(n)%dimname = trim(dimname)
   arr(n)%dimlen = dimlen

   if (verbose) then
      write(6,*)'add_dim: ', trim(dimname), dimlen
   end if

end subroutine add_dim

!-------------------------------------------------------------------------------

integer function get_dimlen (arr, dimname)

   type(info) :: arr(:)
   character*(nf_max_name) :: dimname

   integer :: n

   do n=1,size (arr,1)
      if (arr(n)%dimname == dimname) then
         get_dimlen = arr(n)%dimlen
         return
      end if
   end do

   write(6,*) 'WARNING: get_dimlen: dimname ',trim(dimname), ' not found'
   get_dimlen=-1
    
end function get_dimlen

!-------------------------------------------------------------------------------

character*8 function get_shape (ncid, vardids, nvdims, dimnames)

   integer, intent(in) :: ncid
   integer, intent(in) :: vardids(nf_max_var_dims) ! variable dimension id's
   integer, intent(in) :: nvdims

   character*(nf_max_name), intent(out) :: dimnames(4)

   ! If a time dimension is present we assume that it's the last dimension.
   ! In that case the dimension name is returned in dimnames, but the shape
   ! is determined by the leading spatial dimensions.

   dimnames(:) = ' '

   if (nvdims > 0) call wrap_inq_dimname(ncid, vardids(1), dimnames(1))
   if (nvdims > 1) call wrap_inq_dimname(ncid, vardids(2), dimnames(2))
   if (nvdims > 2) call wrap_inq_dimname(ncid, vardids(3), dimnames(3))
   if (nvdims > 3) call wrap_inq_dimname(ncid, vardids(4), dimnames(4))

   get_shape = 'unknown'

   if ( is_ncoldim (dimnames(1) ) ) then
      if(is_zdim (dimnames(2))) then
         get_shape='nz'
      else
         get_shape='n'
      end if
   else if (     is_ewdim (dimnames(1)) .and. is_nsdim (dimnames(2)) .and. &
      .not. is_zdim (dimnames(3))) then

      get_shape = 'xy'

   else if (is_ewdim (dimnames(1)) .and. is_nsdim (dimnames(2)) .and. &
            is_zdim (dimnames(3))) then

      get_shape = 'xyz'

   else if (is_ewdim (dimnames(1)) .and. is_zdim  (dimnames(2)) .and. &
             is_nsdim (dimnames(3))) then

      get_shape = 'xzy'

   end if

end function get_shape

end module dimensions
