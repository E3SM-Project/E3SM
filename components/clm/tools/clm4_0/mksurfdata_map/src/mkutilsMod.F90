module mkutilsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkutils
!
! !DESCRIPTION:
! General-purpose utilities for mksurfdata_map
!
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
 
   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: slightly_below
   public :: slightly_above
   public :: convert_latlon
!
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!EOP
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slightly_below
!
! !INTERFACE:
logical function slightly_below(a, b, eps)
!
! !DESCRIPTION:
! Returns true if a is slightly below b; false if a is significantly below b or if a is
! greater than or equal to b
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: a
   real(r8), intent(in) :: b

   ! if provided, eps gives the relative error allowed for checking the "slightly"
   ! condition; if not provided, the tolerance defaults to the value given by eps_default
   real(r8), intent(in), optional :: eps
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: l_eps
   real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
!------------------------------------------------------------------------------

   if (present(eps)) then
      l_eps = eps
   else
      l_eps = eps_default
   end if

   if (a < b .and. (b - a)/b < l_eps) then
      slightly_below = .true.
   else
      slightly_below = .false.
   end if

end function slightly_below
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slightly_above
!
! !INTERFACE:
logical function slightly_above(a, b, eps)
!
! !DESCRIPTION:
! Returns true if a is slightly above b; false if a is significantly above b or if a is
! less than or equal to b
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: a
   real(r8), intent(in) :: b

   ! if provided, eps gives the relative error allowed for checking the "slightly"
   ! condition; if not provided, the tolerance defaults to the value given by eps_default
   real(r8), intent(in), optional :: eps
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: l_eps
   real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
!------------------------------------------------------------------------------

   if (present(eps)) then
      l_eps = eps
   else
      l_eps = eps_default
   end if

   if (a > b .and. (a - b)/b < l_eps) then
      slightly_above = .true.
   else
      slightly_above = .false.
   end if

end function slightly_above
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_latlon
!
! !INTERFACE:
  subroutine convert_latlon(ncid, varname, data)
!
! !DESCRIPTION:
! Convert a latitude or longitude variable from its units in the input file to degrees E /
! degrees N. Currently, this just handles conversions from radians to degrees.
!
! Assumes that the longitude / latitude variable has already been read from file, into
! the variable given by 'data'. ncid & varname give the file ID and variable name from
! which this variable was read (needed to obtain the variable's units).
!
! !USES:
    use mkncdio
    use shr_const_mod, only : SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
    integer         , intent(in)   :: ncid      ! ID of open netcdf file
    character(len=*), intent(in)   :: varname   ! name of lat or lon variable that was read into 'data'
    real(r8)        , intent(inout):: data(:)   ! latitude or longitude data
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ier                             ! error return code
    integer :: varid                           ! netCDF variable id
    integer :: units_len                       ! length of units attribute on file
    character(len=256) :: units                ! units attribute
    character(len= 32) :: subname = 'convert_latlon'
!-----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, varname, varid), subname)
    ier = nf_inq_attlen(ncid, varid, 'units', units_len)

    ! Only do the following processing if there is no error; if ier /= NF_NOERR, that
    ! probably means there isn't a units attribute -- in that case, assume units are
    ! degrees and need no conversion
    if (ier == NF_NOERR) then
       if (units_len > len(units)) then
          write(6,*) trim(subname), ' ERROR: units variable not long enough to hold attributue'
          call abort()
       end if

       call check_ret(nf_get_att_text(ncid, varid, 'units', units), subname)
       
       if (units(1:7) == 'radians') then
          ! convert from radians to degrees
          data(:) = data(:) * 180._r8 / SHR_CONST_PI
       end if
    end if

  end subroutine convert_latlon
!------------------------------------------------------------------------------


end module mkutilsMod
