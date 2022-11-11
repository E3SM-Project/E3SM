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
   public :: normalize_classes_by_gcell  ! renormalize array so values are given as % of total grid cell area
   public :: remove_small_cover          ! remove too-small cover types from pft or column-level arrays
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
! !IROUTINE: normalize_classes_by_gcell
!
! !INTERFACE:
subroutine normalize_classes_by_gcell(classes_pct_tot, sums, classes_pct_gcell)
!
! !DESCRIPTION:
! Renormalizes an array (gcell x class) so that values are given as % of total grid cell area
!
! Specifically: Given (1) an array specifying the % cover of different classes, as a % of
! some total ('classes_pct_tot'), and (2) a vector giving these totals ('sums'), expressed
! as % of grid cell area: Returns an array ('classes_pct_gcell') of the same
! dimensionality as classes_pct_tot, where the values now give the % cover of each class
! as a % of total grid cell area.
!
! The size of 'sums' should match the size of the first dimension in 'classes_pct_tot' and
! 'classes_pct_gcell'
!
! For example, if classes_pct_tot(n,i) gives the % of the urban area in grid cell n that is
! in urban class #i, and sums(n) gives the % of grid cell n that is urban, then
! classes_pct_gcell(n,i) will give the % of the total area of grid cell n that is in urban
! class #i.
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: classes_pct_tot(:,:)   ! % cover of classes as % of total
   real(r8), intent(in) :: sums(:)                ! totals, as % of grid cell
   real(r8), intent(out):: classes_pct_gcell(:,:) ! % cover of classes as % of grid cell
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer :: n, n_max

   character(len=*), parameter :: subname = "normalize_classes_by_gcell"
!------------------------------------------------------------------------------

   ! Error-check inputs
   
   n_max = size(sums)
   if (size(classes_pct_tot, 1)   /= n_max .or. &
       size(classes_pct_gcell, 1) /= n_max) then
      write(6,*) subname//' ERROR: array size mismatch'
      write(6,*) 'size(sums)                 = ', n_max
      write(6,*) 'size(classes_pct_tot, 1)   = ', size(classes_pct_tot, 1)
      write(6,*) 'size(classes_pct_gcell, 1) = ', size(classes_pct_gcell, 1)
      call abort()
   end if

   if (size(classes_pct_tot, 2) /= size(classes_pct_gcell, 2)) then
      write(6,*) subname//' ERROR: array size mismatch'
      write(6,*) 'size(classes_pct_tot, 2)   = ', size(classes_pct_tot, 2)
      write(6,*) 'size(classes_pct_gcell, 2) = ', size(classes_pct_gcell, 2)
      call abort()
   end if
   
   ! Do the work

   do n = 1, n_max
      classes_pct_gcell(n,:) = classes_pct_tot(n,:) * (sums(n)/100._r8)
   end do
end subroutine normalize_classes_by_gcell
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: remove_small_cover
!
! !INTERFACE:
subroutine remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small, suma)
!
! !DESCRIPTION:
! Remove any small PFTs (or columns), and adjust pct_gcell and pct_lunit and renormalize
! other elements of the pct_lunit array as necessary. Also returns the number of small
! PFTs/columns found.
!
! Note that suma is just used to determine what the final cover will eventually be once we
! renormalize all landunits to sum to exactly 100%. This is generally unimportant for the
! purposes of this routine, unless suma is substantially different from 100% (e.g., less
! than 50%).
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  real(r8), intent(inout) :: pct_gcell     ! % of the landunit on the grid cell
  real(r8), intent(inout) :: pct_lunit(:)  ! % of each PFT or column on the landunit (adds to 100% before and after this routine)
  integer , intent(out)   :: nsmall        ! number of small (but non-zero) PFTs/columns found
  real(r8), intent(in)    :: too_small     ! threshold for considering a PFT/column too small (% of gridcell)
  real(r8), intent(in), optional :: suma   ! current total % cover of ALL landunits on the gridcell (if not given, assumed to be 100%)
!
! !CALLED FROM:
! subroutine normalizencheck_landuse in mksurfdat.F90
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: my_suma                               ! local version of suma
  logical,  dimension(size(pct_lunit)) :: is_small  ! whether each cover is considered too small (but not 0)
  logical,  dimension(size(pct_lunit)) :: is_zero   ! whether each cover is exactly 0
  integer  :: i                                     ! counter
  real(r8) :: correction                            ! correction needed to pct_gcell

  character(len=*), parameter :: subname = 'remove_small_cover'
!-----------------------------------------------------------------------

  if (present(suma)) then
     if (suma <= 0._r8) then
        write(6,*) 'suma must be > 0; suma = ', suma
        call abort()
     end if
     my_suma = suma
  else
     my_suma = 100._r8
  end if

  ! Determine which array elements are zero, and which are too-small
  is_zero(:)  = (pct_lunit(:) * pct_gcell == 0._r8)
  is_small(:) = ((pct_lunit(:) * pct_gcell / my_suma) < too_small .and. .not. is_zero(:))

  nsmall = count(is_small(:))

  if (nsmall > 0) then

     if (all(is_zero(:) .or. is_small(:))) then
        ! If all PFTs are either 0 or small, then set pct_gcell to 0, but don't touch
        ! pct_lunit(:)
        ! (We do NOT set pct_lunit to all 0 in this case, because we need to maintain
        ! sum(pct_lunit) = 100%)
        pct_gcell = 0._r8

     else
        ! If there are some big PFTs, then we need to adjust pct_lunit as well as
        ! pct_gcell (setting pct_lunit to 0 for the small elements, and renormalizing the
        ! others)

        correction = 0._r8
        do i = 1, size(pct_lunit)
           if (is_small(i)) then
              ! Note that the parenthesized term in the following is the grid cell-level
              ! cover of this point (ignoring the difference between suma and 100, which
              ! is the right thing to do here, since the correction of this difference is
              ! done later, outside this routine)
              correction = correction + (pct_gcell * pct_lunit(i) / 100._r8)
              pct_lunit(i) = 0._r8
           end if
        end do

        pct_gcell = pct_gcell - correction
        pct_lunit(:) = pct_lunit(:) * 100._r8 / sum(pct_lunit(:))
     end if
  end if

end subroutine remove_small_cover


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
