module reweightMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: reweightMod
!
! !DESCRIPTION:
! Handles modifications and error-checks related to changing subgrid weights
!
! ----- Requirements for subgrid weights that are enforced here -----
!
! (These requirements are checked in checkWeights/weightsOkay)
!
! Note: in the following, 'active' refers to a pft, column, landunit or grid cell over
! which computations are performed, and 'inactive' refers to a pft, column or landunit
! where computations are NOT performed (grid cells are always active).
! 
! (1) For all columns, landunits and grid cells, the sum of all subgrid weights of its
!     children (or grandchildren, etc.) is equal to 1. For example:
!     - For all columns, the sum of all pft weights on the column equals 1
!     - For all landunits, the sum of all col weights on the landunit equals 1
!     - For all grid cells, the sum of all pft weights on the grid cell equals 1
!     - etc.
! 
! (2) For all ACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
!     its ACTIVE children (or grandchildren, etc.) is equal to 1. For example:
!     - For all active columns, the sum of all pft weights on the column equals 1 when
!       just considering active pfts
!     - For all active landunits, the sum of all col weights on the landunit equals 1 when
!       just considering active cols
!     - For ALL grid cells, the sum of all pft weights on the grid cell equals 1 when
!       just considering active pfts -- note that all grid cells are considered active!
!     - etc.
!
! (3) For all INACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
!     its ACTIVE children, grandchildren, etc. are equal to either 0 or 1. For example:
!     - For all inactive columns, the sum of all pft weights on the column equals either 0
!       or 1 when just considering active pfts
!     - For all inactive landunits, the sum of all col weights on the landunit equals
!       either 0 or 1 when just considering active cols
!     - etc.
!
! Another way of stating (2) and (3) is that the sum of weights of all ACTIVE pfts, cols
! or landunits on their parent/grandparent/etc. is always equal to either 0 or 1 -- and
! must be equal to 1 if this parent/grandparent, etc. is itself active.
!
! Note that, together, conditions (1) and (2) imply that any pft, col or landunit whose
! weight on the grid cell is non-zero must be active. In addition, these conditions imply
! that any pft whose weight on the column is non-zero must be active if the column is
! active (and similarly for any pft on an active landunit, and any col on an active
! landunit).
!
!
! ----- Implications of these requirements for computing subgrid averages -----
!
! The preferred way to average from, say, pft to col is:
!    colval(c) = 0
!    do p = pfti(c), pftf(c)
!       if (active(p)) colval(c) = colval(c) + pftval(p) * wtcol(p)
! (where wtcol(p) is the weight of the pft on the column)
! If column c is active, then the above conditions guarantee that the pwtcol values
! included in the above sum will sum to 1. If column c is inactive, then the above
! conditions guarantee that the pwtcol values included in the above sum will sum to
! either 1 or 0; if they sum to 0, then colval(c) will remain 0.
!
! Another acceptable method is the following; this method accommodates some unknown
! fraction of pftval's being set to spval, and leaves colval set at spval if there are no
! valid pft values:
!    colval(c) = spval
!    sumwt(c) = 0
!    do p = pfti(c), pftf(c)
!       if (active(p) .and. wtcol(p) /= 0) then
!          if (pftval(p) /= spval) then
!             if (sumwt(c) == 0) colval(c) = 0
!             colval(c) = colval(c) + pftval(p) * wtcol(p)
!             sumwt(c) = sumwt(c) + wtcol(p)
!          end if
!       end if
!    end do
!    if (sumwt(c) /= 0) then
!       colval(c) = colval(c) / sumwt(c)
!    end if
! Note that here we check the condition (active(p) .and. wtcol(p) /= 0). We need to
! include a check for wtcol(p) /= 0 because we don't want to set colval(c) = 0 for zero-
! weight pfts in this line:
!             if (sumwt(c) == 0) colval(c) = 0
! And we include a check for active(p) because we don't want to assume that pftval(p) has
! been set to spval for inactive pfts -- we want to allow for the possibility that
! pftval(p) will be NaN for inactive pfts.
!
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only : endrun
  use clm_varctl, only : iulog
!
! PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: reweightWrapup  ! do modifications and error-checks after modifying subgrid weights
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: setActive      ! set 'active' flags at pft, column & landunit level
  private :: is_active_p    ! determine whether the given pft is active
  private :: is_active_c    ! determine whether the given column is active
  private :: is_active_l    ! determine whether the given landunit is active
  private :: checkWeights   ! check subgrid weights
  private :: weightsOkay    ! determine if sum of weights satisfies requirements laid out above
!
!EOP
!-----------------------------------------------------------------------

contains


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: reweightWrapup
!
! !INTERFACE:
  subroutine reweightWrapup(nc)
!
! !DESCRIPTION:
! Do additional modifications and error-checks that should be done after modifying subgrid
! weights
!
! This should be called whenever any weights change (e.g., pft weights on the column,
! landunit weights on the grid cell, etc.).
!
! !USES:
    use filterMod, only : setFilters
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nc       ! clump index
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
!------------------------------------------------------------------------

    call setActive(nc)
    call checkWeights(nc, active_only=.false.)
    call checkWeights(nc, active_only=.true.)
    call setFilters(nc)

  end subroutine reweightWrapup
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: setActive
!
! !INTERFACE:
  subroutine setActive( nc )
!
! !DESCRIPTION:
! Set 'active' flags at the pft, column and landunit level
! (note that grid cells are always active)
!
! This should be called whenever any weights change (e.g., pft weights on the column,
! landunit weights on the grid cell, etc.).
!
! Ensures that we don't have any active pft on an inactive column, or an active column on an
! inactive landunit (since these conditions could lead to garbage data)
!
! !USES:
    use clmtype
    use decompMod , only : get_clump_bounds
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nc       ! clump index
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
    integer :: l,c,p       ! loop counters
    logical :: error_found ! true if we find an error

    character(len=*), parameter :: subname = 'setActive'
!------------------------------------------------------------------------

    ! Determine clump boundaries

    call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    error_found = .false.

    do l = begl,endl
       lptr%active(l) = is_active_l(l)
    end do

    do c = begc,endc
       l = cptr%landunit(c)
       cptr%active(c) = is_active_c(c)
       if (cptr%active(c) .and. .not. lptr%active(l)) then
          write(iulog,*) trim(subname),' ERROR: active column found on inactive landunit', &
                         'at c = ', c, ', l = ', l
          error_found = .true. 
       end if
    end do

    do p = begp,endp
       c = pptr%column(p)
       pptr%active(p) = is_active_p(p)
       if (pptr%active(p) .and. .not. cptr%active(c)) then
          write(iulog,*) trim(subname),' ERROR: active pft found on inactive column', &
                         'at p = ', p, ', c = ', c
          error_found = .true. 
       end if
    end do

    if (error_found) then
       call endrun()
    end if
      
  end subroutine setActive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_active_p
!
! !INTERFACE:
  logical function is_active_p(p)
!
! !DESCRIPTION:
! Determine whether the given pft is active
!
! !USES:
    use clmtype
    use clm_varcon, only : istice_mec
    use domainMod , only : ldomain
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p   ! pft index
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
    integer :: l  ! landunit index
    integer :: g  ! grid cell index
!------------------------------------------------------------------------
    
    l = clm3%g%l%c%p%landunit(p)
    g = clm3%g%l%c%p%gridcell(p)
    
    is_active_p = .false.

    if (clm3%g%l%c%p%wtgcell(p) > 0) is_active_p = .true.

    ! always run over ice_mec landunits within the glcmask, because this is where glc
    ! might need input from virtual (0-weight) landunits
    if (clm3%g%l%itype(l) == istice_mec .and. ldomain%glcmask(g) == 1) is_active_p = .true.

  end function is_active_p
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_active_c
!
! !INTERFACE:
  logical function is_active_c(c)
!
! !DESCRIPTION:
! Determine whether the given column is active
!
! !USES:
    use clmtype
    use clm_varcon, only : istice_mec
    use domainMod , only : ldomain
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: c   ! column index
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
    integer :: l  ! landunit index
    integer :: g  ! grid cell index
!------------------------------------------------------------------------
    
    l = clm3%g%l%c%landunit(c)
    g = clm3%g%l%c%gridcell(c)

    is_active_c = .false.

    if (clm3%g%l%c%wtgcell(c) > 0) is_active_c = .true.

    ! always run over ice_mec landunits within the glcmask, because this is where glc
    ! might need input from virtual (0-weight) landunits
    if (clm3%g%l%itype(l) == istice_mec .and. ldomain%glcmask(g) == 1) is_active_c = .true.
    
  end function is_active_c
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_active_l
!
! !INTERFACE:
  logical function is_active_l(l)
!
! !DESCRIPTION:
! Determine whether the given landunit is active
!
! !USES:
    use clmtype
    use clm_varcon, only : istice_mec
    use domainMod , only : ldomain
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: l   ! landunit index
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
    integer :: g  ! grid cell index
!------------------------------------------------------------------------
    
    g = clm3%g%l%gridcell(l)

    is_active_l = .false.
    
    if (clm3%g%l%wtgcell(l) > 0) is_active_l = .true.

    ! always run over ice_mec landunits within the glcmask, because this is where glc
    ! might need input from virtual (0-weight) landunits
    if (clm3%g%l%itype(l) == istice_mec .and. ldomain%glcmask(g) == 1) is_active_l = .true.

  end function is_active_l
!------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: checkWeights
!
! !INTERFACE:
  subroutine checkWeights (nc, active_only)
!
! !DESCRIPTION:
! Check subgrid weights.
!
! This routine operates in two different modes, depending on the value of active_only. If
! active_only is true, then we check the sum of weights of the ACTIVE children,
! grandchildren, etc. of a given point. If active_only is false, then we check the sum of
! weights of ALL children, grandchildren, etc. of a given point. 
!
! Normally this routine will be called twice: once with active_only=false, and once with
! active_only=true.
!
! !USES
    use clmtype, only : clm3, gridcell_type, landunit_type, &
                        column_type, pft_type
    use decompMod , only : get_clump_bounds

! !ARGUMENTS
    implicit none
    integer, intent(in) :: nc          ! clump index
    logical, intent(in) :: active_only ! true => check sum of weights just of ACTIVE children, grandchildren, etc.
!
! !REVISION HISTORY:
! 2011/10 Created by Zack Subin
! 2012/02 Modified by Bill Sacks for new subgrid weighting convention
! 2012/06 Modified by Bill Sacks for use of 'active' flags
!
!
! !LOCAL VARIABLES:
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
    integer :: g,l,c,p     ! loop counters
    real(r8), allocatable :: sumwtcol(:), sumwtlunit(:), sumwtgcell(:)
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
    logical :: error_found                ! true if we find an error

    character(len=*), parameter :: subname = 'checkWeights'
!EOP
!------------------------------------------------------------------------------

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    call get_clump_bounds(nc,begg,endg,begl,endl,begc,endc,begp,endp)

    allocate(sumwtcol(begc:endc), sumwtlunit(begl:endl), sumwtgcell(begg:endg))

    error_found = .false.

    ! Check PFT-level weights
    sumwtcol(:) = 0._r8
    sumwtlunit(:) = 0._r8
    sumwtgcell(:) = 0._r8

    do p = begp,endp
       c = pptr%column(p)
       l = pptr%landunit(p)
       g = pptr%gridcell(p)

       if ((active_only .and. pptr%active(p)) .or. .not. active_only) then 
          sumwtcol(c) = sumwtcol(c) + pptr%wtcol(p)
          sumwtlunit(l) = sumwtlunit(l) + pptr%wtlunit(p)
          sumwtgcell(g) = sumwtgcell(g) + pptr%wtgcell(p)
       end if
    end do

    do c = begc,endc
       if (.not. weightsOkay(sumwtcol(c), active_only, cptr%active(c))) then
          write(iulog,*) trim(subname),' ERROR: at c = ',c,'total PFT weight is ',sumwtcol(c), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do l = begl,endl
       if (.not. weightsOkay(sumwtlunit(l), active_only, lptr%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total PFT weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do g = begg,endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total PFT weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check col-level weights
    sumwtlunit(:) = 0._r8
    sumwtgcell(:) = 0._r8

    do c = begc,endc
       l = cptr%landunit(c)
       g = cptr%gridcell(c)

       if ((active_only .and. cptr%active(c)) .or. .not. active_only) then
          sumwtlunit(l) = sumwtlunit(l) + cptr%wtlunit(c)
          sumwtgcell(g) = sumwtgcell(g) + cptr%wtgcell(c)
       end if
    end do

    do l = begl,endl
       if (.not. weightsOkay(sumwtlunit(l), active_only, lptr%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total col weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do
    
    do g = begg,endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total col weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check landunit-level weights
    sumwtgcell(:) = 0._r8

    do l = begl,endl
       g = lptr%gridcell(l)
       if ((active_only .and. lptr%active(l)) .or. .not. active_only) then
          sumwtgcell(g) = sumwtgcell(g) + lptr%wtgcell(l)
       end if
    end do

    do g = begg,endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total lunit weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    deallocate(sumwtcol, sumwtlunit, sumwtgcell)

    if (error_found) then
       call endrun()
    end if

    ! Success

  end subroutine checkWeights
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: weightsOkay
!
! !INTERFACE:
  logical function weightsOkay(sumwts, active_weights_only, i_am_active)
! !DESCRIPTION:
! Determine if sumwts (the sum of weights of children, grandchildren or
! great-grandchildren of a column, landunit or grid cell) satisfies the requirements laid
! out above.
!
! The way this is determined depends on the values of two other variables:
! - active_weights_only: does sumwts just include weights of active children,
!   grandchildren or great-grandchilden? (alternative is that it includes weights of ALL
!   children, grandchildren or great-grandchildren)
! - i_am_active: true if the column, landunit or grid cell of interest is active
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: sumwts              ! sum of weights of children, grandchildren or great-grandchildren
    logical , intent(in) :: active_weights_only ! true if sumwts just includes active children, etc.
    logical , intent(in) :: i_am_active         ! true if the current point is active
!
! !REVISION HISTORY:
! Created by Bill Sacks
!
!EOP
!
! LOCAL VARIABLES:
    logical :: weights_equal_1
    real(r8), parameter :: tolerance = 1.e-7_r8  ! tolerance for checking whether weights sum to 1
!------------------------------------------------------------------------

    weights_equal_1 = (abs(sumwts - 1._r8) <= tolerance)

    if (active_weights_only) then
       if (i_am_active) then        ! condition (2) above
          weightsOkay = weights_equal_1
       else                         ! condition (3) above
          weightsOkay = (sumwts == 0._r8 .or. weights_equal_1)
       end if
    else                            ! condition (1) above
       ! (note that i_am_active is irrelevant in this case)
       weightsOkay = weights_equal_1
    end if

  end function weightsOkay
!------------------------------------------------------------------------

end module reweightMod
