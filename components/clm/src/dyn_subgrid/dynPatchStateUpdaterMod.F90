module dynPatchStateUpdaterMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Class for adjusting patch-level (aboveground) state variables due to transient patch
  ! areas.
  !
  ! In each time step, the object should be set up with:
  !
  !    - call patch_state_updater%set_old_weights (before dyn subgrid weight updates)
  !
  !    - call patch_state_updater%set_new_weights (after dyn subgrid weight updates)
  !
  ! Then it can be used to update each state variable with a call to:
  !
  !    - call patch_state_updater%update_patch_state
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use decompMod            , only : bounds_type, BOUNDS_LEVEL_PROC
  use VegetationType       , only : veg_pp
  use ColumnType           , only : col_pp
  use clm_varpar           , only : mxpft
  use abortutils           , only : endrun
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: patch_state_updater_type

  type patch_state_updater_type
     private
     real(r8), allocatable :: pwtgcell_old(:) ! old patch weights on the gridcell
     real(r8), allocatable :: pwtgcell_new(:) ! new patch weights on the gridcell

     real(r8), allocatable :: cwtgcell_old(:) ! old column weights on the gridcell

     ! (pwtgcell_new - pwtgcell_old) from last call to set_new_weights
     real(r8), allocatable :: dwt(:)

     ! (pwtgcell_old / pwtgcell_new) from last call to set_new_weights; only valid for
     ! growing patches
     real(r8), allocatable :: growing_old_fraction(:)

     ! (dwt / pwtgcell_new) from last call to set_new_weights; only valid for growing
     ! patches
     real(r8), allocatable :: growing_new_fraction(:)

   contains
     ! Public routines
     procedure, public :: set_old_weights     ! set weights before dyn subgrid updates
     procedure, public :: set_new_weights     ! set weights after dyn subgrid updates

     ! Update a patch-level state variable and compute associated fluxes based on changing
     ! patch areas
     procedure, public :: update_patch_state

     ! Update a patch-level state variable and compute associated fluxes based on
     ! changing patch areas, with flux out partitioned into two fluxes. Partitioning is
     ! based on pft type.
     procedure, public :: update_patch_state_partition_flux_by_type

     ! returns a patch-level logical array that is true wherever the patch weight was zero
     ! prior to weight updates
     procedure, public :: old_weight_was_zero

     ! returns a patch-level logical array that is true wherever the patch grew in this
     ! time step
     procedure, public :: patch_grew

     ! returns a patch-level logical array that is true wherever a patch newly has
     ! non-zero weight in this time step
     procedure, public :: patch_initiating

  end type patch_state_updater_type

  interface patch_state_updater_type
     module procedure constructor
  end interface patch_state_updater_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(bounds) result(this)
    !
    ! !DESCRIPTION:
    ! Initialize a patch_state_updater_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(patch_state_updater_type) :: this  ! function result
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, errMsg(sourcefile, __LINE__))

    begp = bounds%begp
    endp = bounds%endp
    begc = bounds%begc
    endc = bounds%endc

    allocate(this%pwtgcell_old(begp:endp))
    this%pwtgcell_old(:) = nan
    allocate(this%pwtgcell_new(begp:endp))
    this%pwtgcell_new(:) = nan
    allocate(this%cwtgcell_old(begc:endc))
    this%cwtgcell_old(:) = nan
    allocate(this%dwt(begp:endp))
    this%dwt(:) = nan
    allocate(this%growing_old_fraction(begp:endp))
    this%growing_old_fraction(:) = nan
    allocate(this%growing_new_fraction(begp:endp))
    this%growing_new_fraction(:) = nan

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine set_old_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights before dyn subgrid updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p
    integer :: c

    character(len=*), parameter :: subname = 'set_old_weights'
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       c = veg_pp%column(p)
       this%pwtgcell_old(p) = veg_pp%wtgcell(p)
       this%cwtgcell_old(c) = col_pp%wtgcell(c)
    end do

  end subroutine set_old_weights

  !-----------------------------------------------------------------------
  subroutine set_new_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights after dyn subgrid updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p

    character(len=*), parameter :: subname = 'set_new_weights'
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       this%pwtgcell_new(p) = veg_pp%wtgcell(p)
       this%dwt(p)          = this%pwtgcell_new(p) - this%pwtgcell_old(p)
       if (this%dwt(p) > 0._r8) then
          this%growing_old_fraction(p) = this%pwtgcell_old(p)/ this%pwtgcell_new(p)
          this%growing_new_fraction(p) = this%dwt(p)         / this%pwtgcell_new(p)
       else
          ! These values are unused in this case, but set them to something reasonable for
          ! safety. (We could set them to NaN, but that requires a more expensive
          ! subroutine call, using the shr_infnan_mod infrastructure.)
          this%growing_old_fraction(p) = 1._r8
          this%growing_new_fraction(p) = 0._r8
       end if
    end do

  end subroutine set_new_weights

  !-----------------------------------------------------------------------
  subroutine update_patch_state(this, bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       var, flux_out_col_area, flux_out_grc_area, &
       seed, seed_addition)
    !
    ! !DESCRIPTION:
    ! Update a patch-level state variable and compute associated fluxes based on changing
    ! patch areas
    !
    ! For growing patches, this subroutine adjusts var. For shrinking patches, this
    ! subroutine accumulates flux in flux_out_col_area.
    !
    ! Changes are only made within the given filter. Note that this filter should include
    ! inactive as well as active patches, so that it includes patches that just became
    ! inactive.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer, intent(in) :: filterp_with_inactive(:) ! patch filter that includes inactive points
    real(r8), intent(inout) :: var( bounds%begp: ) ! patch-level state variable

    ! Accumulated flux from shrinking areas, expressed as mass per unit area COLUMN, using
    ! the OLD column weight. (The use of the old column weight is appropriate if the
    ! fluxes from shrinking patches are applied to column-level state variables BEFORE
    ! doing column-level state adjustments via the column_state_updater). For shrinking
    ! areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_col_area( bounds%begp: )

    ! Accumulated flux from shrinking areas, expressed as mass per unit area GRIDCELL. For
    ! shrinking areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_grc_area( bounds%begp: )

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed( bounds%begp: )

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition( bounds%begp: )
    !
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, c

    character(len=*), parameter :: subname = 'update_patch_state'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(var) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    if (present(flux_out_col_area)) then
       SHR_ASSERT_ALL((ubound(flux_out_col_area) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if

    if (present(flux_out_grc_area)) then
       SHR_ASSERT_ALL((ubound(flux_out_grc_area) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if

    if (present(seed)) then
       SHR_ASSERT_ALL((ubound(seed) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if

    if (present(seed_addition)) then
       if (.not. present(seed)) then
          call endrun(subname//' ERROR: seed_addition can only be provided if seed is provided')
       end if
       SHR_ASSERT_ALL((ubound(seed_addition) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    end if

    do fp = 1, num_filterp_with_inactive
       p = filterp_with_inactive(fp)
       c = veg_pp%column(p)

       if (this%dwt(p) > 0._r8) then
          var(p) = var(p) * this%growing_old_fraction(p)
          if (present(seed)) then
             var(p) = var(p) + seed(p) * this%growing_new_fraction(p)
             if (present(seed_addition)) then
                seed_addition(p) = seed_addition(p) + seed(p) * this%dwt(p)
             end if
          end if

       else if (this%dwt(p) < 0._r8) then
          if (present(flux_out_grc_area)) then
             flux_out_grc_area(p) = flux_out_grc_area(p) + var(p) * this%dwt(p)
          end if
          if (present(flux_out_col_area)) then
             ! No need to check for divide by 0 here: If dwt < 0 then we must have had
             ! cwtgcell_old > 0.
             flux_out_col_area(p) = flux_out_col_area(p) + &
                  var(p) * (this%dwt(p) / this%cwtgcell_old(c))
          end if
       end if
    end do

  end subroutine update_patch_state

  !-----------------------------------------------------------------------
  subroutine update_patch_state_partition_flux_by_type(this, bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       flux1_fraction_by_pft_type, &
       var, flux1_out, flux2_out, &
       seed, seed_addition)
    !
    ! !DESCRIPTION:
    ! Update a patch-level state variable and compute associated fluxes based on
    ! changing patch areas, with flux out partitioned into two fluxes. Partitioning is
    ! based on pft type.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer, intent(in) :: filterp_with_inactive(:) ! patch filter that includes inactive points
    real(r8), intent(in) :: flux1_fraction_by_pft_type( 0: ) ! fraction of flux that goes into flux1_out, indexed by pft type
    real(r8), intent(inout) :: var( bounds%begp: ) ! patch-level state variable

    ! Accumulated fluxes from shrinking areas. For shrinking areas, these are given as
    ! NEGATIVE quantities. Even though these are patch-level arrays, they are expressed
    ! as mass per unit area GRIDCELL (so these are equivalent to the flux_out_grc_area
    ! argument in the main update_patch_state routine).
    real(r8), intent(inout) :: flux1_out( bounds%begp: )
    real(r8), intent(inout) :: flux2_out( bounds%begp: )

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed( bounds%begp: )

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p
    real(r8) :: total_flux_out(bounds%begp:bounds%endp)
    real(r8) :: my_flux1_fraction

    character(len=*), parameter :: subname = 'update_patch_state_partition_flux_by_type'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(flux1_fraction_by_pft_type) == (/mxpft/)), errMsg(sourcefile, __LINE__))

    total_flux_out(bounds%begp:bounds%endp) = 0._r8
    call this%update_patch_state(bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       var, flux_out_grc_area = total_flux_out, &
       seed = seed, seed_addition = seed_addition)

    do fp = 1, num_filterp_with_inactive
       p = filterp_with_inactive(fp)
       my_flux1_fraction = flux1_fraction_by_pft_type(veg_pp%itype(p))
       flux1_out(p) = flux1_out(p) + total_flux_out(p) * my_flux1_fraction
       flux2_out(p) = flux2_out(p) + total_flux_out(p) * (1._r8 - my_flux1_fraction)
    end do

  end subroutine update_patch_state_partition_flux_by_type


  !-----------------------------------------------------------------------
  function old_weight_was_zero(this, bounds)
    !
    ! !DESCRIPTION:
    ! Returns a patch-level logical array that is true wherever the patch weight was zero
    ! prior to weight updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    logical :: old_weight_was_zero(bounds%begp:bounds%endp)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: p

    character(len=*), parameter :: subname = 'old_weight_was_zero'
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       old_weight_was_zero(p) = (this%pwtgcell_old(p) == 0._r8)
    end do

  end function old_weight_was_zero

  !-----------------------------------------------------------------------
  function patch_grew(this, bounds)
    !
    ! !DESCRIPTION:
    ! Returns a patch-level logical array that is true wherever the patch grew in this
    ! time step
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    logical :: patch_grew(bounds%begp:bounds%endp)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: p

    character(len=*), parameter :: subname = 'patch_grew'
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       patch_grew(p) = (this%dwt(p) > 0._r8)
    end do

  end function patch_grew

  !-----------------------------------------------------------------------
  function patch_initiating(this, bounds)
    !
    ! !DESCRIPTION:
    ! Returns a patch-level logical array wherever the patch is initiating - i.e., growing
    ! from zero area to non-zero area - in this time step
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(patch_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    logical :: patch_initiating(bounds%begp:bounds%endp)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: p

    character(len=*), parameter :: subname = 'patch_initiating'
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       patch_initiating(p) = ( &
            this%pwtgcell_old(p) == 0._r8 .and. &
            this%pwtgcell_new(p) > 0._r8)
    end do

  end function patch_initiating


end module dynPatchStateUpdaterMod
