module reweightMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Top level driver for things that happen when subgrid weights are changed. This is in
  ! a separate module from subgridWeightsMod in order to keep subgridWeightsMod lower-
  ! level - and particularly to break its dependency on filterMod.
  !
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! PUBLIC TYPES:
  implicit none
  save

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: reweight_wrapup               ! do modifications and error-checks after modifying subgrid weights
  
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine reweight_wrapup(bounds)
    !
    ! !DESCRIPTION:
    ! Do additional modifications and error-checks that should be done after modifying subgrid
    ! weights
    !
    ! This should be called whenever any weights change (e.g., pft weights on the column,
    ! landunit weights on the grid cell, etc.).
    !
    ! !USES:
    use filterMod         , only : setFilters
    use subgridWeightsMod , only : set_active, check_weights
    use decompMod         , only : bounds_type, BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! clump bounds
    !------------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_CLUMP, errMsg(__FILE__, __LINE__))

    call set_active(bounds)
    call check_weights(bounds, active_only=.false.)
    call check_weights(bounds, active_only=.true.)
    call setFilters(bounds)

  end subroutine reweight_wrapup

end module reweightMod
