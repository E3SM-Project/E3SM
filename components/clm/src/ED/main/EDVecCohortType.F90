module EDVecCohortType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! cohortype.  mimics CLM vector subgrid types.  For now this holds ED data that is
  ! necessary in the rest of CLM
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type 
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  public
  !
  type, public :: EDVecCohort_type
     integer , pointer :: gridcell(:) !index into gridcell level quantities
   contains
     procedure, public :: Init
  end type EDVecCohort_type

  type(EDVecCohort_type), public :: coh
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(EDVecCohort_type) :: this
    type(bounds_type), intent(in) :: bounds
    !------------------------------------------------------------------------

    ! FIX(SPM,032414) pull this out and put in own ED source

    allocate(this%gridcell(bounds%begCohort:bounds%endCohort))

  end subroutine Init

end module EDVecCohortType
