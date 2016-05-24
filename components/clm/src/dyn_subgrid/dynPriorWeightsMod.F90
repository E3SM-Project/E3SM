module dynPriorWeightsMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type and associated methods for working with prior subgrid weights
  ! (i.e., before the weight updates of this time step)
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type, BOUNDS_LEVEL_PROC
  use ColumnType     , only : col
  use PatchType      , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: prior_weights_type

  type prior_weights_type
     ! Components are public for ease-of-use and efficiency. However, these components
     ! should be treated as read-only!
     real(r8), allocatable, public :: pwtcol(:)     ! prior pft weight on the column
     logical , allocatable, public :: cactive(:)    ! prior col%active flags
   contains
     procedure :: set_prior_weights      ! set prior weights to current weights
  end type prior_weights_type

  interface prior_weights_type
     module procedure constructor   ! initialize a prior_weights_type object
  end interface prior_weights_type

contains
  
  ! ======================================================================
  ! Constructors
  ! ======================================================================

  ! ----------------------------------------------------------------------
  type(prior_weights_type) function constructor(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize a prior_weights_type object
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds   ! processor bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'prior_weights_type constructor'
    ! ----------------------------------------------------------------------
     
    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')
    
    allocate(constructor%pwtcol(bounds%begp:bounds%endp))
    allocate(constructor%cactive(bounds%begc:bounds%endc))
  end function constructor


  ! ======================================================================
  ! Public methods
  ! ======================================================================
  
  ! ----------------------------------------------------------------------
  subroutine set_prior_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set prior weights to current weights
    !
    ! !ARGUMENTS:
    class(prior_weights_type) , intent(inout) :: this   ! this object
    type(bounds_type)         , intent(in)    :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer :: p, c   ! patch & col indices
    ! ----------------------------------------------------------------------
    
    do p = bounds%begp, bounds%endp
       this%pwtcol(p) = pft%wtcol(p)
    end do

    do c = bounds%begc, bounds%endc
       this%cactive(c) = col%active(c)
    end do
  end subroutine set_prior_weights

end module dynPriorWeightsMod
