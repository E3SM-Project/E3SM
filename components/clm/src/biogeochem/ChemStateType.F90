module ChemStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  
  implicit none
  save
  private
  !----------------------------------------------------
  ! column chemical state variables structure
  !----------------------------------------------------
  type, public :: chemstate_type

     real(r8), pointer :: soil_pH(:,:)    ! soil pH (-nlevsno+1:nlevgrnd)
    
  contains
    procedure, public  :: Init         
    procedure, private :: InitAllocate   
  end type chemstate_type
  
  contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(chemstate_type)         :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init
  
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use clm_varpar            , only : nlevsoi
    !
    ! !ARGUMENTS:
    class(chemstate_type)            :: this
    type(bounds_type), intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------
    begc = bounds%begc;
    endc = bounds%endc
    lbj  = 1;
    ubj  = nlevsoi
    
    allocate(this%soil_pH(begc:endc, lbj:ubj))
    
  end subroutine InitAllocate
end module ChemStateType
