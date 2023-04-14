module FATESFireBase

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for FATES fire data object
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use FireDataBaseType                   , only : fire_base_type

  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_base_type

  type, abstract, extends(fire_base_type) :: fates_fire_base_type
    
    private
    ! !PRIVATE MEMBER DATA:

    contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure(GetLight24_interface),    public, deferred :: GetLight24     ! Return the 24-hour averaged lightning data
      procedure(GetGDP_interface),        public, deferred :: GetGDP         ! Return the global gdp data
      procedure(InitAccBuffer_interface), public, deferred :: InitAccBuffer  ! Initialize accumulation processes
      procedure(InitAccVars_interface),   public, deferred :: InitAccVars    ! Initialize accumulation variables
      procedure(UpdateAccVars_interface), public, deferred :: UpdateAccVars  ! Update/extract accumulations vars
      procedure(need_lightning_and_popdens_interface), public, deferred :: &
                                                  need_lightning_and_popdens ! Returns true if need lightning & popdens
      
  end type fates_fire_base_type

  !-----------------------

  abstract interface
  !-----------------------------------------------------------------------

  function need_lightning_and_popdens_interface(this) result(need_lightning_and_popdens)
    !
    ! !DESCRIPTION:
    ! Returns true if need lightning and popdens, false otherwise
    !
    ! USES
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type), intent(in) :: this
    logical :: need_lightning_and_popdens  ! function result
    !-----------------------------------------------------------------------
  end function need_lightning_and_popdens_interface
  
  !------------------------------------------------------------------------
  function GetLight24_interface( this ) result(lnfm24)
    !
    ! !DESCRIPTION: Get the 24-hour averaged lightning data
    ! !USES
    use shr_kind_mod   , only: r8 => shr_kind_r8
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    real(r8), pointer :: lnfm24(:)
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
  end function GetLight24_interface
  
  !------------------------------------------------------------------------
  function GetGDP_interface( this ) result(gdp)
    !
    ! !DESCRIPTION: Get the global gross domestic product data
    ! !USES
    use shr_kind_mod   , only: r8 => shr_kind_r8
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    real(r8), pointer :: gdp(:)
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
  end function GetGDP_interface

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer_interface (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize the accumulation buffers
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine InitAccBuffer_interface

  !-----------------------------------------------------------------------
  subroutine InitAccVars_interface(this, bounds)
    !
    ! !DESCRIPTION:
    !      Initialize the accumulation variables
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine InitAccVars_interface

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars_interface (this, bounds)
    !
    ! !DESCRIPTION:
    ! Update accumulation variables
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine UpdateAccVars_interface

  end interface

end module FATESFireBase
