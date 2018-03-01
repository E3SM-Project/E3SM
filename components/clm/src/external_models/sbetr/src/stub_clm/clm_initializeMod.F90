module clm_initializeMod

  implicit none
  save
  public

  contains

  subroutine initialize(bounds)
    !
    ! !DESCRIPTION:
    ! CLM initialization
  use decompMod       , only : bounds_type
  use clm_instMod     , only : clm_inst
  implicit none
  type(bounds_type), intent(in) :: bounds

  call clm_inst(bounds)

  end subroutine initialize
end module clm_initializeMod
