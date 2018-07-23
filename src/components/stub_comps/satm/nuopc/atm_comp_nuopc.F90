module atm_comp_nuopc

  ! This is only needed here to satisfy the current cime build requirements
  public :: ATMSetServices

  contains

  subroutine ATMSetServices(gcomp, rc)
      use ESMF, only : ESMF_GridComp
      implicit none
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc
    end subroutine ATMSetServices


end module atm_comp_nuopc
