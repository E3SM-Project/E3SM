module lnd_comp_nuopc

  ! This is only needed here to satisfy the current cime build requirements
  public :: LNDSetServices

  contains

  subroutine LNDSetServices(gcomp, rc)
      use ESMF, only : ESMF_GridComp
      implicit none
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc
    end subroutine LNDSetServices


end module lnd_comp_nuopc
