module ice_comp_nuopc

  ! This is only needed here to satisfy the current cime build requirements
  public :: ICESetServices

  contains

  subroutine ICESetServices(gcomp, rc)
      use ESMF, only : ESMF_GridComp
      implicit none
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc
    end subroutine ICESetServices

end module ice_comp_nuopc
