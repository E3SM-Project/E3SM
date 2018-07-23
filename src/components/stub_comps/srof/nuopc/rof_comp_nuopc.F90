module rof_comp_nuopc

  ! This is only needed here to satisfy the current cime build requirements
  public :: ROFSetServices

  contains

  subroutine ROFSetServices(gcomp, rc)
    use ESMF, only : ESMF_GridComp
    implicit none
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
  end subroutine ROFSetServices


end module rof_comp_nuopc
