module wav_comp_nuopc  ! This is only needed here to satisfy the current cime build requirements
  public :: WAVSetServices
  contains
    subroutine WAVSetServices(gcomp, rc)
      use ESMF, only : ESMF_GridComp
      implicit none
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc
    end subroutine WAVSetServices
end module wav_comp_nuopc
