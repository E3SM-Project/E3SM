module shoc_intr
!-------------------------------------------------------------------------------
! Dummy module to override src/physics/cam/shoc_intr.F90
!-------------------------------------------------------------------------------
public :: mmf_micro_p3_utils_init
contains
!===============================================================================
subroutine shoc_readnl(nlfile)
  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  return
end subroutine shoc_readnl
!===============================================================================
end module shoc_intr
