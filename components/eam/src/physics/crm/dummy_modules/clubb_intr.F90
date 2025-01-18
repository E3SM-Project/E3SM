module clubb_intr
!-------------------------------------------------------------------------------
! Dummy module to override src/physics/cam/clubb_intr.F90
!-------------------------------------------------------------------------------
use shr_kind_mod,  only: r8=>shr_kind_r8
public :: clubb_implements_cnst
public :: clubb_init_cnst
public :: clubb_readnl
contains
!===============================================================================
function clubb_implements_cnst(name)
  ! Return true if specified constituent is implemented
  character(len=*), intent(in) :: name ! constituent name
  logical :: clubb_implements_cnst     ! return value
  clubb_implements_cnst = .false.
end function clubb_implements_cnst
!===============================================================================
subroutine clubb_init_cnst(name, q, gcid)
  ! Initialize the state if clubb_do_adv
  character(len=*), intent(in)  :: name     ! constituent name
  real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
  integer,          intent(in)  :: gcid(:)  ! global column id
  return
end subroutine clubb_init_cnst
!===============================================================================
subroutine clubb_readnl(nlfile)
  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  return
end subroutine clubb_readnl
!===============================================================================
end module clubb_intr
