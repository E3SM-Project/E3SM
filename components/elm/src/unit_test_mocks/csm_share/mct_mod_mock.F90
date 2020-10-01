module mct_mod

  ! This is a mock of mct_mod, which only includes the bare minimum needed to build CLM
  ! unit tests

  implicit none

  public :: mct_gsMap
  public :: mct_gsMap_OP

  type mct_gsMap
     ! Empty, dummy type
  end type mct_gsMap

contains

  subroutine mct_gsMap_OP(GSMap, PEno, Points)
    ! do-nothing routine that simply matches the signature of mct_gsMap_OP
    type(mct_gsMap), intent(in) :: GSMap
    integer, intent(in) :: PEno
    integer,dimension(:),pointer :: Points
  end subroutine mct_gsMap_OP

end module mct_mod
