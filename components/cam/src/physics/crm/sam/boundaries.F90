module boundaries_mod
  use periodic_mod
  use task_util_mod
  implicit none

contains

  subroutine boundaries(ncrms,flag)
    use grid, only: dompi
    implicit none
    integer, intent(in) :: ncrms,flag

    call periodic(ncrms,flag)

  end subroutine boundaries
end module boundaries_mod
