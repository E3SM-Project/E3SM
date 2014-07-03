module abortutils

  use shr_sys_mod, only: endrun => shr_sys_abort

  implicit none
  private
  save

  public :: endrun

end module abortutils
