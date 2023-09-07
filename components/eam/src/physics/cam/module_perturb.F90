module module_perturb
  use phys_debug_util, only: phys_debug_col
  use cam_instance,    only: inst_index

  implicit none
  
  integer :: kprnt = 65
contains

  function icolprnt(lchnk) result(icol)
    integer, intent(in)  :: lchnk
    !local
    integer :: icol
    icol = 0
    icol = phys_debug_col(lchnk)
  end function icolprnt
end module module_perturb
