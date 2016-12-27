#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module solver_init_mod_base
!
! 
!  Each model provides a solver_init_mod.F90 which can replace these functions
!  with model specific initialization, if needed.
!
!  implicit none
  private

  public :: solver_init2


contains


  subroutine solver_init2( elem , deriv )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    !do nothing
  end subroutine solver_init2


end module solver_init_mod_base
