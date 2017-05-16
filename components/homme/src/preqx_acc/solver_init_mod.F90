#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module solver_init_mod
  !OVERWRITING: solver_init2
  use solver_init_mod_base, only: 
  use dimensions_mod, only: nelemd
  implicit none
  private

  public :: solver_init2


contains


  subroutine solver_init2( elem , deriv )
    use element_mod, only: element_t
    use element_state, only: state_qdp, derived_vn0, derived_divdp, derived_divdp_proj
    use derivative_mod, only: derivative_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    integer :: ie
    !$omp barrier
    !$omp master

    !$acc enter data pcreate(state_Qdp,derived_vn0,derived_divdp,derived_divdp_proj)
    !$acc enter data pcopyin(elem(1:nelemd),deriv)
    do ie = 1 , nelemd
      !$acc enter data pcopyin(elem(ie)%desc)
      !$acc enter data pcopyin(elem(ie)%desc%putmapP,elem(ie)%desc%getmapP,elem(ie)%desc%reverse)
    enddo

    !$omp end master
    !$omp barrier
  end subroutine solver_init2


end module solver_init_mod
