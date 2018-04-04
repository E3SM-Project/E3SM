#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  model specific initialization code, called from prim_init2
!
module model_init_mod
  use element_mod,    only: element_t
  use hybrid_mod,         only: hybrid_t
  use derivative_mod, only: derivative_t
  use time_mod,       only: timelevel_t
  use hybvcoord_mod,  only: hvcoord_t
  use dimensions_mod, only: np, nelemd, nlev

  implicit none

  public :: model_init2

contains


  subroutine model_init2( elem , hybrid, deriv ,hvcoord,tl,nets,nete)
    use element_state, only: state_qdp, derived_vn0, derived_divdp, derived_divdp_proj, deriv_dvv, hvcoord_dp0, derived_dpdiss_ave, derived_dp, derived_omega_p, derived_eta_dot_dpdn
    use dimensions_mod, only: nelemd

    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(hybrid_t)    , intent(in) :: hybrid
    type(derivative_t), intent(in) :: deriv
    type (hvcoord_t)  , intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    integer                        :: nets,nete

    integer :: ie

    !$omp barrier
    !$omp master

    allocate(deriv_dvv(np,np))
    allocate(hvcoord_dp0(nlev))
    deriv_dvv = deriv%dvv
    hvcoord_dp0 = hvcoord%dp0
  
    !$acc enter data pcreate(state_Qdp,derived_vn0,derived_divdp,derived_divdp_proj,derived_dpdiss_ave,derived_dp,derived_omega_p,derived_eta_dot_dpdn)
    !$acc enter data pcopyin(elem(1:nelemd),deriv,deriv_dvv,hvcoord_dp0)
    do ie = 1 , nelemd
      !$acc enter data pcopyin(elem(ie)%desc)
      !$acc enter data pcopyin(elem(ie)%desc%putmapP,elem(ie)%desc%getmapP,elem(ie)%desc%reverse)
    enddo

    !$omp end master
    !$omp barrier
  end subroutine model_init2


end module model_init_mod
