
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module arch_switch_mod
#if USE_OPENACC
  use prim_advection_openacc_mod, only: prim_advec_init1, prim_advec_init2, deriv, Prim_Advec_Tracers_remap
  use openacc_utils_mod, only: arch_init2
#else
  use prim_advection_mod, only: prim_advec_init1, prim_advec_init2, deriv, Prim_Advec_Tracers_remap
#endif
  implicit none

contains

#if (! USE_OPENACC)
  subroutine arch_init2( elem , deriv )
    use element_mod, only: element_t, state_qdp, derived_vn0, derived_divdp, derived_divdp_proj
    use derivative_mod, only: derivative_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    !CPU case, nothing to do
  end subroutine arch_init2
#endif

end module arch_switch_mod

