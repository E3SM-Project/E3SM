#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
  use viscosity_base, only: biharmonic_wk, compute_zeta_C0, compute_div_C0, compute_zeta_C0_contra, compute_div_C0_contra, make_c0, make_c0_vector, neighbor_minmax
  use viscosity_base, only: biharmonic_wk_scalar, neighbor_minmax_start,neighbor_minmax_finish, biharmonic_wk_dp3d, smooth_phis
  implicit none
end module viscosity_mod
