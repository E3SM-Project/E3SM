#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
  use viscosity_mod_base, only: biharmonic_wk, compute_zeta_C0, compute_div_C0, &
    compute_zeta_C0_contra, compute_div_C0_contra, make_c0, neighbor_minmax, &
    make_c0_vector
  implicit none
end module viscosity_mod
