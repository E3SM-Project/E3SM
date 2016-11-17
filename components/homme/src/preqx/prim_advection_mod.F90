#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_mod
  use prim_advection_mod_base, only: vertical_remap, Prim_Advec_Tracers_remap, prim_advec_init1, prim_advec_init2, &
                                     deriv
  implicit none
end module prim_advection_mod
