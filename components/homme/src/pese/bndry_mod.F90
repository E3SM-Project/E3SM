#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use bndry_mod_base, only: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation, bndry_exchangeS, bndry_exchangeS_start, bndry_exchangeS_finish, sort_neighbor_buffer_mapping
  implicit none
end module bndry_mod
