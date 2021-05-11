
#pragma once

#include "samxx_const.h"
#include "vars.h"

void bound_exchange(real4d &f, int dimz, int i_1, int i_2, int j_1, int j_2, int id);
void bound_exchange(real5d &f, int offL, int dimz, int i_1, int i_2, int j_1, int j_2, int id);

YAKL_INLINE int constexpr _IDX(int const l1, int const u1, int const i1, int const l2, int const u2, 
                               int const i2, int const l3, int const u3, int const i3, int const l4, 
                               int const u4, int const i4) {
  return ( ((i4)-(l4))*((u3)-(l3)+1)*((u2)-(l2)+1)*((u1)-(l1)+1) + 
           ((i3)-(l3))              *((u2)-(l2)+1)*((u1)-(l1)+1) + 
           ((i2)-(l2))                            *((u1)-(l1)+1) + 
           ((i1)-(l1)) );
}

