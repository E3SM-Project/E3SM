#include "advect_mom.h"

void advect_mom() {

  if (!docolumn) {
    advect2_mom_xy();
    advect2_mom_z();
  }

}
