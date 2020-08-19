
#include "diffuse_mom.h"

void diffuse_mom() {
  if (RUN3D) {
    diffuse_mom3D(sgs_field_diag);
  } else {
    diffuse_mom2D(sgs_field_diag);
  }
}

