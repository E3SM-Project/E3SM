
#include "setparm.h"
#include "samxx_utils.h"

extern "C" void setparm() {
  doprecip  = true;
  dosgs     = true;
#ifdef MMF_MOMENTUM_FEEDBACK
  dosurface = false;
#else
  dosurface = true;
#endif
  dodamping = true;
  dt        = crm_dt;
  dx        = crm_dx;
  dy        = crm_dy;
  docloud   = true;
  rank      = 0;   // in MMF model, rank = 0

  if (RUN2D) {
    dy=dx;
  }

  if (RUN2D && YES3D) {
    std::cout << "Error: 2D run and YES3D is set to 1. Exitting...";
    exit(-1);
  }
  if (RUN3D && !YES3D) {
    std::cout << "Error: 3D run and YES3D is set to 0. Exitting...";
    exit(-1);
  }
  if (ny == 1) {
    dy=dx;
  }

  dtn = dt;

  // Turbulence scheme options
  if (is_same_str(turbulence_scheme, "smag") == 0) {
    dosmagor   = true;
    advect_sgs = false;
  }
  if (is_same_str(turbulence_scheme, "shoc") == 0) {
    dosmagor   = false;
    advect_sgs = true;
  }

  if (rank==0) {
    masterproc = true;
  } else {
    masterproc = false;
  }
}


