#include "task_init.h"

void task_init() {
  if(YES3D != 1 && YES3D != 0) {
    std::cout << "YES3D is not 1 or 0. STOP";
    exit(-1);
  }

  if(YES3D == 1 && ny_gl < 4) {
    std::cout <<"ny_gl is too small for a 3D case.STOP";
    exit(-1);
  }

  if(YES3D == 0 && ny_gl != 1) {
    std::cout << "ny_gl should be 1 for a 2D case. STOP";
    exit(-1);
  }

  if (nsubdomains == 1) {
    rank =0;
    dompi = false;
  }
}


