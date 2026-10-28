#ifndef E3SM_EMULATOR_CAPI_TYPES_HPP
#define E3SM_EMULATOR_CAPI_TYPES_HPP

#include "grid_types.hpp"

extern "C" {

struct GridDesc_Ctype {
  int grid_type;
  int nx;
  int ny;
  int mpi_comm;
  int* col_gids;
  double* lat;
  double* lon;
  double* area;
};

struct ExternalField_Ctype {
  const char* label;
  int extent;
};

}
#endif
