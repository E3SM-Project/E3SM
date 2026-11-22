#ifndef E3SM_EMULATOR_CAPI_TYPES_HPP
#define E3SM_EMULATOR_CAPI_TYPES_HPP

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

void* harness_create();
void harness_destroy(void* handle);
void harness_register_field(void* handle, const char* name, int extent,
                            int dtype, void* data);

}
#endif
