/**
 * @file emulator_c_api.hpp
 * @brief Define structs and functions to hide all implementation details from Fortran API via opaque pointers
*/
#ifndef EMULATOR_C_API
#define EMULATOR_C_API

#include "coupler_api/coupler_types.hpp"

extern "C" {

/// Opaque handle type in C/Fortran:
/// actually points to an EmulatorComp in C++.
void* emulator_create(const char* kind,
                      const EmulatorCreateConfig* cfg);

void  emulator_set_grid_data(void* handle,
                             const EmulatorGridDesc* grid);

void  emulator_setup_coupling(void* handle,
                              CouplingDesc* cpl);

void emulator_init_coupling_indices(void* handle, const char* import_fields, const char* export_fields);

void  emulator_init(void* handle);
void  emulator_run(void* handle, int dt);
void  emulator_finalize(void* handle);
void  emulator_print_info(void* handle);

/**
 * @brief Destroy an emulator instance created by emulator_create.
 *
 * Removes the instance from the EmulatorRegistry, which drops the
 * shared_ptr reference and runs the C++ destructor.  Call this after
 * emulator_finalize (or instead of it when error-aborting).
 *
 * @param handle Opaque pointer previously returned by emulator_create.
 */
void  emulator_destroy(void* handle);

} // extern "C"
#endif
