#ifndef SCREAM_HOMME_INTERFACE_HPP
#define SCREAM_HOMME_INTERFACE_HPP

#include "share/grid/abstract_grid.hpp"

#include "Hommexx_Session.hpp"
#include "Types.hpp"
#include "Context.hpp"
#include "mpi/Comm.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_assert.hpp"

#include <mpi.h>
#include <string>

namespace Homme {
extern "C" {
  void initialize_dp3d_from_ps_c ();
} // extern "C"
}

namespace scream {

extern "C"
{

// These are F90 routines, called from c

// Status queries
bool is_parallel_inited_f90 ();
bool is_params_inited_f90 ();
bool is_geometry_inited_f90 ();
bool is_data_structures_inited_f90 ();
bool is_model_inited_f90 ();
bool is_hommexx_functors_inited_f90 ();

// Generic setup
void init_parallel_f90 (const int& f_comm);
void init_params_f90 (const char*& fname);
void init_grids_f90 (const int*& pg_types, const int num_pg_types);
void cleanup_grid_init_data_f90 ();
void finalize_geometry_f90 ();

// Prim init/run/finalize
void prim_init_data_structures_f90 ();
void prim_complete_init1_phase_f90 ();
void prim_set_hvcoords_f90 (const double& ps0,
                            Homme::CF90Ptr& hyai_ptr,
                            Homme::CF90Ptr& hybi_ptr,
                            Homme::CF90Ptr& hyam_ptr,
                            Homme::CF90Ptr& hybm_ptr);
void prim_init_model_f90 ();
void prim_run_f90 (const int nsplit_iteration);
void prim_finalize_f90 ();

// Grids specs
int get_nlev_f90 ();
int get_np_f90 ();
int get_num_local_columns_f90 (const int pgN);
int get_num_global_columns_f90 (const int pgN);
int get_num_local_elems_f90 ();
int get_num_global_elems_f90 ();
void get_dyn_grid_data_f90 (AbstractGrid::gid_type* const& dg_gids,
                            AbstractGrid::gid_type* const& cg_gids,
                            int* const& elgp,
                            double* const& lat, double* const& lon);
void get_phys_grid_data_f90 (const int& pg_type,
                             AbstractGrid::gid_type* const& gids,
                             double* const& lat, double* const& lon, double* const& area);
int get_homme_nsplit_f90 (const int& atm_dt);

// Parmaters getters/setters
int get_homme_int_param_f90(const char** name);
double get_homme_real_param_f90 (const char** name);
bool get_homme_bool_param_f90(const char** name);
void set_homme_int_param_f90(const char** name, const int& value);
void set_homme_real_param_f90 (const char** name, const double& value);
void set_homme_bool_param_f90(const char** name, const bool& value);
void set_homme_log_file_name_f90(const char** fname);

} // extern "C"

inline void init_hommexx (const ekat::Comm& comm) {

  if (!::Homme::Session::m_inited) {

    // First, set mpi comm
    auto& c = ::Homme::Context::singleton();
    c.create_if_not_there<::Homme::Comm>();
    c.get<::Homme::Comm>().reset_mpi_comm(comm.mpi_comm());

    // Inform Homme that kokkos is handled externally
    ::Homme::Session::m_handle_kokkos = false;

    // Make Homme throw rather than abort. In Homme, abort causes finalization of Kokkos,
    // which is bad, since scream still has outstanding views.
    ::Homme::Session::m_throw_instead_of_abort = true;

    // Now let Homme initialize itself
    // Note: this does pretty much nothing, simply prints homme's config settings
    ::Homme::initialize_hommexx_session ();
  }
}

template<typename T>
T get_homme_param (const std::string& name);
template<typename T>
void set_homme_param (const std::string& name, const T& value);

template<>
inline
int get_homme_param<int> (const std::string& name) {
  const char* c_name = name.c_str();
  return get_homme_int_param_f90(&c_name);
}

template<>
inline
double get_homme_param<double> (const std::string& name) {
  const char* c_name = name.c_str();
  return get_homme_real_param_f90(&c_name);
}

template<>
inline
bool get_homme_param<bool> (const std::string& name) {
  const char* c_name = name.c_str();
  return get_homme_bool_param_f90(&c_name);
}

template<>
inline
void set_homme_param<int> (const std::string& name, const int& value) {
  const char* c_name = name.c_str();
  set_homme_int_param_f90(&c_name,value);
}

template<>
inline
void set_homme_param<double> (const std::string& name, const double& value) {
  const char* c_name = name.c_str();
  set_homme_real_param_f90(&c_name,value);
}

template<>
inline
void set_homme_param<bool> (const std::string& name, const bool& value) {
  const char* c_name = name.c_str();
  set_homme_bool_param_f90(&c_name,value);
}

} // namespace scream

#endif // SCREAM_HOMME_INTERFACE_HPP
