#ifndef SCREAM_HOMME_INTERFACE_HPP
#define SCREAM_HOMME_INTERFACE_HPP

#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_assert.hpp"
#include <mpi.h>
#include <string>

namespace scream {

extern "C"
{

// These are F90 routines, called from c
void init_homme1_f90 (const int& f_comm);
void init_homme2_f90 ();
bool was_init_homme1_called_f90 ();
bool was_init_homme2_called_f90 ();
void run_homme_f90 (const double& dt);
void finalize_homme_f90 ();

int get_num_owned_columns_f90 ();
void get_elem_cols_gids_f90 (long* const& gids);
void get_unique_cols_f90 (long* const& gids, int* const& elgp);
int get_homme_int_param_value_f90(const char** name);
double get_homme_real_param_value_f90 (const char** name);
bool get_homme_bool_param_value_f90(const char** name);

} // extern "C"

inline void init_homme1 (const Comm& comm) {
  if (!was_init_homme1_called_f90()) {
    // Make sure MPI is inited
    int flag;
    MPI_Initialized(&flag);
    scream_require_msg(flag==1, "Error! MPI does not seem to be initialized yet.\n");

    // Initialize homme, via F90 calls
    auto f_comm = MPI_Comm_c2f(comm.mpi_comm()); 
    init_homme1_f90(f_comm);
  }
}

template<typename T>
T get_homme_param_value (const std::string& name);

template<>
inline
int get_homme_param_value<int> (const std::string& name) {
  const char* c_name = name.c_str();
  return get_homme_int_param_value_f90(&c_name);
}

template<>
inline
double get_homme_param_value<double> (const std::string& name) {
  const char* c_name = name.c_str();
  return get_homme_real_param_value_f90(&c_name);
}

} // namespace scream

extern "C" {

// These are F90-C interoperable routines in Homme.
// In Homme, they are just used from F90, but they could be useful
// in scream too (even from C++), so we declare them here
namespace Homme {
void init_boundary_exchanges_c ();
}

} // extern "C"

#endif // SCREAM_HOMME_INTERFACE_HPP
