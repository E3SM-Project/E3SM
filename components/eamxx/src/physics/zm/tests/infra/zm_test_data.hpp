#ifndef SCREAM_ZM_FUNCTIONS_F90_HPP
#define SCREAM_ZM_FUNCTIONS_F90_HPP

// #include "physics/zm/gw_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

//
// Bridge functions to call fortran version of ZM functions from C++
//

namespace scream {
namespace zm {

struct zm_data_find_mse_max : public PhysicsTestData {
  // Inputs
  Int  pcols;
  Int  ncol;
  Int  pver;
  Int  num_msg;
  bool pergro_active;

  Int  *msemax_top_k;
  Real *temperature;
  Real *zmid;
  Real *sp_humidity;

  // Inputs/Outputs
  Int  *msemax_klev;
  Real *mse_max_val;

  // Constructor
  zm_data_find_mse_max(Int pcols_, Int ncol_, Int pver_, Int num_msg_, bool pergro_active_)
    : PhysicsTestData(
        // dims: group by type and shape.
        {
          {pcols_, pver_},           // (Real, 2D)
          {pcols_},                  // (Real, 1D)
          {pcols_}                   // (Int, 1D)
        },
        // reals: group pointers by shape, in order of appearance above
        {
          {&temperature, &zmid, &sp_humidity}, // (pcols, pver)
          {&mse_max_val}                       // (pcols)
        },
        // ints: group pointers by shape, in order of appearance above
        {
          {&msemax_top_k, &msemax_klev}        // (pcols)
        },
        // bools: if you have 1D bool arrays, list here (none in this case)
        {}
      ),
      pcols(pcols_), ncol(ncol_), pver(pver_), num_msg(num_msg_), pergro_active(pergro_active_)
  {}

  PTD_STD_DEF(zm_data_find_mse_max, 5, pcols, ncol, pver, num_msg, pergro_active);

};

// Glue functions to call fortran from from C++ with the Data struct
void zm_find_mse_max(zm_data_find_mse_max& d);

extern "C" { // _f function decls
}

}  // namespace zm
}  // namespace scream

#endif
