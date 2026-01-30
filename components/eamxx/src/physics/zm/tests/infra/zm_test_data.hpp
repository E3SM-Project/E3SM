#ifndef SCREAM_ZM_FUNCTIONS_F90_HPP
#define SCREAM_ZM_FUNCTIONS_F90_HPP

#include "share/physics/physics_test_data.hpp"
#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

// Bridge functions to call fortran version of ZM functions from C++

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

struct IentropyData : public PhysicsTestData {
  // Inputs
  Real s, p, qt, tfg;

  // Outputs
  Real t, qst;

  IentropyData(Real s_, Real p_, Real qt_, Real tfg_, Real t_, Real qst_) :
    PhysicsTestData({}, {}),
    s(s_), p(p_), qt(qt_), tfg(tfg_), t(t_), qst(qst_)
  {}

  PTD_STD_DEF(IentropyData, 6, s, p, qt, tfg, t, qst);
};

struct EntropyData : public PhysicsTestData {
  // Inputs
  Real tk, p, qtot;

  // Outputs
  Real entropy;

  EntropyData(Real tk_, Real p_, Real qtot_, Real entropy_) :
    PhysicsTestData({}, {}),
    tk(tk_), p(p_), qtot(qtot_), entropy(entropy_)
  {}

  PTD_STD_DEF(EntropyData, 4, tk, p, qtot, entropy);
};

// Glue functions for host test data. We can call either fortran or CXX with this data (_f -> fortran)
void zm_find_mse_max(zm_data_find_mse_max& d);
void ientropy_f(IentropyData& d);
void ientropy(IentropyData& d);
void entropy_f(EntropyData& d);
void entropy(EntropyData& d);
// End glue function decls

}  // namespace zm
}  // namespace scream

#endif
