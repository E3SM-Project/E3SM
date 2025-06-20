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

// // The Data struct is special; it is used to do ZM initialization, which
// // must be called before any ZM function.
// struct ZMInit : public PhysicsTestData {
//   // Inputs
//   Int pver;

//   ZMInit( Int pver_ ) :
//     PhysicsTestData({
//       {pver_ + 1}
//     },
//     // {
//     //   {&cref},
//     //   {&alpha}
//     // }
//     ),
//     pver(pver_)
//   {
//     // Assert valid init data?
//     assert(ktop <= pver);
//   }

//   PTD_STD_DEF(ZMInit, 11, pver );
// };


// struct zm_data_common : public PhysicsTestData {
//   // Inputs
//   Int ncol, pver;

//   // Inputs/Outputs
//   // Real *tau;

//   // Outputs
//   // Real *gwut, *utgw, *vtgw;

//   zm_data_common(Int ncol_, Int pver_) :
//     PhysicsTestData({
//       {ncol_},
//       {ncol_, init_.pver},
//     //   {ncol_, 2*init_.pgwv},
//     //   {ncol_, 2*init_.pgwv, init_.pver + 1},
//     //   {ncol_, init_.pver, 2*ngwv_},
//     //   {ncol_}
//     },
//     {
//     //   {&lat, &xv, &yv},
//     //   {&dpm, &rdpm, &ubm, &t, &nm, &utgw, &vtgw},
//     //   {&c},
//     //   {&tau},
//     //   {&gwut}
//     // },
//     // {
//     //   {&tend_level}
//     }),
//     ncol(ncol_), pver(pver_)
//   {}

//   PTD_STD_DEF_INIT(zm_data_find_mse_max, 2, ncol, pver );
// };

/*
zm_find_mse_max_c( 
pcols, 
ncol, 
pver, 
num_msg, 
msemax_top_k, 
pergro_active, 
temperature, 
zmid, 
sp_humidity, 
zm_const, 
zm_param, 
msemax_klev, 
mse_max_val )
*/

struct zm_data_find_mse_max : public PhysicsTestData {
  // Inputs
  Int  pcols;
  Int  pver;
  Int  ncol;
  Int  num_msg;
  bool pergro_active;

  Int  *msemax_top_k;
  Real *temperature;
  Real *zmid;
  Real *sp_humidity;

  // Inputs/Outputs
  Int  *msemax_klev;
  Real *mse_max_val;

  // zm_data_find_mse_max( Int  pcols_, 
  //                       Int  pver_,
  //                       Int  ncol_,
  //                       Int  num_msg_,
  //                       bool pergro_active_ ) :
  //   PhysicsTestData(
  //     {
  //       {pcols_},           // dimension group 1: pcols
  //       {pcols_, pver_},    // dimension group 2: pcols x pver
  //       {pcols_}            // dimension group 3: pcols
  //     },
  //     {
  //       {&mse_max_val},                     // Real arrays with dims {pcols}
  //       {&temperature, &zmid, &sp_humidity} // Real arrays with dims {pcols, pver}
  //     },
  //     {
  //       {&msemax_top_k, &msemax_klev}       // Int arrays with dims {pcols}
  //     }
  //   ),
  //   pcols(pcols_), ncol(ncol_), pver(pver_), num_msg(num_msg_), pergro_active(pergro_active_)
  // {}

  // Constructor
  zm_data_find_mse_max(Int pcols_, Int pver_, Int ncol_, Int num_msg_, bool pergro_active_)
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
      pcols(pcols_), pver(pver_), ncol(ncol_), num_msg(num_msg_), pergro_active(pergro_active_)
  {}

  PTD_STD_DEF(zm_data_find_mse_max, 5, pcols, pver, ncol, num_msg, pergro_active);

};

// Glue functions to call fortran from from C++ with the Data struct
void zm_find_mse_max(zm_data_find_mse_max& d);

extern "C" { // _f function decls
}

// void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d);
// extern "C" { // _f function decls
// }

}  // namespace zm
}  // namespace scream

#endif
