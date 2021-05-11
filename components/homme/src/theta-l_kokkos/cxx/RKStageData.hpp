#ifndef HOMMEXX_RK_STAGE_DATA_HPP
#define HOMMEXX_RK_STAGE_DATA_HPP

#include <Types.hpp>

namespace Homme {

struct RKStageData {

  RKStageData () = default;

  RKStageData (const int nm1_in, const int n0_in, const int np1_in, const int n0_qdp_in,
               const Real dt_in, const Real eta_ave_w_in,
               const Real scale1_in = 1.0, const Real scale2_in = 1.0, const Real scale3_in = 1.0)
   : nm1    (nm1_in)
   , n0     (n0_in)
   , np1    (np1_in)
   , n0_qdp (n0_qdp_in)
   , dt     (dt_in)
   , eta_ave_w (eta_ave_w_in)
   , scale1 (scale1_in)
   , scale2 (scale2_in)
   , scale3 (scale3_in)
  {}

  int     nm1;
  int     n0;
  int     np1;
  int     n0_qdp;

  Real    dt;
  Real    eta_ave_w;

  Real    scale1;
  Real    scale2;
  Real    scale3;
};

} // namespace Homme

#endif // HOMMEXX_RK_STAGE_DATA_HPP
