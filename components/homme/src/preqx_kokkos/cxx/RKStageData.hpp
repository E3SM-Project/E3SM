#ifndef HOMMEXX_RK_STAGE_DATA_HPP
#define HOMMEXX_RK_STAGE_DATA_HPP

#include <Types.hpp>

namespace Homme {

struct RKStageData {

  RKStageData () = default;

  RKStageData (const int nm1_in, const int n0_in, const int np1_in, const int n0_qdp_in,
               const Real dt_in, const Real eta_ave_w_in)
   : nm1    (nm1_in)
   , n0     (n0_in)
   , np1    (np1_in)
   , n0_qdp (n0_qdp_in)
   , dt     (dt_in)
   , eta_ave_w (eta_ave_w_in)
  {}

  int     nm1;
  int     n0;
  int     np1;
  int     n0_qdp;

  Real    dt;
  Real    eta_ave_w;
};

} // namespace Homme

#endif // HOMMEXX_RK_STAGE_DATA_HPP
