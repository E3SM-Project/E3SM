/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYBRID_V_COORD_HPP
#define HOMMEXX_HYBRID_V_COORD_HPP

#include "Types.hpp"

namespace Homme
{

struct HybridVCoord
{
  HybridVCoord () = default;

  // This method should only be called from the host
  void init(const Real ps0_in,
            CRCPtr hybrid_am_ptr,
            CRCPtr hybrid_ai_ptr,
            CRCPtr hybrid_bm_ptr,
            CRCPtr hybrid_bi_ptr);

  void random_init(int seed);
  void compute_deltas ();

  Real ps0;
  Real hybrid_ai0;

  // hybrid ai
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> hybrid_ai;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_ai_delta;

  // hybrid bi
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> hybrid_bi;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_bi_delta;

  ExecViewManaged<Scalar[NUM_LEV]> dp0;
};

} // namespace Homme

#endif // HOMMEXX_HYBRID_V_COORD_HPP
