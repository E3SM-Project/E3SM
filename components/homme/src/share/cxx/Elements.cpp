/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Elements.hpp"

#include "utilities/SubviewUtils.hpp"

#include <limits>

namespace Homme {

void Elements::init(const int num_elems, const bool consthv) {
  m_num_elems = num_elems;

  m_geometry.init(num_elems,consthv);
  m_state.init(num_elems);
  m_derived.init(num_elems);
  m_forcing.init(num_elems);

  m_inited = true;
}

//test for tensor hv is needed
void Elements::random_init(const int num_elems, const int seed, const Real max_pressure) {
  m_num_elems = num_elems;

  m_geometry.random_init(num_elems,seed);
  m_state.random_init(num_elems,seed,max_pressure);

  Real dp3d_min = std::numeric_limits<Real>::max();
  for (int ie = 0; ie < m_num_elems; ++ie) {
    // Because this constraint is difficult to satisfy for all of the tensors,
    // incrementally generate the view
    for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
      for (int igp = 0; igp < NP; ++igp) {
        for (int jgp = 0; jgp < NP; ++jgp) {
          ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp3d =
              Homme::subview(m_state.m_dp3d, ie, tl, igp, jgp);
          auto h_dp3d = Kokkos::create_mirror_view(pt_dp3d);
          Kokkos::deep_copy(h_dp3d,pt_dp3d);
          for (int ilev=0; ilev<NUM_LEV; ++ilev) {
            for (int iv=0; iv<VECTOR_SIZE; ++iv) {
              dp3d_min = std::min(dp3d_min,h_dp3d(ilev)[iv]);
            }
          }
        }
      }
    }
  }

  m_derived.random_init(num_elems,seed,dp3d_min);

  m_inited = true;
}

} // namespace Homme
