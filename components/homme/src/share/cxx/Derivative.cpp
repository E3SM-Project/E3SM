/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Derivative.hpp"
#include "utilities/TestUtils.hpp"

#include <random>

namespace Homme {

Derivative::Derivative()
    : m_dvv("dvv")
{
  // Nothing to be done here
}

void Derivative::init(CF90Ptr &dvv_ptr) {
  ExecViewManaged<Real[NP][NP]>::HostMirror dvv_host =
      Kokkos::create_mirror_view(m_dvv);

  int k_dvv = 0;
  for (int igp = 0; igp < NP; ++igp) {
    for (int jgp = 0; jgp < NP; ++jgp, ++k_dvv) {
      dvv_host(igp, jgp) = dvv_ptr[k_dvv];
    }
  }

  Kokkos::deep_copy(m_dvv, dvv_host);
}

void Derivative::random_init() {
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> random_dist(16.0, 8192.0);
  genRandArray(m_dvv, engine, random_dist);
}

void Derivative::dvv(Real *dvv_ptr) {
  ExecViewManaged<Real[NP][NP]>::HostMirror dvv_f90(dvv_ptr);
  Kokkos::deep_copy(dvv_f90, m_dvv);
}

} // namespace Homme
