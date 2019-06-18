/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ReferenceElement.hpp"
#include "utilities/TestUtils.hpp"

#include <random>

namespace Homme {

ReferenceElement::ReferenceElement()
  : m_deriv  ("deriv")
  , m_mass   ("mass")
  , m_inited (false)
{
  // Nothing to be done here
}

void ReferenceElement::init(CF90Ptr &deriv_ptr, CF90Ptr& mass_ptr) {
  init_deriv(deriv_ptr);
  init_mass(mass_ptr);

  m_inited  = true;
}

void ReferenceElement::init_deriv(CF90Ptr &deriv_ptr) {
  ExecViewManaged<Real[NP][NP]>::HostMirror deriv_host =
      Kokkos::create_mirror_view(m_deriv);

  int k_deriv = 0;
  for (int igp = 0; igp < NP; ++igp) {
    for (int jgp = 0; jgp < NP; ++jgp, ++k_deriv) {
      deriv_host(igp, jgp) = deriv_ptr[k_deriv];
    }
  }

  Kokkos::deep_copy(m_deriv, deriv_host);
}

void ReferenceElement::init_mass(CF90Ptr &mass_ptr) {
  ExecViewManaged<Real[NP][NP]>::HostMirror mass_host =
      Kokkos::create_mirror_view(m_mass);

  int k_mass = 0;
  for (int igp = 0; igp < NP; ++igp) {
    for (int jgp = 0; jgp < NP; ++jgp, ++k_mass) {
      mass_host(igp, jgp) = mass_ptr[k_mass];
    }
  }

  Kokkos::deep_copy(m_mass, mass_host);
}

void ReferenceElement::random_init(const int seed) {
  std::mt19937_64 engine(seed);
  using urd = std::uniform_real_distribution<Real>;

  genRandArray(m_deriv, engine, urd(16.0,8192.0));
  genRandArray(m_mass, engine, urd(0.1,1000));

  m_inited = true;
}

} // namespace Homme
