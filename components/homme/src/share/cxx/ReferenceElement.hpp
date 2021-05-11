/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_REFERENCE_ELEMENT_HPP
#define HOMMEXX_REFERENCE_ELEMENT_HPP

#include "Types.hpp"

namespace Homme {

 /*
  *  ReferenceElement data
  *
  *  This small class contains some basic data on the reference element
  *  for the Spectral Element Method. In particular:
  *    - m_deriv: contains the (non-symmetric) pseudospectral derivative matrix D,
  *      with D(i,j) being the derivative of basis function j at node/quad_point i.
  *    - m_mass: contains the (symmetric) mass matrix M
  */

class ReferenceElement {
public:
  using deriv_type = ExecViewManaged<const Real[NP][NP]>;
  using mass_type  = ExecViewManaged<const Real[NP][NP]>;

  ReferenceElement();

  void init(CF90Ptr& deriv, CF90Ptr& mass);
  void init_deriv(CF90Ptr& deriv);
  void init_mass(CF90Ptr& mass);

  void random_init(const int seed);

  KOKKOS_INLINE_FUNCTION
  ExecViewManaged<const Real[NP][NP]> get_deriv() const { return m_deriv; }

  KOKKOS_INLINE_FUNCTION
  ExecViewManaged<const Real[NP][NP]> get_mass() const { return m_mass; }

  bool inited () const { return m_inited; }

private:

  // Non-const, cause we have to fill the views at some point
  deriv_type::non_const_type  m_deriv;
  mass_type::non_const_type   m_mass;

  bool m_inited;
};

} // namespace Homme

#endif // HOMMEXX_REFERENCE_ELEMENT_HPP
