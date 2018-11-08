/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_DERIVATIVE_HPP
#define HOMMEXX_DERIVATIVE_HPP

#include "Types.hpp"

namespace Homme {

class Derivative {
public:
  Derivative();

  void init(CF90Ptr &dvv);

  void random_init();

  void dvv(Real *dvv);

  KOKKOS_INLINE_FUNCTION
  ExecViewUnmanaged<const Real[NP][NP]> get_dvv() const { return m_dvv; }

private:
  ExecViewManaged<Real[NP][NP]> m_dvv;
};

} // namespace Homme

#endif // HOMMEXX_DERIVATIVE_HPP
