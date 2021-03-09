#include "ElementsForcing.hpp"

namespace Homme {

void ElementsForcing::init (const int num_elems) {
  m_num_elems = num_elems;

  m_fm = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("F_Momentum",    m_num_elems);
  m_ft = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]>("F_Temperature", m_num_elems);
}

} // namespace Homme
