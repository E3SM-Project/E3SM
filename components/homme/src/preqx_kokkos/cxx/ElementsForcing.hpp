#ifndef HOMMEXX_ELEMENTS_FORCING_HPP
#define HOMMEXX_ELEMENTS_FORCING_HPP

#include <Types.hpp>

namespace Homme {

class ElementsForcing {
public:

  // Per Element Forcings
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> m_fm;  // Momentum (? units are wrong in apply_cam_forcing...) forcing
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>    m_ft;  // Temperature forcing

  ElementsForcing() = default;

  void init(const int num_elems);

private:
  int m_num_elems;
};

} // namespace Homme

#endif // HOMMEXX_ELEMENTS_FORCING_HPP
