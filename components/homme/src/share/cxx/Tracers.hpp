/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_TRACERS_HPP
#define HOMMEXX_TRACERS_HPP

#include "Types.hpp"

namespace Homme {

struct Tracers {
  Tracers() : m_inited(false) {}
  Tracers(const int num_elems, const int num_tracers);

  void init (const int num_elems, const int num_tracers);

  void randomize(const int seed);

  void pull_qdp(CF90Ptr &state_qdp);
  void push_qdp(F90Ptr &state_qdp) const;

  KOKKOS_INLINE_FUNCTION
  int num_tracers() const {
    return nt;
  }

  int num_elems () const {
    return ne;
  }

  bool inited () const { return m_inited; }

  ExecViewManaged<Scalar*[Q_NUM_TIME_LEVELS][QSIZE_D][NP][NP][NUM_LEV]> qdp;
  ExecViewManaged<Scalar*[QSIZE_D][NP][NP][NUM_LEV]>                    qtens_biharmonic; // Also doubles as just qtens.
  ExecViewManaged<Scalar*[QSIZE_D][2][NUM_LEV]>                         qlim;
  ExecViewManaged<Scalar**[NP][NP][NUM_LEV]>                    Q;
  ExecViewManaged<Scalar**[NP][NP][NUM_LEV]>                    fq;

private:
  int nt;
  int ne;
  bool m_inited;
};

} // namespace Homme

#endif // HOMMEXX_TRACERS_HPP
