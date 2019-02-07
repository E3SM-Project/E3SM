/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_HPP
#define HOMMEXX_ELEMENTS_HPP

#include "Types.hpp"

#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsBuffers.hpp"
#include "ElementsDerivedState.hpp"
#include "ElementsForcing.hpp"

namespace Homme {

class HybridVCoord;

/*
 *  A class to store all the elements stuff.
 *
 *  This class contains views that are element-dependent.
 *  The views are grouped into four sub structures:
 *    - ElementsGeometry:     2d views that depend only on geometry (and quadrature)
 *    - ElementsState:        time-level dependent views of the dynamics states (u,v,t,dp,ps_v)
 *    - ElementsDerivedState: diagnostics variables and storage variables for dynamics subcycling
 *    - ElementsBuffers:      buffers to be used for temporary variables inside kernels
 */

class Elements {
public:

  ElementsGeometry      m_geometry;
  ElementsState         m_state;
  ElementsDerivedState  m_derived;
  ElementsBuffers       m_buffers;
  ElementsForcing       m_forcing;

  Elements () : m_inited(false) {}

  int num_elems () const { return m_num_elems; }

  void init (int num_elems, bool consthv);
  void random_init (int num_elems, Real max_pressure = 1.0);

  bool inited () const { return m_inited; }

protected:
  int  m_num_elems;
  bool m_inited;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_HPP
