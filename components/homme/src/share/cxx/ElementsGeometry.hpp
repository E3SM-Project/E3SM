/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_GEOMETRY_HPP
#define HOMMEXX_ELEMENTS_GEOMETRY_HPP

#include "Types.hpp"

namespace Homme {

/*
 *  A structure to hold views depending on the geometry/quadrature
 *
 *  This structure contains 2d views that depend only on the mesh 
 *  (and possibly on quadrature). We also include the geopotential
 *  and the tensor viscosity, which also depend on some simulation
 *  parameter (e.g., how much the geopotential must be smoothed,
 *  or the hyperviscosity scaling).
 */
class ElementsGeometry {
public:
  // Coriolis term
  ExecViewManaged<Real * [NP][NP]> m_fcor;

  // Quadrature weights and metric tensor
  ExecViewManaged<Real * [NP][NP]>        m_spheremp;
  ExecViewManaged<Real * [NP][NP]>        m_rspheremp;
  ExecViewManaged<Real * [2][2][NP][NP]>  m_metinv;
  ExecViewManaged<Real * [NP][NP]>        m_metdet;
  ExecViewManaged<Real * [2][2][NP][NP]>  m_tensorvisc;
  ExecViewManaged<Real * [2][3][NP][NP]>  m_vec_sph2cart;

  // Prescribed surface geopotential height at eta = 1
  ExecViewManaged<Real *    [NP][NP]> m_phis;
  ExecViewManaged<Real * [2][NP][NP]> m_gradphis;

  // D (map for covariant coordinates) and D^{-1}
  ExecViewManaged<Real * [2][2][NP][NP]> m_d;
  ExecViewManaged<Real * [2][2][NP][NP]> m_dinv;

  ElementsGeometry() : m_num_elems(0) {}

  void init (const int num_elems, const bool consthv, const bool alloc_gradphis);

  void randomize (const int seed);

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

  // Fill the exec space views with data coming from F90 pointers
  void set_elem_data (const int ie,
                      CF90Ptr& D, CF90Ptr& Dinv, CF90Ptr& fcor,
                      CF90Ptr& spheremp, CF90Ptr& rspheremp,
                      CF90Ptr& metdet, CF90Ptr& metinv,
                      CF90Ptr& tensorvisc,
                      CF90Ptr& vec_sph2cart, const bool consthv);

  void set_phis (const int ie, CF90Ptr& phis);
private:
  bool m_consthv;
  int m_num_elems;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_GEOMETRY_HPP
