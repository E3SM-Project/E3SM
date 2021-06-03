/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_GLLFVREMAP_HPP
#define HOMMEXX_GLLFVREMAP_HPP

#include "Types.hpp"
#include <memory>

namespace Homme {

class FunctorsBuffersManager;
class SimulationParams;
class GllFvRemapImpl;

class GllFvRemap {
public:
  GllFvRemap();
  GllFvRemap(const GllFvRemap &) = delete;
  GllFvRemap &operator=(const GllFvRemap &) = delete;

  ~GllFvRemap();

  void reset(const SimulationParams& params);

  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  typedef ExecViewUnmanaged<Real**>   Phys1T; // ie, col
  typedef ExecViewUnmanaged<Real***>  Phys2T; // ie, col, lev
  typedef ExecViewUnmanaged<Real****> Phys3T; // ie, col, idx, lev
  typedef Phys2T::const_type CPhys2T;
  typedef Phys3T::const_type CPhys3T;

  void init_data(const int nf, const int nf_max, bool theta_hydrostatic_mode,
                 const Real* fv_metdet, const Real* g2f_remapd, const Real* f2g_remapd,
                 const Real* D_f, const Real* Dinv_f);

  void run_dyn_to_fv_phys(const int time_idx,
                          // ps,phis(ie,col)
                          const Phys1T& ps, const Phys1T& phis,
                          // T,omega(ie,col,lev)
                          const Phys2T& T, const Phys2T& omega,
                          // uv(ie, col, 0 or 1, lev)
                          const Phys3T& uv, 
                          // q(ie,col,idx,lev)
                          const Phys3T& q);
  void run_fv_phys_to_dyn(const int time_idx, const Real dt,
                          const CPhys2T& T, const CPhys3T& uv, const CPhys3T& q);
  void run_fv_phys_to_dyn_dss();

private:
  std::unique_ptr<GllFvRemapImpl> m_impl;
};

extern "C" void
init_gllfvremap_c(int nelemd, int np, int nf, int nf_max, int ftype,
                  const bool theta_hydrostatic_mode,
                  CF90Ptr fv_metdet, CF90Ptr g2f_remapd,
                  CF90Ptr f2g_remapd, CF90Ptr D_f, CF90Ptr Dinv_f);

} // Namespace Homme

#endif // HOMMEXX_GLLFVREMAP_HPP
