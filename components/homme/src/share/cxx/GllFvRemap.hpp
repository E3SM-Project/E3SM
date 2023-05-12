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

  // Need this only if the object was created before the various other objects
  // were initialized.
  void setup();
  
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

  // The following three routines provide dynamics-physics coupling for the
  // atmosphere model.
  //   Remap dynamics state to physics state.
  void run_dyn_to_fv_phys(const int time_idx,
                          // ps,phis(ie,col)
                          const Phys1T& ps, const Phys1T& phis,
                          // T,omega(ie,col,lev)
                          const Phys2T& T, const Phys2T& omega,
                          // uv(ie, col, 0 or 1, lev)
                          const Phys3T& uv, 
                          // q(ie,col,idx,lev)
                          const Phys3T& q,
                          // Optionally return dp
                          const Phys2T* dp = nullptr);
  //   Remap physics state and tendencies to dynamics state and tendencies.
  void run_fv_phys_to_dyn(const int time_idx, const CPhys2T& T, const CPhys3T& uv,
                          const CPhys3T& q);
  //   DSS the remapped dynamics tendencies and state. Call this after
  //   run_fv_phys_to_dyn if the dynamics-physics coupler does not already
  //   provide it.
  void run_fv_phys_to_dyn_dss();

  // Remap nq tracers, with nq <= qsize. This is a convenience routine for use
  // in, e.g., model initialization; it is not a core part of the API for the
  // dynamics-physics coupler. time_idx is used to access dp-related data.
  void remap_tracer_dyn_to_fv_phys(const int time_idx, const int nq,
                                   // q_dyn,fv(ie,col,idx,lev), idx = 0:nq-1
                                   //   dyn: col = 0:np^2-1
                                   //    fv: col = 0:nf^2-1
                                   const CPhys3T& q_dyn, const Phys3T& q_fv);

private:
  std::unique_ptr<GllFvRemapImpl> m_impl;
};

extern "C" void
init_gllfvremap_c(int nelemd, int np, int nf, int nf_max,
                  const bool theta_hydrostatic_mode,
                  CF90Ptr fv_metdet, CF90Ptr g2f_remapd,
                  CF90Ptr f2g_remapd, CF90Ptr D_f, CF90Ptr Dinv_f);

} // Namespace Homme

#endif // HOMMEXX_GLLFVREMAP_HPP
