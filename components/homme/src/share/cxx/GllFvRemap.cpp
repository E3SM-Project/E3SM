/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "GllFvRemap.hpp"
#include "GllFvRemapImpl.hpp"
#include "Context.hpp"
#include "ErrorDefs.hpp"
#include "profiling.hpp"

#include <assert.h>
#include <type_traits>

namespace Homme {

void init_gllfvremap_c (int nelemd, int np, int nf, int nf_max,
                        bool theta_hydrostatic_mode,
                        CF90Ptr fv_metdet, CF90Ptr g2f_remapd,
                        CF90Ptr f2g_remapd, CF90Ptr D_f, CF90Ptr Dinv_f) {
  auto& c = Context::singleton();
  auto& g = c.get<GllFvRemap>();
  g.init_data(nf, nf_max, theta_hydrostatic_mode, fv_metdet, g2f_remapd,
              f2g_remapd, D_f, Dinv_f);
}

GllFvRemap::GllFvRemap () {
  m_impl.reset(new GllFvRemapImpl());
}

void GllFvRemap::setup () {
  m_impl->setup();
}

void GllFvRemap::reset (const SimulationParams& params) {
  m_impl->reset(params);
}

GllFvRemap::~GllFvRemap () = default;

int GllFvRemap::requested_buffer_size () const {
  return m_impl->requested_buffer_size();
}

void GllFvRemap::init_buffers (const FunctorsBuffersManager& fbm) {
  m_impl->init_buffers(fbm);
}

void GllFvRemap::init_boundary_exchanges () {
  m_impl->init_boundary_exchanges();
}

void GllFvRemap
::init_data (const int nf, const int nf_max, bool theta_hydrostatic_mode,
             const Real* fv_metdet, const Real* g2f_remapd, const Real* f2g_remapd,
             const Real* D_f, const Real* Dinv_f) {
  m_impl->init_data(nf, nf_max, theta_hydrostatic_mode, fv_metdet,
                    g2f_remapd, f2g_remapd, D_f, Dinv_f);
}

void GllFvRemap
::run_dyn_to_fv_phys (const int time_idx, const Phys1T& ps, const Phys1T& phis,
                      const Phys2T& T, const Phys2T& omega, const Phys3T& uv,
                      const Phys3T& q, const Phys2T* dp) {
  m_impl->run_dyn_to_fv_phys(time_idx, ps, phis, T, omega, uv, q, dp);
}

void GllFvRemap
::run_fv_phys_to_dyn (const int time_idx, const CPhys2T& T, const CPhys3T& uv,
                      const CPhys3T& q) {
  m_impl->run_fv_phys_to_dyn(time_idx, T, uv, q);
}

void GllFvRemap::run_fv_phys_to_dyn_dss () { m_impl->run_fv_phys_to_dyn_dss(); }

void GllFvRemap
::remap_tracer_dyn_to_fv_phys (const int time_idx, const int nq,
                               const CPhys3T& q_dyn, const Phys3T& q_fv) {
  m_impl->remap_tracer_dyn_to_fv_phys(time_idx, nq, q_dyn, q_fv);
}

} // Namespace Homme

