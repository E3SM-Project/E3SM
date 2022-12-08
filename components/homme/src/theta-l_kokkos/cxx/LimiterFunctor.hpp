/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_LIMITER_FUNCTOR_HPP
#define HOMMEXX_LIMITER_FUNCTOR_HPP

#include "Types.hpp"
#include "Elements.hpp"
#include "ColumnOps.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "ReferenceElement.hpp"
#include "SimulationParams.hpp"
#include "kokkos_utils.hpp"

#include "utilities/SubviewUtils.hpp"
#include "utilities/ViewUtils.hpp"

#include "profiling.hpp"
#include "ErrorDefs.hpp"

#include <assert.h>

namespace Homme {

struct LimiterFunctor {

  struct Buffers {
    static constexpr int num_3d_scalar_mid_buf = 1;
    ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   buffer1;
  };

  int                 m_np1;
  const int           m_num_elems;
  const bool          m_theta_hydrostatic_mode;

  HybridVCoord        m_hvcoord;
  ElementsState       m_state;
  Buffers             m_buffers;
  ElementsGeometry    m_geometry;
  
  double              m_dp3d_thresh;
  double              m_vtheta_thresh;           

  struct TagDp3dLimiter {};

  // Policies
#ifndef NDEBUG
  template<typename Tag>
  using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace,Kokkos::LaunchBounds<512,1>,Tag>;
#else
  template<typename Tag>
  using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace,Tag>;
#endif

  TeamPolicyType<TagDp3dLimiter>   m_policy_dp3d_lim;

  TeamUtils<ExecSpace> m_tu;

  LimiterFunctor(const Elements &elements, const HybridVCoord &hvcoord, const SimulationParams& params)
      : m_num_elems(elements.num_elems())
      , m_theta_hydrostatic_mode(params.theta_hydrostatic_mode)
      , m_hvcoord(hvcoord)
      , m_state(elements.m_state)
      , m_geometry(elements.m_geometry)
      , m_policy_dp3d_lim (Homme::get_default_team_policy<ExecSpace,TagDp3dLimiter>(m_num_elems))
      , m_tu(m_policy_dp3d_lim)
      , m_dp3d_thresh(params.dp3d_thresh)
      , m_vtheta_thresh(params.vtheta_thresh)
  {
    m_np1 = -1;
  }

  int requested_buffer_size () const {
    // Ask the buffers manager to allocate enough buffers
    const int nslots = m_tu.get_num_ws_slots();

    int num_scalar_mid_buf = Buffers::num_3d_scalar_mid_buf;

    return num_scalar_mid_buf  *NP*NP*NUM_LEV  *VECTOR_SIZE*nslots ;
  }

  void init_buffers (const FunctorsBuffersManager& fbm) {
    Errors::runtime_check(fbm.allocated_size()>=requested_buffer_size(), "Error! Buffers size not sufficient.\n");

    Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
    const int nslots = m_tu.get_num_ws_slots();

    // Midpoints scalars
    m_buffers.buffer1       = decltype(m_buffers.buffer1)(mem,nslots);
    mem += m_buffers.buffer1.size();

    assert ((reinterpret_cast<Real*>(mem) - fbm.get_memory())==requested_buffer_size());
  }

  void run (const int& tl)
  {
    profiling_resume();

    GPTLstart("caar limiter");
    m_np1 = tl;
    Kokkos::parallel_for("caar loop dp3d limiter", m_policy_dp3d_lim, *this);
    Kokkos::fence();
    GPTLstop("caar limiter");

    profiling_pause();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagDp3dLimiter&, const TeamMember &team) const {
    KernelVariables kv(team, m_tu);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      const auto& spheremp = m_geometry.m_spheremp(kv.ie,igp,jgp);

      // Check if the minimum dp3d in this column is blow a certain threshold
      auto dp = Homme::subview(m_state.m_dp3d,kv.ie,m_np1,igp,jgp);
      auto& dp0 = m_hvcoord.dp0;
      auto diff = Homme::subview(m_buffers.buffer1,kv.team_idx,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        diff(ilev) = (dp(ilev) - m_dp3d_thresh*dp0(ilev))*spheremp;
      });

      kv.team_barrier();

      Real min_diff = Kokkos::reduction_identity<Real>::min();
      auto diff_as_real = Homme::viewAsReal(diff);
      auto dp_as_real   = Homme::viewAsReal(dp);
      auto dp0_as_real  = Homme::viewAsReal(dp0);
      Kokkos::Min<Real,ExecSpace> reducer(min_diff);
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                              [&](const int k,Real& result) {
#ifndef HOMMEXX_BFB_TESTING
        if(diff_as_real(k) < 0){
          printf("WARNING:CAAR: dp3d too small. k=%d, dp3d(k)=%f, dp0=%f \n",
           k+1,dp_as_real(k),dp0_as_real(k));
        }
#endif
        result = result<=diff_as_real(k) ? result : diff_as_real(k);
      }, reducer);

      auto vtheta_dp = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);

      if (min_diff<0) {
        // Compute vtheta = vtheta_dp/dp
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                             [&](const int ilev) {
          vtheta_dp(ilev) /= dp(ilev);
        });

        // Gotta apply vertical mixing, to prevent levels from getting too thin.
        Real mass = 0.0;
        ColumnOps::column_reduction<NUM_PHYSICAL_LEV>(kv.team,diff,mass);

        if (mass<0) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                               [&](const int ilev) {
            diff(ilev) *= -1.0;
          });
        }

        kv.team_barrier();

        // This loop must be done over physical levels, unless we implement
        // masks, like it has been done in the E3SM/scream project
        Real mass_new = 0.0;
        Dispatch<>::parallel_reduce(kv.team,
                                    Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                                    [&](const int k, Real& accum) {
          auto& val = diff_as_real(k);
          val = (val<0 ? 0.0 : val);
          accum += val;
        }, mass_new);

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                             [&](const int ilev) {
          if (mass_new>0) {
            diff(ilev) *= fabs(mass)/mass_new;
          }
          if (mass<0) {
            diff(ilev) *= -1.0;
          }

          dp(ilev) = diff(ilev)/spheremp + m_dp3d_thresh*dp0(ilev);
          vtheta_dp(ilev) *= dp(ilev);
        });
      } //end of min_diff < 0

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        // Check if vtheta is too low
        // Note: another place where scream's masks could help
        for (int ivec=0; ivec<VECTOR_SIZE; ++ivec) {
          if ( (vtheta_dp(ilev)[ivec] - m_vtheta_thresh*dp(ilev)[ivec]) < 0) {
#ifndef HOMMEXX_BFB_TESTING
             printf("WARNING:CAAR: k=%d,theta(k)=%f<%f=th_thresh, applying limiter \n",
               ilev*VECTOR_SIZE+ivec+1,vtheta_dp(ilev)[ivec]/dp(ilev)[ivec],m_vtheta_thresh);
#endif
             vtheta_dp(ilev)[ivec]=m_vtheta_thresh*dp(ilev)[ivec];
          }
        }
      });
    });
    kv.team_barrier();
  }

};

} // Namespace Homme

#endif // HOMMEXX_LIMITER_FUNCTOR_HPP
