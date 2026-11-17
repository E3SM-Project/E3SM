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
#include <limits>

namespace Homme {

struct LimiterFunctor {

  struct Buffers {
    static constexpr int num_3d_scalar_mid_buf = 1;
    ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   buffer1;
  };

  int                 m_np1;
  int                 m_nold;
  int                 m_nstep = -1;
  Real                m_step_start_seconds = -1;
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
#if defined(KOKKOS_ENABLE_CUDA) && !defined(NDEBUG)
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
      , m_dp3d_thresh(params.dp3d_thresh)
      , m_vtheta_thresh(params.vtheta_thresh)
      , m_policy_dp3d_lim (Homme::get_default_team_policy<ExecSpace,TagDp3dLimiter>(m_num_elems))
      , m_tu(m_policy_dp3d_lim)
  {
    m_np1 = -1;
    m_nold = -1;
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

  // tl:     time level to limit.
  // tl_old: time level holding the state before the update that produced tl.
  //         It is used to rebuild theta in the layers that get limited.
  //         If no such state is available, pass tl_old=tl.
  void run (const int& tl, const int& tl_old)
  {
    GPTLstart("caar limiter");
    m_np1 = tl;
    m_nold = tl_old;
    Kokkos::parallel_for("caar loop dp3d limiter", m_policy_dp3d_lim, *this);
    Kokkos::fence();
    GPTLstop("caar limiter");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagDp3dLimiter&, const TeamMember &team) const {
    KernelVariables kv(team, m_tu);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      const bool has_latlon = m_geometry.m_sphere_latlon.size() != 0;
      constexpr Real rad_to_deg = 57.2957795130823208768;
      const Real lat = has_latlon
        ? m_geometry.m_sphere_latlon(kv.ie,igp,jgp,0)*rad_to_deg
        : std::numeric_limits<Real>::quiet_NaN();
      const Real lon = has_latlon
        ? m_geometry.m_sphere_latlon(kv.ie,igp,jgp,1)*rad_to_deg
        : std::numeric_limits<Real>::quiet_NaN();

      const auto& spheremp = m_geometry.m_spheremp(kv.ie,igp,jgp);

      // Check if the minimum dp3d in this column is blow a certain threshold
      auto dp = Homme::subview(m_state.m_dp3d,kv.ie,m_np1,igp,jgp);
      auto& dp0 = m_hvcoord.dp0;
      auto diff = Homme::subview(m_buffers.buffer1,kv.team_idx,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        diff(ilev) = (dp(ilev) - m_dp3d_thresh*dp0(ilev))*spheremp;
      });

      Real min_diff = Kokkos::reduction_identity<Real>::min();
      auto diff_as_real = Homme::viewAsReal(diff);
      auto dp_as_real   = Homme::viewAsReal(dp);
      auto dp0_as_real  = Homme::viewAsReal(dp0);
      Kokkos::Min<Real,ExecSpace> reducer(min_diff);
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                              [&](const int k,Real& result) {
#ifndef HOMMEXX_DISABLE_LIMITER_WARNINGS
        if(diff_as_real(k) < 0){
          Kokkos::printf("WARNING:CAAR: dp3d too small. k=%d, dp3d(k)=%f, dp0=%f, "
                         "lat_deg=%.8f, lon_deg=%.8f, dyn_step=%d, step_start_seconds=%.9f \n",
                         k+1,dp_as_real(k),dp0_as_real(k),lat,lon,m_nstep,m_step_start_seconds);
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

        // In the layers that fell below the threshold, dp can be close to zero
        // or negative, so vtheta_dp/dp is unreliable. Rebuild theta there from
        // the state before the update:
        //    theta = theta_old + (vtheta_dp - theta_old*dp)/max(dp_old,dp_ref)
        // with dp_ref the hybrid coordinate reference thickness for this column.
        // Then clamp it to the range spanned by the nearest layers above and
        // below that did not need limiting (widened by theta_tol).
        Kokkos::single(Kokkos::PerThread(kv.team), [&] () {
          constexpr Real theta_tol = 0.05;
          auto vtheta_as_real     = Homme::viewAsReal(vtheta_dp);
          auto dp_old_as_real     = Homme::viewAsReal(Homme::subview(m_state.m_dp3d,kv.ie,m_nold,igp,jgp));
          auto vtheta_old_as_real = Homme::viewAsReal(Homme::subview(m_state.m_vtheta_dp,kv.ie,m_nold,igp,jgp));
          auto dai_as_real        = Homme::viewAsReal(m_hvcoord.hybrid_ai_delta);
          auto dbi_as_real        = Homme::viewAsReal(m_hvcoord.hybrid_bi_delta);
          const bool have_old = (m_nold != m_np1);

          Real ps = m_hvcoord.hybrid_ai0*m_hvcoord.ps0;
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            ps += dp_as_real(k);
          }

          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            if (diff_as_real(k) < 0) {
              const Real th_raw = vtheta_as_real(k);
              Real th_old = th_raw;
              Real th = th_raw;
              if (have_old && dp_old_as_real(k) > 0) {
                th_old = vtheta_old_as_real(k)/dp_old_as_real(k);
                const Real dp_ref = dai_as_real(k)*m_hvcoord.ps0 + dbi_as_real(k)*ps;
                const Real denom = dp_old_as_real(k) > dp_ref ? dp_old_as_real(k) : dp_ref;
                if (dp_as_real(k) != 0) {
                  th = th_old + dp_as_real(k)*(th_raw - th_old)/denom;
                } else {
                  th = th_old;
                }
              }

              // Nearest layers above (ka) and below (kb) that were not limited
              int ka = -1;
              int kb = -1;
              for (int j=k-1; j>=0; --j) {
                if (diff_as_real(j) >= 0) { ka = j; break; }
              }
              for (int j=k+1; j<NUM_PHYSICAL_LEV; ++j) {
                if (diff_as_real(j) >= 0) { kb = j; break; }
              }
              if (ka>=0 || kb>=0) {
                Real th_lo, th_hi;
                if (ka>=0 && kb>=0) {
                  th_lo = vtheta_as_real(ka) < vtheta_as_real(kb) ? vtheta_as_real(ka) : vtheta_as_real(kb);
                  th_hi = vtheta_as_real(ka) > vtheta_as_real(kb) ? vtheta_as_real(ka) : vtheta_as_real(kb);
                } else if (ka>=0) {
                  th_lo = vtheta_as_real(ka);
                  th_hi = vtheta_as_real(ka);
                } else {
                  th_lo = vtheta_as_real(kb);
                  th_hi = vtheta_as_real(kb);
                }
                th_lo *= (1.0-theta_tol);
                th_hi *= (1.0+theta_tol);
                if (th < th_lo) th = th_lo;
                if (th > th_hi) th = th_hi;
              }
#ifndef HOMMEXX_DISABLE_LIMITER_WARNINGS
              Kokkos::printf("WARNING:CAAR: theta in limited layer. k=%d, theta_raw=%f, theta_old=%f, theta_new=%f, "
                             "lat_deg=%.8f, lon_deg=%.8f, dyn_step=%d, step_start_seconds=%.9f \n",
                             k+1,th_raw,th_old,th,lat,lon,m_nstep,m_step_start_seconds);
#endif
              vtheta_as_real(k) = th;
            }
          }
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
      } // end of min_diff < 0

      // Check if vtheta is too low. If so, take theta from the nearest layers
      // above and below that pass the check (the smaller of the two), instead
      // of resetting it to vtheta_thresh. Only if no layer in the column passes
      // the check do we fall back to vtheta_thresh.
      // The neighbour search is serial, so first check in parallel whether any
      // level needs it; usually none does.
      auto vtheta_dp_as_real = Homme::viewAsReal(vtheta_dp);
      int num_low = 0;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                              [&](const int k, int& count) {
        if ( (vtheta_dp_as_real(k) - m_vtheta_thresh*dp_as_real(k)) < 0) {
          ++count;
        }
      }, num_low);

      if (num_low>0) {
        Kokkos::single(Kokkos::PerThread(kv.team), [&] () {
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            if ( (vtheta_dp_as_real(k) - m_vtheta_thresh*dp_as_real(k)) < 0) {
              int ka = -1;
              int kb = -1;
              for (int j=k-1; j>=0; --j) {
                if ( (vtheta_dp_as_real(j) - m_vtheta_thresh*dp_as_real(j)) >= 0) { ka = j; break; }
              }
              for (int j=k+1; j<NUM_PHYSICAL_LEV; ++j) {
                if ( (vtheta_dp_as_real(j) - m_vtheta_thresh*dp_as_real(j)) >= 0) { kb = j; break; }
              }
              Real th = m_vtheta_thresh;
              if (ka>=0 && kb>=0) {
                const Real th_a = vtheta_dp_as_real(ka)/dp_as_real(ka);
                const Real th_b = vtheta_dp_as_real(kb)/dp_as_real(kb);
                th = th_a < th_b ? th_a : th_b;
              } else if (ka>=0) {
                th = vtheta_dp_as_real(ka)/dp_as_real(ka);
              } else if (kb>=0) {
                th = vtheta_dp_as_real(kb)/dp_as_real(kb);
              }
#ifndef HOMMEXX_DISABLE_LIMITER_WARNINGS
              Kokkos::printf("WARNING:CAAR: k=%d,theta(k)=%f<%f=th_thresh, applying limiter, theta_new=%f, "
                             "lat_deg=%.8f, lon_deg=%.8f, dyn_step=%d, step_start_seconds=%.9f \n",
                             k+1,vtheta_dp_as_real(k)/dp_as_real(k),m_vtheta_thresh,th,
                             lat,lon,m_nstep,m_step_start_seconds);
#endif
              vtheta_dp_as_real(k) = th*dp_as_real(k);
            }
          }
        });
      } // end of num_low > 0
    });
    kv.team_barrier();
  }

};

} // Namespace Homme

#endif // HOMMEXX_LIMITER_FUNCTOR_HPP
