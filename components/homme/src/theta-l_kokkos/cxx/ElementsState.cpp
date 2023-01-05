/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ElementsState.hpp"
#include "ElementOps.hpp"
#include "EquationOfState.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "HybridVCoord.hpp"
#include "ElementsGeometry.hpp"
#include "Context.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/Comm.hpp"

#include <limits>
#include <random>
#include <assert.h>

namespace Homme {

void RefStates::init(const int num_elems) {
  dp_ref = decltype(dp_ref)("dp_ref",num_elems);
  phi_i_ref = decltype(phi_i_ref)("phi_i_ref",num_elems);
  theta_ref = decltype(theta_ref)("theta_ref",num_elems);

  m_num_elems = num_elems;

  m_policy = get_default_team_policy<ExecSpace>(num_elems);
  m_tu     = TeamUtils<ExecSpace>(m_policy);
}

void ElementsState::init(const int num_elems) {
  m_num_elems = num_elems;

  m_v         = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV  ]>("Horizontal velocity", num_elems);
  m_w_i       = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]>("Vertical velocity at interfaces", num_elems);
  m_vtheta_dp = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]>("Virtual potential temperature", num_elems);
  m_phinh_i   = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]>("Geopotential at interfaces", num_elems);
  m_dp3d      = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]>("Delta p at levels", num_elems);

  m_ps_v = ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]>("PS_V", num_elems);

  m_ref_states.init(num_elems);

  m_policy = get_default_team_policy<ExecSpace>(m_num_elems*NUM_TIME_LEVELS);
  m_tu     = TeamUtils<ExecSpace>(m_policy);
}

void ElementsState::randomize(const int seed) {
  randomize(seed,1.0);
}

void ElementsState::randomize(const int seed, const Real max_pressure) {
  randomize(seed,max_pressure,max_pressure/100,0.0);
}

void ElementsState::randomize(const int seed,
                              const Real max_pressure,
                              const Real ps0, const Real hyai0,
                              const ExecViewUnmanaged<const Real*[NP][NP]>& phis) {
  randomize(seed,max_pressure,ps0,hyai0);

  // Re-do phinh so it satisfies phinh_i(bottom)=phis

  // Sanity check
  assert(phis.extent_int(0)==m_num_elems);

  std::mt19937_64 engine(seed);

  // Note: to avoid errors in the equation of state, we need phi to be increasing.
  //       Rather than using a constraint (which may call the function many times,
  //       we simply ask that there are no duplicates, then we sort it later.
  auto sort_and_chek = [](const ExecViewManaged<Real[NUM_PHYSICAL_LEV]>::HostMirror v)->bool {
    Real* start = reinterpret_cast<Real*>(v.data());
    Real* end   = reinterpret_cast<Real*>(v.data()) + NUM_PHYSICAL_LEV;
    std::sort(start,end);
    std::reverse(start,end);
    auto it = std::unique(start,end);
    return it==end;
  };

  auto h_phis = Kokkos::create_mirror_view(phis);
  Kokkos::deep_copy(h_phis,phis);
  for (int ie=0; ie<m_num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++ jgp) {
        const Real phis_ij = h_phis(ie,igp,jgp);
        // Ensure generated values are larger than phis
        std::uniform_real_distribution<Real> random_dist(1.001*phis_ij,100.0*phis_ij);
        for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
          // Get column
          auto phi_col = Homme::viewAsReal(Homme::subview(m_phinh_i,ie,itl,igp,jgp));

          // Stuff phis at the bottom
          Kokkos::deep_copy(Kokkos::subview(phi_col,NUM_PHYSICAL_LEV),phis_ij);

          // Generate values except at bottom
          ExecViewUnmanaged<Real[NUM_PHYSICAL_LEV]> phi_no_bottom(phi_col.data());
          genRandArray(phi_no_bottom,engine,random_dist,sort_and_chek);
        }
      }
    }
  }
}

void ElementsState::randomize(const int seed,
                              const Real max_pressure,
                              const Real ps0,
                              const Real hyai0) {
  // Check elements were inited
  assert (m_num_elems>0);

  // Check data makes sense
  assert (max_pressure>ps0);
  assert (ps0>0);
  assert (hyai0>=0);

  // Arbitrary minimum value to generate
  constexpr const Real min_value = 0.015625;

  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);
  std::uniform_real_distribution<Real> pdf_vtheta_dp(100.0, 1000.0);

  genRandArray(m_v,         engine, random_dist);
  genRandArray(m_w_i,       engine, random_dist);
  genRandArray(m_vtheta_dp, engine, pdf_vtheta_dp);
  // Note: to avoid errors in the equation of state, we need phi to be increasing.
  //       Rather than using a constraint (which may call the function many times,
  //       we simply ask that there are no duplicates, then we sort it later.
  auto sort_and_chek = [](const ExecViewManaged<Scalar[NUM_LEV_P]>::HostMirror v)->bool {
    Real* start = reinterpret_cast<Real*>(v.data());
    Real* end   = reinterpret_cast<Real*>(v.data()) + NUM_LEV_P*VECTOR_SIZE;
    std::sort(start,end);
    std::reverse(start,end);
    auto it = std::unique(start,end);
    return it==end;
  };
  for (int ie=0; ie<m_num_elems; ++ie) {
    for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++ jgp) {
          genRandArray(Homme::subview(m_phinh_i,ie,itl,igp,jgp),engine,random_dist,sort_and_chek);
        }
      }
    }
  }

  // This ensures the pressure in a single column is monotonically increasing
  // and has fixed upper and lower values
  const auto make_pressure_partition = [=](
      ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp) {

    auto h_pt_dp = Kokkos::create_mirror_view(pt_dp);
    Kokkos::deep_copy(h_pt_dp,pt_dp);
    Real* data     = reinterpret_cast<Real*>(h_pt_dp.data());
    Real* data_end = data + NUM_PHYSICAL_LEV;

    Real p[NUM_INTERFACE_LEV];
    Real* p_start = &p[0];
    Real* p_end   = p_start+NUM_INTERFACE_LEV;

    for (int i=0; i<NUM_PHYSICAL_LEV; ++i) {
      p[i+1] = data[i];
    }
    p[0] = ps0*hyai0;
    p[NUM_INTERFACE_LEV-1] = max_pressure;

    // Put in monotonic order
    std::sort(p_start, p_end);

    // Check for no repetitions
    if (std::unique(p_start,p_end)!=p_end) {
      return false;
    }

    // Compute dp from p (we assume p(last interface)=max_pressure)
    for (int i=0; i<NUM_PHYSICAL_LEV; ++i) {
      data[i] = p[i+1]-p[i];
    }

    // Check that dp>=dp_min
    const Real min_dp = std::numeric_limits<Real>::epsilon()*1000;
    for (auto it=data; it!=data_end; ++it) {
      if (*it < min_dp) {
        return false;
      }
    }

    // Fill remainder of last vector pack with quiet nan's
    Real* alloc_end = data+NUM_LEV*VECTOR_SIZE;
    for (auto it=data_end; it!=alloc_end; ++it) {
      *it = std::numeric_limits<Real>::quiet_NaN();
    }

    Kokkos::deep_copy(pt_dp,h_pt_dp);

    return true;
  };

  std::uniform_real_distribution<Real> pressure_pdf(min_value, max_pressure);

  for (int ie = 0; ie < m_num_elems; ++ie) {
    // Because this constraint is difficult to satisfy for all of the tensors,
    // incrementally generate the view
    for (int igp = 0; igp < NP; ++igp) {
      for (int jgp = 0; jgp < NP; ++jgp) {
        for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
          ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp3d =
              Homme::subview(m_dp3d, ie, tl, igp, jgp);
          do {
            genRandArray(pt_dp3d, engine, pressure_pdf);
          } while (make_pressure_partition(pt_dp3d)==false);
        }
      }
    }
  }

  // Generate ps_v so that it is equal to sum(dp3d).
  HybridVCoord hvcoord;
  hvcoord.ps0 = ps0;
  hvcoord.hybrid_ai0 = hyai0;
  hvcoord.m_inited = true;
  auto dp = m_dp3d;
  auto ps = m_ps_v;
  const auto tu = m_tu;
  Kokkos::parallel_for(m_policy, KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team, tu);
    const int ie = kv.ie / NUM_TIME_LEVELS;
    const int tl = kv.ie % NUM_TIME_LEVELS;
    hvcoord.compute_ps_ref_from_dp(kv,Homme::subview(dp,ie,tl),
                                      Homme::subview(ps,ie,tl));
  });
  Kokkos::fence();
}

void ElementsState::pull_from_f90_pointers (CF90Ptr& state_v,         CF90Ptr& state_w_i,
                                            CF90Ptr& state_vtheta_dp, CF90Ptr& state_phinh_i,
                                            CF90Ptr& state_dp3d,      CF90Ptr& state_ps_v) {
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ][2][NP][NP]> state_v_f90         (state_v,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_w_i_f90       (state_w_i,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_vtheta_dp_f90 (state_vtheta_dp,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_phinh_i_f90   (state_phinh_i,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_dp3d_f90      (state_dp3d,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS]                      [NP][NP]> ps_v_f90            (state_ps_v,m_num_elems);

  sync_to_device(state_v_f90,         m_v);
  sync_to_device(state_w_i_f90,       m_w_i);
  sync_to_device(state_vtheta_dp_f90, m_vtheta_dp);
  sync_to_device(state_phinh_i_f90,   m_phinh_i);
  sync_to_device(state_dp3d_f90,      m_dp3d);

  // F90 ptrs to arrays (np,np,num_time_levels,nelemd) can be stuffed directly in an unmanaged view
  // with scalar Real*[NUM_TIME_LEVELS][NP][NP] (with runtime dimension nelemd)

  auto ps_v_host = Kokkos::create_mirror_view(m_ps_v);
  Kokkos::deep_copy(ps_v_host,ps_v_f90);
  Kokkos::deep_copy(m_ps_v,ps_v_host);
}

void ElementsState::push_to_f90_pointers (F90Ptr& state_v, F90Ptr& state_w_i, F90Ptr& state_vtheta_dp,
                                          F90Ptr& state_phinh_i, F90Ptr& state_dp3d) const {
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ][2][NP][NP]> state_v_f90         (state_v,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_w_i_f90       (state_w_i,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_vtheta_dp_f90 (state_vtheta_dp,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_phinh_i_f90   (state_phinh_i,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_dp3d_f90      (state_dp3d,m_num_elems);

  sync_to_host(m_v,         state_v_f90);
  sync_to_host(m_w_i,       state_w_i_f90);
  sync_to_host(m_vtheta_dp, state_vtheta_dp_f90);
  sync_to_host(m_phinh_i,   state_phinh_i_f90);
  sync_to_host(m_dp3d,      state_dp3d_f90);
}

static bool all_good_elems (const ElementsState& s, const int tlvl) {
  using Kokkos::ALL;
  using Kokkos::parallel_for;

  const int nelem = s.num_elems();
  const int nplev = NUM_PHYSICAL_LEV;

  const auto vtheta_dp = Kokkos::subview(s.m_vtheta_dp, ALL, tlvl, ALL, ALL, ALL);
  const auto dp3d = Kokkos::subview(s.m_dp3d, ALL, tlvl, ALL, ALL, ALL);
  const auto phinh_i = Kokkos::subview(s.m_phinh_i, ALL, tlvl, ALL, ALL, ALL);

  // Write nerr = 1 if there is a problem; else do nothing.
  const auto check = KOKKOS_LAMBDA (const TeamMember& team, int& nerr) {
    const auto ie = team.league_rank();
    const auto g = [&] (const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      const Real* const v = reinterpret_cast<const Real*>(&vtheta_dp(ie,igp,jgp,0));
      const Real* const p = reinterpret_cast<const Real*>(&     dp3d(ie,igp,jgp,0));
      const Real* const h = reinterpret_cast<const Real*>(&  phinh_i(ie,igp,jgp,0));
      // Write races but doesn't matter since any single nerr = 1 write is
      // sufficient to conclude there is a problem.
      const auto f1 = [&] (const int k) { if (v[k] < 0 || isnan(v[k])) nerr = 1; };
      const auto f2 = [&] (const int k) { if (p[k] < 0 || isnan(p[k])) nerr = 1; };
      // The isnan checks on k and k+1 are redundant except for the last k, but
      // it's probably the fastest way to check all interface values.
      const auto f3 = [&] (const int k) { if (h[k] < h[k+1] || isnan(h[k]) || isnan(h[k+1])) nerr = 1; };
      const auto tvr = Kokkos::ThreadVectorRange(team, nplev);
      parallel_for(tvr, f1);
      parallel_for(tvr, f2);
      parallel_for(tvr, f3);
    };
    parallel_for(Kokkos::TeamThreadRange(team, NP*NP), g);
  };
  int nerr;
  parallel_reduce(get_default_team_policy<ExecSpace>(nelem), check, nerr);

  return nerr == 0;
}

void check_print_abort_on_bad_elems (const std::string& label, const int tlvl,
                                     const int error_code) {
  const auto& s = Context::singleton().get<ElementsState>();

  // On device and, thus, efficient.
  if (all_good_elems(s, tlvl)) return;

  // Now that we know there is an error, we can do the rest inefficiently.
  const auto& geometry = Context::singleton().get<ElementsGeometry>();
  const int nelem = s.num_elems();
  const int nplev = NUM_PHYSICAL_LEV, ntl = NUM_TIME_LEVELS, vecsz = VECTOR_SIZE;
  const auto& comm = Context::singleton().get<Connectivity>().get_comm();

  const auto vtheta_dp_h = Kokkos::create_mirror_view(s.m_vtheta_dp);
  Kokkos::deep_copy(vtheta_dp_h, s.m_vtheta_dp);
  const auto dp3d_h = Kokkos::create_mirror_view(s.m_dp3d);
  Kokkos::deep_copy(dp3d_h, s.m_dp3d);
  const auto phinh_i_h = Kokkos::create_mirror_view(s.m_phinh_i);
  Kokkos::deep_copy(phinh_i_h, s.m_phinh_i);
  const auto sphere_latlon = Kokkos::create_mirror_view(geometry.m_sphere_latlon);
  Kokkos::deep_copy(sphere_latlon, geometry.m_sphere_latlon);

  HostView<Real*****>
    vtheta_dp(reinterpret_cast<Real*>(&vtheta_dp_h(0,0,0,0,0)), nelem, ntl, NP, NP, NUM_LEV  *vecsz),
    dp3d     (reinterpret_cast<Real*>(&dp3d_h     (0,0,0,0,0)), nelem, ntl, NP, NP, NUM_LEV  *vecsz),
    phinh_i  (reinterpret_cast<Real*>(&phinh_i_h  (0,0,0,0,0)), nelem, ntl, NP, NP, NUM_LEV_P*vecsz);

  bool first = true;
  FILE* fid = nullptr;
  std::string filename;
  for (int ie = 0; ie < nelem; ++ie) {
    for (int gi = 0; gi < NP; ++gi)
      for (int gj = 0; gj < NP; ++gj) {
        int k_bad = -1;
        bool v = true, d = true, p = true;
        for (int k = 0; k < nplev; ++k) {
          v = std::isnan(vtheta_dp(ie,tlvl,gi,gj,k)) || vtheta_dp(ie,tlvl,gi,gj,k) < 0;
          d = std::isnan(dp3d(ie,tlvl,gi,gj,k)) || dp3d(ie,tlvl,gi,gj,k) < 0;
          p = (std::isnan(phinh_i(ie,tlvl,gi,gj,k)) || std::isnan(phinh_i(ie,tlvl,gi,gj,k+1)) ||
               phinh_i(ie,tlvl,gi,gj,k) < phinh_i(ie,tlvl,gi,gj,k+1));
          if (v || d || p) {
            k_bad = k;
            break;
          }
        }
        if (k_bad >= 0) {
          if (first) {
            filename = (std::string("hommexx.errlog.") +
                        std::to_string(comm.size()) + "." +
                        std::to_string(comm.rank()));
            fid = fopen(filename.c_str(), "w");
            fprintf(fid, "label: %s\ntime-level %d\n", label.c_str(), tlvl);
            first = false;
          }
          fprintf(fid, "lat %22.15e lon %22.15e\n",
                  sphere_latlon(ie,gi,gj,0), sphere_latlon(ie,gi,gj,1));
          fprintf(fid, "ie %d igll %d jgll %d lev %d:", ie, gi, gj, k_bad);
          if (v) fprintf(fid, " bad vtheta_dp");
          if (d) fprintf(fid, " bad dp3d");
          if (p) fprintf(fid, " bad dphi");
          fprintf(fid, "\n");
          fprintf(fid, "level                   dphi                   dp3d              vtheta_dp\n");
          for (int k = 0; k < nplev; ++k)
            fprintf(fid, "%5d %22.15e %22.15e %22.15e\n",
                    k, phinh_i(ie,tlvl,gi,gj,k+1) - phinh_i(ie,tlvl,gi,gj,k),
                    dp3d(ie,tlvl,gi,gj,k), vtheta_dp(ie,tlvl,gi,gj,k));
        }
      }
  }
  if (fid) fclose(fid);

  Errors::runtime_abort(std::string("Bad dphi, dp3d, or vtheta_dp; label: '") +
                        label + "'; see " + filename,
                        error_code < 0 ? Errors::err_bad_column_value : error_code);
}

} // namespace Homme
