/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

/* This is the second trajectory method for semi-Lagrangian transport; the first
   is in ComposeTransportImplTrajectory.cpp.

   Usage

   To use it, add the following setting to the Homme namelist:
       semi_lagrange_trajectory_nsubstep = N ! where N > 0
   A value of 0 (default) means Homme will use the original method.

   Another option makes the method use more than two velocity snapshots per time
   step:
       semi_lagrange_trajectory_nvelocity = N ! where N > 2
   A value <= 2 other than -1 means Homme will use the standard two velocity
   snapshots (time step end points). -1 (default) triggers an internal
   calculation based on nsubstep. This option has no effect if nsubstep=0
   (default).

   Summary

   This method provides multiple benefits over the original method, depending on
   configuration:
     * Supports (much) longer time steps than the original method.
     * Maximizes flexibility in specifying the various atmosphere time steps.
     * Greater accuracy than the original method for time steps the original
       method can handle.
     * Extreme accuracy in hypothetical niche applications.

   Method overview

   Recall that semi-Lagrangian tracer transport has six phases:
     1. At time step n+1, for each GLL point on the Eulerian grid, compute a
   trajectory backward in time to the departure point at time n. This step is
   independent of number of tracers.
     2. Simultaneously reconstruct vertically Lagrangian levels at time
   n+1. This step is independent of number of tracers.
     3. In each level, for each level-midpoint departure point, interpolate tracer
   mixing ratios at time n to the point. These mixing ratios are then the new
   ratios at time n+1 on the Eulerian grid.
     4. Optionally apply hyperviscosity.
     5. Apply the Communication-Efficient Density Reconstructor (CEDR).
     6. Vertically remap tracers at time n+1 from the reconstructed vertically
   Lagrangian levels to the vertically Eulerian grid.

   Trajectory methods implement phases 1 and 2.

   The key capability of this enhanced trajectory method (ETM) is to be able to
   take multiple substeps when computing the departure points. The original
   method cannot substep. Each substep has second-order accuracy, so the overall
   method is always second-order accurate. But as the number of substeps
   increases, so does accuracy.

   A second capability of the ETM is to use more velocity snapshots than just
   the tracer time step end-point snapshots. For example, if there are two
   substeps, the method can use three velocity snapshots: beginning, middle,
   end. The first substep uses (middle, end), and the second uses (beginning,
   middle). (This might be the opposite of what you expected, but recall that
   the trajectory is computed backward in time.)

   At a software level, a third capability is to use an arbitrarily large
   element halo when computing departure points. The original method is limited
   to two halo layers. Extra halo layers do not increase the cost of search for
   a fixed time step because of the ordering of elements in the layer.

   Speedup comes from the fact that taking a longer time step means phases 3-6,
   the most expensive phases, run less often. Phases 1 and 2 also run less often
   but take more time per run, summing to about the same cost over a fixed time
   duration T as the original method.

   Algorithm outline

     This method works principally in the eta coordinate. eta is constant in a
   level on the Eulerian grid.
     Let time t1 > t0 and consider a trajectory substep from t1 to t0.
     Terminology: An arrival point is the time-t1 point of a trajectory. A
   departure point is the time-t0 point.
     For each interface node at times t0 and t1, compute eta_dot dp/deta
   (calc_eta_dot_dpdn).
     For each midpoint node at times t0 and t1, compute
         eta_dot = eta_dot dp/deta/(A_eta p0 + B_eta ps)
   (calc_etadotmid_from_etadotdpdnint).
     Use eta_dot, the horizontal velocity data, and the update formula described
   in the comment for calc_nodal_velocities to compute the velocity term in the
   update formula at the Eulerian vertical-midpoint nodes. Call the result V_upd
   (calc_vel_horiz_formula_node_ref_mid, calc_eta_dot_formula_node_ref_mid).
     In general, the trajectory arrival point at t1 is not on the grid, but it
   is in the first substep. If it is not on the grid, interpolate V_upd to the
   arrival points to produce V_upd' (interp_v_update). A detail here is we
   should actually think of the original velocity data as being interpolated,
   and then V_upd' is computed from the interpolated data. But the formula to
   compute V_upd is linear in these data, so we can defer interpolation to the
   end and do it just once.
     Integrate the trajectory backward from t1 at the arrival point to t0 at the
   departure point using V_upd or V_upd' (update_dep_points).
     The algorithm up to this point can be substepped, running multiple times to
   go backward from t1 to t0 in multiple steps.
     After substepping is finished, there is one final part.
     So far we have computed departure points corresponding to Eulerian-grid
   arrival points. But now we need to account for the fact that the levels are
   vertically Lagrangian ("floating"). The arrival points we actually need are
   those on the floating levels at time t1 corresponding to Eulerian levels at
   time t0. This is implemented in
       describe interp_departure_points_to_floating_level_midpoints.
     We use the following notation: yi = I[y(x)](xi) is an interpolant
   constructed from y(x) and evaluated at xi.
     On input, we have departure-point level-midpoint eta, eta_dep_mid, and the
   corresponding horizontal position, p_dep_mid. We also know eta level
   interfaces and midpoints on the reference grid, eta_ref_int and eta_ref_mid,
   and the top and bottom boundary values of eta, eta(0) and eta(1).
     First, reconstruct Lagrangian levels at t1 on the arrival column:
       eta_arr_int = I[eta_ref_mid([eta(0),eta_dep_mid,eta(1)])](eta_ref_int).
   Second, compute the Lagrangian level midpoints at t1 on the arrival column:
       eta_arr_mid = I[eta_ref_mid([eta(0),eta_dep_mid,eta(1)])](eta_ref_mid).
   Third, compute the departure horizontal points corresponding to the arrival
   Lagrangian level midpoints:
       p_dep_mid(eta_arr_mid) = I[p_dep_mid(eta_ref_mid)](eta_arr_mid).
     These calcualtions provide us with the final results: p_dep_mid is used in
   phase 3, and eta_arr_int is used in phase 6.
 */

#include "ComposeTransportImplEnhancedTrajectoryImpl.hpp"

namespace Homme {

using CTI = ComposeTransportImpl;

Real& CTI::VelocityRecord::t_vel (const int i) {
  assert(_nvel > 2);
  assert(i >= 0 and i < _nvel);
  return _t_vel[i];
}

int& CTI::VelocityRecord::obs_slots (const int n, const int k) {
  assert(_nvel > 2);
  assert(n >= 0 and n <= _dtf);
  assert(k >= 0 and k <= 1);
  return _obs_slots[2*n+k];
}

Real& CTI::VelocityRecord::obs_wts (const int n, const int k) {
  assert(_nvel > 2);
  assert(n >= 0 and n <= _dtf);
  assert(k >= 0 and k <= 1);
  return _obs_wts[2*n+k];  
}

int& CTI::VelocityRecord::run_step (const int i) {
  assert(_nvel > 2);
  assert(i >= 0 and i <= _nsub);
  return _run_step[i];
}

Real CTI::VelocityRecord::t_vel (const int i) const {
  return const_cast<VelocityRecord*>(this)->t_vel(i);
}

int CTI::VelocityRecord::obs_slots (const int n, const int k) const {
  return const_cast<VelocityRecord*>(this)->obs_slots(n,k);
}

Real CTI::VelocityRecord::obs_wts (const int n, const int k) const {
  return const_cast<VelocityRecord*>(this)->obs_wts(n,k);
}

int CTI::VelocityRecord::run_step (const int i) const {
  return const_cast<VelocityRecord*>(this)->run_step(i);
}

void CTI::VelocityRecord
::init (const int dtf, const int drf_param, const int nsub, const int nvel_param) {
  const int
    drf = drf_param == 0 ? 1 : drf_param,
    navail = dtf/drf + 1,
    nvel = std::min(nvel_param == -1 ?
                    (2 + (nsub-1)/2) :  // default value
                    nvel_param,
                    std::min(nsub+1,    // can't use more than this
                             navail));  // this is the max available

  _dtf = dtf; _drf = drf; _nsub = nsub; _nvel = nvel;

  // nsub <= 1: No substepping.
  // nvel <= 2: Save velocity only at endpoints, as always occurs.
  if (nsub <= 1 or nvel <= 2) {
    _nvel = 2;
    return;
  }

  _t_vel.resize(nvel);
  _obs_slots.resize(2*dtf); _obs_wts.resize(2*dtf);
  _run_step.resize(nsub+1);

  // Times at which velocity data are available.
  std::vector<int> t_avail(navail); {
    t_avail[0] = 0;
    int i = 1;
    for (int n = 0; n < dtf; ++n) {
      if ((n+1) % drf != 0) continue;
      t_avail[i] = n+1;
      i = i + 1;
    }
    assert(i == navail);
    assert(t_avail[navail-1] == dtf);
  }

  // Times to which we associate velocity data.
  for (int n = 0; n < nvel; ++n) {
    t_vel(n) = ((n*dtf) % (nvel-1) == 0 ?
                /**/ (n*dtf) / (nvel-1) :
                Real (n*dtf) / (nvel-1));
    assert(t_vel(n) >= 0 and t_vel(n) <= dtf);
    assert(n == 0 or t_vel(n) > t_vel(n-1));
  }

  // Build the tables mapping n in 0:dtf-1 to velocity slots to accumulate into.
  for (int n = 0; n < dtf; ++n) {
    for (int k = 0; k < 2; ++k) {
      obs_slots(n,k) = -1;
      obs_wts(n,k) = 0;
    }
    if (n == dtf-1) continue;
    if ((n+1) % drf != 0) continue;
    const int time = n+1;
    int iav = -1;
    for (int i = 1; i <= navail-1; ++i)
      if (time == t_avail[i]) {
        iav = i;
        break;
      }
    assert(iav > 0 and iav < navail-1);
    for (int i = 1; i < nvel-1; ++i) {
      if (t_avail[iav-1] < t_vel(i) and time > t_vel(i)) {
        obs_slots(n,0) = i;
        obs_wts(n,0) = ((t_vel(i) - t_avail[iav-1]) /
                        (t_avail[iav] - t_avail[iav-1]));
      }
      if (time <= t_vel(i) and t_avail[iav+1] > t_vel(i)) {
        obs_slots(n,1) = i;
        obs_wts(n,1) = ((t_avail[iav+1] - t_vel(i)) /
                        (t_avail[iav+1] - t_avail[iav]));
      }
    }
  }

  // Build table mapping n to interval to use. The trajectories go backward in
  // time, and this table reflects that.
  run_step(0) = nvel-1;
  run_step(nsub) = 1;
  for (int n = 1; n < nsub; ++n) {
    const auto time = Real((nsub-n)*dtf)/nsub;
    int ifnd = -1;
    for (int i = 0; i < nvel-1; ++i)
      if (t_vel(i) <= time and time <= t_vel(i+1)) {
        ifnd = i;
        break;
      }
    assert(ifnd >= 0 and ifnd < nvel-1);
    run_step(n) = ifnd + 1;
  }
}

namespace {

// Set dep_points_all to level-midpoint arrival points.
void init_dep_points (const CTI& c, const cti::DeparturePoints& dep_pts) {
  const auto independent_time_steps = c.m_data.independent_time_steps;
  const auto& sphere_cart = c.m_geometry.m_sphere_cart;
  const auto& etai = c.m_hvcoord.etai;
  assert(not independent_time_steps or dep_pts.extent_int(4) == 4);
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < 3; ++d)
      dep_pts(ie,lev,i,j,d) = sphere_cart(ie,i,j,d);
    if (independent_time_steps)
      dep_pts(ie,lev,i,j,3) = etai(lev);
  };
  c.launch_ie_physlev_ij(f);
}

void update_dep_points (
  const CTI& c, const Real dtsub, const cti::DeparturePoints& vdep,
  const cti::DeparturePoints& dep_pts)
{
  const auto independent_time_steps = c.m_data.independent_time_steps;
  const auto is_sphere = c.m_data.geometry_type == 0;
  const auto scale_factor = c.m_geometry.m_scale_factor;
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    // Update horizontal position.
    Real p[3];
    for (int d = 0; d < 3; ++d)
      p[d] = dep_pts(ie,lev,i,j,d) - dtsub*vdep(ie,lev,i,j,d)/scale_factor;
    if (is_sphere) {
      const auto norm = std::sqrt(square(p[0]) + square(p[1]) + square(p[2]));
      for (int d = 0; d < 3; ++d)
        p[d] /= norm;
    }
    for (int d = 0; d < 3; ++d)
      dep_pts(ie,lev,i,j,d) = p[d];
    if (independent_time_steps) {
      // Update vertical position.
      dep_pts(ie,lev,i,j,3) -= dtsub*vdep(ie,lev,i,j,3);
    }
  };
  c.launch_ie_physlev_ij(f);
}

// Interpolate eta_dot at interfaces. The support data are not midpoint data,
// though; rather, they're interface data collected at different horizontal
// points. Thus, to be clear, this is not midpoint-to-interface interpolation of
// eta_dot.
void interp_etadot_at_interfaces (const CTI& c, const cti::DeparturePoints& vdep) {
  assert(vdep.extent_int(4) == 5);
  const auto& etai = c.m_hvcoord.etai;
  const CRNV<NUM_PHYSICAL_LEV> etam(cti::cpack2real(c.m_hvcoord.etam));
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    if (lev == 0) return;
    const auto a = (etai(lev) - etam(lev-1)) / (etam(lev) - etam(lev-1));
    // Safe to write to this slot in parallel b/c only this slot reads from it.
    vdep(ie,lev,i,j,3) = ((1-a)*vdep(ie,lev-1,i,j,4) +
                          (  a)*vdep(ie,lev  ,i,j,3));
  };
  c.launch_ie_physlev_ij(f);
}

/* Evaluate a formula to provide an estimate of nodal velocities that are use to
   create a 2nd-order update to the trajectory. The fundamental formula for the
   update in position p from arrival point p1 to departure point p0 is
       p0 = p1 - dt/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
   Here we compute the velocity estimate at the nodes:
       1/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
*/
template <typename Snapshots>
void calc_nodal_velocities (
  const CTI& c, const Real dtsub, const Snapshots& snaps,
  const cti::DeparturePoints& vnode)
{
  using Kokkos::ALL;
  const auto& d = c.m_data;
  const auto& h = c.m_hvcoord;
  const auto& sphere_ops = c.m_sphere_ops;
  const auto& vec_sph2cart = c.m_geometry.m_vec_sph2cart;
  const bool independent_time_steps = d.independent_time_steps;
  const auto ps0 = h.ps0;
  const auto hyai0 = h.hybrid_ai0;
  const auto& hybi = h.hybrid_bi_packed;
  const auto& hydai = h.hybrid_ai_delta;
  const auto& hydbi = h.hybrid_bi_delta;
  const auto& db_deta_i = d.db_deta_i;
  const auto& hyetai = h.etai;
  const auto& hyetam = h.etam;
  const auto& hydetai = d.hydetai;
  const auto& buf1a = d.buf1o[0]; const auto& buf1b = d.buf1o[1];
  const auto& buf1c = d.buf1o[2]; const auto& buf1d = d.buf1o[3];
  const auto& buf2a = d.buf2 [0]; const auto& buf2b = d.buf2 [1];
  const auto& buf2c = d.buf2 [2]; const auto& buf2d = d.buf2 [3];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto  wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto  wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto vwrk1 = Homme::subview(buf2a, kv.team_idx);
    const auto vwrk2 = Homme::subview(buf2b, kv.team_idx);
    CSelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                           Homme::subview(buf1d, kv.team_idx)};
    {
      SelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                            Homme::subview(buf1d, kv.team_idx)};
      if (independent_time_steps) {
        calc_eta_dot_ref(kv, sphere_ops, snaps,
                         ps0, hyai0, hybi,
                         hydai, hydbi, hydetai, db_deta_i,
                         wrk1, wrk2, vwrk1,
                         eta_dot);
      } else {
        for (int t = 0; t < 2; ++t) {
          const auto& ed = eta_dot[t];
          const auto f = [&] (const int i, const int j, const int k) {
            ed(i,j,k) = 0;
          };
          cti::loop_ijk<cti::num_lev_pack>(kv, f);
        }
      }
    }
    // Collect the horizontal nodal velocities.
    auto* vm1 = Homme::subview(buf2c, kv.team_idx).data();
    auto* vm2 = Homme::subview(buf2d, kv.team_idx).data();
    CS2elNlev vsph[] = {CS2elNlev(vm1), CS2elNlev(vm2)};
    {
      S2elNlev vsph[] = {S2elNlev(vm1), S2elNlev(vm2)};
      const auto e = snaps.get_element(ie);
      for (int t = 0; t < 2; ++t) {
        const auto& v = vsph[t];
        for (int d = 0; d < 2; ++d) {
          const auto f = [&] (const int i, const int j, const int k) {
            v(d,i,j,k) = e.combine_v(t, d, i, j, k);
          };
          cti::loop_ijk<cti::num_lev_pack>(kv, f);
        }
      }
    }
    kv.team_barrier();
    // Given the vertical and horizontal nodal velocities at time endpoints,
    // evaluate the velocity estimate formula, providing the final horizontal
    // and vertical velocity estimates at midpoint nodes.
    const auto vnode_ie = Kokkos::subview(vnode, ie, ALL,ALL,ALL,ALL);
    const auto vec_sph2cart_ie = Homme::subview(vec_sph2cart, ie);
    calc_vel_horiz_formula_node_ref_mid(kv, sphere_ops,
                                        hyetam, vec_sph2cart_ie,
                                        dtsub, vsph, eta_dot,
                                        wrk1, vwrk1, vwrk2,
                                        vnode_ie);
    if (independent_time_steps) {
      kv.team_barrier();
      calc_eta_dot_formula_node_ref_int(kv, sphere_ops,
                                        hyetai, hyetam,
                                        dtsub, vsph, eta_dot,
                                        wrk1, vwrk1,
                                        vnode_ie);
    }
  };
  Kokkos::parallel_for(c.m_tp_ne, f);
}

// Determine the departure points corresponding to the vertically Lagrangian
// grid's arrival midpoints, where the floating levels are those that evolve
// over the course of the full tracer time step. Also compute divdp, which holds
// the floating levels' dp values for later use in vertical remap.
int interp_departure_points_to_floating_level_midpoints (const CTI& c, const int np1) {
  using Kokkos::ALL;
  const int nlev = NUM_PHYSICAL_LEV, nlevp = nlev+1;
  const auto is_sphere = c.m_data.geometry_type == 0;
  const auto& d = c.m_data;
  const auto& h = c.m_hvcoord;
  const auto ps0 = h.ps0;
  const auto hyai0 = h.hybrid_ai0;
  const auto& hybi = h.hybrid_bi;
  const auto& hyetai = h.etai;
  const CRNV<NUM_PHYSICAL_LEV> hyetam(cti::cpack2real(h.etam));
  const CRNV<NUM_PHYSICAL_LEV> detai(cti::cpack2real(d.hydetai));
  const auto deta_tol = d.deta_tol;
  const auto& dep_pts = d.dep_pts;
  const auto& dp3d = c.m_state.m_dp3d;
  const auto& buf1a = d.buf1e[0]; const auto& buf1b = d.buf1e[1];
  const auto& buf1c = d.buf1e[2]; const auto& buf1d = d.buf1e[3];
  const auto& buf2a = d.buf2[0];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team, int& limcnt) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto wrk3 = Homme::subview(buf1c, kv.team_idx);
    const auto wrk4 = Homme::subview(buf1d, kv.team_idx);
    const auto vwrk = Homme::subview(buf2a, kv.team_idx);
    // Reconstruct Lagrangian levels at t1 on arrival column.
    const auto eta = p2rel(wrk3.data(), nlev);
    const auto f = [&] (const int i, const int j, const int k) {
      eta(i,j,k) = dep_pts(ie,k,i,j,3);
    };
    cti::loop_ijk<cti::num_phys_lev>(kv, f);
    kv.team_barrier();
    limcnt += limit_etai(kv, nlev,
                         hyetai, detai, deta_tol,
                         p2rel(wrk1.data(), nlev), p2rel(wrk2.data(), nlev),
                         eta);
    kv.team_barrier();
    {
      // Compute Lagrangian level interfaces at t1 on arrival column.
      const auto etai_arr = p2rel(wrk4.data(), nlevp);
      // eta_arr_int = I[eta_ref_int(eta_dep_int)](eta_ref_int)
      eta_interp_eta(kv, nlev, hyetai,
                     nlevp-2, 1, eta, hyetai,
                     p2rel(wrk1.data(), nlev+1), RnV(cti::pack2real(wrk2), nlev+1),
                     nlevp-2, 1, hyetai, etai_arr);
      const auto f = [&] (const int i, const int j) {
        etai_arr(i,j,0) = hyetai(0);
        etai_arr(i,j,nlev) = hyetai(nlev);
      };
      c.loop_ij(kv, f);
      // Compute divdp.
      const ExecViewUnmanaged<Real[NP][NP]> ps(cti::pack2real(vwrk));
      calc_ps(kv, nlev,
              ps0, hyai0,
              Homme::subview(dp3d, ie, np1),
              ps);
      kv.team_barrier();
      eta_to_dp(kv, nlev,
                ps0, hybi, hyetai,
                ps, etai_arr,
                p2rel(wrk2.data(), nlev+1),
                RelnV(cti::pack2real(Homme::subview(c.m_derived.m_divdp, ie)),
                      NP, NP, NUM_LEV*VECTOR_SIZE));
      kv.team_barrier();
    }
    // Compute Lagrangian level midpoints at t1 on arrival column.
    const auto etam_arr = p2rel(wrk4.data(), nlev);
    // eta_arr_mid = I[eta_ref_int(eta_dep_int)](eta_ref_mid)
    eta_interp_eta(kv, nlev, hyetai,
                   nlevp-2, 1, eta, hyetai,
                   p2rel(wrk1.data(), nlev+1), RnV(cti::pack2real(wrk2), nlev+1),
                   nlev, 0, hyetam, etam_arr);
    kv.team_barrier();
    // Compute departure horizontal points corresponding to arrival
    // Lagrangian level midpoints:
    //     p_dep_mid(eta_arr_mid) = I[p_dep_mid(eta_ref_mid)](eta_arr_mid)
    {
      const RelnV dpts_in(cti::pack2real(vwrk), NP, NP, nlev);
      const RelnV dpts_out(dpts_in.data() + NP*NP*nlev, NP, NP, nlev);
      for (int d = 0; d < 3; ++d) {
        const auto f = [&] (const int i, const int j, const int k) {
          dpts_in(i,j,k) = dep_pts(ie,k,i,j,d);
        };
        c.loop_ijk<cti::num_phys_lev>(kv, f);
        kv.team_barrier();
        eta_interp_horiz(kv, nlev,
                         hyetai,
                         hyetam, dpts_in,
                         RnV(cti::pack2real(wrk2), nlev+2), p2rel(wrk1.data(), nlev+2),
                         etam_arr, dpts_out);
        kv.team_barrier();
        const auto g = [&] (const int i, const int j, const int k) {
          dep_pts(ie,k,i,j,d) = dpts_out(i,j,k);
        };
        c.loop_ijk<cti::num_phys_lev>(kv, g);
        kv.team_barrier();
      }
      if (is_sphere) {
        // Normalize.
        const auto h = [&] (const int i, const int j, const int k) {
          Real norm = 0;
          for (int d = 0; d < 3; ++d) norm += square(dep_pts(ie,k,i,j,d));
          norm = std::sqrt(norm);
          for (int d = 0; d < 3; ++d) dep_pts(ie,k,i,j,d) /= norm;
        };
        c.loop_ijk<cti::num_phys_lev>(kv, h);
      }
    }
  };
  int limcnt = 0;
  Kokkos::parallel_reduce(c.m_tp_ne, f, limcnt);
  return limcnt;
}

void dss_vnode (const CTI& c, const cti::DeparturePoints& vnode) {
  const int ndim = c.m_data.independent_time_steps ? 4 : 3;
  const auto& spheremp = c.m_geometry.m_spheremp;
  const auto& rspheremp = c.m_geometry.m_rspheremp;
  const auto& vp = c.m_tracers.qtens_biharmonic;
  const ExecViewUnmanaged<Real**[NP][NP][NUM_LEV*VECTOR_SIZE]>
    v(cti::pack2real(vp), vp.extent_int(0), vp.extent_int(1));
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < ndim; ++d)
      v(ie,d,i,j,lev) = vnode(ie,lev,i,j,d)*spheremp(ie,i,j)*rspheremp(ie,i,j);
  };
  c.launch_ie_physlev_ij(f);
  Kokkos::fence();
  const auto be = c.m_v_dss_be[c.m_data.independent_time_steps ? 1 : 0];
  be->exchange();
  Kokkos::fence();
  const auto g = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < ndim; ++d)
      vnode(ie,lev,i,j,d) = v(ie,d,i,j,lev);
  };
  c.launch_ie_physlev_ij(g);
}

} // namespace anon

// For limit_etam.
void ComposeTransportImpl
::setup_enhanced_trajectory (const SimulationParams& params, const int num_elems) {
  assert(params.dt_tracer_factor > 0);

  const auto etai = cmvdc(m_hvcoord.etai);
  const Real deta_ave = (etai(num_phys_lev) - etai(0)) / num_phys_lev;
  m_data.deta_tol = 10*std::numeric_limits<Real>::epsilon()*deta_ave;

  // diff(etai)
  m_data.hydetai = decltype(m_data.hydetai)("hydetai");
  const auto detai_pack = Kokkos::create_mirror_view(m_data.hydetai);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> detai(pack2real(detai_pack));
  for (int k = 0; k < num_phys_lev; ++k)
    detai(k) = etai(k+1) - etai(k);
  Kokkos::deep_copy(m_data.hydetai, detai_pack);

  const auto etamp = cmvdc(m_hvcoord.etam);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> etam(pack2real(etamp));
  
  // hydetam_ref
  m_data.hydetam_ref = decltype(m_data.hydetam_ref)("hydetam_ref");
  {
    const auto m = Kokkos::create_mirror_view(m_data.hydetam_ref);
    const int nlev = num_phys_lev;
    m(0) = etam(0) - etai(0);
    for (int k = 1; k < nlev; ++k) m(k) = etam(k) - etam(k-1);
    m(nlev) = etai(nlev) - etam(nlev-1);
    Kokkos::deep_copy(m_data.hydetam_ref, m);
  }

  // etam
  homme::compose::set_hvcoord(etai.data(), etam.data());

  // Initialization for semi_lagrange_trajectory_nvelocity > 2.
  m_data.vrec = std::make_shared<VelocityRecord>(
    params.dt_tracer_factor, params.dt_remap_factor, m_data.trajectory_nsubstep,
    m_data.trajectory_nvelocity);
  const auto nv = m_data.vrec->nvel();
  if (nv > 2) {
    m_data.dp_extra_snapshots = DpSnaps( "dp_extra_snapshots", num_elems, nv-2);
    m_data.vel_extra_snapshots = VSnaps("vel_extra_snapshots", num_elems, nv-2);
  }

  // B_eta at interfaces
  m_data.db_deta_i = decltype(m_data.db_deta_i)("db_deta_i");
  {
    const auto m_p = Kokkos::create_mirror_view(m_data.db_deta_i);
    HostViewUnmanaged<Real[NUM_INTERFACE_LEV]> m(pack2real(m_p));
    m(0) = m(NUM_INTERFACE_LEV-1) = 0; // unused
    const auto hybi = cmvdc(m_hvcoord.hybrid_bi);
    for (int k = 2; k < NUM_PHYSICAL_LEV-1; ++k)
      m(k) = approx_derivative(etai(k-1), etai(k), etai(k+1),
                               hybi(k-1), hybi(k), hybi(k+1));
    Kokkos::deep_copy(m_data.db_deta_i, m_p);
  }
}

void ComposeTransportImpl::observe_velocity (const TimeLevel& tl, const int step) {
  // Optionally observe the dynamics velocity snapshot available at this step.

  if (not m_data.vrec or m_data.vrec->nvel() == 2) return;

  const auto& v = *m_data.vrec;
  assert((step + 1) % v.drf() == 0);

  if (step + 1 == v.drf()) {
    // This is either (1) the first vertical remap step in the tracer step or
    // (2) the first dynamics step and we're running vertically Eulerian. In
    // either case, zero the quantities accumulated over the step.
    Kokkos::deep_copy(m_data.dp_extra_snapshots, 0);
    Kokkos::deep_copy(m_data.vel_extra_snapshots, 0);
  }

  const auto& state = Context::singleton().get<ElementsState>();
  const auto np1 = tl.np1;
  for (int t = 0; t < 2; ++t) {
    const auto slot = v.obs_slots(step, t);
    if (slot == -1) continue;
    assert(slot > 0 and slot < v.nvel()-1);
    const auto wt = v.obs_wts(step, t);
    const auto& dp_snap = m_data.dp_extra_snapshots;
    const auto& v_snap = m_data.vel_extra_snapshots;
    const auto& dp3d = state.m_dp3d;
    const auto& v = state.m_v;
    const auto f = KOKKOS_LAMBDA (const int idx) {
      int ie, igp, jgp, ilev;
      idx_ie_ij_nlev<num_lev_pack>(idx, ie, igp, jgp, ilev);
      dp_snap(ie,slot-1,igp,jgp,ilev) += wt * dp3d(ie,np1,igp,jgp,ilev);
      for (int d = 0; d < 2; ++d)
        v_snap(ie,slot-1,d,igp,jgp,ilev) += wt * v(ie,np1,d,igp,jgp,ilev);
    };
    launch_ie_ij_nlev<num_lev_pack>(f);
  }
}

void ComposeTransportImpl
::calc_enhanced_trajectory (const int nstep, const int np1, const Real dt) {
  GPTLstart("compose_calc_enhanced_trajectory");

  const auto& dep_pts = m_data.dep_pts;
  const auto& vnode = m_data.vnode;
  const auto& vdep = m_data.vdep;

  init_dep_points(*this, dep_pts);

  const int nelemd = m_data.nelemd;
  const Real dtsub = dt / m_data.trajectory_nsubstep;
  const int nsubstep = m_data.trajectory_nsubstep;
  for (int step = 0; step < nsubstep; ++step) {
    {
      Kokkos::fence();
      GPTLstart("compose_vnode");
      const CVSnap v1(m_derived.m_vstar.data(), nelemd, 1);
      const CDpSnap dp1(m_derived.m_dp.data(), nelemd, 1);
      const auto& v2 = m_state.m_v;
      const auto& dp2 = m_state.m_dp3d;
      if (m_data.vrec->nvel() == 2) {
        const Real alpha[] = {Real(nsubstep-step-1)/nsubstep,
                              Real(nsubstep-step  )/nsubstep};
        EndpointSnapshots snaps(alpha, dp1, v1, 0, dp2, v2, np1);
        calc_nodal_velocities(*this, dtsub, snaps, vnode);
      } else {
        ManySnapshots snaps(*m_data.vrec,
                            dp1, v1, 0, dp2, v2, np1,
                            m_data.dp_extra_snapshots, m_data.vel_extra_snapshots,
                            nsubstep, step);
        calc_nodal_velocities(*this, dtsub, snaps, vnode);
      }
        
      Kokkos::fence();
      GPTLstop("compose_vnode");
    }

    GPTLstart("compose_v_bexchv");
    dss_vnode(*this, vnode);
    Kokkos::fence();
    GPTLstop("compose_v_bexchv");

    if (step == 0) {
      update_dep_points(*this, dtsub, vnode, dep_pts);
    } else {
      GPTLstart("compose_vdep");
      homme::compose::interp_v_update(step, dtsub);
      Kokkos::fence();
      GPTLstop("compose_vdep");

      interp_etadot_at_interfaces(*this, vdep);

      update_dep_points(*this, dtsub, vdep, dep_pts);
    }
  }
  Kokkos::fence();

  if (m_data.independent_time_steps) {
    GPTLstart("compose_floating_dep_pts");
    const int limcnt =
      interp_departure_points_to_floating_level_midpoints(*this, np1);
    Kokkos::fence();
    GPTLstop("compose_floating_dep_pts");
    if (m_data.diagnostics & 1) {
      const auto& c = Context::singleton();
      const auto& comm = c.get<Comm>();
      int glimcnt;
      MPI_Allreduce(&limcnt, &glimcnt, 1, MPI_INT, MPI_SUM, comm.mpi_comm());
      if (glimcnt > 0 and comm.root())
        printf("COMPOSE> nstep %10d limiter_active_count %10d\n",
               nstep, glimcnt);
    }
  }

  GPTLstop("compose_calc_enhanced_trajectory");
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
