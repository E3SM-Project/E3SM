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
   update formula at the Eulerian vertical-midpoint nodes. Call the result V
   (calc_vel_horiz_formula_node_ref_mid, calc_eta_dot_formula_node_ref_mid).
     In general, the trajectory arrival point at t1 is not on the grid, but it
   is in the first substep. If it is not on the grid, interpolate V to the
   arrival points to produce V_dep (calc_v_departure). A detail here is we
   should actually think of the original velocity data as being interpolated,
   and then V_dep is computed from the interpolated data. But the formula to
   compute V is linear in these data, so we can defer interpolation to the end
   and do it just once.
     Integrate the trajectory backward from t1 at the arrival point to t0 at the
   departure point using V_dep (update_dep_points).
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
namespace {

// Set dep_points_all to level-midpoint arrival points.
void init_dep_points (const CTI& c, const cti::DeparturePoints& dep_pts) {
  const auto independent_time_steps = c.m_data.independent_time_steps;
  const auto& sphere_cart = c.m_geometry.m_sphere_cart;
  const CRNV<NUM_PHYSICAL_LEV> hyetam(cti::cpack2real(c.m_hvcoord.etam));
  assert(not independent_time_steps or dep_pts.extent_int(4) == 4);
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < 3; ++d)
      dep_pts(ie,lev,i,j,d) = sphere_cart(ie,i,j,d);
    if (independent_time_steps)
      dep_pts(ie,lev,i,j,3) = hyetam(lev);
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

/* Evaluate a formula to provide an estimate of nodal velocities that are use to
   create a 2nd-order update to the trajectory. The fundamental formula for the
   update in position p from arrival point p1 to departure point p0 is
       p0 = p1 - dt/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
   Here we compute the velocity estimate at the nodes:
       1/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
*/
void calc_nodal_velocities (
  const CTI& c, const Real dtsub, const Real halpha[2],
  const cti::CVSlot& v1, const cti::CDpSlot& dp1, const int idx1,
  const cti::CVSlot& v2, const cti::CDpSlot& dp2, const int idx2,
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
  const auto& hyetam = h.etam;
  const auto& hyetai = h.etai;
  const auto& hydetai = d.hydetai;
  const auto& buf1a = d.buf1o[0]; const auto& buf1b = d.buf1o[1];
  const auto& buf1c = d.buf1o[2]; const auto& buf1d = d.buf1o[3];
  const auto& buf2a = d.buf2 [0]; const auto& buf2b = d.buf2 [1];
  const auto& buf2c = d.buf2 [2]; const auto& buf2d = d.buf2 [3];
  const auto alpha0 = halpha[0], alpha1 = halpha[1];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto  wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto  wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto vwrk1 = Homme::subview(buf2a, kv.team_idx);
    const auto vwrk2 = Homme::subview(buf2b, kv.team_idx);
    const auto v1_ie = Homme::subview(v1, ie, idx1);
    const auto v2_ie = Homme::subview(v2, ie, idx2);
    const Real alpha[] = {alpha0, alpha1};
    CSelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                           Homme::subview(buf1d, kv.team_idx)};
    {
      SelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                            Homme::subview(buf1d, kv.team_idx)};
      if (independent_time_steps) {
        const auto dp1_ie = Homme::subview(dp1, ie, idx1);
        const auto dp2_ie = Homme::subview(dp2, ie, idx2);
        calc_eta_dot_ref_mid(kv, sphere_ops,
                             ps0, hyai0, hybi, hydai, hydbi, hydetai,
                             alpha, v1_ie, dp1_ie, v2_ie, dp2_ie,
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
    // Collect the horizontal nodal velocities. v1,2 are on Eulerian levels. v1
    // is from time t1 < t2.
    auto* vm1 = Homme::subview(buf2c, kv.team_idx).data();
    auto* vm2 = Homme::subview(buf2d, kv.team_idx).data();
    CS2elNlev vsph[] = {CS2elNlev(vm1), CS2elNlev(vm2)};
    {
      S2elNlev vsph[] = {S2elNlev(vm1), S2elNlev(vm2)};
      for (int t = 0; t < 2; ++t) {
        const auto& v = vsph[t];
        for (int d = 0; d < 2; ++d) {
          const auto f = [&] (const int i, const int j, const int k) {
            v(d,i,j,k) = ((1 - alpha[t])*v1_ie(d,i,j,k) +
                          /**/ alpha[t] *v2_ie(d,i,j,k));
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
      calc_eta_dot_formula_node_ref_mid(kv, sphere_ops,
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
void interp_departure_points_to_floating_level_midpoints (const CTI& c, const int np1) {
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
  const auto& detam_ref = d.hydetam_ref;
  const auto deta_tol = d.deta_tol;
  const auto& dep_pts = d.dep_pts;
  const auto& dp3d = c.m_state.m_dp3d;
  const auto& buf1a = d.buf1e[0]; const auto& buf1b = d.buf1e[1];
  const auto& buf1c = d.buf1e[2]; const auto& buf1d = d.buf1e[3];
  const auto& buf2a = d.buf2[0];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto wrk3 = Homme::subview(buf1c, kv.team_idx);
    const auto wrk4 = Homme::subview(buf1d, kv.team_idx);
    const auto vwrk = Homme::subview(buf2a, kv.team_idx);
    // Reconstruct Lagrangian levels at t1 on arrival column:
    //     eta_arr_int = I[eta_ref_mid([eta(0),eta_dep_mid,eta(1)])](eta_ref_int)
    const auto etam = p2rel(wrk3.data(), nlev);
    const auto f = [&] (const int i, const int j, const int k) {
      etam(i,j,k) = dep_pts(ie,k,i,j,3);
    };
    cti::loop_ijk<cti::num_phys_lev>(kv, f);
    kv.team_barrier();
    limit_etam(kv, nlev,
               hyetai, detam_ref, deta_tol,
               p2rel(wrk1.data(), nlevp), p2rel(wrk2.data(), nlevp),
               etam);
    kv.team_barrier();
    {
      // Compute eta_arr_int.
      const auto etai_arr = p2rel(wrk4.data(), nlevp);
      eta_interp_eta(kv, nlev,
                     hyetai,
                     etam, hyetam,
                     p2rel(wrk1.data(), nlev+2), RnV(cti::pack2real(wrk2), nlev+2),
                     nlevp-2, hyetai, etai_arr, 1);
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
    // Compute Lagrangian level midpoints at t1 on arrival column:
    //     eta_arr_mid = I[eta_ref_mid([eta(0),eta_dep_mid,eta(1)])](eta_ref_mid)
    const auto etam_arr = p2rel(wrk4.data(), nlev);
    eta_interp_eta(kv, nlev,
                   hyetai,
                   etam, hyetam,
                   p2rel(wrk1.data(), nlev+2), RnV(cti::pack2real(wrk2), nlev+2),
                   nlev, hyetam, etam_arr);
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
  Kokkos::parallel_for(c.m_tp_ne, f);
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
void ComposeTransportImpl::setup_enhanced_trajectory () {
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
  
  // hydetam_ref.
  m_data.hydetam_ref = decltype(m_data.hydetam_ref)("hydetam_ref");
  const auto m = Kokkos::create_mirror_view(m_data.hydetam_ref);
  const int nlev = num_phys_lev;
  m(0) = etam(0) - etai(0);
  for (int k = 1; k < nlev; ++k) m(k) = etam(k) - etam(k-1);
  m(nlev) = etai(nlev) - etam(nlev-1);
  Kokkos::deep_copy(m_data.hydetam_ref, m);

  // etam
  homme::compose::set_hvcoord(etai(0), etai(num_phys_lev), etam.data());
}

void ComposeTransportImpl::calc_enhanced_trajectory (const int np1, const Real dt) {
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
      const Real alpha[] = {Real(nsubstep-step-1)/nsubstep,
                            Real(nsubstep-step  )/nsubstep};
      const CVSlot v1(m_derived.m_vstar.data(), nelemd, 1);
      const CDpSlot dp1(m_derived.m_dp.data(), nelemd, 1);
      const auto& v2 = m_state.m_v;
      const auto& dp2 = m_state.m_dp3d;
      calc_nodal_velocities(*this, dtsub, alpha,
                            v1, dp1, 0, v2, dp2, np1,
                            vnode);
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
      homme::compose::calc_v_departure(step, dtsub);
      Kokkos::fence();
      GPTLstop("compose_vdep");

      update_dep_points(*this, dtsub, vdep, dep_pts);
    }
  }
  Kokkos::fence();

  if (m_data.independent_time_steps) {
    GPTLstart("compose_floating_dep_pts");
    interp_departure_points_to_floating_level_midpoints(*this, np1);
    Kokkos::fence();
    GPTLstop("compose_floating_dep_pts");
  }

  GPTLstop("compose_calc_enhanced_trajectory");
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
