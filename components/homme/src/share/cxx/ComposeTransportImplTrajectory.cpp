/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"
#include "PhysicalConstants.hpp"

#include "compose_test.hpp"

namespace Homme {
using cti = ComposeTransportImpl;

KOKKOS_FUNCTION
static void ugradv_sphere (
  const SphereOperators& sphere_ops, const KernelVariables& kv,
  const typename ViewConst<ExecViewUnmanaged<Real[2][3][NP][NP]> >::type& vec_sphere2cart,
  // velocity, latlon
  const typename ViewConst<ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> >::type& u,
  const typename ViewConst<ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> >::type& v,
  const ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>& v_cart,
  const ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]>& ugradv_cart,
  // [u dot grad] v, latlon
  const ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]>& ugradv)
{
  for (int d_cart = 0; d_cart < 3; ++d_cart) {
    const auto f1 = [&] (const int i, const int j, const int k) {
      v_cart(i,j,k) = (vec_sphere2cart(0,d_cart,i,j) * v(0,i,j,k) +
                       vec_sphere2cart(1,d_cart,i,j) * v(1,i,j,k));      
    };
    cti::loop_ijk<NUM_LEV>(kv, f1);
    kv.team_barrier();

    sphere_ops.gradient_sphere<NUM_LEV>(kv, v_cart, ugradv_cart);

    const auto f2 = [&] (const int i, const int j, const int k) {
      if (d_cart == 0) ugradv(0,i,j,k) = ugradv(1,i,j,k) = 0;
      for (int d_latlon = 0; d_latlon < 2; ++d_latlon)
        ugradv(d_latlon,i,j,k) +=
          vec_sphere2cart(d_latlon,d_cart,i,j)*
          (u(0,i,j,k) * ugradv_cart(0,i,j,k) + u(1,i,j,k) * ugradv_cart(1,i,j,k));
    };
    cti::loop_ijk<NUM_LEV>(kv, f2);
  }
}

typedef typename ViewConst<ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]> >::type CSNlev;
typedef typename ViewConst<ExecViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]> >::type CRNlev;
typedef typename ViewConst<ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV_P]> >::type CSNlevp;
typedef typename ViewConst<ExecViewUnmanaged<Real[NP][NP][NUM_LEV_P*VECTOR_SIZE]> >::type CRNlevp;
typedef typename ViewConst<ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> >::type CS2Nlev;
typedef ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]> SNlev;
typedef ExecViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]> RNlev;
typedef ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV_P]> SNlevp;
typedef ExecViewUnmanaged<Real[NP][NP][NUM_LEV_P*VECTOR_SIZE]> RNlevp;
typedef ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> S2Nlev;
typedef ExecViewUnmanaged<Real[2][NP][NP][NUM_LEV*VECTOR_SIZE]> R2Nlev;
typedef ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV_P]> S2Nlevp;

/* Form a 3rd-degree Lagrange polynomial over (x(k-1:k+1), y(k-1:k+1)) and set
   yi(k) to its derivative at x(k). yps(:,:,0) is not written.
 */
KOKKOS_FUNCTION static void approx_derivative (
  const KernelVariables& kv, const CSNlevp& xs, const CSNlevp& ys,
  const SNlev& yps) // yps(:,:,0) is undefined
{
  CRNlevp x(cti::cpack2real(xs));
  CRNlevp y(cti::cpack2real(ys));
  RNlev yp(cti::pack2real(yps));
  const auto f = [&] (const int i, const int j, const int k) {
    if (k == 0) return;
    const auto& xkm1 = x(i,j,k-1);
    const auto& xk   = x(i,j,k  ); // also the interpolation point
    const auto& xkp1 = x(i,j,k+1);
    yp(i,j,k) = (y(i,j,k-1)*((         1 /(xkm1 - xk  ))*((xk - xkp1)/(xkm1 - xkp1))) +
                 y(i,j,k  )*((         1 /(xk   - xkm1))*((xk - xkp1)/(xk   - xkp1)) +
                             ((xk - xkm1)/(xk   - xkm1))*(         1 /(xk   - xkp1))) +
                 y(i,j,k+1)*(((xk - xkm1)/(xkp1 - xkm1))*(         1 /(xkp1 - xk  ))));
  };
  cti::loop_ijk<cti::num_phys_lev>(kv, f);
}

// Pad by an amount ~ smallest level to keep the computed dp > 0.
void ComposeTransportImpl::set_dp_tol () {
  const auto dp0h = cmvdc(m_hvcoord.dp0);
  const Real* dp0 = pack2real(dp0h);
  Real min_dp0 = dp0[0];
  for (int i = 1; i < num_phys_lev; ++i) min_dp0 = std::min(min_dp0, dp0[i]);
  m_data.dp_tol = 10*std::numeric_limits<Real>::epsilon()*min_dp0;
}

KOKKOS_FUNCTION static void
calc_p (const KernelVariables& kv, const Real& ps0, const Real& hybrid_ai0,
        const CSNlev& dpp, const SNlevp& pp) {
  CRNlev dp(cti::cpack2real(dpp));
  RNlevp p(cti::pack2real(pp));
  const auto f = [&] (const int i, const int j) {
    p(i,j,0) = hybrid_ai0*ps0;
    for (int k = 0; k < cti::num_phys_lev; ++k)
      p(i,j,k+1) = p(i,j,k) + dp(i,j,k);
  };
  cti::loop_ij(kv, f);
}

// Move mass around in a column as needed to make dp nonnegative.
KOKKOS_FUNCTION static void
reconstruct_and_limit_dp (const KernelVariables& kv, const CSNlev& dprefp, const Real& dt,
                          const Real& dp_tol, const CSNlevp& eta_dot_dpdn_p,
                          const SNlev& dpreconp) {
  const CRNlev dpref(cti::cpack2real(dprefp));
  const CRNlevp eta_dot_dpdn(cti::cpack2real(eta_dot_dpdn_p));
  const RNlev dprecon(cti::pack2real(dpreconp));
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV);
  const auto f = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g1 = [&] (const int k, Kokkos::Real2& sums) {
      // Store the full update rather than reconstructing eta_dot_dpdn. See
      // comment in sl_vertically_remap_tracers for more.
      dprecon(i,j,k) = dpref(i,j,k) + dt*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k));
      if (dprecon(i,j,k) < dp_tol) {
        sums.v[0] += dprecon(i,j,k) - dp_tol;
        dprecon(i,j,k) = dp_tol;
      } else {
        sums.v[1] += dprecon(i,j,k) - dp_tol;
      }
    };
    Kokkos::Real2 sums;
    Dispatch<>::parallel_reduce(kv.team, tvr, g1, sums);
    kv.team_barrier();
    const Real nmass = sums.v[0];
    if (nmass == 0) return;
    // Compensate for clipping.
    const Real wsum = sums.v[1];
    const auto g2 = [&] (const int k) {
      const Real w = dprecon(i,j,k) - dp_tol;
      dprecon(i,j,k) += nmass*(w/wsum);
    };
    Kokkos::parallel_for(tvr, g2);
  };
  Kokkos::parallel_for(ttr, f);
}

/*
  Reconstruct the vertically Lagrangian levels, thus permitting the dynamics
  vertical remap time step to be shorter than the tracer time step.
    Recall
      p(eta,ps) = A(eta) p0 + B(eta) ps
      => dp/dt = p_eta deta/dt + p_ps dps/dt
               = (A_eta p0 + B_eta ps) deta/dt + B(eta) dps/dt
  dp3d already accounts for B(eta) dps/dt, so it does not appear in what
  follows. We use eta as the vertical coordinate in these calculations. But
  here, to focus on the key ideas, consider a p = (x,z) system with
  velocities v = (u,w). Let 0, h, 1, suffixes denote start, middle, and end
  of the time step. The algorithm is as follows.
    The trajectory is computed forward in time, so the departure 0 points
  are on the grid. We need to compute z(x0,t1), the floating height at
  horizontal grid point x0 at time t1. We start by computing z1 = z(x1,t1),
  the floating height at the non-grid point x1 at time t1 or, in other
  words, the non-grid arrival point:
      z1 = z0 + dt/2 (w(p0,t0) + w(p1,t1)) + O(dt^3)
      w(p1,t1) = w(p0,t1) + grad w(p0,t1) (p1 - p0) + O(dt^2)
      x1 - x0 = dt u(p0,t0) + O(dt^2)
      z1 - z0 = dt w(p0,t0) + O(dt^2)
      z1 = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
                      dt (w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0))) + O(dt^3)  (*)
  Now we compute z(x0,t1). First, we need
      x0 - x1 = -dt u(p0,t0) + O(dt^2)
  and
      z_x(x1,t1) = z0_x + dt w_x(p0,t1) + O(dt^2)
                 = dt w_x(p0,t1) + O(dt^2).
      z_xx(x1,t1) = z0_xx + dt w_xx(p0,t1) + O(dt^2)
                  = dt w_xx(p0,t1) + O(dt^2).
  Then we expand z(x0,t1) in Taylor series and substitute:
      z(x0,t1) = z(x1,t1) + z_x(x1,t1) (x0 - x1) + z_xx(x1,t1) (x0 - x1)^2
                 + O(|x0 - x1|^3)
               = z(x1,t1) - dt z_x(x1,t1) u(p0,t0)
                 + dt^2 z_xx(x1,t1) u(p0,t0)^2
                 + O(dt^3) + O(|x0 - x1|^3)
               = z(x1,t1) - dt^2 w_x(p0,t1) u(p0,t0)
                 + dt^3 w_xx(x0,t1) u(p0,t0)^2
                 + O(dt^3) + O(|x0 - x1|^3).
  Gather all O(dt^3) terms, with O(|x0 - x1|^p) = O(dt^p):
      z(x0,t1) = z(x1,t1) - dt^2 w_x(p0,t1) u(p0,t0) + O(dt^3).
  Now substitute (*) for z(x1,t1):
      z(x0,t1) = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
                            dt (w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0)))
                 - dt^2 w_x(p0,t1) u(p0,t0) + O(dt^3)
               = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
                            dt (-w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0)))
                 + O(dt^3)
  This is locally accurate to O(dt^3) and so globally 2nd-order
  accurate. Notably, compared with (*), this formula differs only in a
  sign. Note also that a straightforward first-order accurate formula is
      z(x0,t1) = z0 + dt w(p0,th) + O(dt^2).
 */
KOKKOS_FUNCTION static void calc_vertically_lagrangian_levels (
  const SphereOperators& sphere_ops, const KernelVariables& kv,
  const Real& ps0, const Real& hybrid_ai0, const ExecViewUnmanaged<Scalar[NUM_LEV_P]>& hybrid_bi,
  const CSNlev& dp3d, const CSNlev& dp, const CS2Nlev& vn0, const CS2Nlev& vstar,
  const Real& dt, const Real& dp_tol,
  const SNlevp& wrk1a, const SNlevp& wrk1b, const SNlevp& wrk1c, const S2Nlevp& wrk2,
  const SNlev& dprecon)
{
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;
  using Kokkos::ThreadVectorRange;

  assert(hybrid_bi(0)[0] == 0);

  const auto ttr = TeamThreadRange(kv.team, NP*NP);
  const auto tvr = ThreadVectorRange(kv.team, NUM_LEV);
  
  // Reconstruct an approximation to endpoint eta_dot_dpdn on Eulerian levels.
  const auto& divdp = wrk1a;
  const auto& vdp = wrk2;
  const SNlevp* eta_dot_dpdn[] = {&wrk1b, &wrk1c};
  for (int t = 0; t < 2; ++t) {
    const auto& edd = *eta_dot_dpdn[t];
    if (t == 0) {
      const auto f = [&] (const int i, const int j, const int kp) {
        for (int d = 0; d < 2; ++d)
          vdp(d,i,j,kp) = vstar(d,i,j,kp)*dp(i,j,kp);
      };
      cti::loop_ijk<cti::num_lev_pack>(kv, f);
    } else {
      const auto f = [&] (const int i, const int j, const int kp) {
        for (int d = 0; d < 2; ++d)
          vdp(d,i,j,kp) = vn0(d,i,j,kp)*dp3d(i,j,kp);
      };
      cti::loop_ijk<cti::num_lev_pack>(kv, f);
    }

    sphere_ops.divergence_sphere(kv, vdp, divdp);

    RNlevp edds(cti::pack2real(edd)), divdps(cti::pack2real(divdp));
    const auto f = [&] (const int idx) {
      const int i = idx / NP, j = idx % NP;
      const auto r = [&] (const int k, Real& dps, const bool final) {
        assert(k != 0 || dps == 0);
        if (final) edds(i,j,k) = dps;
        dps += divdps(i,j,k);
      };
      Dispatch<>::parallel_scan(kv.team, cti::num_phys_lev, r);
      const int kend = cti::num_phys_lev - 1;
      const Real dps = edds(i,j,kend) + divdps(i,j,kend);
      assert(hybrid_bi(0)[0] == 0);
      const auto s = [&] (const int kp) {
        edd(i,j,kp) = hybrid_bi(kp)*dps - edd(i,j,kp);
        if (kp == 0) edd(i,j,kp)[0] = 0;
      };
      parallel_for(tvr, s);
      assert(edds(i,j,0) == 0);
      const int bottom = cti::num_phys_lev;
      edds(i,j,bottom) = 0;
    };
    parallel_for(ttr, f);
  }
  
  // Use p0 as the reference coordinate system. p0 differs from p1 by B(eta)
  // (ps1 - ps0); dp3d already accounts for this term w.r.t. derived%dp. Recall
  //     eta_dot_dpdn = p_eta eta_dot = (A_eta p0 + B_eta ps) deta/dt,
  // except that in the code eta_dot_dpdn is actually dp deta/dt rather than
  // dp/deta deta/dt. eta_dot_dpdn is the motion of a pressure level excluding
  // its motion due to dps/dt.
  const auto& pref = wrk1a;
  calc_p(kv, ps0, hybrid_ai0, dp, pref);
    
  // Gradient of eta_dot_dpdn = p_eta deta/dt at final time
  // w.r.t. horizontal sphere coords.
  const auto& grad = S2Nlev(wrk2.data());
  sphere_ops.gradient_sphere(kv, *eta_dot_dpdn[1], grad);
  
  // Gradient of eta_dot_dpdn = p_eta deta/dt at final time w.r.t. p at initial
  // time.
  const auto& ptp0 = dprecon;
  approx_derivative(kv, pref, *eta_dot_dpdn[1], ptp0);

  {
    const auto& edd = *eta_dot_dpdn[0];
    const R2Nlev vstars(cti::pack2real(vstar));
    const auto f_v = [&] (const int i, const int j, const int kp) {
      // Horizontal velocity at initial time.
      Scalar v[2];
      if (kp == 0) {
        for (int d = 0; d < 2; ++d)
          for (int s = 1; s < cti::packn; ++s)
            v[d][s] = 0.5*(vstars(d,i,j,s-1) + vstars(d,i,j,s));
      } else {
        const int os = kp*cti::packn;
        for (int d = 0; d < 2; ++d)
          for (int s = 0; s < cti::packn; ++s)
            v[d][s] = 0.5*(vstars(d,i,j,os+s-1) + vstars(d,i,j,os+s));
      }
      // Reconstruct eta_dot_dpdn over the time interval. Boundary points are
      // always 0.
#define SL_ADVECTION_TRAJ_OLD
#ifdef SL_ADVECTION_TRAJ_OLD
      // Reconstruct departure level coordinate at final time.
      edd(i,j,kp) = (pref(i,j,kp) +
                     0.5*dt*(edd(i,j,kp) + (*eta_dot_dpdn[1])(i,j,kp) +
                             dt*(ptp0(i,j,kp)*edd(i,j,kp)
                                 - grad(0,i,j,kp)*v[0] - grad(1,i,j,kp)*v[1])));
      edd(i,j,kp) = (edd(i,j,kp) - pref(i,j,kp))/dt;
#else
      edd(i,j,kp) = (0.5*(edd(i,j,kp) + (*eta_dot_dpdn[1])(i,j,kp) +
                          dt*(ptp0(i,j,kp)*edd(i,j,kp)
                              - grad(0,i,j,kp)*v[0] - grad(1,i,j,kp)*v[1])));
#endif
    };
    cti::loop_ijk<cti::num_lev_pack>(kv, f_v);
    if (static_cast<int>(cti::num_lev_pack) ==
        static_cast<int>(cti::max_num_lev_pack)) {
      // Re-zero eta_dot_dpdn at bottom.
      RNlevp edds(cti::pack2real(edd));
      const auto f = [&] (const int idx) {
        const int i = idx / NP, j = idx % NP;
        const int bottom = cti::num_phys_lev;
        edds(i,j,bottom) = 0;
      };
      parallel_for(ttr, f);
    }
  }

  reconstruct_and_limit_dp(kv, dp3d, dt, dp_tol, *eta_dot_dpdn[0], dprecon);
}

/* Calculate the trajectory at second order using Taylor series expansion. Also
   DSS the vertical velocity data if running the 3D algorithm.

   Derivation:
       p is position, v velocity
       p0 = p1 - dt/2 (v(p0,t0) + v(p1,t1)) + O(dt^3)
       O((p1 - p0)^2) = O(dt^2)
       p0 - p1 = -dt v(p1,t1) + O(dt^2)
       v(p0,t0) = v(p1,t0) + grad v(p1,t0) (p0 - p1) + O((p0 - p1)^2)
                = v(p1,t0) + grad v(p1,t0) (p0 - p1) + O(dt^2)
                = v(p1,t0) - dt grad v(p1,t0) v(p1,t1) + O(dt^2)
       p0 = p1 - dt/2 (v(p0,t0) + v(p1,t1)) + O(dt^3)
          = p1 - dt/2 (v(p1,t0) - dt grad v(p1,t0) v(p1,t1) + v(p1,t1)) + O(dt^3)
          = p1 - dt/2 (v(p1,t0) + v(p1,t1) - dt grad v(p1,t0) v(p1,t1)) + O(dt^3)
   In the code, v(p1,t0) = vstar, v(p1,t1) is vn0.
 */
void ComposeTransportImpl::calc_trajectory (const int np1, const Real dt) {
  GPTLstart("compose_calc_trajectory");
  const auto sphere_ops = m_sphere_ops;
  const auto geo = m_geometry;
  const auto m_vec_sph2cart = geo.m_vec_sph2cart;
  const auto m_vstar = m_derived.m_vstar;
  const auto tu_ne = m_tu_ne;
  { // Calculate midpoint velocity.
    const auto buf1a = m_data.buf1[0]; const auto buf1b = m_data.buf1[1];
    const auto buf1c = m_data.buf1[2]; const auto buf2a = m_data.buf2[0];
    const auto buf2b = m_data.buf2[1];
    const auto m_spheremp = geo.m_spheremp;
    const auto m_rspheremp = geo.m_rspheremp;
    const auto m_v = m_state.m_v;
    const auto m_vn0 = m_derived.m_vn0;
    const auto independent_time_steps = m_data.independent_time_steps;
    const Real dp_tol = m_data.dp_tol;
    const Real ps0 = m_hvcoord.ps0;
    const Real hybrid_ai0 = m_hvcoord.hybrid_ai0;
    const auto hybrid_bi = m_hvcoord.hybrid_bi_packed;
    const auto m_dp3d = m_state.m_dp3d;
    const auto m_dp = m_derived.m_dp;
    const auto m_divdp = m_derived.m_divdp;
    if (m_data.independent_time_steps) {
      GPTLstart("compose_3d_levels");
      const auto copy_v = KOKKOS_LAMBDA (const int idx) {
        int ie, lev, i, j;
        cti::idx_ie_packlev_ij(idx, ie, lev, i, j);
        for (int d = 0; d < 2; ++d)
          m_vn0(ie,d,i,j,lev) = m_v(ie,np1,d,i,j,lev);
      };
      launch_ie_packlev_ij(copy_v);
      Kokkos::fence();
      const auto calc_dprecon = KOKKOS_LAMBDA (const MT& team) {
        KernelVariables kv(team, tu_ne);
        const auto ie = kv.ie;
        const auto vn0 = Homme::subview(m_vn0, ie);
        const auto vstar = Homme::subview(m_vstar, ie);
        const auto dprecon = Homme::subview(m_divdp, ie);
        calc_vertically_lagrangian_levels(
          sphere_ops, kv, ps0, hybrid_ai0, hybrid_bi,
          Homme::subview(m_dp3d, ie, np1), Homme::subview(m_dp, ie),
          vn0, vstar, dt, dp_tol,
          Homme::subview(buf1a, kv.team_idx), Homme::subview(buf1b, kv.team_idx),
          Homme::subview(buf1c, kv.team_idx), Homme::subview(buf2a, kv.team_idx),
          dprecon);
      };
      Kokkos::parallel_for(m_tp_ne, calc_dprecon);
      Kokkos::fence();
      remap_v(m_dp3d, np1, m_divdp, m_vn0);
      const auto sphere = KOKKOS_LAMBDA (const MT& team) {
        KernelVariables kv(team, tu_ne);
        const auto ie = kv.ie;
        const auto spheremp = Homme::subview(m_spheremp, ie);
        const auto rspheremp = Homme::subview(m_rspheremp, ie);
        const auto f = [&] (const int i, const int j, const int kp) {
          const auto dprecon = Homme::subview(m_divdp, ie);
          dprecon(i,j,kp) = dprecon(i,j,kp)*spheremp(i,j)*rspheremp(i,j);
        };
        cti::loop_ijk<num_lev_pack>(kv, f);
      };
      Kokkos::fence();
      Kokkos::parallel_for(m_tp_ne, sphere);
      Kokkos::fence();
      GPTLstop("compose_3d_levels");
    }
    GPTLstart("compose_v_bexchv");
    const auto calc_midpoint_velocity = KOKKOS_LAMBDA (const MT& team) {
      KernelVariables kv(team, tu_ne);
      const auto ie = kv.ie;
      const auto vn0 = (independent_time_steps ?
                        Homme::subview(m_vn0, ie) :
                        Homme::subview(m_v, ie, np1));
      const auto vstar = Homme::subview(m_vstar, ie);
      const auto spheremp = Homme::subview(m_spheremp, ie);
      const auto rspheremp = Homme::subview(m_rspheremp, ie);
      const auto ugradv = S2Nlev(Homme::subview(buf2b, kv.team_idx).data());
      ugradv_sphere(sphere_ops, kv, Homme::subview(m_vec_sph2cart, ie), vn0, vstar,
                    SNlev(Homme::subview(buf1a, kv.team_idx).data()),
                    S2Nlev(Homme::subview(buf2a, kv.team_idx).data()),
                    ugradv);
      // Write the midpoint velocity to vstar.
      const auto f = [&] (const int i, const int j, const int k) {
        for (int d = 0; d < 2; ++d)
          vstar(d,i,j,k) = (((vn0(d,i,j,k) + vstar(d,i,j,k))/2 - dt*ugradv(d,i,j,k)/2)*
                            spheremp(i,j)*rspheremp(i,j));
      };
      cti::loop_ijk<num_lev_pack>(kv, f);
    };
    Kokkos::parallel_for(m_tp_ne, calc_midpoint_velocity);
  }
  { // DSS velocity.
    Kokkos::fence();
    const auto be = m_v_dss_be[m_data.independent_time_steps ? 1 : 0];
    be->exchange();
    Kokkos::fence();
  }
  GPTLstop("compose_v_bexchv");
  { // Calculate departure point.
    GPTLstart("compose_v2x");
    const int packn = this->packn;
    const int num_phys_lev = this->num_phys_lev;
    const auto m_sphere_cart = geo.m_sphere_cart;
    const auto scale_factor = geo.m_scale_factor;
    const auto m_dep_pts = m_data.dep_pts;
    const auto calc_departure_point = KOKKOS_LAMBDA (const MT& team) {
      KernelVariables kv(team, tu_ne);
      const auto ie = kv.ie;
      const auto vstar = Homme::subview(m_vstar, ie);
      const auto vec_sphere2cart = Homme::subview(m_vec_sph2cart, ie);
      const auto sphere_cart = Homme::subview(m_sphere_cart, ie);
      const auto dep_pts = Homme::subview(m_dep_pts, ie);
      const auto f = [&] (const int i, const int j, const int k) {
        // dp = p1 - dt v/scale_factor
        Scalar dp[3];
        for (int d = 0; d < 3; ++d) {
          const auto vel_cart = (vec_sphere2cart(0,d,i,j)*vstar(0,i,j,k) +
                                 vec_sphere2cart(1,d,i,j)*vstar(1,i,j,k));
          dp[d] = sphere_cart(i,j,d) - dt*vel_cart/scale_factor;
        }
        const auto r2 = square(dp[0]) + square(dp[1]) + square(dp[2]);
        // Pack -> scalar storage.
        const auto os = packn*k;
        for (int s = 0; s < packn; ++s) {
          const auto oss = os + s;
          if (num_phys_lev % packn != 0 && // compile out this conditional when possible
              oss >= num_phys_lev) break;
          // No vec call for sqrt.
          const auto r = std::sqrt(r2[s]);
          for (int d = 0; d < 3; ++d)
            dep_pts(oss,i,j,d) = dp[d][s]/r;
        }
      };
      cti::loop_ijk<num_lev_pack>(kv, f);
    };
    Kokkos::parallel_for(m_tp_ne, calc_departure_point);
    Kokkos::fence();
    GPTLstop("compose_v2x");
  }
  GPTLstop("compose_calc_trajectory");
}

static int test_approx_derivative () {
  const Real a = 1.5, b = -0.7, c = 0.2;
  int nerr = 0;
  const auto policy = Homme::get_default_team_policy<ExecSpace>(1);
  ExecView<Scalar[NP][NP][NUM_LEV_P]> xp("xp"), yp("yp");
  ExecView<Real[NP][NP][NUM_LEV_P*VECTOR_SIZE]>
    xs(cti::pack2real(xp)), ys(cti::pack2real(yp));
  { // Fill xs and ys with manufactured coordinates and function.
    const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
      KernelVariables kv(team);
      const auto f = [&] (const int i, const int j, const int k) {
        const Real x = 1.7*(k + 1e-1*i*k + 1e-2*j*k*k)/cti::num_phys_lev;
        xs(i,j,k) = x;
        ys(i,j,k) = (a*x + b)*x + c;
      };
      cti::loop_ijk<cti::num_phys_lev+1>(kv, f);
    };
    Kokkos::parallel_for(policy, f);
  }
  ExecView<Scalar[NP][NP][NUM_LEV]> yip("yp");
  { // Run approx_derivative.
    const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
      KernelVariables kv(team);
      approx_derivative(kv, xp, yp, yip);
    };
    Kokkos::fence();
    Kokkos::parallel_for(policy, f);
  }
  { // Check answer.
    Kokkos::fence();
    const auto xsh = cti::cmvdc(xs);
    ExecView<Real[NP][NP][NUM_LEV*VECTOR_SIZE]> yis(cti::pack2real(yip));
    const auto yish = cti::cmvdc(yis);
    for (int i = 0; i < cti::np; ++i)
      for (int j = 0; j < cti::np; ++j)
        for (int k = 1; // k = 0 is not written
             k < cti::num_phys_lev; ++k) {
          const Real x = xsh(i,j,k);
          const Real ypp = 2*a*x + b;
          const Real err = std::abs(yish(i,j,k) - ypp);
          if (err > 1e3*std::numeric_limits<Real>::epsilon()) {
            ++nerr;
            printf("%2d %d %d %1.2f %6.2f %6.2f %9.2e\n",
                   i, j, k, x, ypp, yish(i,j,k), err);
          }
        }
  }
  if (nerr) printf("test_approx_derivative failed %d\n", nerr);
  return nerr;
}

static int test_calc_p (const HybridVCoord& hvcoord) {
  int nerr = 0;
  ExecView<Scalar[NP][NP][NUM_LEV]> dpp("dpp");
  ExecView<Real[NP][NP][NUM_LEV*VECTOR_SIZE]> dp(cti::pack2real(dpp));
  const auto dph = Kokkos::create_mirror_view(dp);
  for (int i = 0; i < cti::np; ++i)
    for (int j = 0; j < cti::np; ++j)
      for (int k = 0; k < cti::num_phys_lev; ++k)
        dph(i,j,k) = k + 3*i + 7*j;
  Kokkos::deep_copy(dp, dph);
  const auto hybrid_ai = cti::cmvdc(hvcoord.hybrid_ai_packed);
  const Real hybrid_ai0 = hybrid_ai(0)[0];
  const Real ps0 = hvcoord.ps0;
  ExecView<Scalar[NP][NP][NUM_LEV_P]> pp("pp");
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    calc_p(kv, hybrid_ai0, ps0, dpp, pp);
  };
  const auto policy = Homme::get_default_team_policy<ExecSpace>(1);
  Kokkos::parallel_for(policy, f);
  Kokkos::fence();
  ExecView<Real[NP][NP][NUM_LEV_P*VECTOR_SIZE]> p(cti::pack2real(pp));
  const auto ph = cti::cmvdc(p);
  for (int i = 0; i < cti::np; ++i)
    for (int j = 0; j < cti::np; ++j) {
      Real accum = hvcoord.hybrid_ai0*hvcoord.ps0;
      for (int k = 0; k < cti::num_phys_lev; ++k)
        accum += dph(i,j,k);
      if (ph(i,j,cti::num_phys_lev) != accum) ++nerr;
    }
  if (nerr) printf("test_calc_p failed %d\n", nerr);
  return nerr;
}

int test_reconstruct_and_limit_dp () {
  const Real eps = std::numeric_limits<Real>::epsilon(),
    dt = 1800, p0 = 1e5, dp_tol = (p0/72)*eps, tol = 1e3*eps;
  ExecView<Scalar[NP][NP][NUM_LEV]> dprefp("dprefp"), dp1limp("dp1limp");
  HostView<Real[NP][NP][NUM_PHYSICAL_LEV]> dp1ul("dp1ul");
  ExecView<Scalar[NP][NP][NUM_LEV_P]> eta_dot_dpdn_p("eta_dot_dpdn_p");
  RNlev dpref(cti::pack2real(dprefp)), dp1lim(cti::pack2real(dp1limp));
  RNlevp eta_dot_dpdn(cti::pack2real(eta_dot_dpdn_p));
  const auto dprefh = Kokkos::create_mirror_view(dpref);
  { // fill dpref, eta_dot_dpdn, and the resulting unlimited dp
    const auto edd = Kokkos::create_mirror_view(eta_dot_dpdn);
    for (int i = 0; i < cti::np; ++i)
      for (int j = 0; j < cti::np; ++j) {
        Real fac = 0;
        for (int k = 0; k < cti::np; ++k) {
          dprefh(i,j,k) = k+1 + 0.1*(i+1) + 0.01*j*j;
          fac += dprefh(i,j,k);
          edd(i,j,k) = (k % 2 == 0 ? -1 : 1)*0.1*(cti::num_phys_lev - k - 1);
        }
        fac = p0/fac;
        for (int k = 0; k < cti::np; ++k) dprefh(i,j,k) *= fac;
        for (int k = 0; k < cti::np; ++k) edd(i,j,k) *= fac;
        edd(i,j,0) = edd(i,j,cti::num_phys_lev) = 0;
      }
    for (int i = 0; i < cti::np; ++i)
      for (int j = 0; j < cti::np; ++j)
        for (int k = 0; k < cti::np; ++k)
          dp1ul(i,j,k) = dprefh(i,j,k) + dt*(edd(i,j,k+1) - edd(i,j,k));
    Kokkos::deep_copy(dpref, dprefh);
    Kokkos::deep_copy(eta_dot_dpdn, edd);
  }
  { // call function
    const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
      KernelVariables kv(team);
      reconstruct_and_limit_dp(kv, dprefp, dt, dp_tol, eta_dot_dpdn_p, dp1limp);
    };
    const auto policy = Homme::get_default_team_policy<ExecSpace>(1);
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
  }
  int nerr = 0;
  { // tests
    const auto dp1l = cti::cmvdc(dp1lim);
    for (int i = 0; i < cti::np; ++i)
      for (int j = 0; j < cti::np; ++j) {
        // mass conservation
        Real mass_ref = 0, mass_ul = 0, mass_lim = 0;
        for (int k = 0; k < cti::num_phys_lev; ++k) mass_ref += dprefh(i,j,k);
        for (int k = 0; k < cti::num_phys_lev; ++k) mass_ul += dp1ul(i,j,k);
        for (int k = 0; k < cti::num_phys_lev; ++k) mass_lim += dp1l(i,j,k);
        if (std::abs(mass_ul  - mass_ref) > tol*mass_ref) ++nerr;
        if (std::abs(mass_lim - mass_ref) > tol*mass_ref) ++nerr;
        // limiter needs to be active
        bool found = false;
        for (int k = 0; k < cti::num_phys_lev; ++k)
          if (dp1ul(i,j,k) < 0) found = true;
        if ( ! found) ++nerr;
        // limiter succeeded
        for (int k = 0; k < cti::num_phys_lev; ++k)
          if (dp1l(i,j,k) < 0) ++nerr;
      }
  }
  return nerr;
}

int ComposeTransportImpl::run_trajectory_unit_tests () {
  return (test_approx_derivative() + test_calc_p(m_hvcoord) +
          test_reconstruct_and_limit_dp());
}

ComposeTransport::TestDepView::HostMirror ComposeTransportImpl::
test_trajectory (Real t0, Real t1, const bool independent_time_steps) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;

  const bool its_save = m_data.independent_time_steps;
  m_data.independent_time_steps = independent_time_steps;

  const auto vstar = create_mirror_view(m_derived.m_vstar);
  const auto vn0 = create_mirror_view(m_derived.m_vn0);
  const auto v = create_mirror_view(m_state.m_v);
  const auto dp3d = create_mirror_view(m_state.m_dp3d);
  const auto dp = create_mirror_view(m_derived.m_dp);
  const auto pll = cmvdc(m_geometry.m_sphere_latlon);
  const auto np1 = 0;
  const int packn = this->packn;
  const compose::test::NonDivergentWindField wf;
  // On host b/c trig isn't BFB between host and device.
  const auto f = [&] (const int ie, const int lev, const int i, const int j) {
    Real latlon[] = {pll(ie,i,j,0), pll(ie,i,j,1)};
    compose::test::offset_latlon(num_phys_lev, lev, latlon[0], latlon[1]);
    Real uv[2];
    wf.eval(t0, latlon, uv);
    const int p = lev/packn, s = lev%packn;
    for (int d = 0; d < 2; ++d) vstar(ie,d,i,j,p)[s] = uv[d];
    wf.eval(t1, latlon, uv);
    for (int d = 0; d < 2; ++d) v(ie,np1,d,i,j,p)[s] = uv[d];
    for (int t = 0; t < NUM_TIME_LEVELS; ++t) dp3d(ie,t,i,j,p)[s] = 1;
    dp(ie,i,j,p)[s] = 1;
  };
  loop_host_ie_plev_ij(f);
  deep_copy(m_derived.m_vstar, vstar);
  deep_copy(m_state.m_v, v);
  deep_copy(m_state.m_dp3d, dp3d);
  deep_copy(m_derived.m_dp, dp);

  calc_trajectory(np1, t1 - t0);
  Kokkos::fence();

  m_data.independent_time_steps = its_save;
  const auto deph = cti::cmvdc(m_data.dep_pts);
  return deph;
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
