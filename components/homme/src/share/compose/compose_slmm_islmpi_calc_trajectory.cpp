/* This file contains the driver routine interp_v_update. It interpolates
   reference-grid velocities to off-grid arrival points. This situation occurs
   when the enhanced trajectory method is run with nsubstep >= 2. In the first
   substep, the arrival points are on-grid; thus, interp_v_update is not
   needed. Subsequent substeps start with off-grid arrival points. These are the
   departure points from the previous substep.
     interp_v_update follows the same computational and communication patterns
   as the main 'step' routine (compose_slmm_islmpi_step.cpp). Both
   interp_v_update and 'step' are essentially computing interpolants. In the
   case of 'step', the mixing ratios q are interpolated; in the case of
   interp_v_update, (xdot, ydot, zdot, etadot) are interpolated, where (xdot,
   ydot, zdot) is the horizontal velocity on the sphere and etadot is the
   vertical velocity. In this sense, interp_v_update can be understood as a
   specialization of 'step' to this case.
     However, there is one bit of complexity in interp_v_update that requires
   further explanation.
     Let the horizontal position (velocity) be abbreviated as h(dot) := (x(dot),
   y(dot), z(dot)), where h(dot) means the statement applies to both h and
   hdot. h(dot) lives on vertical midpoints. In contrast, the eta(dot) part of
   the trajectory lives on vertical interfaces. A previous version of this
   method interpolated etadot to midpoints, then interpolated the vertical part
   of departure points from midpoints to interfaces. This caused issues in
   certain chemistry parameterizations around the tropopause. The current
   version keeps eta and etadot at interfaces, which adds complexity since
   (h(dot), eta(dot)) are at a mix of midpoints and interfaces. This change
   removes the dependence of the reconstruction of floating Lagrangian levels on
   the vertical grid in some (but, unavoidably, not all) terms of the
   approximation. Simply stated, staying on interfaces as much as possible,
   where eta(dot) naturally lives, is better than moving to midpoints, where
   eta(dot) has to be interpolated.
     The high-level key ideas are as follows. Let _ref denote quantities on the
   reference grid; _arr (for "arrival point"), quantities on the vertically
   Lagrangian grid.
     Interpolate (hdot_ref, etadot_ref) at (h_arr, midpoint(eta_arr)) to get
   hdot_arr. In this step, eta_arr is interpolated to midpoints. But this step
   does not compute etadot_arr.
     Next, we need to interpolate for etadot_arr. The simplest procedure would
   be to compute extra full trajectories at interfaces. But this would be
   expensive. Instead, we reuse the midpoint trajectories as follows.
     At each interface k, interpolate for etadot_k_ref at (h_(k-1)_arr,
   eta_k_arr) and (h_k_arr, eta_k_arr). Then combine these two values to get the
   final etadot_k_arr. The key here is that only interface etadot_ref values are
   used in computing these steps; a derived midpoint etadot_arr is never used.
     This calculation is a bit subtle, but it's set up such that it adds very
   little cost to the much simpler calculation that uses eta(dot) at
   midpoints. It adds two additional entries to the communicated arrays and a
   small amount of extra computation.
     Again, the conceptually simplest and most accurate approach would double
   the number of trajectories to compute. The approach above is one
   approximation to this. Other approximations are possible, but the one above
   seems to use the least computation and communication among the available
   approximations.
     The following notes provide precise implementation details about this
   calculation, outlining the code in the rest of this file.

     Level arrangement. x is horizontal position. A suffixed 'd' means 'dot'.
          0 i etai(0), etaid(0) = 0
          0 m x or xdot (= xd)
          1 i e or ed
          1 m x(d)
          2 i e(d)
          2 m x(d)
       nlev i
     Algorithm:
       x[k] is the horizontal position at midpoint k
       analyze_dep_points for x[k]
       send/recv p[k] = (x[k], e[k], e[k+1])
         where e is at interfaces
       compute (rx,ry) using x[k]
       interp3d for xd[k]   at (x[k], (e[k] + e[k+1])/2)
         using eta_ref_mid
       interp3d for ed[k]   at (x[k], e[k])
         using eta_ref_int if k > 0, else 0
       interp3d for ed[k+1] at (x[k], e[k+1])
         using eta_ref_int if k+1 < nlev, else 0
       send/recv pd[k] = (xd[k], ed[k], ed[k+1])
       vdep[k,0:3] = pd[k].xd[k]
       vdep[k,3] = linterp(eta_ref_mid[k-1:k+1],
                           (pd[k-1].ed[k], pd[k].ed[k]),
                           eta_ref_int[k])
         if k > 0, else 0
 */

#include "compose_slmm_islmpi.hpp"
#include "compose_slmm_islmpi_interpolate.hpp"
#include "compose_slmm_islmpi_buf.hpp"

namespace homme {
namespace islmpi {

template <typename T> using CA4 = ko::View<T****, ko::LayoutRight, ko::HostSpace>;

// Interpolate at eta interfaces, rather than midpoints as below.
template <Int np, typename EtaT, typename VnodeT> SLMM_KF
Real interpolate_at_vertical_interfaces (
  const Int nlev, const EtaT& etai, const VnodeT& vnode, const Int src_lid,
  const Int lev, const Real rx[np], const Real ry[np], Real etai_dep)
{
  if (etai_dep < etai[0   ]) etai_dep = etai[0   ];
  if (etai_dep > etai[nlev]) etai_dep = etai[nlev];

  // Search for the eta interface values that support etai_dep.
  Int lev_dep = lev;
  if (etai_dep != etai(lev)) {
    if (etai_dep < etai(lev)) {
      for (lev_dep = lev-1; lev_dep >= 0; --lev_dep)
        if (etai_dep >= etai(lev_dep))
          break;
    } else {
      for (lev_dep = lev; lev_dep < nlev; ++lev_dep)
        if (etai_dep < etai(lev_dep+1))
          break;
    }
  }
  slmm_kernel_assert(lev_dep >= 0 and lev_dep <= nlev);
  if (lev_dep == nlev) return 0; // etai_dep == etai_end => eta_dot = 0
  const Real a = (etai_dep - etai(lev_dep)) / (etai(lev_dep+1) - etai(lev_dep));
  // Linear interp coefficients.
  const Real alpha[] = {1-a, a};
  const int dim = 3;
  Real etai_dot = 0;
  for (int i = 0; i < 2; ++i) {
    if (lev_dep+i == 0 or lev_dep+i == nlev)
      continue; // vel_nodes[:] = 0
    slmm_kernel_assert(lev_dep+i >= 0 and lev_dep+i < nlev);
    Real vel_nodes[np*np];
    for (int k = 0; k < np*np; ++k)
      vel_nodes[k] = vnode(src_lid,lev_dep+i,k,dim);
    etai_dot += alpha[i]*calc_q_tgt(rx, ry, vel_nodes);
  }
  return etai_dot;
}

template <Int np, typename EtaT, typename VnodeT> SLMM_KF void
interpolate_vertical (const Int nlev, const EtaT& etai, const EtaT& etam,
                      const VnodeT& vnode, const Int src_lid, const Int lev,
                      const Real etai_lev, const Real etai_levp1,
                      const Real rx[np], const Real ry[np],
                      Real* const v_tgt) {
  Real eta_mid_dep = (etai_lev + etai_levp1)/2;
  if (eta_mid_dep < etai[0   ]) eta_mid_dep = etai[0   ];
  if (eta_mid_dep > etai[nlev]) eta_mid_dep = etai[nlev];
  
  // Search for the eta midpoint values that support the departure point's eta
  // value.
  Int lev_dep = lev;
  if (eta_mid_dep != etam(lev)) {
    if (eta_mid_dep < etam(lev)) {
      for (lev_dep = lev-1; lev_dep >= 0; --lev_dep)
        if (eta_mid_dep >= etam(lev_dep))
          break;
    } else {
      for (lev_dep = lev; lev_dep < nlev-1; ++lev_dep)
        if (eta_mid_dep < etam(lev_dep+1))
          break;
    }
  }
  slmm_kernel_assert(lev_dep >= -1 and lev_dep < nlev);
  slmm_kernel_assert(lev_dep == -1 or eta_mid_dep >= etam(lev_dep));
  Real a;
  if (lev_dep == -1) {
    lev_dep = 0;
    a = 0;
  } else if (lev_dep == nlev-1) {
    a = 0;
  } else {
    a = ((eta_mid_dep - etam(lev_dep)) /
         (etam(lev_dep+1) - etam(lev_dep)));
  }
  // Linear interp coefficients.
  const Real alpha[] = {1-a, a};

  for (int d = 0; d < 4; ++d)
    v_tgt[d] = 0;
  for (int i = 0; i < 2; ++i) {
    if (alpha[i] == 0) continue;
    const int ndim = 3;
    for (int d = 0; d < ndim; ++d) {
      Real vel_nodes[np*np];
      for (int k = 0; k < np*np; ++k)
        vel_nodes[k] = vnode(src_lid,lev_dep+i,k,d);
      v_tgt[d] += alpha[i]*calc_q_tgt(rx, ry, vel_nodes);
    }
  }

  v_tgt[3] = interpolate_at_vertical_interfaces<np>(
    nlev, etai, vnode, src_lid, lev, rx, ry, etai_lev);
  v_tgt[4] = interpolate_at_vertical_interfaces<np>(
    nlev, etai, vnode, src_lid, lev, rx, ry, etai_levp1);
}

template <Int np, typename VnodeT, typename MT>
void calc_v (const IslMpi<MT>& cm, const VnodeT& vnode,
             const Int src_lid, const Int lev,
             const Real* const dep, Real* const v_tgt) {
  // Horizontal interpolation.
  Real rx[np], ry[np]; {
    Real ref_coord[2];
    const auto& m = cm.advecter->local_mesh(src_lid);
    cm.advecter->s2r().calc_sphere_to_ref(src_lid, m, dep,
                                          ref_coord[0], ref_coord[1]);
    interpolate<MT>(cm.advecter->alg(), ref_coord, rx, ry);
  }

  if (not cm.traj_3d) {
    for (int d = 0; d < cm.dep_points_ndim; ++d) {
      Real vel_nodes[np*np];
      for (int k = 0; k < np*np; ++k)
        vel_nodes[k] = vnode(src_lid,lev,k,d);
      v_tgt[d] = calc_q_tgt(rx, ry, vel_nodes);
    }
    return;
  }

  // Vertical Interpolation.
  slmm_kernel_assert(cm.dep_points_ndim == 4);
  interpolate_vertical<np>(cm.nlev, cm.etai, cm.etam, vnode,
                           src_lid, lev, dep[3], dep[4], rx, ry, v_tgt);
}

template <typename MT>
struct CalcVData {
  typedef slmm::Advecter<MT> Adv;
  const typename Adv::LocalMeshesD local_meshes;
  const typename Adv::Alg::Enum interp_alg;
  const slmm::SphereToRef<typename MT::DES> s2r;
  const bool traj_3d;
  const int dep_points_ndim;
  const int nlev;
  const typename IslMpi<MT>::template ArrayD<Real*> etai, etam;

  CalcVData (const IslMpi<MT>& cm)
    : local_meshes(cm.advecter->local_meshes()),
      interp_alg(cm.advecter->alg()),
      s2r(cm.advecter->s2r()),
      traj_3d(cm.traj_3d),
      dep_points_ndim(cm.dep_points_ndim),
      nlev(cm.nlev),
      etai(cm.etai), etam(cm.etam)
  {}
};

template <Int np, typename VnodeT, typename MT> SLMM_KF
void calc_v (const CalcVData<MT>& cvd, const VnodeT& vnode,
             const Int src_lid, const Int lev,
             const Real* const dep, Real* const v_tgt) {
  // Horizontal interpolation.
  Real rx[np], ry[np]; {
    Real ref_coord[2];
    const auto& m = cvd.local_meshes(src_lid);
    cvd.s2r.calc_sphere_to_ref(src_lid, m, dep,
                               ref_coord[0], ref_coord[1]);
    interpolate<MT>(cvd.interp_alg, ref_coord, rx, ry);
  }

  if (not cvd.traj_3d) {
    for (int d = 0; d < cvd.dep_points_ndim; ++d) {
      Real vel_nodes[np*np];
      for (int k = 0; k < np*np; ++k)
        vel_nodes[k] = vnode(src_lid,lev,k,d);
      v_tgt[d] = calc_q_tgt(rx, ry, vel_nodes);
    }
    return;
  }

  // Vertical Interpolation.
  slmm_kernel_assert(cvd.dep_points_ndim == 4);
  interpolate_vertical<np>(cvd.nlev, cvd.etai, cvd.etam, vnode,
                           src_lid, lev, dep[3], dep[4], rx, ry, v_tgt);
}

template <int np, typename VnodeT, typename MT>
void traj_calc_rmt_next_step (IslMpi<MT>& cm, const VnodeT& vnode) {
  calc_rmt_q_pass1(cm, true);
  const auto xsz = cm.traj_msg_sz;
  const auto& rmt_xs = cm.rmt_xs;
  const auto& sendbuf = cm.sendbuf;
  const auto& recvbuf = cm.recvbuf;
  CalcVData<MT> cvd(cm);
#ifdef COMPOSE_PORT
  ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, cm.nrmt_xs),
                   COMPOSE_LAMBDA (const Int it)
#else
# ifdef COMPOSE_HORIZ_OPENMP
#  pragma omp for
# endif
  for (Int it = 0; it < cm.nrmt_xs; ++it)
#endif
  {
    const Int
      ri = rmt_xs(5*it), lid = rmt_xs(5*it + 1), lev = rmt_xs(5*it + 2),
      xos = rmt_xs(5*it + 3), vos = xsz*rmt_xs(5*it + 4);
    const auto&& xs = recvbuf(ri);
    auto&& v = sendbuf(ri);
    calc_v<np>(cvd, vnode, lid, lev, &xs(xos), &v(vos));
  }
#ifdef COMPOSE_PORT
  );
#endif
}

template <int np, typename VnodeT, typename VdepT, typename MT>
void traj_calc_own_next_step (IslMpi<MT>& cm, const DepPoints<MT>& dep_points,
                              const VnodeT& vnode, const VdepT& vdep) {
  const auto xsz = cm.traj_msg_sz;
  const auto ndim = cm.dep_points_ndim;
  const auto etai_end = cm.etai_end;
#ifdef COMPOSE_PORT
  const auto& ed_d = cm.ed_d;
  const auto& own_dep_list = cm.own_dep_list;
  CalcVData<MT> cvd(cm);
  const auto f = COMPOSE_LAMBDA (const Int& it) {
    const Int tci = own_dep_list(it,0);
    const Int tgt_lev = own_dep_list(it,1);
    const Int tgt_k = own_dep_list(it,2);
    const auto& ed = ed_d(tci);
    const Int slid = ed.nbrs(ed.src(tgt_lev, tgt_k)).lid_on_rank;
    Real dep[5], v_tgt[5];
    for (Int d = 0; d < ndim; ++d)
      dep[d] = dep_points(tci,tgt_lev,tgt_k,d);
    dep[ndim] = (tgt_lev+1 == cvd.nlev ?
                 etai_end :
                 dep_points(tci,tgt_lev+1,tgt_k,ndim-1));
    calc_v<np>(cvd, vnode, slid, tgt_lev, dep, v_tgt);
    for (int d = 0; d < xsz; ++d)
      vdep(tci,tgt_lev,tgt_k,d) = v_tgt[d];
  };
  ko::parallel_for(
    ko::RangePolicy<typename MT::DES>(0, cm.own_dep_list_len), f);
#else
  const int tid = get_tid();
  for (Int tci = 0; tci < cm.nelemd; ++tci) {
    auto& ed = cm.ed_d(tci);
    const Int ned = ed.own.n();
#ifdef COMPOSE_HORIZ_OPENMP
#   pragma omp for
#endif
    for (Int idx = 0; idx < ned; ++idx) {
      const auto& e = ed.own(idx);
      const Int slid = ed.nbrs(ed.src(e.lev, e.k)).lid_on_rank;
      Real dep[5], v_tgt[5];
      for (Int d = 0; d < ndim; ++d)
        dep[d] = dep_points(tci,e.lev,e.k,d);
      dep[ndim] = (e.lev+1 == cm.nlev ?
                   etai_end :
                   dep_points(tci,e.lev+1,e.k,ndim-1));
      calc_v<np>(cm, vnode, slid, e.lev, dep, v_tgt);
      for (int d = 0; d < xsz; ++d)
        vdep(tci,e.lev,e.k,d) = v_tgt[d];
    }
  }
#endif
}

template <typename VdepT, typename MT>
void traj_copy_next_step (IslMpi<MT>& cm, const VdepT& vdep) {
#ifndef NDEBUG
  const auto myrank = cm.p->rank();
#endif
  const auto xsz = cm.traj_msg_sz;
#ifdef COMPOSE_PORT
  const auto& mylid_with_comm = cm.mylid_with_comm_d;
  const auto& ed_d = cm.ed_d;
  const auto& recvbufs = cm.recvbuf;
  const Int nlid = cm.mylid_with_comm_h.size();
  const Int nlev = cm.nlev, np2 = cm.np2;
  const auto f = COMPOSE_LAMBDA (const Int& it) {
    const Int tci = mylid_with_comm(it/(np2*nlev));
    const Int rmt_id = it % (np2*nlev);
    auto& ed = ed_d(tci);
    if (rmt_id >= ed.rmt.size()) return;
    const auto& e = ed.rmt(rmt_id);
    slmm_kernel_assert(ed.nbrs(ed.src(e.lev, e.k)).rank != myrank);
    const Int ri = ed.nbrs(ed.src(e.lev, e.k)).rank_idx;
    const auto&& recvbuf = recvbufs(ri);
    for (int d = 0; d < xsz; ++d)
      vdep(tci,e.lev,e.k,d) = recvbuf(e.q_ptr + d);
  };
  ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, nlid*np2*nlev), f);
#else
  const int tid = get_tid();
  for (Int ptr = cm.mylid_with_comm_tid_ptr_h(tid),
           end = cm.mylid_with_comm_tid_ptr_h(tid+1);
       ptr < end; ++ptr) {
    const Int tci = cm.mylid_with_comm_d(ptr);
    auto& ed = cm.ed_d(tci);
    for (const auto& e: ed.rmt) {
      slmm_assert(ed.nbrs(ed.src(e.lev, e.k)).rank != myrank);
      const Int ri = ed.nbrs(ed.src(e.lev, e.k)).rank_idx;
      const auto&& recvbuf = cm.recvbuf(ri);
      for (int d = 0; d < xsz; ++d)
        vdep(tci,e.lev,e.k,d) = recvbuf(e.q_ptr + d);
    }
  }
#endif
}

// vnode and vdep are indexed as (ie,lev,k,dim), On entry, vnode contains nodal
// velocity data. These data are used to provide updates at departure points for
// both own and remote departure points, writing to vdep. dim = 0:2 is for the
// 3D Cartesian representation of the horizontal velocity; dim = 3 is for
// eta_dot.
template <typename MT> void
interp_v_update (IslMpi<MT>& cm, const Int nets, const Int nete,
                 const Int step, const Real dtsub,
                 Real* dep_points_r, const Real* vnode_r, Real* vdep_r)
{
  const int np = 4;

  slmm_assert(cm.np == np);
  slmm_assert((cm.traj_3d and cm.dep_points_ndim == 4) or
              (not cm.traj_3d and cm.dep_points_ndim == 3));
#ifdef COMPOSE_PORT
  slmm_assert(nets == 0 and nete+1 == cm.nelemd);
#endif

  // If step = 0, the departure points are at the nodes and no interpolation is
  // needed. interp_v_update should not have been called; rather, the calling
  // routine should use vnode instead of vdep in subsequent calculations.
  slmm_assert(step > 0);

  const auto ndim = cm.dep_points_ndim;

#ifdef COMPOSE_PORT
  const auto& vnode = cm.tracer_arrays->vnode;
  const auto& vdep  = cm.tracer_arrays->vdep;
#else
  CA4<const Real> vnode(vnode_r, cm.nelemd, cm.nlev, cm.np2, ndim);
  CA4<      Real> vdep (vdep_r , cm.nelemd, cm.nlev, cm.np2, ndim+1);
#endif
  slmm_assert(vnode.extent_int(3) == ndim);
  slmm_assert(vdep .extent_int(3) == ndim+1);

#ifdef COMPOSE_PORT
  const auto& dep_points = cm.tracer_arrays->dep_points;
#else
  DepPointsH<MT> dep_points(dep_points_r, cm.nelemd, cm.nlev, cm.np2, ndim);
#endif
  slmm_assert(dep_points.extent_int(3) == cm.dep_points_ndim);

  // See comments in homme::islmpi::step for details. Each substep follows
  // essentially the same pattern.
  if (cm.mylid_with_comm_tid_ptr_h.capacity() == 0)
    init_mylid_with_comm_threaded(cm, nets, nete);
  setup_irecv(cm);
  analyze_dep_points(cm, nets, nete, dep_points);
  pack_dep_points_sendbuf_pass1(cm, true /* trajectory */);
  pack_dep_points_sendbuf_pass2(cm, dep_points, true /* trajectory */);
  isend(cm);
  recv_and_wait_on_send(cm);
  traj_calc_rmt_next_step<np>(cm, vnode);
  Kokkos::fence();
  isend(cm, true /* want_req */, true /* skip_if_empty */);
  setup_irecv(cm, true /* skip_if_empty */);
  traj_calc_own_next_step<np>(cm, dep_points, vnode, vdep);
  recv(cm, true /* skip_if_empty */);
  traj_copy_next_step(cm, vdep);
  wait_on_send(cm, true /* skip_if_empty */);
}

template void interp_v_update(
  IslMpi<ko::MachineTraits>&, const Int, const Int, const Int, const Real,
  Real*, const Real*, Real*);

} // namespace islmpi
} // namespace homme
