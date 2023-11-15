#include "compose_slmm_islmpi.hpp"
#include "compose_slmm_islmpi_interpolate.hpp"
#include "compose_slmm_islmpi_buf.hpp"

namespace homme {
namespace islmpi {

template <typename T> using CA4 = ko::View<T****, ko::LayoutRight, ko::HostSpace>;

template <Int np, typename EtamT, typename VnodeT> SLMM_KF void
interpolate_vertical (const Int nlev, const Real etai_beg, const Real etai_end,
                      const EtamT& etam, const VnodeT& vnode,
                      const Int src_lid, const Int lev, const Real eta_dep,
                      const Real rx[np], const Real ry[np], Real* const v_tgt) {
  slmm_kernel_assert(eta_dep > etai_beg && eta_dep < etai_end);
  
  // Search for the eta midpoint values that support the departure point's eta
  // value.
  Int lev_dep = lev;
  if (eta_dep != etam(lev)) {
    if (eta_dep < etam(lev)) {
      for (lev_dep = lev-1; lev_dep >= 0; --lev_dep)
        if (eta_dep >= etam(lev_dep))
          break;
    } else {
      for (lev_dep = lev; lev_dep < nlev-1; ++lev_dep)
        if (eta_dep < etam(lev_dep+1))
          break;
    }
  }
  slmm_kernel_assert(lev_dep >= -1 && lev_dep < nlev);
  slmm_kernel_assert(lev_dep == -1 || eta_dep >= etam(lev_dep));
  Real a;
  bool bdy = false;
  if (lev_dep == -1) {
    lev_dep = 0;
    a = 0;
    bdy = true;
  } else if (lev_dep == nlev-1) {
    a = 0;
    bdy = true;
  } else {
    a = ((eta_dep - etam(lev_dep)) /
         (etam(lev_dep+1) - etam(lev_dep)));
  }
  // Linear interp coefficients.
  const Real alpha[] = {1-a, a};

  for (int d = 0; d < 4; ++d)
    v_tgt[d] = 0;
  for (int i = 0; i < 2; ++i) {
    if (alpha[i] == 0) continue;
    for (int d = 0; d < 4; ++d) {
      Real vel_nodes[np*np];
      for (int k = 0; k < np*np; ++k)
        vel_nodes[k] = vnode(src_lid,lev_dep+i,k,d);
      v_tgt[d] += alpha[i]*calc_q_tgt(rx, ry, vel_nodes);
    }
  }
  // Treat eta_dot specially since eta_dot goes to 0 at the boundaries.
  if (bdy) {
    slmm_kernel_assert(etam(0) > etai_beg);
    slmm_kernel_assert(etam(nlev-1) < etai_end);
    if (lev_dep == 0)
      v_tgt[3] *= (eta_dep - etai_beg)/(etam(0) - etai_beg);
    else
      v_tgt[3] *= (etai_end - eta_dep)/(etai_end - etam(nlev-1));
  }
}

template <Int np, typename VnodeT, typename MT>
void calc_v (const IslMpi<MT>& cm, const VnodeT& vnode,
             const Int src_lid, const Int lev,
             const Real* const dep_point, Real* const v_tgt) {
  // Horizontal interpolation.
  Real rx[np], ry[np]; {
    Real ref_coord[2];
    const auto& m = cm.advecter->local_mesh(src_lid);
    cm.advecter->s2r().calc_sphere_to_ref(src_lid, m, dep_point,
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
  interpolate_vertical<np>(cm.nlev, cm.etai_beg, cm.etai_end, cm.etam, vnode,
                           src_lid, lev, dep_point[3], rx, ry, v_tgt);
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
  const Real etai_beg, etai_end;
  const typename IslMpi<MT>::template ArrayD<Real*> etam;

  CalcVData (const IslMpi<MT>& cm)
    : local_meshes(cm.advecter->local_meshes()),
      interp_alg(cm.advecter->alg()),
      s2r(cm.advecter->s2r()),
      traj_3d(cm.traj_3d),
      dep_points_ndim(cm.dep_points_ndim),
      nlev(cm.nlev),
      etai_beg(cm.etai_beg), etai_end(cm.etai_end),
      etam(cm.etam)
  {}
};

template <Int np, typename VnodeT, typename MT> SLMM_KF
void calc_v (const CalcVData<MT>& cvd, const VnodeT& vnode,
             const Int src_lid, const Int lev,
             const Real* const dep_point, Real* const v_tgt) {
  // Horizontal interpolation.
  Real rx[np], ry[np]; {
    Real ref_coord[2];
    const auto& m = cvd.local_meshes(src_lid);
    cvd.s2r.calc_sphere_to_ref(src_lid, m, dep_point,
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
  interpolate_vertical<np>(cvd.nlev, cvd.etai_beg, cvd.etai_end, cvd.etam, vnode,
                           src_lid, lev, dep_point[3], rx, ry, v_tgt);
}

template <int np, typename VnodeT, typename MT>
void traj_calc_rmt_next_step (IslMpi<MT>& cm, const VnodeT& vnode) {
  calc_rmt_q_pass1(cm, true);
  const auto ndim = cm.dep_points_ndim;
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
      xos = rmt_xs(5*it + 3), vos = ndim*rmt_xs(5*it + 4);
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
  const auto ndim = cm.dep_points_ndim;
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
    Real v_tgt[4];
    calc_v<np>(cvd, vnode, slid, tgt_lev, &dep_points(tci,tgt_lev,tgt_k,0), v_tgt);
    for (int d = 0; d < ndim; ++d)
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
      Real v_tgt[4];
      calc_v<np>(cm, vnode, slid, e.lev, &dep_points(tci,e.lev,e.k,0), v_tgt);
      for (int d = 0; d < ndim; ++d)
        vdep(tci,e.lev,e.k,d) = v_tgt[d];
    }
  }
#endif
}

template <typename VdepT, typename MT>
void traj_copy_next_step (IslMpi<MT>& cm, const VdepT& vdep) {
  const auto myrank = cm.p->rank();
  const auto ndim = cm.dep_points_ndim;
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
    for (int d = 0; d < ndim; ++d)
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
      for (int d = 0; d < ndim; ++d)
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
calc_v_departure (IslMpi<MT>& cm, const Int nets, const Int nete,
                  const Int step, const Real dtsub,
                  Real* dep_points_r, const Real* vnode_r, Real* vdep_r)
{
  const int np = 4;

  slmm_assert(cm.np == np);
  slmm_assert((cm.traj_3d and cm.dep_points_ndim == 4) or
              (not cm.traj_3d and cm.dep_points_ndim == 3));
#ifdef COMPOSE_PORT
  slmm_assert(nets == 0 && nete+1 == cm.nelemd);
#endif

  // If step = 0, the departure points are at the nodes and no interpolation is
  // needed. calc_v_departure should not have been called; rather, the calling
  // routine should use vnode instead of vdep in subsequent calculations.
  slmm_assert(step > 0);

  const auto ndim = cm.dep_points_ndim;

#ifdef COMPOSE_PORT
  const auto& vnode = cm.tracer_arrays->vnode;
  const auto& vdep  = cm.tracer_arrays->vdep;
#else
  CA4<const Real> vnode(vnode_r, cm.nelemd, cm.nlev, cm.np2, ndim);
  CA4<      Real> vdep (vdep_r , cm.nelemd, cm.nlev, cm.np2, ndim);
#endif
  slmm_assert(vnode.extent_int(3) == ndim);
  slmm_assert(vdep .extent_int(3) == ndim);

#ifdef COMPOSE_PORT
  const auto& dep_points = cm.tracer_arrays->dep_points;
#else
  DepPointsH<MT> dep_points(dep_points_r, cm.nelemd, cm.nlev, cm.np2, ndim);
#endif
  slmm_assert(dep_points.extent_int(3) == ndim);

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

template void calc_v_departure(
  IslMpi<ko::MachineTraits>&, const Int, const Int, const Int, const Real,
  Real*, const Real*, Real*);

} // namespace islmpi
} // namespace homme
