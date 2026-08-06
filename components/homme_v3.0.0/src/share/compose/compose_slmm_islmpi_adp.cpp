#include "compose_slmm_islmpi.hpp"
#include "compose_slmm_departure_point.hpp"

#if 0
#ifdef COMPOSE_MIMIC_GPU
# pragma message "COMPOSE_MIMIC_GPU"
#endif
#ifdef COMPOSE_HORIZ_OPENMP
# pragma message "COMPOSE_HORIZ_OPENMP"
#endif
#ifdef COMPOSE_COLUMN_OPENMP
# pragma message "COMPOSE_COLUMN_OPENMP"
#endif
#ifdef COMPOSE_PORT
# pragma message "COMPOSE_PORT"
#endif
#ifdef COMPOSE_PORT_SEPARATE_VIEWS
# pragma message "COMPOSE_PORT_SEPARATE_VIEWS"
#endif
#ifdef COMPOSE_PACK_NOSCAN
# pragma message "COMPOSE_PACK_NOSCAN"
#endif
#endif

namespace slmm {
template <typename ES> SLMM_KIF
int get_nearest_point (const LocalMesh<ES>& m, Real* v, const Int my_ic) {
  nearest_point::calc(m, v);
  return get_src_cell(m, v, my_ic);
}
} // namespace slmm

namespace homme {
namespace islmpi {

template <typename MT>
ko::EnableIfNotOnGpu<MT> throw_on_sci_error (
  const slmm::LocalMesh<typename MT::DES>& mesh,
  const typename IslMpi<MT>::ElemDataD& ed,
  const bool nearest_point_permitted, const DepPoints<MT>& dep_points,
  Int k, Int lev, Int tci)
{
  const auto tgt_idx = mesh.tgt_elem;
  std::stringstream ss;
  ss.precision(17);
  const auto* v = &dep_points(tci,lev,k,0);
  ss << "Departure point is outside of halo:\n"
     << "  nearest point permitted: "
     << nearest_point_permitted
     << "\n  elem LID " << tci
     << " elem GID " << ed.me->gid
     << " (lev, k) (" << lev << ", " << k << ")"
     << " v " << v[0] << " " << v[1] << " " << v[2]
     << "\n  tgt_idx " << tgt_idx
     << " local mesh:\n  " << slmm::to_string(mesh) << "\n";
  slmm_throw_if(true, ss.str());
}

template <typename MT>
SLMM_KF ko::EnableIfOnGpu<MT> throw_on_sci_error (
  const slmm::LocalMesh<typename MT::DES>& mesh,
  const typename IslMpi<MT>::ElemDataD& ed,
  const bool nearest_point_permitted, const DepPoints<MT>& dep_points,
  Int k, Int lev, Int tci)
{
  ko::abort("throw_on_sci_error");
}

// Find where each departure point is.
template <typename MT>
void analyze_dep_points (IslMpi<MT>& cm, const Int& nets, const Int& nete,
                         const DepPoints<MT>& dep_points) {
  const auto& plane = cm.advecter->get_plane();
#ifdef COMPOSE_PORT
  const auto myrank = cm.p->rank();
  const Int np2 = cm.np2, nlev = cm.nlev;
  const auto& own_dep_mask = cm.own_dep_mask;
  cm.bla.zero();
  cm.nx_in_lid.zero();
  cm.nx_in_rank.zero();
  ko::deep_copy(own_dep_mask, 0);
  {
    const Int nearest_point_permitted_lev_bdy =
      cm.advecter->nearest_point_permitted_lev_bdy();
    const auto& local_meshes = cm.advecter->local_meshes();
    const auto& ed_d = cm.ed_d;
    const auto& nx_in_lid = cm.nx_in_lid;
    const auto& bla = cm.bla;
    const auto& nx_in_rank = cm.nx_in_rank;
    const auto f = COMPOSE_LAMBDA (const Int& ki) {
      const Int tci = nets + ki/(nlev*np2);
      const Int   k = (ki/nlev) % np2;
      const Int lev = ki % nlev;
      const auto& mesh = local_meshes(tci);
      const auto tgt_idx = mesh.tgt_elem;
      auto& ed = ed_d(tci);
      Int sci = slmm::get_src_cell(mesh, &dep_points(tci,lev,k,0), tgt_idx);
      if (sci == -1) {
        const bool npp = slmm::Advecter<MT>::nearest_point_permitted(
          nearest_point_permitted_lev_bdy, lev);
        if (npp) sci = slmm::get_nearest_point(mesh, &dep_points(tci,lev,k,0), tgt_idx);
        if (sci == -1) throw_on_sci_error<MT>(mesh, ed, npp, dep_points, k, lev, tci);
      }
      ed.src(lev,k) = sci;
      if (ed.nbrs(sci).rank == myrank)
        own_dep_mask(tci,lev,k) = 1;
      else {
        const auto ri = ed.nbrs(sci).rank_idx;
        const auto lidi = ed.nbrs(sci).lid_on_rank_idx;
        ko::atomic_increment(static_cast<volatile Int*>(&nx_in_lid(ri,lidi)));
        ko::atomic_increment(static_cast<volatile Int*>(&bla(ri,lidi,lev).xptr));
        ko::atomic_increment(static_cast<volatile Int*>(&nx_in_rank(ri)));
      }
    };
    ko::fence();
    ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, (nete - nets + 1)*nlev*np2), f);
  }
  {
    const auto& own_dep_list = cm.own_dep_list;
    const auto f = COMPOSE_LAMBDA (const Int ki, Int& slot, const bool fin) {
      const Int tci = nets + ki/(nlev*np2);
      const Int   k = (ki/nlev) % np2;
      const Int lev = ki % nlev;
      if (own_dep_mask(tci,lev,k)) {
        if (fin) {
          own_dep_list(slot,0) = tci;
          own_dep_list(slot,1) = lev;
          own_dep_list(slot,2) = k;
        }
        ++slot;
      }
    };
    Int slot;
    ko::fence();
    ko::parallel_scan(ko::RangePolicy<typename MT::DES>(0, (nete - nets + 1)*nlev*np2),
                      f, slot);
    cm.own_dep_list_len = slot;
  }
#else // COMPOSE_PORT
  const auto myrank = cm.p->rank();
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  ConstExceptGnu Int np2 = cm.np2, nlev = cm.nlev;
  cm.bla.zero();
  cm.nx_in_lid.zero();
  {
    const Int nearest_point_permitted_lev_bdy =
      cm.advecter->nearest_point_permitted_lev_bdy();
    const auto& local_meshes = cm.advecter->local_meshes();
    const auto& ed_d = cm.ed_d;
    const auto& bla = cm.bla;
    const auto& nx_in_lid = cm.nx_in_lid;
#ifdef COMPOSE_HORIZ_OPENMP
    const auto horiz_openmp = cm.horiz_openmp;
    const auto& ri_lidi_locks = cm.ri_lidi_locks;
#endif
    ko::parallel_for(ko::RangePolicy<typename MT::DES>(nets, nete+1),
                     COMPOSE_LAMBDA (const Int& tci) { ed_d(tci).own.clear(); });
    const auto f = COMPOSE_LAMBDA (const Int& ki) {
      const Int tci = nets + ki/(nlev*np2);
#if 0
      const Int   k = (ki/nlev) % np2;
      const Int lev = ki % nlev;
#else
      const Int   k = ki % np2;
      const Int lev = (ki/np2) % nlev;      
#endif
      const auto& mesh = local_meshes(tci);
      const auto tgt_idx = mesh.tgt_elem;
      auto& ed = ed_d(tci);
      Int sci = slmm::get_src_cell(mesh, &dep_points(tci,lev,k,0), tgt_idx);
      if (sci == -1) {
        const bool npp = slmm::Advecter<MT>::nearest_point_permitted(
          nearest_point_permitted_lev_bdy, lev);
        if (npp) sci = slmm::get_nearest_point(mesh, &dep_points(tci,lev,k,0), tgt_idx);
        if (sci == -1) throw_on_sci_error<MT>(mesh, ed, npp, dep_points, k, lev, tci);
      }
      ed.src(lev,k) = sci;
      if (ed.nbrs(sci).rank == myrank) {
        auto& t = ed.own.atomic_inc_and_return_next();
        t.lev = lev;
        t.k = k;
      } else {
        const auto ri = ed.nbrs(sci).rank_idx;
        const auto lidi = ed.nbrs(sci).lid_on_rank_idx;
#ifdef COMPOSE_HORIZ_OPENMP
        omp_lock_t* lock;
        if (horiz_openmp) {
          lock = &ri_lidi_locks(ri,lidi);
          omp_set_lock(lock);
        }
#endif
        ++nx_in_lid(ri,lidi);
        ++bla(ri,lidi,lev).xptr;
#ifdef COMPOSE_HORIZ_OPENMP
        if (horiz_openmp) omp_unset_lock(lock);
#endif
      }
    };
    ko::fence();
    ko::parallel_for(
      ko::RangePolicy<typename MT::DES>(0, (nete - nets + 1)*nlev*np2), f);
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto& nx_in_rank = cm.nx_in_rank(ri);
    nx_in_rank = 0;
    for (Int i = 0, n = cm.lid_on_rank_h(ri).n(); i < n; ++i)
      nx_in_rank += cm.nx_in_lid(ri,i);
  }
#endif
}

template void
analyze_dep_points(IslMpi<ko::MachineTraits>& cm, const Int& nets,
                   const Int& nete, const DepPoints<ko::MachineTraits>& dep_points);

} // namespace islmpi
} // namespace homme
