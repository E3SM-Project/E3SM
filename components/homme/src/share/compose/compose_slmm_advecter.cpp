#include "compose_slmm_advecter.hpp"

namespace slmm {

template <typename MT>
void Advecter<MT>::fill_nearest_points_if_needed () {
  if (geometry_ == Geometry::Type::plane) {
    // Shift local coords again after halo expansion. make_continuous is
    // idempotent, so there is no harm in calling it twice even when not needed.
    for (Int ie = 0; ie < local_mesh_h_.extent_int(0); ++ie)
      make_continuous(plane_, local_mesh_h_(ie));
  }
  if (nearest_point_permitted_lev_bdy_ >= 0)
    for (Int ie = 0; ie < local_mesh_h_.extent_int(0); ++ie)
      nearest_point::fill_perim(local_mesh_h_(ie));
}

template <typename MT>
void Advecter<MT>
::init_meta_data (const Int nelem_global, const Int* lid2facenum) {
  const auto nelemd = local_mesh_h_.extent_int(0);
  lid2facenum_ = Ints<DES>("Advecter::lid2facenum", nelemd);
  lid2facenum_h_ = ko::create_mirror_view(lid2facenum_);
  std::copy(lid2facenum, lid2facenum + nelemd, lid2facenum_h_.data());
  ko::deep_copy(lid2facenum_, lid2facenum_h_);
  s2r_.init(geometry_, cubed_sphere_map_, nelem_global, lid2facenum_, plane_);
}

template <typename MT>
void Advecter<MT>::check_ref2sphere (const Int ie, const Real* p_homme) {
  const auto& m = local_mesh_host(ie);
  Real ref_coord[2];
  siqk::sqr::Info info;
  SphereToRef<typename MT::HES> s2r;
  s2r.init(geometry_, cubed_sphere_map_, s2r_.nelem_global(), lid2facenum_h_,
           plane_);
  const Real tol = s2r_.tol();
  s2r.calc_sphere_to_ref(ie, m, p_homme, ref_coord[0], ref_coord[1], &info);
  const slmm::Basis basis(4, 0);
  const slmm::GLL gll;
  const Real* x, * wt;
  gll.get_coef(basis, x, wt);
  int fnd[2] = {0};
  Real min[2] = {1,1};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 2; ++j) {
      const Real d = std::abs(ref_coord[j] - x[i]);
      min[j] = std::min(min[j], d);
      if (d < tol)
        fnd[j] = 1;
    }
  if ( ! fnd[0] || ! fnd[1])
    printf("COMPOSE check_ref2sphere: %1.15e %1.15e (%1.2e %1.2e) %d %d\n",
           ref_coord[0], ref_coord[1], min[0], min[1],
           info.success, info.n_iterations);
  if ( ! s2r.check(ie, m))
    printf("COMPOSE SphereToRef::check return false: ie = %d\n", ie);
}

template <typename MT>
ko::EnableIfDiffSpace<MT>
deep_copy (typename Advecter<MT>::LocalMeshesD& d,
           typename Advecter<MT>::LocalMeshesD::HostMirror& m,
           const typename Advecter<MT>::LocalMeshesH& s) {
  const Int nlm = s.extent_int(0);
  d = typename Advecter<MT>::LocalMeshesD("LocalMeshes", nlm);
  m = ko::create_mirror_view(d);
  for (Int i = 0; i < nlm; ++i)
    deep_copy(m(i), s(i));
  ko::deep_copy(d, m);
}

template <typename MT>
ko::EnableIfSameSpace<MT>
deep_copy (typename Advecter<MT>::LocalMeshesD& d,
           typename Advecter<MT>::LocalMeshesD::HostMirror& m,
           const typename Advecter<MT>::LocalMeshesH& s) {
  d = s;
}

template <typename MT>
void Advecter<MT>::sync_to_device() {
  deep_copy<MT>(local_mesh_d_, local_mesh_m_, local_mesh_h_);
}

template class Advecter<>;

} // namespace slmm
