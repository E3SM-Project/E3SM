#ifndef INCLUDE_COMPOSE_SLMM_ADVECTER_HPP
#define INCLUDE_COMPOSE_SLMM_ADVECTER_HPP

#include "compose_slmm.hpp"
#include "compose_slmm_siqk.hpp"
#include "compose_slmm_departure_point.hpp"

#include <limits>
#include <memory>

namespace slmm {

template <typename ES>
using Ints = Kokkos::View<Int*, ES>;

// Wrap call to siqk::sqr::calc_sphere_to_ref. That impl supports only
// cubed_sphere_map=2, and we want to keep it that way. This wrapper supports,
// in addition, cubed_sphere_map=0.
template <typename ES>
struct SphereToRef {
  void init (const Geometry::Type geometry, const Int cubed_sphere_map,
             const Int nelem_global, const typename Ints<ES>::const_type& lid2facenum,
             const Plane& plane) {
    geometry_ = geometry;
    cubed_sphere_map_ = cubed_sphere_map;
    lid2facenum_ = lid2facenum;
    nelem_global_ = nelem_global;
    ne_ = static_cast<Int>(std::round(std::sqrt((nelem_global / 6))));
    slmm_throw_if( ! (cubed_sphere_map_ != 0 || 6*ne_*ne_ == nelem_global),
                  "If cubed_sphere_map = 0, then the mesh must be a "
                  "regular cubed-sphere.");
    plane_ = plane;
    calc_length_scale();
  }

  void set_plane (const Plane& plane) {
    plane_ = plane;
    calc_length_scale();
  }

  Int nelem_global () const { return nelem_global_; }

  Real tol () const {
    return 1e3 * ne_ * std::numeric_limits<Real>::epsilon();
  }

  // See siqk::sqr::calc_sphere_to_ref for docs.
  SLMM_KF void calc_sphere_to_ref (
    const Int& ie, const LocalMesh<ES>& m,
    const Real q[3],
    Real& a, Real& b,
    siqk::sqr::Info* const info = nullptr,
    const Int max_its = 10,
    const Real tol = 1e2*std::numeric_limits<Real>::epsilon()) const
  {
    if (cubed_sphere_map_ == 2) {
      // Tolerance is absolute, so multiply the rel tol by the length scale.
      if (geometry_ == Geometry::Type::sphere)
        siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, m.tgt_elem), q, a, b,
                                      info, max_its, tol*length_scale_);
      else {
        Real q_per[] = {q[0], q[1], q[2]};
        continuous2periodic(plane_, m, q_per);
        siqk::sqr::calc_plane_to_ref(m.p, slice(m.e, m.tgt_elem), q_per, a, b,
                                     info, max_its, tol*length_scale_);
      }
    } else {
      const Int face = lid2facenum_(ie); //assume: ie corresponds to m.tgt_elem.
      map_sphere_coord_to_face_coord(face-1, q[0], q[1], q[2], a, b);
      a = map_face_coord_to_cell_ref_coord(a);
      b = map_face_coord_to_cell_ref_coord(b);
      if (info) { info->success = true; info->n_iterations = 1; }
    }
  }

  bool check (const Int& ie, const LocalMesh<ES>& m) const {
    if (cubed_sphere_map_ != 0) return true;
    const Int face = lid2facenum_(ie); //assume: ie corresponds to m.tgt_elem.
    Real cent[3] = {0};
    const auto cell = slice(m.e, m.tgt_elem);
    for (int i = 0; i < 4; ++i)
      for (int d = 0; d < 3; ++d)
        cent[d] += 0.25*m.p(cell[i], d);
    const Int cf = get_cube_face_idx(cent[0], cent[1], cent[2]) + 1;
    return face == cf;
  }

private:
  Geometry::Type geometry_;
  Int ne_, nelem_global_, cubed_sphere_map_;
  typename Ints<ES>::const_type lid2facenum_;
  Plane plane_;
  Real length_scale_;

  void calc_length_scale () {
    length_scale_ = (geometry_ == Geometry::Type::sphere ?
                     1 :
                     std::max(plane_.Lx, plane_.Ly));
  }

  // Follow the description given in
  //     coordinate_systems_mod::unit_face_based_cube_to_unit_sphere.
  SLMM_KF static Int get_cube_face_idx (const Real& x, const Real& y, const Real& z) {
    const Real ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
    if (ax >= ay) {
      if (ax >= az) return x > 0 ? 0 : 2;
      else return z > 0 ? 5 : 4;
    } else {
      if (ay >= az) return y > 0 ? 1 : 3;
      else return z > 0 ? 5 : 4;
    }
  }

  SLMM_KF static void map_sphere_coord_to_face_coord (
    const Int& face_idx, const Real& x, const Real& y, const Real& z,
    Real& fx, Real& fy)
  {
    static constexpr Real theta_max = 0.25*M_PI;
    Real d;
    switch (face_idx) {
    case  0: d = std::abs(x); fx =  y/d; fy =  z/d; break;
    case  1: d = std::abs(y); fx = -x/d; fy =  z/d; break;
    case  2: d = std::abs(x); fx = -y/d; fy =  z/d; break;
    case  3: d = std::abs(y); fx =  x/d; fy =  z/d; break;
    case  4: d = std::abs(z); fx =  y/d; fy =  x/d; break;
    default: d = std::abs(z); fx =  y/d; fy = -x/d;
    }
    fx = std::atan(fx) / theta_max;
    fy = std::atan(fy) / theta_max;
  }

  SLMM_KF Real map_face_coord_to_cell_ref_coord (Real a) const {
    a = (0.5*(1 + a))*ne_;
    a = 2*(a - std::floor(a)) - 1;
    return a;
  }

  // See comments before 'step' in compose_slmm_islmpi.hpp re: continuous and
  // periodic coordinate values.
  //   Move x periodically, as needed, so that it's closest to the identified
  // source cell in its periodic coordinates.
  SLMM_KF static void
  continuous2periodic (const Plane& p, const LocalMesh<ES>& m, Real* const x) {
    const auto cell = siqk::slice(m.e, m.tgt_elem);
    const auto sx = (m.p(cell[0],0) + m.p(cell[1],0) +
                     m.p(cell[2],0) + m.p(cell[3],0))/4;
    const auto sy = (m.p(cell[0],1) + m.p(cell[1],1) +
                     m.p(cell[2],1) + m.p(cell[3],1))/4;
    Real dist = 1e20, mdx = 0;
    for (const Real dx : {-p.Lx, 0.0, p.Lx}) {
      const Real d = std::abs(x[0] + dx - sx);
      if (d < dist) { dist = d; mdx = dx; }
    }
    Real mdy = 0;
    dist = 1e20;
    for (const Real dy : {-p.Ly, 0.0, p.Ly}) {
      const Real d = std::abs(x[1] + dy - sy);
      if (d < dist) { dist = d; mdy = dy; }
    }
    x[0] += mdx;
    x[1] += mdy;
  }
};

// Advecter has purely mesh-local knowledge, with once exception noted below.
template <typename MT = ko::MachineTraits>
struct Advecter {
  typedef typename MT::HES HES;
  typedef typename MT::DES DES;

  template <typename ES> using LocalMeshes = ko::View<LocalMesh<ES>*, ES>;
  typedef LocalMeshes<HES> LocalMeshesH;
  typedef LocalMeshes<DES> LocalMeshesD;

  typedef std::shared_ptr<Advecter> Ptr;
  typedef std::shared_ptr<const Advecter> ConstPtr;

  struct Alg {
    enum Enum {
      jct,             // Cell-integrated Jacobian-combined transport.
      qos,             // Cell-integrated quadrature-on-sphere transport.
      csl_gll,         // Intended to mimic original Fortran CSL.
      csl_gll_subgrid, // Stable np=4 subgrid bilinear interp.
      csl_gll_exp,     // Stabilized np=4 subgrid reconstruction.
    };
    static Enum convert (Int alg) {
      switch (alg) {
      case 2: case 29:  return jct;
      case 3: case 39:  return qos;
      case 10: case 18: return csl_gll;
      case 11: case 17: return csl_gll_subgrid;
      case 12: case 19: return csl_gll_exp;
      default: slmm_throw_if(true, "transport_alg " << alg << " is invalid.");
      }
    }
    static bool is_cisl (const Enum& alg) { return alg == jct || alg == qos; }
  };

  Advecter (const Int np, const Int nelem, const Int transport_alg,
            const Int cubed_sphere_map, const Geometry::Type geometry,
            const Int nearest_point_permitted_lev_bdy)
    : alg_(Alg::convert(transport_alg)),
      np_(np), np2_(np*np), np4_(np2_*np2_),
      cubed_sphere_map_(geometry == Geometry::Type::plane ? 2 : cubed_sphere_map),
      geometry_(geometry),
      tq_order_(alg_ == Alg::qos ? 14 : 12),
      nearest_point_permitted_lev_bdy_(nearest_point_permitted_lev_bdy)
  {
    slmm_throw_if(cubed_sphere_map == 0 && Alg::is_cisl(alg_),
                  "When cubed_sphere_map = 0, SLMM supports only ISL methods.");
    local_mesh_h_ = LocalMeshesH("local_mesh_h_", nelem);
  }

  void init_plane (Real Sx, Real Sy, Real Lx, Real Ly) {
    slmm_assert(geometry_ == Geometry::Type::plane);
    plane_.Sx = Sx; plane_.Sy = Sy;
    plane_.Lx = Lx; plane_.Ly = Ly;
    s2r_.set_plane(plane_);
  }

  void fill_nearest_points_if_needed();

  // After Advecter is fully initialized, send all the data to device.
  void sync_to_device();

  Int np  () const { return np_ ; }
  Int np2 () const { return np2_; }
  Int np4 () const { return np4_; }
  Int nelem () const { return local_mesh_h_.extent_int(0); }
  Int tq_order () const { return tq_order_; }
  typename Alg::Enum alg () const { return alg_; }
  bool is_cisl () const { return Alg::is_cisl(alg_); }

  Int cubed_sphere_map () const { return cubed_sphere_map_; }
  Geometry::Type geometry () const { return geometry_; }
  bool is_sphere () const { return geometry_ == Geometry::Type::sphere; }
  const Plane& get_plane () const { return plane_; }
  const Ints<DES>& lid2facenum () const { return lid2facenum_; }

  // nelem_global is used only if cubed_sphere_map = 0, to deduce ne in
  // nelem_global = 6 ne^2. That is b/c cubed_sphere_map = 0 is supported in
  // Homme only for regular meshes (not RRM), and ne is essential to using the
  // efficiency it provides.
  void init_meta_data(const Int nelem_global, const Int* lid2facenum);

  template <typename Array3D>
  void init_local_mesh_if_needed(const Int ie, const Array3D& corners,
                                 const Real* p_inside);

  // Check that our ref <-> sphere map agrees with Homme's. p_homme is a GLL
  // point on the sphere. Check that we map it to a GLL ref point.
  void check_ref2sphere(const Int ie, const Real* p_homme);

  const LocalMesh<HES>& local_mesh_host (const Int ie) const {
    slmm_assert(ie < static_cast<Int>(local_mesh_h_.size()));
    return local_mesh_h_(ie);
  }
  LocalMesh<HES>& local_mesh_host (const Int ie) {
    slmm_assert(ie < static_cast<Int>(local_mesh_h_.size()));
    return local_mesh_h_(ie);
  }

  const LocalMesh<DES>& local_mesh (const Int ie) const {
    slmm_assert(ie < static_cast<Int>(local_mesh_d_.size()));
    return local_mesh_d_(ie);
  }
  const LocalMeshesD& local_meshes () const {
    return local_mesh_d_;
  }

  SLMM_KIF static bool nearest_point_permitted (
    const Int& nearest_point_permitted_lev_bdy, const Int& lev)
  {
    return lev <= nearest_point_permitted_lev_bdy;
  }
  bool nearest_point_permitted (const Int& lev) const {
    return nearest_point_permitted(nearest_point_permitted_lev_bdy_, lev);
  }
  Int nearest_point_permitted_lev_bdy () const {
    return nearest_point_permitted_lev_bdy_;
  }

  const SphereToRef<DES>& s2r () const { return s2r_; }

private:
  const typename Alg::Enum alg_;
  const Int np_, np2_, np4_, cubed_sphere_map_;
  Geometry::Type geometry_;
  LocalMeshesH local_mesh_h_;
  LocalMeshesD local_mesh_d_;
  typename LocalMeshesD::HostMirror local_mesh_m_; // handle managed allocs
  // For CISL:
  const Int tq_order_;
  // For recovery from get_src_cell failure:
  Int nearest_point_permitted_lev_bdy_;
  // Meta data obtained at initialization that can be used later.
  Ints<HES> lid2facenum_h_;
  Ints<DES> lid2facenum_;
  SphereToRef<DES> s2r_;
  Plane plane_;
};

template <typename MT>
template <typename Array3D>
void Advecter<MT>
::init_local_mesh_if_needed (const Int ie, const Array3D& corners,
                             const Real* p_inside) {
  slmm_assert(ie < static_cast<Int>(local_mesh_h_.size()));
  if (local_mesh_h_(ie).p.extent_int(0) != 0) return;
  const Int
    nd = 3,
    nvert = corners.extent_int(1),
    ncell = corners.extent_int(2),
    N = nvert*ncell;
  auto& m = local_mesh_h_(ie);
  m.geometry = geometry_;
  typedef LocalMesh<HES> LM;
  m.p = typename LM::RealArray("p", N);
  m.e = typename LM::IntArray("e", ncell, nvert);
  for (Int ci = 0, k = 0; ci < ncell; ++ci)
    for (Int vi = 0; vi < nvert; ++vi, ++k) {
      for (int j = 0; j < nd; ++j)
        m.p(k,j) = corners(j,vi,ci);
      m.e(ci,vi) = k;
    }
  fill_normals(m);
  m.tgt_elem = slmm::get_src_cell(m, p_inside);
  slmm_assert(m.tgt_elem >= 0 && m.tgt_elem < ncell);
  if (geometry_ == Geometry::Type::plane) {
    // Shift local coords here for checks.
    make_continuous(plane_, m);
  }
}

} // namespace slmm

#endif
