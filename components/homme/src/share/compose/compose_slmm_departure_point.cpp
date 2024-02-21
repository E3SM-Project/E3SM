#include "compose_slmm_departure_point.hpp"

namespace slmm {

void fill_normals (LocalMesh<ko::MachineTraits::HES>& m) {
  if (m.geometry == Geometry::Type::plane) {
#ifndef NDEBUG
    {
      // In principle the plane could be embedded in a general orientation in 3D
      // space. In practice, it's sensibly implemented in Homme as || to the x-y
      // plane. Be sure of that here.
      const Int n = nslices(m.p);
      for (Int i = 0; i < n; ++i) assert(m.p(i,2) == 0);
    }
#endif
    // Following the assumption above, be sure the third component is 0 since
    // PlaneGeometry::fill_normals won't write to it.
    const Int n = nslices(m.nml);
    for (Int i = 0; i < n; ++i) m.nml(i,2) = 0;
    fill_normals<siqk::PlaneGeometry>(m);
  } else
    fill_normals<siqk::SphereGeometry>(m);
}

// Modify m.p so that the elements surrounding the target element have
// continuous rather than periodic coordinate values.
void make_continuous (const Plane& p, LocalMesh<ko::MachineTraits::HES>& m) {
  slmm_assert(m.tgt_elem >= 0 && m.tgt_elem < nslices(m.e));
  // The translations below assume each cell has its own four vertex point data.
  slmm_assert(nslices(m.p) == 4*nslices(m.e));
  // Geometric center of cell.
  const auto getctr = [&] (const Int ie, Real ctr[3]) {
    const auto cell = slice(m.e, ie);
    for (Int d = 0; d < 3; ++d) ctr[d] = 0;
    for (Int i = 0; i < 4; ++i)
      for (Int d = 0; d < 3; ++d)
        ctr[d] += m.p(cell[i],d);
    for (Int d = 0; d < 3; ++d) ctr[d] /= 4;
  };
  Real tctr[3];
  getctr(m.tgt_elem, tctr);
  // Distance using cell center.
  const auto dist = [&] (const Int ie, const Real dx, const Real dy) {
    Real ctr[3];
    getctr(ie, ctr);
    Real dist = 0;
    Real trans[3] = {dx, dy, 0};
    for (Int d = 0; d < 3; ++d) dist += siqk::square(ctr[d] + trans[d] - tctr[d]);
    return std::sqrt(dist);
  };
  // Translate cell dim d by delta.
  const auto translate = [&] (const Int ie, const Int d, const Real delta) {
    const auto cell = slice(m.e, ie);
    for (Int i = 0; i < 4; ++i)
      m.p(cell[i],d) += delta;
  };
  const Int ncell = nslices(m.e);
  for (Int ie = 0; ie < ncell; ++ie) {
    if (ie == m.tgt_elem) continue;
    // Find the distance-minimizing translation, including 0 translation.
    Real min_dist = 1e20, mdx = 0, mdy = 0;
    for (const Real dx : {-p.Lx, 0.0, p.Lx})
      for (const Real dy : {-p.Ly, 0.0, p.Ly}) {
        const Real d = dist(ie, dx, dy);
        if (d < min_dist) { min_dist = d; mdx = dx; mdy = dy; }
      }
    // Apply the translation.
    translate(ie, 0, mdx);
    translate(ie, 1, mdy);
  }
}

namespace nearest_point {

Int test_canpoa (const bool sphere) {
  Int nerr = 0;
  using geo = siqk::SphereGeometry;
  Real p0[] = {-0.5,0.5,0}, p1[] = {0.1,0.8,0}, nml[] = {0,0,-1};
  geo::normalize(p0);
  geo::normalize(p1);
  const Real points[][3] = {{-100,1,7}, {0,1,-7}, {11,-1,7}};
  Real correct[3][3];
  geo::copy(correct[0], p0);
  correct[1][0] = 0; correct[1][1] = 1; correct[1][2] = 0;
  geo::copy(correct[2], p1);
  for (Int i = 0; i < 3; ++i) {
    Real nr[3];
    nearest_point::calc_approx_nearest_point_on_arc(p0, p1, nml, points[i], nr, sphere);
    Real d[3];
    geo::axpbyz(1, nr, -1, correct[i], d);
    if (std::sqrt(geo::norm2(d)) > 10*std::numeric_limits<Real>::epsilon())
      ++nerr;
  }
  return nerr;
}

template <typename ES>
Int test_calc (const LocalMesh<ES>& m, const Int& tgt_ic, const Real length_scale) {
  Int nerr = 0;
  // If a point is on the perim, then calc should return it as the nearest
  // point.
  const Int np_perim = len(m.perimp);
  for (Int k = 0; k < np_perim; ++k) {
    const auto p0 = slice(m.p, m.perimp(k));
    const auto p1 = slice(m.p, m.perimp((k+1) % np_perim));
    const auto L = 1e2*length_scale;
    const auto tol = std::numeric_limits<Real>::epsilon()*L;
    Real v[3];
    for (Int d = 0; d < 3; ++d) v[d] = p0[d];
    calc(m, v);
    if (calc_dist(p0, v) > tol) ++nerr;
    Real on[3];
    siqk::SphereGeometry::combine(p0, p1, 0.4, on);
    if (m.is_sphere()) siqk::SphereGeometry::normalize(on);
    for (Int d = 0; d < 3; ++d) v[d] = on[d];
    calc(m, v);
    if (calc_dist(on, v) > tol) ++nerr;
  }
  return nerr;
}

template <typename geo, typename ES>
Int test_fill_perim (const LocalMesh<ES>& m, const Int& tgt_ic) {
  Int nerr = 0;
  const auto tgt_cell = slice(m.e, tgt_ic);
  const Int np_perim = len(m.perimp);
  { // Check that no target-cell point is on the perim.
    Int on = 0;
    for (Int ip = 0; ip < 4; ++ip) {
      const auto p = slice(m.p, tgt_cell[ip]);
      for (Int k = 0; k < np_perim; ++k) {
        const auto perim = slice(m.p, m.perimp(k));
        Real diff = 0;
        for (Int i = 0; i < 3; ++i) diff += std::abs(p[i] - perim[i]);
        if (diff == 0) ++on;
      }
    }
    if (on > 0) ++nerr;
  }
  { // Check that adjacent perimeter points agree with the edge's normal.
    for (Int k = 0; k < np_perim; ++k) {
      const auto p0 = slice(m.p, m.perimp(k));
      const auto p1 = slice(m.p, m.perimp((k+1) % np_perim));
      const auto nml = slice(m.nml, m.perimnml(k));
      const auto dot = geo::dot_c_amb(nml, p1, p0);
      const auto L = calc_dist(p0, p1);
      // Normal orthogonal to edge?
      if (std::abs(dot) > 1e2*std::numeric_limits<Real>::epsilon()/L) ++nerr;
      // Oriented inwards?
      Real outward[3] = {0};
      geo::edge_normal(p1, p0, outward);
      if (geo::dot(nml, outward) >= 0) ++nerr;
    }
  }
  return nerr;
}
} // namespace nearest_point

Int unittest (LocalMesh<ko::MachineTraits::HES>& m, const Int tgt_elem,
              const Real length_scale) {
  Int nerr = 0, ne = 0;
  const Int nc = len(m.e);
  for (Int ic = 0; ic < nc; ++ic) {
    const auto cell = slice(m.e, ic);
    static const Real alphas[] = { 0.01, 0.99, 0, 1 };
    static const int nalphas = sizeof(alphas)/sizeof(*alphas);
    for (Int i = 0; i < nalphas; ++i)
      for (Int j = 0; j < nalphas; ++j) {
        const Real a = alphas[i], oma = 1-a, b = alphas[j], omb = 1-b;
        Real v[3] = {0};
        for (Int d = 0; d < 3; ++d)
          v[d] = (  b*(a*m.p(cell[0], d) + oma*m.p(cell[1], d)) + 
                  omb*(a*m.p(cell[3], d) + oma*m.p(cell[2], d)));
        if (i < 2 && j < 2) {
          if (get_src_cell(m, v, tgt_elem) != ic) ++ne;
        } else {
          if (get_src_cell(m, v, tgt_elem) == -1) ++ne;
        }
      }
  }
  if (ne) pr("slmm::unittest: get_src_cell failed");
  nerr += ne;
  ne = nearest_point::test_canpoa(true);
  if (ne) pr("slmm::unittest: test_canpoa sphere failed");
  nerr += ne;
  ne = nearest_point::test_canpoa(false);
  if (ne) pr("slmm::unittest: test_canpoa plane failed");
  nerr += ne;
  {
    nearest_point::fill_perim(m);
    if (m.is_sphere())
      ne = nearest_point::test_fill_perim<siqk::SphereGeometry>(m, tgt_elem);
    else
      ne = nearest_point::test_fill_perim<siqk::PlaneGeometry>(m, tgt_elem);
    if (ne) pr("slmm::unittest: test_fill_perim failed");
    nerr += ne;
    ne = nearest_point::test_calc(m, tgt_elem, length_scale);
    if (ne) pr("slmm::unittest: test_calc failed");
    nerr += ne;
  }
  {
    ne = siqk::sqr::test::test_sphere_to_ref(m.p, m.e, m.is_sphere());
    if (ne) pr("slmm::unittest: siqk::sqr::test::test_sphere_to_ref failed");
    nerr += ne;
  }
  return nerr;
}

std::string to_string (const LocalMesh<ko::MachineTraits::HES>& m) {
  std::stringstream ss;
  ss.precision(17);
  ss << "(mesh nnode " << nslices(m.p) << " nelem " << nslices(m.e);
  for (Int ie = 0; ie < nslices(m.e); ++ie) {
    ss << "\n  (elem " << ie;
    ss << " (";
    for (Int iv = 0; iv < szslice(m.e); ++iv) ss << " " << m.e(ie,iv);
    ss << ")";
    for (Int iv = 0; iv < szslice(m.e); ++iv) {
      ss << "\n     (p";
      for (Int d = 0; d < 3; ++d)
        ss << " " << m.p(m.e(ie,iv),d);
      ss << ")";
      ss << "\n     (nml";
      for (Int d = 0; d < 3; ++d)
        ss << " " << m.nml(m.en(ie,iv),d);
      ss << ")";
    }
    ss << "))";
  }
  ss << ")";
  return ss.str();
}

} // namespace slmm
