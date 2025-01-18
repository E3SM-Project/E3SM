#ifndef INCLUDE_COMPOSE_SLMM_DEPARTURE_POINT_HPP
#define INCLUDE_COMPOSE_SLMM_DEPARTURE_POINT_HPP

#include "compose.hpp"
#include "compose_slmm.hpp"

namespace slmm {

// Plane dimensions for planar doubly-periodic domain.
struct Plane { Real Sx, Sy, Lx, Ly; };

/* A local mesh is described by the following arrays:
       p: 3 x #nodes, the array of vertices.
       e: max(#verts) x #elems, the array of element base-0 indices.
       nml: 3 x #edges, the array of edge normals.
       en: max(#verts) x #elems, the array of edge-normal base-0 indices.
     e. e indexes p. e(i,j) == -1 in column j indicates that j:end are not used.
     nml. As a mesh is refined, cancellation error makes an edge normal based
   off of an element's vertices increasingly inaccurate. Roughly, if an edge
   subtends angle phi of the sphere, -log10(phi/(2 pi)) digits are lost in the
   edge normal. Therefore, we compute edge normals offline, since in certain
   meshes, they can be computed by an accurate means. E.g., in a cubed-sphere
   mesh, the whole line of a square face can be used to compute the edge
   normal. Furthermore, there are far fewer unique edge normals than edges.
 */
template <typename ES>
struct LocalMesh {
  using RealArray = ko::View<Real*[3], siqk::Layout, ES>;
  using IntArray = ko::View<Int**, siqk::Layout, ES>;
  using Ints = ko::View<Int*, ES>;

  // Local mesh data.
  Geometry::Type geometry;
  RealArray p, nml;
  IntArray e, en;

  // If the nearest-point procedure is active, these data are allocated.
  // mesh.p(perimp(k),:) is the k'th vertex in the mesh's perimeter. Similarly,
  // mesh.nml(perimnml(k),:) is the k'th edge's normal.
  Ints perimp, perimnml;

  // Index of the target element in this local mesh.
  Int tgt_elem;

  SLMM_KIF bool is_sphere () const { return geometry == Geometry::Type::sphere; }
};

// Inward-oriented normal. In practice, we want to form high-quality normals
// using information about the cubed-sphere mesh. This is a low-quality
// brute-force calculation.
template <typename geo>
void fill_normals (LocalMesh<ko::MachineTraits::HES>& m) {
  // Count number of edges.
  Int ne = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1) break; else ++ne;
  // Fill.
  siqk::Idxs::HostMirror en("en", nslices(m.e), szslice(m.e));
  ko::deep_copy(en, -1);
  siqk::Vec3s::HostMirror nml("nml", ne);
  Int ie = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1)
        break;
      else {
        // Somewhat complicated next node index.
        const Int iv_next = (iv+1 == szslice(m.e) ? 0 :
                             (m.e(ip,iv+1) == -1 ? 0 : iv+1));
        geo::edge_normal(slice(m.p, m.e(ip, iv)), slice(m.p, m.e(ip, iv_next)),
                         slice(nml, ie));
        en(ip,iv) = ie;
        ++ie;
      }
  m.en = en;
  m.nml = nml;
}

void fill_normals(LocalMesh<ko::MachineTraits::HES>& m);

// For doubly-periodic planar mesh, make the local mesh patch's vertices
// continuous in space, anchored at m.tgt_elem.
void make_continuous(const Plane& p, LocalMesh<ko::MachineTraits::HES>& m);

template <typename ESD, typename ESS>
void deep_copy (LocalMesh<ESD>& d, const LocalMesh<ESS>& s) {
  siqk::resize_and_copy(d.p, s.p);
  siqk::resize_and_copy(d.nml, s.nml);
  siqk::resize_and_copy(d.e, s.e);
  siqk::resize_and_copy(d.en, s.en);
  d.tgt_elem = s.tgt_elem;
  siqk::resize_and_copy(d.perimp, s.perimp);
  siqk::resize_and_copy(d.perimnml, s.perimnml);
}

// Is v inside, including on, the quad ic?
template <typename ES> SLMM_KIF
bool is_inside (const LocalMesh<ES>& m,
                const Real* v, const Real& atol, const Int& ic) {
  using slmm::slice;
  const auto cell = slice(m.e, ic);
  const auto celln = slice(m.en, ic);
  const Int ne = 4;
  assert(szslice(m.e) == ne);
  assert(szslice(m.en) == ne);
  bool inside = true;
  for (Int ie = 0; ie < ne; ++ie) {
    // We can use SphereGeometry here regardless of plane or sphere because the
    // third component is 0; see comments in fill_normals.
    if (siqk::SphereGeometry::dot_c_amb(slice(m.nml, celln[ie]),
                                        v, slice(m.p, cell[ie]))
        < -atol) {
      inside = false;
      break;
    }
  }
  return inside;
}

// Both cubed_sphere_map=0 and cubed_sphere_map=2 can use this method.
// (cubed_sphere_map=1 is not impl'ed in Homme.)
//   This method is natural for cubed_sphere_map=2, so RRM works automatically.
//   cubed_sphere_map=0 is supported in Homme only for regular cubed-sphere
// meshes, not RRM; in that case, edges are great arcs, so again this impl
// works.
//   In contrast, calc_sphere_to_ref has to be specialized on cubed_sphere_map.
//   In the case of planar geometry, this method again works because the only
// difference is in the edge normals. The SphereGeometry linear algebra below
// works because the third component is 0; see the comments in fill_normals.
template <typename ES> SLMM_KF
int get_src_cell (const LocalMesh<ES>& m, // Local mesh.
                  const Real* v, // 3D Cartesian point.
                  const Int my_ic = -1) { // Target cell in the local mesh.
  using slmm::len;
  const Int nc = len(m.e);
  Real atol = 0;
  for (Int trial = 0; trial < 3; ++trial) {
    if (trial > 0) {
      if (trial == 1) {
        using slmm::slice;
        // If !inside in the first sweep, pad each cell. Recall we're operating
        // on the unit sphere, so we don't have to worry about a radius in the
        // following.
        //   Get a representative edge length.
        const int ic = my_ic == -1 ? 0 : my_ic;
        const auto cell = slice(m.e, ic);
        Real d[3];
        siqk::SphereGeometry::axpbyz( 1, slice(m.p, cell[1]),
                                     -1, slice(m.p, cell[0]),
                                     d);
        const Real L = std::sqrt(siqk::SphereGeometry::norm2(d));
        // We can expect to lose approx. -log10(L) digits due to cancellation in
        // the formation of the edge normal and in dot_c_amb. Multiply by 100
        // for a little extra padding.
        atol = 1e2 * ko::NumericTraits<Real>::epsilon() / L;
      } else {
        // Ok, we really didn't do that very well. We're still failing to find
        // the element. Ramp up the atol even more.
        atol = ko::max(1e3*atol,
                       std::sqrt(ko::NumericTraits<Real>::epsilon()));
      }
    }
    if (my_ic != -1 && is_inside(m, v, atol, my_ic)) return my_ic;
    for (Int ic = 0; ic < nc; ++ic) {
      if (ic == my_ic) continue;
      if (is_inside(m, v, atol, ic)) return ic;
    }
  }
  return -1;
}

namespace nearest_point {
/* Get external segments in preproc step.
   Get approximate nearest point in each segment.
     Project onto plane of segment and normalize.
     If outside of arc, choose nearest endpoint.
     Approx distance as cartesian distance.
   Of the approx dists, take the point associated with smallest.
 */

template <typename V3a, typename V3b>
inline Real calc_dist (const V3a& p0, const V3b& p1) {
  Real len = 0;
  for (Int d = 0; d < 3; ++d) len += siqk::square(p0[d] - p1[d]);
  return std::sqrt(len);    
}

template <typename ES>
inline Real calc_dist (const LocalMesh<ES>& m, const Int& i0, const Int& i1) {
  return calc_dist(slice(m.p, i0), slice(m.p, i1));
}

template <typename ES>
void find_external_edges (const LocalMesh<ES>& m, std::vector<Int>& external_edges) {
  const Int ncell = nslices(m.e);
  // Get the minimum edge length to normalize distance.
  Real min_edge_len = 10;
  for (Int ic = 0; ic < ncell; ++ic) {
    const auto cell = slice(m.e, ic);
    for (Int ie = 0; ie < 4; ++ie)
      min_edge_len = std::min(min_edge_len,
                              calc_dist(m, cell[ie], cell[(ie+1)%4]));
  }
  // If two edges are sufficiently close, they are interpreted as the same and
  // so internal. That leaves just the external edges.
  for (Int ic0 = 0; ic0 < ncell; ++ic0) {
    const auto cell0 = slice(m.e, ic0);
    for (Int ie0 = 0; ie0 < 4; ++ie0) {
      bool fnd_match = false;
      for (Int ic1 = 0; ic1 < ncell; ++ic1) {
        if (ic1 == ic0) continue;
        const auto cell1 = slice(m.e, ic1);
        for (Int ie1 = 0; ie1 < 4; ++ie1) {
          const auto edist = 
            std::min(calc_dist(m, cell0[ie0], cell1[ie1]) +
                     calc_dist(m, cell0[(ie0+1)%4], cell1[(ie1+1)%4]),
                     calc_dist(m, cell0[ie0], cell1[(ie1+1)%4]) +
                     calc_dist(m, cell0[(ie0+1)%4], cell1[ie1]));
          if (edist < 0.01*min_edge_len) {
            fnd_match = true;
            break;
          }
        }
        if (fnd_match) break;
      }
      if ( ! fnd_match) external_edges.push_back(4*ic0 + ie0);
    }
  }
  const Int nee = external_edges.size();
  slmm_throw_if((nslices(m.e) == 9 && nee != 12) ||
                (nslices(m.e) == 8 && nee != 10),
                "The local mesh has " << nslices(m.e)
                << " edges but fill_perim found " << nee
                << " external edges. Mesh:\n"
                << to_string(m));
  // Order edges so the leading point of edge i is the same as the trailing
  // point of edge (i+1) % nee.
  for (Int i = 0; i < nee-2; ++i) {
    const Int ic = external_edges[i] / 4, ie = external_edges[i] % 4;
    const auto icell = slice(m.e, ic);
    Real min_dist = 10;
    Int min_j = -1;
    for (Int j = i+1; j < nee; ++j) {
      const Int jc = external_edges[j] / 4, je = external_edges[j] % 4;
      const auto jcell = slice(m.e, jc);
      // External edges are all oriented CCW, so the leading point of one
      // matches the trailing point of the next.
      const auto d = calc_dist(m, icell[(ie+1)%4], jcell[je]);
      if (d < min_dist) {
        min_dist = d;
        min_j = j;
      }
    }
    if (min_j != i+1)
      std::swap(external_edges[i+1], external_edges[min_j]);
  }
}

template <typename ES>
void fill_perim (LocalMesh<ES>& m) {
  std::vector<Int> external_edges;
  find_external_edges(m, external_edges);
  const Int nee = external_edges.size();
  m.perimp = typename LocalMesh<ES>::Ints("perimp", nee);
  m.perimnml = typename LocalMesh<ES>::Ints("perimnml", nee);
  for (Int k = 0; k < nee; ++k) {
    const auto i = external_edges[k];
    const auto ic = i / 4, ie = i % 4;
    m.perimp(k) = m.e(ic, ie);
    m.perimnml(k) = m.en(ic, ie);
  }
}

template <typename V3> SLMM_KF
void calc_approx_nearest_point_on_arc (
  const V3& p0, const V3& p1, const V3& nml, const Real* point, Real* nearest,
  const bool sphere)
{
  // The approximation is to project point onto the plane of the arc, then to
  // normalize to the arc. If point is on the arc, then nearest = point, as
  // desired.
  using geo = siqk::SphereGeometry;
  if (sphere) {
    const auto dot = geo::dot(point, nml);
    for (Int d = 0; d < 3; ++d) nearest[d] = point[d] - dot*nml[d];
    geo::normalize(nearest);
  } else {
    Real dot = 0;
    for (Int d = 0; d < 3; ++d) dot += (point[d] - p0[d])*nml[d];
    for (Int d = 0; d < 3; ++d) nearest[d] = point[d] - dot*nml[d];
  }
  // If this initial value for nearest is outside of the arc, then 'nearest' is
  // set to the nearer of p0 and p1.
  //   Find the nearest point on line in parameterized coord alpha. alpha is in
  // [0,1] if 'nearest' is on the arc segment.
  Real pdiff[3];
  for (Int d = 0; d < 3; ++d) pdiff[d] = p0[d] - p1[d];
  Real p1mn[3];
  for (Int d = 0; d < 3; ++d) p1mn[d] = nearest[d] - p1[d];
  const Real alpha = geo::dot(p1mn, pdiff) / geo::norm2(pdiff);
  if      (alpha <= 0) for (Int d = 0; d < 3; ++d) nearest[d] = p1[d];
  else if (alpha >= 1) for (Int d = 0; d < 3; ++d) nearest[d] = p0[d];
}

template <typename ES> SLMM_KF
void calc (const LocalMesh<ES>& m, Real* v) {
  using geo = siqk::SphereGeometry;
  ConstExceptGnu Int nedge = m.perimp.size();
  const bool sphere = m.is_sphere();
  const auto canpoa = [&] (const Int& ie, Real* vn) {
    calc_approx_nearest_point_on_arc(slice(m.p, m.perimp(ie)),
                                     slice(m.p, m.perimp((ie+1) % nedge)),
                                     slice(m.nml, m.perimnml(ie)),
                                     v, vn, sphere);
  };
  Real min_dist2 = 1e20;
  Int min_ie = -1;
  for (Int ie = 0; ie < nedge; ++ie) {
    Real vn[3];
    canpoa(ie, vn);
    Real d[3];
    for (Int i = 0; i < 3; ++i) d[i] = vn[i] - v[i];
    const Real dist2 = geo::norm2(d);
    if (dist2 < min_dist2) {
      min_ie = ie;
      min_dist2 = dist2;
    }
  }
  {
    Real vn[3];
    canpoa(min_ie, vn);
    for (Int d = 0; d < 3; ++d) v[d] = vn[d];
  }
}
} // namespace nearest_point

Int unittest(LocalMesh<ko::MachineTraits::HES>& m, const Int tgt_elem,
             const Real length_scale = 1);

std::string to_string(const LocalMesh<ko::MachineTraits::HES>& m);

} // namespace slmm

#endif
