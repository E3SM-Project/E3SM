#include <mpi.h>
#include <cassert>
#include <memory>
#include <vector>

#include "compose_homme.hpp"
#include "compose_test.hpp"

namespace compose {
namespace test {

const char* WindFieldType::inputs[] = {"nondivergent", "divergent"};

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells", "constant"};

std::string format_strings_as_list (const char** strings, const Size n) {
  std::stringstream ss;
  ss << "{";
  for (Size i = 0; i < n-1; ++i) ss << strings[i] << ", ";
  ss << strings[n-1] << "}";
  return ss.str();
}

std::string& tolower (std::string& s) {
  for (auto& c: s)
    c = std::tolower(c);
  return s;
}

// Convert from (u,v), where u is velocity along latitude and v is velocity
// along longitude, to (x,y,z), which is velocity in the global cartesian
// coordinate system. Add a w (local vertical) component to push the position
// (X,Y,Z) back to the unit sphere.
inline void uv2xyz (
  const Real X, const Real Y, const Real Z, // position
  const Real u, const Real v, // velocity in tangent coord system
  Real& x, Real& y, Real& z) // velocity in global coord system
{
  // r should be 1 but will numerically drift, so measure it ...
  const Real r = std::sqrt(X*X + Y*Y + Z*Z);
  // ... and then add a local vertical velocity to project back to the sphere.
  const Real w = (1 - r)/consts::earth_radius_m;
  Real R[9]; // Row major.
  // The local vertical is just the position vector.
  R[2] = X/r; R[5] = Y/r; R[8] = Z/r;
  // The local along-latitude vector.
  R[0] = -Y; R[3] = X; R[6] = 0;
  const Real den = std::sqrt(R[0]*R[0] + R[3]*R[3]);
  R[0] /= den; R[3] /= den;
  // Local vertical x along-latitude.
  R[1] = R[5]*R[6] - R[8]*R[3];
  R[4] = R[8]*R[0] - R[2]*R[6];
  R[7] = R[2]*R[3] - R[5]*R[0];
  // Transform.
  x = R[0]*u + R[1]*v + R[2]*w;
  y = R[3]*u + R[4]*v + R[5]*w;
  z = R[6]*u + R[7]*v + R[8]*w;
}

void InitialCondition
::init (const Shape shape, const Size n, const Real* const lat,
        const Real* const lon, Real* const u) {
  const Real lon1 = 5*(M_PI/6), lat1 = 0, lon2 = -5*(M_PI/6), lat2 = 0;
  Real x1, y1, z1, x2, y2, z2;
  ll2xyz(lat1, lon1, x1, y1, z1);
  ll2xyz(lat2, lon2, x2, y2, z2);
  switch (shape) {
  case XYZTrig: {
    for (Size i = 0; i < n; ++i) {
      Real x, y, z;
      ll2xyz(lat[i], lon[i], x, y, z, 1);
      u[i] = 0.5*(1 + std::sin(3*x)*std::sin(3*y)*std::sin(4*z));
    }
  } break;
  case GaussianHills: {
    for (Size i = 0; i < n; ++i) {
      Real x, y, z;
      ll2xyz(lat[i], lon[i], x, y, z, 1);
      u[i] = GH(x, y, z, x1, y1, z1) + GH(x, y, z, x2, y2, z2);
    }
  } break;
  case CosineBells: {
    const Real r = 0.5, b = 0.1, c = 0.9;
    for (Size i = 0; i < n; ++i) {
      const Real r1 = great_circle_dist(lat[i], lon[i], lat1, lon1);
      Real h = 0;
      if (r1 < r)
        h = CB(r1, r);
      else {
        const Real r2 = great_circle_dist(lat[i], lon[i], lat2, lon2);
        if (r2 < r)
          h = CB(r2, r);
      }
      u[i] = b + c*h;
    }
  } break;
  case SlottedCylinders: {
    const Real b = 0.1, c = 1, R = 1, r = 0.5*R, lon_thr = r/(6*R),
      lat_thr = 5*(r/(12*R));
    for (Size i = 0; i < n; ++i) {
      const Real r1 = great_circle_dist(lat[i], lon[i], lat1, lon1);
      if (r1 <= r) {
        if (std::abs(lon[i] - lon1) >= lon_thr) {
          u[i] = c;
          continue;
        }
        if (std::abs(lon[i] - lon1) < lon_thr && lat[i] - lat1 < -lat_thr) {
          u[i] = c;
          continue;
        }
      }
      const Real r2 = great_circle_dist(lat[i], lon[i], lat2, lon2);
      if (r2 <= r) {
        if (std::abs(lon[i] - lon2) >= lon_thr) {
          u[i] = c;
          continue;
        }
        if (std::abs(lon[i] - lon2) < lon_thr && lat[i] - lat2 > lat_thr) {
          u[i] = c;
          continue;
        }
      }
      u[i] = b;
    }      
  } break;
  case CorrelatedCosineBells: {
    const Real a = -0.8, b = 0.9;
    init(CosineBells, n, lat, lon, u);
    for (Size i = 0; i < n; ++i)
      u[i] = a*square(u[i]) + b;
  } break;
  default: assert(0);
  }
}

struct StandaloneTracersTester {
  typedef std::shared_ptr<StandaloneTracersTester> Ptr;

  struct Plane { Real Sx, Sy, Lx, Ly; };

  const Int np, nlev, qsize, qsize_d, nelemd, testno;
  const bool is_sphere;
  const Plane plane;

  // testno is 0 for the 2D transport test, 1 for 3D. Right now, the 3D test is
  // not actually a valid test; it's used just to create w values for BFB
  // comparison of F90 and C++ trajectory impls. Eventually, it might duplicate
  // test_src/dcmip2012_test1_conv.F90.
  StandaloneTracersTester (Int np_, Int nlev_, Int qsize_, Int qsize_d_,
                           Int nelemd_, Int testno_, Int geometry_,
                           Real Sx, Real Sy, Real Lx, Real Ly)
    : np(np_), nlev(nlev_), qsize(qsize_), qsize_d(qsize_d_), nelemd(nelemd_),
      testno(testno_), is_sphere(geometry_ == 0), plane{Sx, Sy, Lx, Ly}
  {}

  void fill_uniform_density (const Int ie, Real* rdp) const {
    for (int i = 0; i < np*np*nlev; ++i)
      rdp[i] = 1;
  }

  void fill_ics (const Int ie, const Real* rp_elem, const Real* rdp,
                 Real* rqdp) const {
    const homme::FA3<const Real> p_elem(rp_elem, 3, np, np); // (rad, lon, lat)
    const homme::FA3<const Real> dp(rdp, np, np, nlev);
    const homme::FA4<Real> qdp(rqdp, np, np, nlev, qsize_d);
    for (Int q = 0; q < qsize; ++q)
      for (Int k = 0; k < nlev; ++k)
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            Real lat, lon;
            if (is_sphere) {
              lat = p_elem(2,i,j);
              lon = p_elem(1,i,j);
            } else {
              // To get something reasonable quickly, map (x,y) to (lon,lat).
              // rp_elem is still the spherical_polar_t data structure, and
              // plane geometry convention is to fill lon with x, lat with y,
              // and set r to 0.
              assert(p_elem(0,i,j) == 0); // check that the convention still holds
              lon = 2*M_PI*((p_elem(1,i,j) - plane.Sx)/plane.Lx);
              lat = M_PI*((p_elem(2,i,j) - plane.Sy)/plane.Ly - 0.5);
            }
            offset_latlon(nlev, k, lat, lon);
            InitialCondition::init(get_ic(qsize,k,q), 1, &lat, &lon, &qdp(i,j,k,q));
            if (rdp)
              qdp(i,j,k,q) *= dp(i,j,k);
          }
  }

  void fill_v (const Int ie, const Real* rp_elem, const Real t, Real* rv) const {
    const homme::FA3<const Real> p_elem(rp_elem, 3, np, np);
    const homme::FA4<Real> v(rv, np, np, 2, nlev);
    if (is_sphere) {
      NonDivergentWindField f;
      for (Int k = 0; k < nlev; ++k)
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            Real lat = p_elem(2,i,j), lon = p_elem(1,i,j);
            offset_latlon(nlev, k, lat, lon);
            const Real latlon[] = {lat, lon};
            Real uv[2];
            f.eval(t, latlon, uv);
            v(i,j,0,k) = uv[0];
            v(i,j,1,k) = uv[1];
          }
    } else {
      const Real f = 1.0/day2sec(12);
      for (Int k = 0; k < nlev; ++k)
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
#if 1
            v(i,j,0,k) =  2.0*f*plane.Lx;
            v(i,j,1,k) = -1.0*f*plane.Ly;
#else
            v(i,j,0,k) = 0;
            v(i,j,1,k) = 0;
#endif
          }
    }
  }

  // Error data.
  std::vector<Real> l2_num_, l2_den_, wrk_, mass0_, massf_;

  void record_begin () {
#ifdef COMPOSE_HORIZ_OPENMP
#   pragma omp master
#endif
    {
      l2_num_.resize(nlev*qsize*nelemd, 0);
      l2_den_.resize(nlev*qsize*nelemd, 0);
      wrk_.resize(np*np*nlev*qsize*nelemd);
      mass0_.resize(qsize*nelemd, 0);
      massf_.resize(qsize*nelemd, 0);
    }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
  }

  void record (const Int ie, const Real* rp_elem, const Real* rspheremp,
               const Real* rdp, const Real* rqdp) {
    const int tid =
#ifdef COMPOSE_HORIZ_OPENMP
      omp_get_thread_num()
#else
      0
#endif
      ;
    const homme::FA3<const Real> p_elem(rp_elem, 3, np, np);
    const homme::FA2<const Real> spheremp(rspheremp, np, np);
    const homme::FA3<const Real> dp(rdp, np, np, nlev);
    const homme::FA4<const Real> qdp(rqdp, np, np, nlev, qsize_d);
    homme::FA3<Real> l2_num(l2_num_.data(), nelemd, nlev, qsize),
                     l2_den(l2_den_.data(), nelemd, nlev, qsize);
    homme::FA2<Real> mass0(mass0_.data(), nelemd, qsize),
                     massf(massf_.data(), nelemd, qsize);
    Real* const wrk = wrk_.data() + np*np*nlev*qsize*tid;
    fill_ics(ie, rp_elem, nullptr, wrk);
    const homme::FA4<const Real> q0(wrk, np, np, nlev, qsize_d);
    for (Int q = 0; q < qsize; ++q) {
      Real m0 = 0, mf = 0;
      for (Int k = 0; k < nlev; ++k) {
        Real num = 0, den = 0;
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            num += spheremp(i,j)*square(qdp(i,j,k,q)/dp(i,j,k) -
                                        q0(i,j,k,q));
            den += spheremp(i,j)*square(q0(i,j,k,q));
            m0 += spheremp(i,j)*q0(i,j,k,q) /* times rho = 1 */;
            mf += spheremp(i,j)*qdp(i,j,k,q);
          }
        l2_num(ie,k,q) = num;
        l2_den(ie,k,q) = den;
      }
      mass0(ie,q) = m0;
      massf(ie,q) = mf;
    }
  }

  void record_end (const Int fcomm, const Int root, const Int rank,
                   Real* l2_errs, Real* mass_errs) {
#ifdef COMPOSE_HORIZ_OPENMP
#   pragma omp barrier
#endif
    // The way I've written this chunk probably looks funny. But if I put the
    // two pragmas directly adjacent, just in this section of code, Intel 17
    // gives "internal error: assertion failed: find_assoc_pragma: pragma not
    // found (shared/cfe/edgcpfe/il.c, line 23535)".
    const int nr = nlev*qsize;
#ifdef COMPOSE_HORIZ_OPENMP
#   pragma omp master
    {
#endif
      {
        homme::FA2<Real> l2_num(wrk_.data(), nlev, qsize),
                         l2_den(wrk_.data() + nr, nlev, qsize);
        compose_repro_sum(l2_num_.data(), l2_num.data(), nelemd, nr, fcomm);
        compose_repro_sum(l2_den_.data(), l2_den.data(), nelemd, nr, fcomm);
        if (rank == root)
          for (int q = 0, cnt = 0; q < qsize; ++q) {
            printf("COMPOSE>");
            for (int k = 0; k < nlev; ++k) {
              const auto err = std::sqrt(l2_num(k,q)/l2_den(k,q));
              l2_errs[cnt++] = err;
              printf("%23.16e", err);
            }
            printf("\n");
          }
      }
      {
        Real* const mass0 = wrk_.data();
        Real* const massf = mass0 + qsize;
        compose_repro_sum(mass0_.data(), mass0, nelemd, qsize, fcomm);
        compose_repro_sum(massf_.data(), massf, nelemd, qsize, fcomm);
        if (rank == root)
          for (int q = 0; q < qsize; ++q) {
            const auto err = (massf[q] - mass0[q])/mass0[q];
            mass_errs[q] = err;
            printf("COMPOSE> mass0 %8.2e mass re %9.2e\n", mass0[q], err);
          }
      }  
#ifdef COMPOSE_HORIZ_OPENMP
    }
#endif
  }
};

static StandaloneTracersTester::Ptr g_stt;
} // namespace test
} // namespace compose

using namespace compose::test;

extern "C" void compose_unittest () { 
  slmm_unittest();
  cedr_unittest();
}

extern "C" void compose_stt_init (
  Int np, Int nlev, Int qsize, Int qsize_d, Int nelemd, Int testno,
  // Data used only if geometry = 1 (plane); 0 is sphere.
  Int geometry, Real Sx, Real Sy, Real Lx, Real Ly)
{
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  g_stt = std::make_shared<StandaloneTracersTester>(
    np, nlev, qsize, qsize_d, nelemd, testno, geometry, Sx, Sy, Lx, Ly);
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

extern "C" void compose_stt_fill_uniform_density (
  Int ie, Int np1, Real* dp3d, Real* dp)
{
  const size_t os = square(g_stt->np) * g_stt->nlev * (np1 - 1);
  g_stt->fill_uniform_density(ie-1, dp3d + os);
  g_stt->fill_uniform_density(ie-1, dp);
}

extern "C" void compose_stt_fill_ics (
  Int ie, Real* p_elem, Real* dp, Int n0_qdp, Real* qdp)
{
  g_stt->fill_ics(ie-1, reinterpret_cast<const Real*>(p_elem), dp,
                  qdp + square(g_stt->np) * g_stt->nlev * g_stt->qsize_d *
                  (n0_qdp - 1));
}

extern "C" void compose_stt_fill_v (Int ie, Real* p_elem, Real* t, Real* v) {
  g_stt->fill_v(ie-1, reinterpret_cast<const Real*>(p_elem), *t, v);
}

extern "C" void compose_stt_begin_record () {
  g_stt->record_begin();
}

extern "C" void compose_stt_record_q (
  Int ie, Real* p_elem, Real* spheremp, Int np1, Real* dp3d, Int n0_qdp,
  Real* qdp)
{
  g_stt->record(ie-1, reinterpret_cast<const Real*>(p_elem), spheremp,
                dp3d + square(g_stt->np) * g_stt->nlev * (np1 - 1),
                qdp + square(g_stt->np) * g_stt->nlev * g_stt->qsize_d *
                (n0_qdp - 1));
}

extern "C" void compose_stt_clear () {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  g_stt = nullptr;
}

extern "C" void compose_stt_finish (Int fcomm, Int root, Int rank, Real* eval) {
  g_stt->record_end(fcomm, root, rank, eval, eval + g_stt->nlev * g_stt->qsize);
  compose_stt_clear();
}
