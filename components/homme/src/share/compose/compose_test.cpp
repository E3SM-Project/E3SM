#include <mpi.h>
#include <cassert>
#include <memory>
#include <vector>

#include <Kokkos_Core.hpp>

#ifdef _OPENMP
# include <omp.h>
#endif

#include "compose_test.hpp"

#define pr(m) do {                              \
    std::stringstream _ss_;                     \
    _ss_ << m << std::endl;                     \
    std::cerr << _ss_.str();                    \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m)"(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template<typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << ": ";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "\n";
}
#define mprarr(m) siqk::prarr(#m, m.data(), m.size())

namespace {
typedef double Real;
typedef int Int;
typedef Int Size;

template <typename T> inline constexpr T square (const T& x) { return x*x; }
template<typename T> inline constexpr T cube (const T& x) { return x*x*x; }

struct consts {
  static constexpr Real earth_radius_m = 6.376e6;
};

template<typename T> inline T sign (const T& a) { return a >= 0 ? 1 : -1; }

inline Real sec2day (const Real sec) { return sec/(24*3600); }
inline Real day2sec (const Real day) { return day*(24*3600); }

// Output is in radians.
//todo Make a version that lets you pass R = mag(x,y,z).
inline void xyz2ll (const Real x, const Real y, const Real z,
                    Real& lat, Real& lon) {
  const Real r = std::sqrt(square(x) + square(y) + square(z));
  lat = std::asin(z/r);
  lon = std::atan2(y, x);
}

// Input is in radians.
inline void ll2xyz (const Real lat, const Real lon, Real& x, Real& y, Real& z,
                    const Real radius = 1) {
  const Real sinl = std::sin(lat), cosl = std::cos(lat);
  x = radius*std::cos(lon)*cosl;
  y = radius*std::sin(lon)*cosl;
  z = radius*sinl;
}

// Eq after eq 10 in Lauritzen et al test cases paper.
inline Real great_circle_dist (
  const Real lat1, const Real lon1, const Real lat2, const Real lon2,
  const Real R = 1)
{
  Real xA, yA, zA;
  Real xB, yB, zB;
  ll2xyz(lat1, lon1, xA, yA, zA);
  ll2xyz(lat2, lon2, xB, yB, zB);
  Real cp1, cp2, cp3, cpnorm, dotprod;
  cp1 = yA*zB - yB*zA;
	cp2 = xB*zA - xA*zB;
	cp3 = xA*yB - xB*yA;
	cpnorm = std::sqrt(cp1*cp1 + cp2*cp2 + cp3*cp3);
	dotprod = xA*xB + yA*yB + zA*zB;	
	return R * std::atan2(cpnorm, dotprod);
}

inline constexpr Real m2radlat (const Real m)
{ return m/consts::earth_radius_m; }

inline Real m2radlon(const Real lat, const Real m)
{ return m2radlat(m)/std::abs(std::cos(lat)); }

inline constexpr Real deg2rad (const Real v) { return v * (M_PI/180); }
inline constexpr Real rad2deg (const Real v) { return v * (180/M_PI); }

inline Real reldif (const Real a, const Real b, const Real abstol = 0)
{ return std::abs(b - a)/(abstol + std::abs(a)); }

std::string& tolower (std::string& s) {
  for (auto& c: s)
    c = std::tolower(c);
  return s;
}

std::string format_strings_as_list (const char** strings, const Size n) {
  std::stringstream ss;
  ss << "{";
  for (Size i = 0; i < n-1; ++i) ss << strings[i] << ", ";
  ss << strings[n-1] << "}";
  return ss.str();
}

class OdeFn {
  mutable int ne_;
public:
  OdeFn () : ne_(0) {}
  virtual bool eval (const Real t, const Real* const d, Real* const f) const = 0;
  virtual void record (const Real t, const Real* const y) const { ++ne_; }
  int ne () const { return ne_; }
};

// From Lauritzen et al, A standard test case suite for two-dimensional linear
// transport on the sphere, Geosci. Model Dev., 2012.
class InitialCondition {
  static const char* inputs[];

  static inline Real GH (const Real x, const Real y, const Real z,
                         const Real xi, const Real yi, const Real zi) {
    const Real h_max = 0.95, b = 5;
    const Real r2 = (square(x - xi) + square(y - yi) +
                     square(z - zi));
    return h_max*std::exp(-b*r2);
  }

  static inline Real CB (const Real ri, const Real r) {
    const Real h_max = 1;
    return 0.5*h_max*(1 + std::cos(M_PI*ri/r));
  }

public:
  enum Shape : int {
    XYZTrig = 0, GaussianHills, CosineBells, SlottedCylinders,
    CorrelatedCosineBells, nshapes
  };

  static Shape from_string (const std::string& si) {
    std::string s(si);
    tolower(s);
    if (s == inputs[0]) return XYZTrig;
    if (s == inputs[1]) return GaussianHills;
    if (s == inputs[2]) return CosineBells;
    if (s == inputs[3]) return SlottedCylinders;
    if (s == inputs[4]) return CorrelatedCosineBells;
    throw std::runtime_error(si + " is not an initial condition.");
  }

  static const char* to_string (const Shape& shape) {
    switch (shape) {
    case XYZTrig: return inputs[0];
    case GaussianHills: return inputs[1];
    case CosineBells: return inputs[2];
    case SlottedCylinders: return inputs[3];
    case CorrelatedCosineBells: return inputs[4];
    }
    throw std::runtime_error("Should not be here.");
  }

  // Input is in radians.
  static void init (const Shape shape, const Size n, const Real* const lat,
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

  static std::string get_inputs ()
  { return format_strings_as_list(inputs, 7); }
};

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells", "constant"};

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

// Also from Lauritzen et al.
struct NonDivergentWindField : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
    const Real theta = d[0];  // latitude
    const Real lambda = d[1]; // longitude
    const Real
      T = day2sec(12),
      R = consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // u
    f[0] = R/T*(10*square(std::sin(lambda_p))*std::sin(2*theta)*cost +
                2*M_PI*costh);
    // v
    f[1] = 10*R/T*std::sin(2*lambda_p)*costh*cost;
    return true;
  }
};

// Also from Lauritzen et al.
struct DivergentWindField : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
    const Real
      theta = d[0],  // latitude
      lambda = d[1], // longitude
      T = day2sec(12),
      R = consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // u
    f[0] = R/T*(-5*square(std::sin(0.5*lambda_p))*std::sin(2*theta)*
                square(costh)*cost + 2*M_PI*costh);
    // v
    f[1] = 2.5*R/T*std::sin(lambda_p)*cube(costh)*cost;
    return true;
  }
};

struct WindFieldType {
  static const char* inputs[];
public:
  enum Enum { NonDivergentWindField, DivergentWindField };
  static Enum from_string (const std::string& si) {
    std::string s(si);
    tolower(s);
    if (s == inputs[1]) return NonDivergentWindField;
    if (s == inputs[2]) return DivergentWindField;
    throw std::runtime_error(si + " is not an ODE function.");
  }
  static std::string get_inputs ()
  { return format_strings_as_list(inputs, 6); }
};

const char* WindFieldType::inputs[] =
  {"nondivergent", "divergent"};

namespace ko = Kokkos;

// Fortran array wrappers.
template <typename T> using FA2 = ko::View<T**,   ko::LayoutLeft,ko::HostSpace>;
template <typename T> using FA3 = ko::View<T***,  ko::LayoutLeft,ko::HostSpace>;
template <typename T> using FA4 = ko::View<T****, ko::LayoutLeft,ko::HostSpace>;
template <typename T> using FA5 = ko::View<T*****,ko::LayoutLeft,ko::HostSpace>;

struct StandaloneTracersTester {
  typedef std::shared_ptr<StandaloneTracersTester> Ptr;

  const Int np, nlev, qsize, qsize_d, nelemd;

  StandaloneTracersTester (Int np_, Int nlev_, Int qsize_, Int qsize_d_,
                           Int nelemd_)
    : np(np_), nlev(nlev_), qsize(qsize_), qsize_d(qsize_d_), nelemd(nelemd_)
  {}

  void fill_uniform_density (const Int ie, Real* rdp) const {
    for (int i = 0; i < np*np*nlev; ++i)
      rdp[i] = 1;
  }

  InitialCondition::Shape get_ic (const Int k, const Int q) const {
    return static_cast<InitialCondition::Shape>((k*qsize + q) %
                                                InitialCondition::nshapes);
  }

  // Make the winds differ by elevation.
  void offset_latlon (const Int k, Real& lat, Real& lon) const {
    lon += k*(0.4/nlev);
  }

  void fill_ics (const Int ie, const Real* rp_elem, const Real* rdp,
                 Real* rqdp) const {
    const FA3<const Real> p_elem(rp_elem, 3, np, np); // (rad, lon, lat)
    const FA3<const Real> dp(rdp, np, np, nlev);
    const FA4<Real> qdp(rqdp, np, np, nlev, qsize_d);
    for (Int q = 0; q < qsize; ++q)
      for (Int k = 0; k < nlev; ++k)
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            Real lat = p_elem(2,i,j), lon = p_elem(1,i,j);
            offset_latlon(k, lat, lon);
            InitialCondition::init(get_ic(k,q), 1, &lat, &lon, &qdp(i,j,k,q));
            if (rdp)
              qdp(i,j,k,q) *= dp(i,j,k);
          }
  }

  void fill_v (const Int ie, const Real* rp_elem, const Real t, Real* rv) const {
    const FA3<const Real> p_elem(rp_elem, 3, np, np);
    const FA4<Real> v(rv, np, np, 2, nlev);
    NonDivergentWindField f;
    for (Int k = 0; k < nlev; ++k)
      for (Int j = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i) {
          Real lat = p_elem(2,i,j), lon = p_elem(1,i,j);
          offset_latlon(k, lat, lon);
          const Real latlon[] = {lat, lon};
          Real uv[2];
          f.eval(t, latlon, uv);
          v(i,j,0,k) = uv[0];
          v(i,j,1,k) = uv[1];
        }
  }

#if 0
  // Debug output on a single rank.
  std::vector<Real> d_;

  void record_begin () {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
    d_.resize(3 * nelemd * np*np);
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
  }

  void record (const Int ie, const Real* rp_elem, const Real* spheremp,
               const Real* rdp, const Real* rqdp) {
    const FA3<const Real> p_elem(rp_elem, 3, np, np);
    const FA4<const Real> qdp(rqdp, np, np, nlev, qsize_d);
    const Int q = 0, k = 0;
    Real* d = d_.data() + 3*np*np*ie;
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i) {
        d[3*(np*j + i) + 0] = p_elem(2,i,j);
        d[3*(np*j + i) + 1] = p_elem(1,i,j);
        d[3*(np*j + i) + 2] = qdp(i,j,k,q);
      }
  }

  void record_end () {
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
    {
#endif
      FILE* fid = fopen("ct.txt", "w");
      for (size_t i = 0; i < d_.size()/3; ++i)
        fprintf(fid, "%1.15e %1.15e %1.15e\n",
                d_[3*i+0], d_[3*i+1], d_[3*i+2]);
      fclose(fid);
#ifdef HORIZ_OPENMP
    }
#endif
  }
#else
  // Error data.
  std::vector<Real> l2_err_, wrk_, mass0_, massf_;

  void record_begin () {
    int nthr = 1;
#ifdef HORIZ_OPENMP
    nthr = omp_get_num_threads();
# pragma omp master
#endif
    {
      l2_err_.resize(2*nlev*qsize*nthr, 0);
      wrk_.resize(np*np*nlev*qsize*nthr);
      mass0_.resize(qsize*nthr, 0);
      massf_.resize(qsize*nthr, 0);
    }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
  }

  void record (const Int ie, const Real* rp_elem, const Real* rspheremp,
               const Real* rdp, const Real* rqdp) {
    const int tid =
#ifdef HORIZ_OPENMP
      omp_get_thread_num()
#else
      0
#endif
      ;
    const FA3<const Real> p_elem(rp_elem, 3, np, np);
    const FA2<const Real> spheremp(rspheremp, np, np);
    const FA3<const Real> dp(rdp, np, np, nlev);
    const FA4<const Real> qdp(rqdp, np, np, nlev, qsize_d);
    Real* const l2_err = l2_err_.data() + 2*nlev*qsize*tid;
    Real* const wrk = wrk_.data() + np*np*nlev*qsize*tid;
    Real* const mass0 = mass0_.data() + qsize*tid;
    Real* const massf = massf_.data() + qsize*tid;
    fill_ics(ie, rp_elem, nullptr, wrk);
    const FA4<const Real> q0(wrk, np, np, nlev, qsize_d);
    for (Int q = 0; q < qsize; ++q)
      for (Int k = 0; k < nlev; ++k) {
        Real num = 0, den = 0;
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            num += spheremp(i,j)*square(qdp(i,j,k,q)/dp(i,j,k) -
                                        q0(i,j,k,q));
            den += spheremp(i,j)*square(q0(i,j,k,q));
            mass0[q] += spheremp(i,j)*q0(i,j,k,q) /* times rho = 1 */;
            massf[q] += spheremp(i,j)*qdp(i,j,k,q);
          }
        l2_err[2*nlev*q + 2*k    ] += num;
        l2_err[2*nlev*q + 2*k + 1] += den;
      }
  }

  void record_end (MPI_Comm comm, const Int root, const Int rank) {
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
    // The way I've written this chunk probably looks funny. But if I put the
    // two pragmas directly adjacent, just in this section of code, Intel 17
    // gives "internal error: assertion failed: find_assoc_pragma: pragma not
    // found (shared/cfe/edgcpfe/il.c, line 23535)".
    const int nr = 2*nlev*qsize;
#ifdef HORIZ_OPENMP
# pragma omp master
    {
#endif
      const int nthr = wrk_.size()/(np*np*nlev*qsize);
      for (int i = 1; i < nthr; ++i) {
        for (int j = 0; j < nr; ++j)
          l2_err_[j] += l2_err_[nr*i + j];
        for (int j = 0; j < qsize; ++j) {
          mass0_[j] += mass0_[qsize*i + j];
          massf_[j] += massf_[qsize*i + j];
        }
      }
      MPI_Reduce(l2_err_.data(), wrk_.data(), nr, MPI_DOUBLE, MPI_SUM,
                 root, comm);
      if (rank == root)
        for (int q = 0; q < qsize; ++q) {
          printf("COMPOSE>");
          for (int k = 0; k < nlev; ++k)
            printf("%9.2e", std::sqrt(wrk_[2*nlev*q + 2*k] /
                                      wrk_[2*nlev*q + 2*k + 1]));
          printf("\n");
        }
      MPI_Reduce(mass0_.data(), wrk_.data(), qsize, MPI_DOUBLE, MPI_SUM,
                 root, comm);
      MPI_Reduce(massf_.data(), wrk_.data() + qsize, qsize, MPI_DOUBLE, MPI_SUM,
                 root, comm);
      if (rank == root)
        for (int q = 0; q < qsize; ++q) {
          printf("COMPOSE>");
          const Real mass0 = wrk_[q], massf = wrk_[qsize + q];
          printf(" mass0 %8.2e mass re %9.2e\n", mass0, (massf - mass0)/mass0);
        }
#ifdef HORIZ_OPENMP
    }
#endif
  }
#endif
};
 
static StandaloneTracersTester::Ptr g_stt;
} // namespace anon

extern "C" void compose_unittest () {
  slmm_unittest();
}

extern "C" void compose_stt_init (
  Int np, Int nlev, Int qsize, Int qsize_d, Int nelemd)
{
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
    g_stt = std::make_shared<StandaloneTracersTester>(
      np, nlev, qsize, qsize_d, nelemd);
#ifdef HORIZ_OPENMP
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

extern "C" void compose_stt_finish (Int fcomm, Int root, Int rank) {
  g_stt->record_end(MPI_Comm_f2c(fcomm), root, rank);
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  g_stt = nullptr;
}
