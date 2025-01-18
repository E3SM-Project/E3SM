#ifndef INCLUDE_COMPOSE_TEST_HPP
#define INCLUDE_COMPOSE_TEST_HPP

#include <mpi.h>

extern "C"
void compose_repro_sum(const double* send, double* recv,
                       int nlocal, int nfld, int fcomm);

namespace compose {
namespace test {

int slmm_unittest();
int cedr_unittest();
int cedr_unittest(MPI_Comm comm);
int interpolate_unittest();

typedef double Real;
typedef int Int;
typedef Int Size;

namespace ko = Kokkos;

template <typename T> KOKKOS_INLINE_FUNCTION T square (const T& x) { return x*x; }
template <typename T> KOKKOS_INLINE_FUNCTION T cube (const T& x) { return x*x*x; }

struct consts {
  static constexpr Real earth_radius_m = 6.376e6;
};

template<typename T> inline T sign (const T& a) { return a >= 0 ? 1 : -1; }

KOKKOS_INLINE_FUNCTION Real sec2day (const Real sec) { return sec/(24*3600); }
KOKKOS_INLINE_FUNCTION Real day2sec (const Real day) { return day*(24*3600); }

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

std::string& tolower(std::string& s);
std::string format_strings_as_list(const char** strings, const Size n);

class OdeFn {
  mutable int ne_;
public:
  OdeFn () : ne_(0) {}
  KOKKOS_FUNCTION virtual bool eval(const Real t, const Real* const d,
                                    Real* const f) const = 0;
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
    case nshapes: break;
    }
    throw std::runtime_error("Should not be here.");
  }

  // Input is in radians.
  static void init(const Shape shape, const Size n, const Real* const lat,
                   const Real* const lon, Real* const u);

  static std::string get_inputs ()
  { return format_strings_as_list(inputs, 7); }
};

// Also from Lauritzen et al.
struct NonDivergentWindField : public OdeFn {
  NonDivergentWindField () {}
  KOKKOS_FUNCTION
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
  DivergentWindField () {}
  KOKKOS_FUNCTION
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

// Function to add some additional variation to the test: make the winds differ
// by elevation.
KOKKOS_INLINE_FUNCTION void
offset_latlon (const Int nlev, const Int k, Real& lat, Real& lon) {
  lon += k*(0.4/nlev);
}

inline InitialCondition::Shape get_ic (const Int qsize, const Int k, const Int q) {
  return static_cast<InitialCondition::Shape>((k*qsize + q) % InitialCondition::nshapes);
}

} // namespace test
} // namespace compose

#endif
