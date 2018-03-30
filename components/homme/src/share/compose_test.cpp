#include <cassert>

#include <Kokkos_Core.hpp>

#include "compose_test.hpp"

namespace {
typedef double Real;
typedef int Int;
typedef Int Size;

template <typename T> inline constexpr T square (const T& x) { return x*x; }
template<typename T> inline constexpr T cube (const T& x) { return x*x*x; }

struct consts {
  static constexpr Real earth_radius_m = 6.37122e6;
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
  bool xyz_form_;
public:
  OdeFn () : ne_(0), xyz_form_(false) {}
  virtual bool eval (const Real t, const Real* const d, Real* const f) const = 0;
  virtual void record (const Real t, const Real* const y) const { ++ne_; }
  int ne () const { return ne_; }
  void set_xyz_form (const bool use_xyz_form) { xyz_form_ = use_xyz_form; }
  bool use_xyz_form () const { return xyz_form_; }
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
  enum Shape {
    XYZTrig, GaussianHills, CosineBells, SlottedCylinders,
    CorrelatedCosineBells, Constant
  };

  static Shape from_string (const std::string& si) {
    std::string s(si);
    tolower(s);
    if (s == inputs[0]) return XYZTrig;
    if (s == inputs[1]) return GaussianHills;
    if (s == inputs[2]) return CosineBells;
    if (s == inputs[3]) return SlottedCylinders;
    if (s == inputs[4]) return CorrelatedCosineBells;
    if (s == inputs[5]) return Constant;
    throw std::runtime_error(si + " is not an initial condition.");
  }

  static const char* to_string (const Shape& shape) {
    switch (shape) {
    case XYZTrig: return inputs[0];
    case GaussianHills: return inputs[1];
    case CosineBells: return inputs[2];
    case SlottedCylinders: return inputs[3];
    case CorrelatedCosineBells: return inputs[4];
    case Constant: return inputs[5];
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
    case Constant: {
      for (Size i = 0; i < n; ++i)
        u[i] = 0.42;
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
    Real theta, lambda;
    if (use_xyz_form())
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    else {
      theta = d[0];  // latitude
      lambda = d[1]; // longitude
    }
    const Real
      T = day2sec(12),
      R = consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // v
    f[0] = 10*R/T*std::sin(2*lambda_p)*costh*cost;
    // u
    f[1] = R/T*(10*square(std::sin(lambda_p))*std::sin(2*theta)*cost +
                2*M_PI*costh);
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {  
      f[0] = m2radlat(f[0]);
      f[1] = m2radlon(theta, f[1]);
    }
    return true;
  }
};

// Also from Lauritzen et al.
struct DivergentWindField : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
    Real theta, lambda;
    if (use_xyz_form())
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    else {
      theta = d[0];  // latitude
      lambda = d[1]; // longitude
    }
    const Real
      T = day2sec(12),
      R = consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // v
    f[0] = 2.5*R/T*std::sin(lambda_p)*cube(costh)*cost;
    // u
    f[1] = R/T*(-5*square(std::sin(0.5*lambda_p))*std::sin(2*theta)*
                square(costh)*cost + 2*M_PI*costh);
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {  
      f[0] = m2radlat(f[0]);
      f[1] = m2radlon(theta, f[1]);
    }
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
 
struct Input {

};

void run_time_integration (const Input& in) {

}
} // namespace anon

extern "C" void compose_unittest_ () {
  slmm_unittest();
  Input in;
  run_time_integration(in);
}
