#ifndef SCREAM_COSP_TEST_UTILS_HPP
#define SCREAM_COSP_TEST_UTILS_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <string>

/*
 * Helpers shared by the BFB tests of the C++/Kokkos COSP port against the Fortran.
 * The tests always use double, since the Fortran library is double (wp = dp).
 *
 * On CPU builds the C++ must match the Fortran exactly. On GPU builds the C++ runs on
 * device, with fused multiply-adds and the device math library, while the Fortran runs
 * on host, so only close agreement is possible: values must agree to a relative
 * tolerance, and a small fraction of entries may differ by more (a roundoff difference
 * can flip a discrete choice, such as a histogram bin or a retrieval succeeding).
 */

// C bindings to the Fortran routines, see cosp_f90_shims.F90
extern "C" {
void cosp_rng_f(int n, int ndraw, const int* seeds, double* rnd);
void cosp_simulator_optics_f(int dim1, int dim2, int dim3, const double* flag,
                             const double* var1, const double* var2, double* varout);
void cosp_modis_optics_partition_f(int npoints, int ncolumns, int nlev,
                                   const double* cloud_water, const double* cloud_ice,
                                   const double* water_size, const double* ice_size,
                                   const double* tau, double* tau_liq, double* tau_ice);
void cosp_modis_optics_f(int npoints, int ncolumns, int nlev,
                         const double* tau_liq, const double* size_liq,
                         const double* tau_ice, const double* size_ice,
                         double* frac_liq, double* g, double* w0);
void cosp_nir_fits_f(int n, int phase, const double* re, double* g, double* w0);
void cosp_hist2d_f(int npts, const double* var1, const double* var2,
                   int nbin1, const double* bin1, int nbin2, const double* bin2, double* hist);
void cosp_constants_f(double* r_undef, double* tau_bnds, double* pres_bnds, double* hgt_bnds,
                      double* misr_tauchk);
void cosp_misr_subcolumn_f(int npoints, int ncol, int nlev, const double* dtau,
                           const double* zfull, const double* at, const int* sunlit,
                           double* tau_out, double* ztop);
void cosp_isccp_constants_f(double* tauchk);
void cosp_icarus_f(int npoints, int ncol, int nlev, const int* sunlit, const double* dtau,
                   const double* dem, const double* skt, double emsfc_lw, const double* qv,
                   const double* at, const double* pfull, const double* phalf,
                   double* boxtau, double* boxptop, double* boxttop, double* fq,
                   double* totalcldarea);
void cosp_misr_column_f(int npoints, int ncol, const double* ztop, const int* sunlit,
                        const double* tau, double* fq);
void cosp_modis_subcolumn_f(int npoints, int ncol, int nlev, const double* phalf,
                            const double* tau, const double* liqfrac, const double* g,
                            const double* w0, const double* isccp_ctp, int* phase, double* ctp,
                            double* tau_ret, double* size_ret);
void cosp_scops_f(int npoints, int nlev, int ncol, const int* seeds, const double* cc,
                  const double* conv, int overlap, double* frac_out);
void cosp_two_stream_f(int n, const double* tau, const double* g, const double* w0,
                       double* ref, double* tra, double* ref2);
void cosp_modis_column_f(int npoints, int ncol, const int* phase, const double* ctp,
                         const double* tau, const double* size_ret, double* fq);
}

namespace scream {
namespace cosp_test {

using F        = cosp::Functions<double, DefaultDevice>;
using KT       = ekat::KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;

// LayoutLeft host views, matching Fortran memory order
using lview_1d  = Kokkos::View<double*,   Kokkos::LayoutLeft, Kokkos::HostSpace>;
using lview_2d  = Kokkos::View<double**,  Kokkos::LayoutLeft, Kokkos::HostSpace>;
using lview_3d  = Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using liview_1d = Kokkos::View<int*,      Kokkos::LayoutLeft, Kokkos::HostSpace>;
using liview_2d = Kokkos::View<int**,     Kokkos::LayoutLeft, Kokkos::HostSpace>;

// Device views in the layout the port uses
using dview_1d  = Kokkos::View<double*,   DefaultDevice>;
using dview_2d  = Kokkos::View<double**,  DefaultDevice>;
using dview_3d  = Kokkos::View<double***, DefaultDevice>;
using diview_1d = Kokkos::View<int*,      DefaultDevice>;
using diview_2d = Kokkos::View<int**,     DefaultDevice>;

using Policy2 = Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>>;
using Policy3 = Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>>;

// Whether the C++ is required to match the Fortran exactly (CPU builds), see above
constexpr bool exact_bfb = !ekat::OnGpu<ExeSpace>::value;

// Whether a C++ value matches the Fortran one. The GPU tolerance is loose enough for
// the single precision arithmetic in the RNG.
inline bool same (const double a, const double b) {
  if (exact_bfb) return a == b;
  constexpr double rtol = 1e-6, atol = 1e-12;
  const double d = std::abs(a - b);
  return d <= atol || d <= rtol*std::max(std::abs(a), std::abs(b));
}

// Number of entries out of n that may fail same(): none on CPU, 0.1% (at least 1) on GPU.
// Roundoff flips are rare (about 5e-6 of the MODIS histogram entries in a month-long
// Frontier run), so this still catches real bugs that affect a few percent of entries.
inline int allowed_diffs (const size_t n) {
  return exact_bfb ? 0 : std::max(1, static_cast<int>(n/1000));
}

inline std::mt19937_64& engine () {
  static std::mt19937_64 e (20261022);
  return e;
}

inline double uniform (const double lo, const double hi) {
  return std::uniform_real_distribution<double>(lo,hi)(engine());
}

// Fill a host view, setting a fraction of the entries to exactly zero
template <typename V>
void fill_random (const V& v, const double lo, const double hi, const double frac_zero = 0) {
  auto* p = v.data();
  for (size_t i = 0; i < v.size(); ++i) {
    p[i] = uniform(0,1) < frac_zero ? 0.0 : uniform(lo,hi);
  }
}

// Copy a LayoutLeft host view to a (LayoutRight) device view of the same extents
template <typename DV, typename HV>
DV to_device (const HV& h, const std::string& name) {
  DV d;
  if constexpr (DV::rank == 1) d = DV(name, h.extent(0));
  if constexpr (DV::rank == 2) d = DV(name, h.extent(0), h.extent(1));
  if constexpr (DV::rank == 3) d = DV(name, h.extent(0), h.extent(1), h.extent(2));
  auto dh = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(dh, h);
  Kokkos::deep_copy(d, dh);
  return d;
}

// Count entries of a device view that don't match a LayoutLeft host view, see same()
template <typename DV, typename HV>
int count_diffs (const DV& d, const HV& h) {
  auto dh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d);
  int ndiff = 0;
  if constexpr (HV::rank == 1) {
    for (size_t i = 0; i < h.extent(0); ++i)
      if (!same(dh(i), h(i))) ++ndiff;
  } else if constexpr (HV::rank == 2) {
    for (size_t i = 0; i < h.extent(0); ++i)
      for (size_t j = 0; j < h.extent(1); ++j)
        if (!same(dh(i,j), h(i,j))) ++ndiff;
  } else {
    for (size_t i = 0; i < h.extent(0); ++i)
      for (size_t j = 0; j < h.extent(1); ++j)
        for (size_t k = 0; k < h.extent(2); ++k)
          if (!same(dh(i,j,k), h(i,j,k))) ++ndiff;
  }
  return ndiff;
}

// Check that a device view matches the Fortran result in a host view, see above
#define COSP_REQUIRE_MATCH(d, h) REQUIRE(count_diffs(d, h) <= allowed_diffs(h.size()))

// Rough atmospheric profiles for npoints columns of nlev levels, ordered top to bottom
// as COSP expects: mid-level height (m) and temperature (K)
inline void make_profiles (const lview_2d& zfull, const lview_2d& at) {
  const int npoints = zfull.extent(0), nlev = zfull.extent(1);
  for (int j = 0; j < npoints; ++j) {
    const double zsfc = uniform(0,1) < 0.3 ? uniform(0,3000) : 0.0;
    const double tsfc = uniform(240,305);
    for (int k = 0; k < nlev; ++k) {
      // Levels get thinner toward the surface, roughly like a hybrid grid
      const double s = double(nlev-1-k)/(nlev-1);  // 1 at the top, 0 at the bottom
      zfull(j,k) = zsfc + 40.0 + 40000.0*s*s;
      at(j,k)    = std::max(tsfc - 6.5e-3*(zfull(j,k)-zsfc), 200.0);
    }
  }
}

// Thermodynamic profiles for npoints columns of nlev levels, top to bottom: interface and
// mid-level pressure (Pa), temperature (K), specific humidity (kg/kg) and skin temperature
inline void make_thermo (const lview_2d& phalf, const lview_2d& pfull, const lview_2d& at,
                         const lview_2d& qv, const lview_1d& skt) {
  const int npoints = pfull.extent(0), nlev = pfull.extent(1);
  for (int j = 0; j < npoints; ++j) {
    const double psfc = uniform(0,1) < 0.2 ? uniform(55000,85000) : uniform(95000,103500);
    const double tsfc = uniform(235,305);
    const double ttrop = uniform(190,225);
    const double ptop = 225.0;
    for (int k = 0; k <= nlev; ++k) {
      // Spacing that is finer near the surface and in the upper troposphere
      const double s = double(k)/nlev;
      phalf(j,k) = ptop + (psfc - ptop)*(0.6*s + 0.4*s*s*s);
    }
    for (int k = 0; k < nlev; ++k) {
      const double p = 0.5*(phalf(j,k) + phalf(j,k+1));
      pfull(j,k) = p;
      const double t_trop = tsfc*std::pow(p/psfc, 0.19);
      // Isothermal tropopause layer, warming above about 100 hPa
      at(j,k) = std::max(t_trop, ttrop + (p < 10000 ? 20.0*(1 - p/10000) : 0.0));
      qv(j,k) = 0.8*3.8e-3*std::exp(17.27*(at(j,k)-273.15)/(at(j,k)-35.85))*(1e5/p)*uniform(0.1,1.0);
      qv(j,k) = std::min(qv(j,k), 0.03);
    }
    skt(j) = tsfc + uniform(-5,5);
  }
}

// Subcolumn optical depths with a mix of clear, thin and thick cloud
inline void make_dtau (const lview_3d& dtau) {
  const int npoints = dtau.extent(0), ncol = dtau.extent(1), nlev = dtau.extent(2);
  for (int j = 0; j < npoints; ++j) {
    const double cld_prob = uniform(0,1) < 0.15 ? 0.0 : uniform(0,0.4);
    const double scale    = std::pow(10.0, uniform(-2.5,1.5));
    for (int i = 0; i < ncol; ++i) {
      for (int k = 0; k < nlev; ++k) {
        dtau(j,i,k) = uniform(0,1) < cld_prob ? scale*uniform(0,1) : 0.0;
      }
    }
  }
}

} // namespace cosp_test
} // namespace scream

#endif // SCREAM_COSP_TEST_UTILS_HPP
