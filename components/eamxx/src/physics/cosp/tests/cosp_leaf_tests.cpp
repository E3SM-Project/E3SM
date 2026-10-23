#include <catch2/catch.hpp>

#include "physics/cosp/cosp_cxx_functions.hpp"

#include <random>
#include <vector>

/*
 * BFB tests of the leaf routines of the C++/Kokkos COSP port against the Fortran COSP
 * library. Both sides get identical inputs; results must match exactly.
 *
 * The tests always use double, since the Fortran library is double (wp = dp).
 */

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
}

namespace {

using namespace scream;

using F   = cosp::Functions<double, DefaultDevice>;
using KT  = ekat::KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;

// LayoutLeft host views, matching Fortran memory order
using lview_1d = Kokkos::View<double*,   Kokkos::LayoutLeft, Kokkos::HostSpace>;
using lview_2d = Kokkos::View<double**,  Kokkos::LayoutLeft, Kokkos::HostSpace>;
using lview_3d = Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

// Device views in the layout the port uses
using dview_1d = Kokkos::View<double*,   DefaultDevice>;
using dview_2d = Kokkos::View<double**,  DefaultDevice>;
using dview_3d = Kokkos::View<double***, DefaultDevice>;

using Policy3 = Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>>;

std::mt19937_64 engine (20261022);

double uniform (const double lo, const double hi) {
  return std::uniform_real_distribution<double>(lo,hi)(engine);
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
  DV d(name, h.extent(0), h.extent(1), h.extent(2));
  auto dh = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(dh, h);
  Kokkos::deep_copy(d, dh);
  return d;
}

// Count entries of a device view that differ from a LayoutLeft host view
template <typename DV, typename HV>
int count_diffs (const DV& d, const HV& h) {
  auto dh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d);
  int ndiff = 0;
  for (size_t i = 0; i < h.extent(0); ++i)
    for (size_t j = 0; j < h.extent(1); ++j)
      for (size_t k = 0; k < h.extent(2); ++k)
        if (dh(i,j,k) != h(i,j,k)) ++ndiff;
  return ndiff;
}

constexpr int npoints = 37;
constexpr int ncols   = 10;
constexpr int nlev    = 72;

} // anonymous namespace

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_rng_bfb", "[cosp]")
{
  const int n = 200, ndraw = 1500;

  // Seeds as cosp_c2f.F90 makes them (1 to 100001), plus edge cases
  std::vector<int> seeds(n);
  for (int i = 0; i < n; ++i) seeds[i] = 1 + static_cast<int>(uniform(0,100000));
  seeds[0] = 0; seeds[1] = 1; seeds[2] = 100001; seeds[3] = 2147483647; seeds[4] = -12345;

  lview_2d rnd_f("rnd_f", n, ndraw);
  cosp_rng_f(n, ndraw, seeds.data(), rnd_f.data());

  Kokkos::View<int*, DefaultDevice> seeds_d("seeds", n);
  Kokkos::deep_copy(seeds_d, Kokkos::View<int*, Kokkos::HostSpace>(seeds.data(), n));
  dview_2d rnd_d("rnd_d", n, ndraw);
  Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,n), KOKKOS_LAMBDA (const int i) {
    F::RngState s = seeds_d(i);
    for (int k = 0; k < ndraw; ++k) rnd_d(i,k) = F::get_rng(s);
  });
  auto rnd_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rnd_d);

  int ndiff = 0;
  double rmin = 1, rmax = 0;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < ndraw; ++k) {
      if (rnd_h(i,k) != rnd_f(i,k)) ++ndiff;
      rmin = std::min(rmin, rnd_h(i,k));
      rmax = std::max(rmax, rnd_h(i,k));
    }
  }
  REQUIRE(ndiff == 0);
  REQUIRE(rmin >= 0.0);
  REQUIRE(rmax < 1.0);
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_simulator_optics_bfb", "[cosp]")
{
  lview_3d flag("flag", npoints, ncols, nlev), out_f("out_f", npoints, ncols, nlev);
  lview_2d var1("var1", npoints, nlev), var2("var2", npoints, nlev);
  auto* fp = flag.data();
  for (size_t i = 0; i < flag.size(); ++i) fp[i] = static_cast<int>(uniform(0,3)); // 0, 1, 2
  fill_random(var1, 0, 50);
  fill_random(var2, 0, 50);

  cosp_simulator_optics_f(npoints, ncols, nlev, flag.data(), var1.data(), var2.data(), out_f.data());

  auto flag_d = to_device<dview_3d>(flag, "flag");
  dview_2d var1_d("var1", npoints, nlev), var2_d("var2", npoints, nlev);
  {
    auto h1 = Kokkos::create_mirror_view(var1_d); Kokkos::deep_copy(h1, var1); Kokkos::deep_copy(var1_d, h1);
    auto h2 = Kokkos::create_mirror_view(var2_d); Kokkos::deep_copy(h2, var2); Kokkos::deep_copy(var2_d, h2);
  }
  dview_3d out_d("out_d", npoints, ncols, nlev);
  Kokkos::parallel_for(Policy3({0,0,0},{npoints,ncols,nlev}),
                       KOKKOS_LAMBDA (const int i, const int j, const int k) {
    out_d(i,j,k) = F::simulator_optics(flag_d(i,j,k), var1_d(i,k), var2_d(i,k));
  });

  REQUIRE(count_diffs(out_d, out_f) == 0);
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_modis_optics_partition_bfb", "[cosp]")
{
  lview_3d cw("cw", npoints, ncols, nlev), ci("ci", npoints, ncols, nlev),
           ws("ws", npoints, ncols, nlev), is("is", npoints, ncols, nlev),
           tau("tau", npoints, ncols, nlev),
           tl_f("tl_f", npoints, ncols, nlev), ti_f("ti_f", npoints, ncols, nlev);
  fill_random(cw,  0, 2e-3, 0.3);
  fill_random(ci,  0, 1e-3, 0.3);
  fill_random(ws,  2e-6, 40e-6);
  fill_random(is,  5e-6, 120e-6);
  fill_random(tau, 0, 80, 0.2);

  cosp_modis_optics_partition_f(npoints, ncols, nlev, cw.data(), ci.data(), ws.data(), is.data(),
                                tau.data(), tl_f.data(), ti_f.data());

  auto cw_d  = to_device<dview_3d>(cw,  "cw");
  auto ci_d  = to_device<dview_3d>(ci,  "ci");
  auto ws_d  = to_device<dview_3d>(ws,  "ws");
  auto is_d  = to_device<dview_3d>(is,  "is");
  auto tau_d = to_device<dview_3d>(tau, "tau");
  dview_3d tl_d("tl_d", npoints, ncols, nlev), ti_d("ti_d", npoints, ncols, nlev);
  Kokkos::parallel_for(Policy3({0,0,0},{npoints,ncols,nlev}),
                       KOKKOS_LAMBDA (const int i, const int j, const int k) {
    F::modis_optics_partition(cw_d(i,j,k), ci_d(i,j,k), ws_d(i,j,k), is_d(i,j,k), tau_d(i,j,k),
                              tl_d(i,j,k), ti_d(i,j,k));
  });

  REQUIRE(count_diffs(tl_d, tl_f) == 0);
  REQUIRE(count_diffs(ti_d, ti_f) == 0);
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_nir_fits_bfb", "[cosp]")
{
  // Sweep sizes across all the branch points (min/max clamps, and 7 microns for liquid)
  const int n = 2001;
  lview_1d re("re", n);
  for (int i = 0; i < n; ++i) re(i) = 120.0*i/(n-1);
  re(1) = 4.0; re(2) = 5.0; re(3) = 7.0; re(4) = 30.0; re(5) = 90.0;

  for (const int phase : {F::C::phaseIsLiquid, F::C::phaseIsIce}) {
    lview_1d g_f("g_f", n), w0_f("w0_f", n);
    cosp_nir_fits_f(n, phase, re.data(), g_f.data(), w0_f.data());

    dview_1d re_d("re_d", n), g_d("g_d", n), w0_d("w0_d", n);
    auto re_dh = Kokkos::create_mirror_view(re_d);
    Kokkos::deep_copy(re_dh, re);
    Kokkos::deep_copy(re_d, re_dh);
    Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,n), KOKKOS_LAMBDA (const int i) {
      g_d(i)  = F::get_g_nir  (phase, re_d(i));
      w0_d(i) = F::get_ssa_nir(phase, re_d(i));
    });
    auto g_h  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), g_d);
    auto w0_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), w0_d);

    int ndiff = 0;
    for (int i = 0; i < n; ++i) {
      if (g_h(i)  != g_f(i))  ++ndiff;
      if (w0_h(i) != w0_f(i)) ++ndiff;
    }
    REQUIRE(ndiff == 0);
  }
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_modis_optics_bfb", "[cosp]")
{
  lview_3d tl("tl", npoints, ncols, nlev), sl("sl", npoints, ncols, nlev),
           ti("ti", npoints, ncols, nlev), si("si", npoints, ncols, nlev),
           fl_f("fl_f", npoints, ncols, nlev), g_f("g_f", npoints, ncols, nlev),
           w0_f("w0_f", npoints, ncols, nlev);
  fill_random(tl, 0, 40, 0.4);
  fill_random(ti, 0, 20, 0.4);
  fill_random(sl, 0, 40, 0.2);   // microns; zero where there is no liquid
  fill_random(si, 0, 120, 0.2);

  cosp_modis_optics_f(npoints, ncols, nlev, tl.data(), sl.data(), ti.data(), si.data(),
                      fl_f.data(), g_f.data(), w0_f.data());

  auto tl_d = to_device<dview_3d>(tl, "tl");
  auto sl_d = to_device<dview_3d>(sl, "sl");
  auto ti_d = to_device<dview_3d>(ti, "ti");
  auto si_d = to_device<dview_3d>(si, "si");
  dview_3d fl_d("fl_d", npoints, ncols, nlev), g_d("g_d", npoints, ncols, nlev),
           w0_d("w0_d", npoints, ncols, nlev);
  Kokkos::parallel_for(Policy3({0,0,0},{npoints,ncols,nlev}),
                       KOKKOS_LAMBDA (const int i, const int j, const int k) {
    F::modis_optics(tl_d(i,j,k), sl_d(i,j,k), ti_d(i,j,k), si_d(i,j,k),
                    fl_d(i,j,k), g_d(i,j,k), w0_d(i,j,k));
  });

  REQUIRE(count_diffs(fl_d, fl_f) == 0);
  REQUIRE(count_diffs(g_d,  g_f)  == 0);
  REQUIRE(count_diffs(w0_d, w0_f) == 0);
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_hist2d_bfb", "[cosp]")
{
  // MODIS-like tau and CTP bins; some samples sit exactly on bin edges, some outside
  const int nbin1 = 7, nbin2 = 7, npts = 10, ntrials = 500;
  const std::vector<double> b1 = {0.3, 1.3, 3.6, 9.4, 23.0, 60.0, 100000.0, 1e7};
  const std::vector<double> b2 = {100000.0, 80000.0, 68000.0, 56000.0, 44000.0, 31000.0, 18000.0, 0.0};
  // Edges must increase for hist2D to count anything, so reverse the pressure bins
  std::vector<double> b2r(b2.rbegin(), b2.rend());

  lview_1d bin1("bin1", nbin1+1), bin2("bin2", nbin2+1);
  for (int i = 0; i <= nbin1; ++i) { bin1(i) = b1[i]; bin2(i) = b2r[i]; }

  lview_2d v1("v1", npts, ntrials), v2("v2", npts, ntrials);
  for (int t = 0; t < ntrials; ++t) {
    for (int n = 0; n < npts; ++n) {
      const double u = uniform(0,1);
      v1(n,t) = u < 0.15 ? b1[static_cast<int>(uniform(0,nbin1+1))] : uniform(-1, 150);
      v2(n,t) = u > 0.85 ? b2r[static_cast<int>(uniform(0,nbin2+1))] : uniform(-1000, 105000);
    }
  }

  lview_3d hist_f("hist_f", nbin1, nbin2, ntrials);
  lview_2d h("h", nbin1, nbin2);
  for (int t = 0; t < ntrials; ++t) {
    cosp_hist2d_f(npts, &v1(0,t), &v2(0,t), nbin1, bin1.data(), nbin2, bin2.data(), h.data());
    for (int i = 0; i < nbin1; ++i)
      for (int k = 0; k < nbin2; ++k)
        hist_f(i,k,t) = h(i,k);
  }

  dview_2d v1_d("v1_d", npts, ntrials), v2_d("v2_d", npts, ntrials);
  dview_1d bin1_d("bin1_d", nbin1+1), bin2_d("bin2_d", nbin2+1);
  {
    auto h = Kokkos::create_mirror_view(v1_d); Kokkos::deep_copy(h, v1); Kokkos::deep_copy(v1_d, h);
    auto g = Kokkos::create_mirror_view(v2_d); Kokkos::deep_copy(g, v2); Kokkos::deep_copy(v2_d, g);
    auto b = Kokkos::create_mirror_view(bin1_d); Kokkos::deep_copy(b, bin1); Kokkos::deep_copy(bin1_d, b);
    auto c = Kokkos::create_mirror_view(bin2_d); Kokkos::deep_copy(c, bin2); Kokkos::deep_copy(bin2_d, c);
  }
  dview_3d hist_d("hist_d", nbin1, nbin2, ntrials);
  Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,ntrials), KOKKOS_LAMBDA (const int t) {
    auto var1 = Kokkos::subview(v1_d, Kokkos::ALL, t);
    auto var2 = Kokkos::subview(v2_d, Kokkos::ALL, t);
    auto hist = Kokkos::subview(hist_d, Kokkos::ALL, Kokkos::ALL, t);
    F::hist2d(npts, var1, var2, nbin1, bin1_d, nbin2, bin2_d, hist);
  });

  REQUIRE(count_diffs(hist_d, hist_f) == 0);

  // Make sure the test actually exercises the histogram
  double total = 0;
  for (size_t i = 0; i < hist_f.size(); ++i) total += hist_f.data()[i];
  REQUIRE(total > 0);
}
