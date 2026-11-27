#include <catch2/catch.hpp>

#include "physics/cosp/tests/cosp_test_utils.hpp"
#include "physics/cosp/cosp_functions.hpp"

#include <cmath>

/*
 * BFB test of the whole C++ COSP call (cosp_main) against the Fortran one EAMxx uses
 * (cosp_c2f_run). The Fortran can only be initialized once per process, so the number of
 * subcolumns is set at compile time, and each value gets its own test executable.
 */

#ifndef COSP_TEST_NSUBCOL
#define COSP_TEST_NSUBCOL 1
#endif

using namespace scream;
using namespace scream::cosp_test;

namespace {

// Inputs of one COSP call, in the Fortran (LayoutLeft) order
struct CospInputs {
  int ncol, nlev;
  lview_1d sunlit, skt;
  lview_2d T_mid, p_mid, p_int, z_mid, qv, qc, qi, cldfrac, reff_qc, reff_qi, dtau067, dtau105;

  CospInputs (const int ncol_, const int nlev_)
   : ncol(ncol_), nlev(nlev_)
   , sunlit("sunlit", ncol), skt("skt", ncol)
   , T_mid("T_mid", ncol, nlev), p_mid("p_mid", ncol, nlev), p_int("p_int", ncol, nlev+1)
   , z_mid("z_mid", ncol, nlev), qv("qv", ncol, nlev), qc("qc", ncol, nlev)
   , qi("qi", ncol, nlev), cldfrac("cldfrac", ncol, nlev), reff_qc("reff_qc", ncol, nlev)
   , reff_qi("reff_qi", ncol, nlev), dtau067("dtau067", ncol, nlev)
   , dtau105("dtau105", ncol, nlev)
  {
    make_thermo(p_int, p_mid, T_mid, qv, skt);
    for (int j = 0; j < ncol; ++j) {
      sunlit(j) = uniform(0,1) < 0.75 ? 1 : 0;
      const double psfc = p_int(j,nlev), zsfc = uniform(0,1) < 0.3 ? uniform(0,3000) : 0.0;
      const double cld_prob = uniform(0,1) < 0.1 ? 0.0 : uniform(0,0.5);
      const double scale    = std::pow(10.0, uniform(-2,1.5));
      for (int k = 0; k < nlev; ++k) {
        z_mid(j,k) = zsfc + 7000.0*std::log(psfc/p_mid(j,k));
        const bool cloudy = uniform(0,1) < cld_prob;
        cldfrac(j,k) = cloudy ? (uniform(0,1) < 0.3 ? 1.0 : uniform(0.01,1)) : 0.0;
        // Liquid only when warm enough, ice only when cold enough, mixed in between
        const double t = T_mid(j,k);
        const double liq = !cloudy || t < 240 ? 0.0 : uniform(1e-6,1e-3);
        const double ice = !cloudy || t > 270 ? 0.0 : uniform(1e-7,5e-4);
        qc(j,k) = liq;
        qi(j,k) = ice;
        reff_qc(j,k) = uniform(2,35);
        reff_qi(j,k) = uniform(3,120);
        dtau067(j,k) = cloudy && (liq > 0 || ice > 0) ? scale*uniform(0,1) : 0.0;
        dtau105(j,k) = dtau067(j,k)*uniform(0.3,0.7);
      }
    }
  }
};

struct CospOutputs {
  lview_1d isccp_cldtot;
  lview_3d isccp_ctptau, modis_ctptau, misr_cthtau;
  CospOutputs (const int ncol, const int ntau, const int nctp, const int ncth)
   : isccp_cldtot("isccp_cldtot", ncol), isccp_ctptau("isccp_ctptau", ncol, ntau, nctp)
   , modis_ctptau("modis_ctptau", ncol, ntau, nctp), misr_cthtau("misr_cthtau", ncol, ntau, ncth)
  {}
};

void run_fortran (const CospInputs& in, const double emsfc_lw, const CospOutputs& out) {
  using C = F::C;
  cosp_c2f_run(in.ncol, COSP_TEST_NSUBCOL, in.nlev, C::ntau, C::npres, C::nhgt, emsfc_lw, in.sunlit.data(),
               in.skt.data(), in.T_mid.data(), in.p_mid.data(), in.p_int.data(),
               in.z_mid.data(), in.qv.data(), in.qc.data(), in.qi.data(), in.cldfrac.data(),
               in.reff_qc.data(), in.reff_qi.data(), in.dtau067.data(), in.dtau105.data(),
               out.isccp_cldtot.data(), out.isccp_ctptau.data(), out.modis_ctptau.data(),
               out.misr_cthtau.data());
}

// Run the C++ and compare with the Fortran outputs
int run_cxx_and_compare (const CospInputs& in, const double emsfc_lw, const CospOutputs& out,
                         const F::Workspace& ws) {
  using C = F::C;
  const int ncol = in.ncol, nlev = in.nlev;
  diview_1d sunlit("sunlit", ncol);
  {
    auto h = Kokkos::create_mirror_view(sunlit);
    for (int j = 0; j < ncol; ++j) h(j) = static_cast<int>(in.sunlit(j));
    Kokkos::deep_copy(sunlit, h);
  }
  dview_1d cldtot("cldtot", ncol);
  dview_3d isccp("isccp", ncol, C::ntau, C::npres), modis("modis", ncol, C::ntau, C::npres),
           misr("misr", ncol, C::ntau, C::nhgt);
  F::cosp_main(ncol, COSP_TEST_NSUBCOL, nlev, emsfc_lw, sunlit, to_device<dview_1d>(in.skt, "skt"),
               to_device<dview_2d>(in.T_mid, "T_mid"), to_device<dview_2d>(in.p_mid, "p_mid"),
               to_device<dview_2d>(in.p_int, "p_int"), to_device<dview_2d>(in.z_mid, "z_mid"),
               to_device<dview_2d>(in.qv, "qv"), to_device<dview_2d>(in.qc, "qc"),
               to_device<dview_2d>(in.qi, "qi"), to_device<dview_2d>(in.cldfrac, "cldfrac"),
               to_device<dview_2d>(in.reff_qc, "reff_qc"),
               to_device<dview_2d>(in.reff_qi, "reff_qi"),
               to_device<dview_2d>(in.dtau067, "dtau067"),
               to_device<dview_2d>(in.dtau105, "dtau105"), cldtot, isccp, modis, misr, ws);
  const int n1 = count_diffs(cldtot, out.isccp_cldtot);
  const int n2 = count_diffs(isccp,  out.isccp_ctptau);
  const int n3 = count_diffs(modis,  out.modis_ctptau);
  const int n4 = count_diffs(misr,   out.misr_cthtau);
  INFO("diffs: isccp_cldtot=" << n1 << " isccp_ctptau=" << n2 << " modis_ctptau=" << n3
       << " misr_cthtau=" << n4);
  // Returns the number of diffs beyond what is allowed, see allowed_diffs
  const int excess = std::max(0, n1 - allowed_diffs(out.isccp_cldtot.size())) +
                     std::max(0, n2 - allowed_diffs(out.isccp_ctptau.size())) +
                     std::max(0, n3 - allowed_diffs(out.modis_ctptau.size())) +
                     std::max(0, n4 - allowed_diffs(out.misr_cthtau.size()));
  CHECK(excess == 0);
  return excess;
}

// Sum of the histogram entries over the sunlit columns, to make sure they are not empty
double sunlit_total (const CospInputs& in, const lview_3d& fq) {
  double s = 0;
  for (int j = 0; j < in.ncol; ++j) {
    if (in.sunlit(j) != 1) continue;
    for (size_t a = 0; a < fq.extent(1); ++a)
      for (size_t b = 0; b < fq.extent(2); ++b)
        s += fq(j,a,b);
  }
  return s;
}

} // anonymous namespace

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_main_bfb", "[cosp]")
{
  using C = F::C;
  const int ncol = 256, nlev = 72;
  const double emsfc_lw = 0.99;

  const int nsubcol = COSP_TEST_NSUBCOL;
  INFO("nsubcol=" << nsubcol);
  cosp_c2f_init(ncol, nsubcol, nlev);
  F::Workspace ws(ncol, nsubcol, nlev);

  // cosp_c2f_init can only run once (cosp_init allocates module arrays), so this is one
  // test case without sections
  {
    for (int rep = 0; rep < 3; ++rep) {
      CospInputs in(ncol, nlev);
      CospOutputs out(ncol, C::ntau, C::npres, C::nhgt);
      run_fortran(in, emsfc_lw, out);
      REQUIRE(run_cxx_and_compare(in, emsfc_lw, out, ws) == 0);
      CHECK(sunlit_total(in, out.isccp_ctptau) > 0);
      CHECK(sunlit_total(in, out.modis_ctptau) > 0);
      CHECK(sunlit_total(in, out.misr_cthtau)  > 0);
      if (nsubcol > 1) {
        // Partly cloudy columns give histogram entries other than 0 and 100
        int npartial = 0;
        for (int j = 0; j < ncol; ++j)
          for (int a = 0; a < C::ntau; ++a)
            for (int b = 0; b < C::npres; ++b) {
              const double v = out.modis_ctptau(j,a,b);
              if (in.sunlit(j) == 1 && v > 0 && v < 100) ++npartial;
            }
        CHECK(npartial > 0);
      }
    }
  }

  {
    // One bad value anywhere switches simulators off for all columns (cosp_errorCheck)
    for (const int which : {0, 1, 2, 3}) {
      CospInputs in(ncol, nlev);
      const int j = ncol/2, k = nlev/2;
      if (which == 0) in.qv(j,k) = -1e-12;       // ISCCP off
      if (which == 1) in.z_mid(j,nlev-1) = -400; // MISR off
      if (which == 2) in.dtau067(j,k) = -1;      // all off (cldfrac > 0 needed)
      if (which == 3) in.T_mid(j,k) = -1;        // ISCCP and MISR off
      if (which == 2) in.cldfrac(j,k) = 0.5;
      CospOutputs out(ncol, C::ntau, C::npres, C::nhgt);
      run_fortran(in, emsfc_lw, out);
      INFO("which=" << which);
      REQUIRE(run_cxx_and_compare(in, emsfc_lw, out, ws) == 0);

      // Check that the Fortran did switch the expected simulators off, everywhere
      int jsun = 0;
      while (in.sunlit(jsun) != 1) ++jsun;
      const bool isccp_off = which != 1, misr_off = which >= 1, modis_off = which == 2;
      CHECK((out.isccp_cldtot(jsun)      == C::R_UNDEF) == isccp_off);
      CHECK((out.misr_cthtau(jsun,0,0)   == C::R_UNDEF) == misr_off);
      CHECK((out.modis_ctptau(jsun,0,0)  == C::R_UNDEF) == modis_off);
    }
  }

  cosp_c2f_final();
}
