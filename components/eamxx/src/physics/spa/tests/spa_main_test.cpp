#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "physics/spa/spa_functions.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_pack.hpp"

#include <random>

namespace {

using namespace scream;
using namespace spa;

template <typename S>
using view_1d = typename KokkosTypes<DefaultDevice>::template view_1d<S>;
template <typename S>
using view_2d = typename KokkosTypes<DefaultDevice>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<DefaultDevice>::template view_3d<S>;

using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;
using Spack = SPAFunc::Spack;

using RPDF = std::uniform_real_distribution<Real>;
using IPDF = std::uniform_int_distribution<int>;

// Helper Functions
std::pair<Real,Real> compute_min_max(const view_1d<Real>::HostMirror& input, const int start, const int end);

template<typename Engine, typename PDF>
void randomize (SPAFunc::SPAData& data, Engine& engine, PDF&& pdf);

struct SPADataHost {
  using view_2d_host = typename view_2d<Real>::HostMirror;
  using view_3d_host = typename view_3d<Real>::HostMirror;

  view_2d_host   ccn3;
  view_3d_host   aer_g_sw;
  view_3d_host   aer_ssa_sw;
  view_3d_host   aer_tau_sw;
  view_3d_host   aer_tau_lw;

  SPADataHost (const SPAFunc::SPAData& spa_out) {
    ccn3       = Kokkos::create_mirror_view (ekat::scalarize(spa_out.CCN3));
    aer_g_sw   = Kokkos::create_mirror_view (ekat::scalarize(spa_out.AER_G_SW));
    aer_ssa_sw = Kokkos::create_mirror_view (ekat::scalarize(spa_out.AER_SSA_SW));
    aer_tau_sw = Kokkos::create_mirror_view (ekat::scalarize(spa_out.AER_TAU_SW));
    aer_tau_lw = Kokkos::create_mirror_view (ekat::scalarize(spa_out.AER_TAU_LW));
  }

  void copy_from_dev (const SPAFunc::SPAData& spa_out) {
    Kokkos::deep_copy (ccn3, ekat::scalarize(spa_out.CCN3));
    Kokkos::deep_copy (aer_g_sw, ekat::scalarize(spa_out.AER_G_SW));
    Kokkos::deep_copy (aer_ssa_sw, ekat::scalarize(spa_out.AER_SSA_SW));
    Kokkos::deep_copy (aer_tau_sw, ekat::scalarize(spa_out.AER_TAU_SW));
    Kokkos::deep_copy (aer_tau_lw, ekat::scalarize(spa_out.AER_TAU_LW));
  }
};

TEST_CASE("spa_main")
{
  auto engine = setup_random_test ();

 /* ====================================================================
  *                  Test Setup, create structures, etc.
  * ====================================================================*/

  int ncols    = IPDF(1,10)(engine);
  int nlevs    = 3;//IPDF(30,70)(engine);
  int nswbands = IPDF(10,20)(engine);
  int nlwbands = IPDF(10,20)(engine);

  std::cout << " +---------------------------+\n";
  std::cout << " |  SPA Interpolation sizes  |\n";
  std::cout << " +---------------------------+\n";
  std::cout << "    ncols   : " << ncols << "\n";
  std::cout << "    nlevs   : " << nlevs << "\n";
  std::cout << "    nswbands: " << nswbands << "\n";
  std::cout << "    nlwbands: " << nlwbands << "\n";
  std::cout << " +---------------------------+\n\n";

  // Set up the set of SPA structures needed to run the test
  // Recall: SPAInput has padded values (for stable vertical interpolation)
  SPAFunc::SPATimeState spa_time_state;
  SPAFunc::SPAInput     spa_beg(ncols, nlevs+2, nswbands, nlwbands);
  SPAFunc::SPAInput     spa_end(ncols, nlevs+2, nswbands, nlwbands);
  SPAFunc::SPAInput     spa_tmp(ncols, nlevs+2, nswbands, nlwbands);
  SPAFunc::SPAOutput    spa_out(ncols, nlevs,   nswbands, nlwbands);

  // Create local host copies of all relevant data:
  SPADataHost data_beg_h(spa_beg.data);
  SPADataHost data_end_h(spa_end.data);
  SPADataHost data_tmp_h(spa_tmp.data);
  SPADataHost data_out_h(spa_out);
  auto ps_beg_h = Kokkos::create_mirror_view(spa_beg.PS);
  auto ps_end_h = Kokkos::create_mirror_view(spa_end.PS);
  auto ps_tmp_h = Kokkos::create_mirror_view(spa_tmp.PS);

  // First initialize the time_state: set for January run
  util::TimeStamp t_beg(1900,1,1,0,0,0);
  util::TimeStamp t_end(1900,2,1,0,0,0);
  spa_time_state.current_month = t_beg.get_month();
  spa_time_state.t_beg_month = t_beg.frac_of_year_in_days();
  spa_time_state.days_this_month = util::days_in_month(t_beg.get_year(),t_beg.get_month());
  spa_time_state.inited = true;

  // Generate random beg/end data
  randomize(spa_beg.data,engine,RPDF(1.0,10.0));
  randomize(spa_end.data,engine,RPDF(1.0,10.0));
  ekat::genRandArray(spa_beg.PS,engine,RPDF(9e4,1.1e5));
  ekat::genRandArray(spa_end.PS,engine,RPDF(9e4,1.1e5));

  // These won't change between tests, so deep copy to host right away.
  data_beg_h.copy_from_dev(spa_beg.data);
  data_end_h.copy_from_dev(spa_end.data);
  Kokkos::deep_copy(ps_beg_h,spa_beg.PS);
  Kokkos::deep_copy(ps_end_h,spa_end.PS);

  // ======================================================== //
  //                Test elemental interpolation              //
  // ======================================================== //

  std::cout << "  -> elemental linear interp\n";
  Spack x0, x1;
  ekat::genRandArray(&x0,1,engine,RPDF(-10.0,10.0));
  ekat::genRandArray(&x1,1,engine,RPDF(-10.0,10.0));

  REQUIRE ( (SPAFunc::linear_interp(x0,x1,0)==x0).all() );
  REQUIRE ( (SPAFunc::linear_interp(x0,x1,1.0)==x1).all() );

  std::cout << "  -> elemental linear interp .............................. OK!\n";

  // ======================================================== //
  //                Test time interpolation                   //
  // ======================================================== //

  // 1. If t=t_beg, we should get back spa_beg
  std::cout << "  -> time interp, t=t_beg\n";

  spa_time_state.t_now = t_beg.frac_of_year_in_days();
  SPAFunc::perform_time_interpolation(spa_time_state,spa_beg,spa_end,spa_tmp);
  data_tmp_h.copy_from_dev(spa_tmp.data);
  Kokkos::deep_copy(ps_tmp_h,spa_tmp.PS);

  for (int i=0; i<ncols; ++i) {
    REQUIRE (ps_tmp_h(i) == ps_beg_h(i));
    for (int k=0; k<nlevs+2; ++k) {
      REQUIRE (data_tmp_h.ccn3(i,k) == data_beg_h.ccn3(i,k));
      for (int n=0; n<nswbands; ++n) {
        REQUIRE (data_tmp_h.aer_g_sw(i,n,k) == data_beg_h.aer_g_sw(i,n,k));
        REQUIRE (data_tmp_h.aer_ssa_sw(i,n,k) == data_beg_h.aer_ssa_sw(i,n,k));
        REQUIRE (data_tmp_h.aer_tau_sw(i,n,k) == data_beg_h.aer_tau_sw(i,n,k));
      }
      for (int n=0; n<nlwbands; ++n) {
        REQUIRE (data_tmp_h.aer_tau_lw(i,n,k) == data_beg_h.aer_tau_lw(i,n,k));
      }
    }
  }
  std::cout << "  -> time interp, t=t_beg ................................. OK!\n";

  // 2. If t=t_end, we should get back spa_end
  std::cout << "  -> time interp, t=t_end\n";

  spa_time_state.t_now = t_end.frac_of_year_in_days();
  SPAFunc::perform_time_interpolation(spa_time_state,spa_beg,spa_end,spa_tmp);
  data_tmp_h.copy_from_dev(spa_tmp.data);
  Kokkos::deep_copy(ps_tmp_h,spa_tmp.PS);

  for (int i=0; i<ncols; ++i) {
    REQUIRE (ps_tmp_h(i) == ps_end_h(i));
    for (int k=0; k<nlevs+2; ++k) {
      REQUIRE (data_tmp_h.ccn3(i,k) == data_end_h.ccn3(i,k));
      for (int n=0; n<nswbands; ++n) {
        REQUIRE (data_tmp_h.aer_g_sw(i,n,k) == data_end_h.aer_g_sw(i,n,k));
        REQUIRE (data_tmp_h.aer_ssa_sw(i,n,k) == data_end_h.aer_ssa_sw(i,n,k));
        REQUIRE (data_tmp_h.aer_tau_sw(i,n,k) == data_end_h.aer_tau_sw(i,n,k));
      }
      for (int n=0; n<nlwbands; ++n) {
        REQUIRE (data_tmp_h.aer_tau_lw(i,n,k) == data_end_h.aer_tau_lw(i,n,k));
      }
    }
  }
  std::cout << "  -> time interp, t=t_end ................................. OK!\n";

  // 3. Since time interpolation is convex, the output is bounded by the input
  std::cout << "  -> time interp, t_beg<t<t_end\n";

  RPDF pdf(0.001,0.999);
  auto dt = (t_end-t_beg) * pdf(engine);

  auto t_mid = t_beg + static_cast<int>(round(dt));
  spa_time_state.t_now = t_mid.frac_of_year_in_days();
  SPAFunc::perform_time_interpolation(spa_time_state,spa_beg,spa_end,spa_tmp);
  data_tmp_h.copy_from_dev(spa_tmp.data);
  Kokkos::deep_copy(ps_tmp_h,spa_tmp.PS);

  auto max_min = [&] (const Real v0, const Real v1) -> std::pair<Real,Real> {
    using pair_t = std::pair<Real,Real>;
    return (v0>v1 ? pair_t{v0,v1} : pair_t{v1,v0});
  };
  auto check_time_bounds = [&] (const Real v, const std::pair<Real,Real>& bounds) {
    REQUIRE ( v <= bounds.first  );
    REQUIRE ( v >= bounds.second );
  };
  for (int i=0; i<ncols; ++i) {
    REQUIRE (ps_tmp_h(i) >= std::min(ps_beg_h(i),ps_end_h(i)));
    REQUIRE (ps_tmp_h(i) <= std::max(ps_beg_h(i),ps_end_h(i)));
    for (int k=0; k<nlevs+2; ++k) {
      auto ccn3_bounds = max_min (data_end_h.ccn3(i,k),data_beg_h.ccn3(i,k));
      check_time_bounds (data_tmp_h.ccn3(i,k),ccn3_bounds);

      for (int n=0; n<nswbands; ++n) {
        auto aer_g_sw_bounds = max_min(data_end_h.aer_g_sw(i,n,k),
                                       data_beg_h.aer_g_sw(i,n,k));
        check_time_bounds (data_tmp_h.aer_g_sw(i,n,k), aer_g_sw_bounds);

        auto aer_ssa_sw_bounds = max_min(data_end_h.aer_ssa_sw(i,n,k),
                                         data_beg_h.aer_ssa_sw(i,n,k));
        check_time_bounds (data_tmp_h.aer_ssa_sw(i,n,k), aer_ssa_sw_bounds);

        auto aer_tau_sw_bounds = max_min(data_end_h.aer_tau_sw(i,n,k),
                                         data_beg_h.aer_tau_sw(i,n,k));
        check_time_bounds (data_tmp_h.aer_tau_sw(i,n,k), aer_tau_sw_bounds);
      }
      for (int n=0; n<nlwbands; ++n) {
        auto aer_tau_lw_bounds = max_min(data_end_h.aer_tau_lw(i,n,k),
                                         data_beg_h.aer_tau_lw(i,n,k));
        check_time_bounds (data_tmp_h.aer_tau_lw(i,n,k), aer_tau_lw_bounds);
      }
    }
  }
  std::cout << "  -> time interp, t_beg<t<t_end ........................... OK!\n";

  using col_type = view_1d<Real>::HostMirror;

  // ======================================================== //
  //                Test vert interpolation                   //
  // ======================================================== //

  // Create a target pressure profile to interpolate onto.
  // NOTE: set a slightly higher surface pressure than the p_src,
  //       to ensure spa vert interpolation extrapolates correctly.
  auto npacks_src = ekat::PackInfo<Spack::n>::num_packs(nlevs+2);
  auto npacks_tgt = ekat::PackInfo<Spack::n>::num_packs(nlevs);
  auto p_src = view_2d<Spack>("",ncols,npacks_src);
  auto p_tgt = view_2d<Spack>("",ncols,npacks_tgt);
  auto p_tgt_h = Kokkos::create_mirror_view(ekat::scalarize(p_tgt));
  auto p_src_h = Kokkos::create_mirror_view(ekat::scalarize(p_src));

  // Make a monotonic tgt pressure profile
  ekat::genRandArray(p_tgt,engine,RPDF(1e3,1e5));
  Kokkos::deep_copy(p_tgt_h, ekat::scalarize(p_tgt));
  for (int i=0; i<ncols; ++i) {
    auto col = ekat::subview(p_tgt_h,i);
    std::sort(col.data(),col.data()+nlevs);
    // For safety, set padding back to NaN's, in case we accidentally use it
    for (int k=nlevs; k<npacks_tgt*Spack::n; ++k) {
      col(k) = ekat::ScalarTraits<Real>::invalid();
    }
  }
  Kokkos::deep_copy(ekat::scalarize(p_tgt),p_tgt_h);

  // 1. If p_src==p_tgt, we should get back the input data
  std::cout << "  -> vert interp, p_tgt = p_src\n";

  for (int i=0; i<ncols; ++i) {
    p_src_h(i,0) = 0;
    for (int k=0; k<nlevs; ++k) {
      if (k<nlevs/2) {
        p_src_h(i,k+1) = p_tgt_h(i,k);
      } else {
        p_src_h(i,k+1) = p_tgt_h(i,k);
      }
    }
    p_src_h(i,nlevs+1) = 10*p_src_h(i,nlevs);
  }
  Kokkos::deep_copy(ekat::scalarize(p_src),p_src_h);

  // Run vertical interpolation
  SPAFunc::perform_vertical_interpolation(p_src,p_tgt,spa_beg.data,spa_out);
  data_out_h.copy_from_dev(spa_out);

  for (int i=0; i<ncols; ++i) {
    for (int k=0; k<nlevs; ++k) {
      REQUIRE (data_out_h.ccn3(i,k) == data_beg_h.ccn3(i,k+1));
      for (int n=0; n<nswbands; ++n) {
        REQUIRE (data_out_h.aer_g_sw(i,n,k) == data_beg_h.aer_g_sw(i,n,k+1));
        REQUIRE (data_out_h.aer_ssa_sw(i,n,k) == data_beg_h.aer_ssa_sw(i,n,k+1));
        REQUIRE (data_out_h.aer_tau_sw(i,n,k) == data_beg_h.aer_tau_sw(i,n,k+1));
      }
      for (int n=0; n<nlwbands; ++n) {
        REQUIRE (data_out_h.aer_tau_lw(i,n,k) == data_beg_h.aer_tau_lw(i,n,k+1));
      }
    }
  }
  std::cout << "  -> vert interp, p_tgt = p_src ........................... OK!\n";

  // 2. If p_src!=p_tgt, for each column verify that each level of output
  //    is bounded by max/min of input over column.
  //    Note: make p_tgt<p_src toward TOM, and p_tgt>p_src toward surface
  std::cout << "  -> vert interp, p_tgt!=p_src and extrapolation needed\n";

  for (int i=0; i<ncols; ++i) {
    p_src_h(i,0) = 0;
    for (int k=1; k<=nlevs; ++k) {
      if (k<nlevs/2) {
        p_src_h(i,k) = 0.9*p_tgt_h(i,k-1);
      } else {
        p_src_h(i,k) = 1.1*p_tgt_h(i,k-1);
      }
    }
    p_src_h(i,nlevs+1) = 10*p_src_h(i,nlevs);
  }
  Kokkos::deep_copy(ekat::scalarize(p_src),p_src_h);

  // Run vertical interpolation
  SPAFunc::perform_vertical_interpolation(p_src,p_tgt,spa_beg.data,spa_out);
  data_out_h.copy_from_dev(spa_out);

  // Helper lambda: verify max/min of output bounded by max/min of data
  auto check_bounds = [&](const col_type& data,const col_type& output) {
    auto data_minmax   = compute_min_max(data,0,nlevs+2); // Recall: data is padded
    auto output_minmax = compute_min_max(output,0,nlevs);
    REQUIRE( data_minmax.first  <= output_minmax.first );
    REQUIRE( data_minmax.second >= output_minmax.second );
  };
  auto sv = [&](const view_3d<Real>::HostMirror& v, const int col, const int band) -> col_type {
    return ekat::subview(v,col,band);
  };

  for (int i=0; i<ncols; ++i) {
    const auto& ccn3_in  = ekat::subview(data_beg_h.ccn3,i);
    const auto& ccn3_out = ekat::subview(data_out_h.ccn3,i);
    check_bounds( ccn3_in, ccn3_out );
    for (int n=0;n<nswbands;n++) {
      check_bounds (sv(data_beg_h.aer_g_sw,i,n),  sv(data_out_h.aer_g_sw,i,n));
      check_bounds (sv(data_beg_h.aer_ssa_sw,i,n),sv(data_out_h.aer_ssa_sw,i,n));
      check_bounds (sv(data_beg_h.aer_tau_sw,i,n),sv(data_out_h.aer_tau_sw,i,n));
    }
    for (int n=0;n<nlwbands;n++) {
      check_bounds (sv(data_beg_h.aer_tau_lw,i,n),sv(data_out_h.aer_tau_lw,i,n));
    }
  }
  std::cout << "  -> vert interp, p_tgt!=p_src and extrapolation needed ... OK!\n\n";
}

// Compute min/max of input over [start,end) indices
std::pair<Real,Real> compute_min_max(const view_1d<Real>::HostMirror& input, const int start, const int end)
{
  std::pair<Real,Real> result;
  // Initialize min and max with first entry in input
  result.first = input[start];
  result.second = input[start];
  // Now go through the rest of the view
  for (int k = start+1; k<end; k++) {
    result.first  = std::min(result.first,input[k]);
    result.second = std::max(result.second,input[k]);
  }
  return result;
}

template<typename Engine,typename PDF>
void randomize (SPAFunc::SPAData& data, Engine& engine, PDF&& pdf) {

  // Need to make sure padding is set correctly in each column:
  //   col(0) = 0
  //   col(end) = col(end-1)
  auto ccn3 = Kokkos::create_mirror_view(ekat::scalarize(data.CCN3));
  auto aer_g_sw = Kokkos::create_mirror_view(ekat::scalarize(data.AER_G_SW));
  auto aer_ssa_sw = Kokkos::create_mirror_view(ekat::scalarize(data.AER_SSA_SW));
  auto aer_tau_sw = Kokkos::create_mirror_view(ekat::scalarize(data.AER_TAU_SW));
  auto aer_tau_lw = Kokkos::create_mirror_view(ekat::scalarize(data.AER_TAU_LW));

  // NOTE: the deep_copy calls in 'genRandArray' are no-ops for host views
  ekat::genRandArray(ccn3,engine,pdf);
  ekat::genRandArray(aer_g_sw,engine,pdf);
  ekat::genRandArray(aer_ssa_sw,engine,pdf);
  ekat::genRandArray(aer_tau_sw,engine,pdf);
  ekat::genRandArray(aer_tau_lw,engine,pdf);

  const int nlevs = data.nlevs-2;
  for (int i=0; i<data.ncols; ++i) {
    ccn3(i,0) = 0;
    ccn3(i,nlevs+1) = ccn3(i,nlevs);
    for (int n=0; n<data.nswbands; ++n) {
      aer_g_sw(i,n,0) = 0;
      aer_g_sw(i,n,nlevs+1) = aer_g_sw(i,n,nlevs);
      aer_ssa_sw(i,n,0) = 0;
      aer_ssa_sw(i,n,nlevs+1) = aer_ssa_sw(i,n,nlevs);
      aer_tau_sw(i,n,0) = 0;
      aer_tau_sw(i,n,nlevs+1) = aer_tau_sw(i,n,nlevs);
    }
    for (int n=0; n<data.nlwbands; ++n) {
      aer_tau_lw(i,n,0) = 0;
      aer_tau_lw(i,n,nlevs+1) = aer_tau_lw(i,n,nlevs);
    }
  }

  Kokkos::deep_copy(ekat::scalarize(data.CCN3),ccn3);
  Kokkos::deep_copy(ekat::scalarize(data.AER_G_SW),aer_g_sw);
  Kokkos::deep_copy(ekat::scalarize(data.AER_SSA_SW),aer_ssa_sw);
  Kokkos::deep_copy(ekat::scalarize(data.AER_TAU_SW),aer_tau_sw);
  Kokkos::deep_copy(ekat::scalarize(data.AER_TAU_LW),aer_tau_lw);
}

} // namespace

