#include <catch2/catch.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>
#include <string>

#include "share/core/eamxx_types.hpp"

#include <mam4xx/mam4.hpp>
#include <yaml-cpp/yaml.h>

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

namespace scream {
namespace impl {

mam4::mo_photo::PhotoTableData read_photo_table(
    const std::string& rsf_file, const std::string& xs_long_file);

}  // namespace impl
}  // namespace scream

namespace {

using Real = scream::Real;
using HostSpace = Kokkos::HostSpace;
using HostView1D = mam4::DeviceType::view_1d<Real>::host_mirror_type;
using HostView5D = mam4::DeviceType::view<Real*****>::host_mirror_type;

using Device     = scream::DefaultDevice;
using ExecSpace  = Device::execution_space;
using KT         = ekat::KokkosTypes<ExecSpace>;
using view_1d    = typename KT::template view_1d<Real>;
using view_2d    = typename KT::template view_2d<Real>;
using view_3d    = typename KT::template view_3d<Real>;
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using MemberType = TeamPolicy::member_type;

inline bool nearly_equal(const Real a, const Real b,
                         const Real rtol = 1e-8,
                         const Real atol = 1e-14) {
  return std::abs(a - b) <= atol + rtol * std::abs(b);
}

std::vector<Real> read_real_vector(const YAML::Node& node) {
  std::vector<Real> vals;
  vals.reserve(node.size());
  for (std::size_t i = 0; i < node.size(); ++i) {
    vals.push_back(node[i].as<Real>());
  }
  return vals;
}

std::vector<int> read_int_vector(const YAML::Node& node) {
  std::vector<int> vals;
  vals.reserve(node.size());
  for (std::size_t i = 0; i < node.size(); ++i) {
    vals.push_back(node[i].as<int>());
  }
  return vals;
}

}  // namespace

TEST_CASE("mam_photo_table_yaml_reference_regression",
          "[mam4][photo][kokkos]") {
  if constexpr (mam4::nlev != 72) return;
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);
  struct ScorpioGuard {
    explicit ScorpioGuard(const ekat::Comm& comm) : comm_(comm) {
      scorpio::init_subsystem(comm_);
    }
    ~ScorpioGuard() {
      scorpio::finalize_subsystem();
    }
    const ekat::Comm& comm_;
  } scorpio_guard(comm);

  const std::string rsf_file =
      std::string(SCREAM_DATA_DIR) + "/mam4xx/photolysis/RSF_GT200nm_v3.0_c080811.nc";
  const std::string xs_long_file =
      std::string(SCREAM_DATA_DIR) + "/mam4xx/photolysis/temp_prs_GT200nm_JPL10_c130206.nc";
  const std::string input_yaml_file = "jlong_input_ts_355.yaml";

  const auto photo_table = scream::impl::read_photo_table(rsf_file, xs_long_file);
  const YAML::Node root = YAML::LoadFile(input_yaml_file);
  REQUIRE(root["input"]);
  REQUIRE(root["input"]["fixed"]);
  const auto fixed = root["input"]["fixed"];

  REQUIRE(photo_table.nw > 0);
  REQUIRE(photo_table.numj == 1);

  auto sza_h       = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.sza);
  auto del_sza_h   = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.del_sza);
  auto alb_h       = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.alb);
  auto del_alb_h   = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.del_alb);
  auto colo3_h     = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.colo3);
  auto o3rat_h     = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.o3rat);
  auto del_o3rat_h = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.del_o3rat);
  auto press_h     = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.press);
  auto prs_h       = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.prs);
  auto dprs_h      = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.dprs);
  auto rsf_tab_h   = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.rsf_tab);
  auto xsqy_h      = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.xsqy);
  auto etfphot_h   = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.etfphot);
  auto lng_indexer_h = Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.lng_indexer);
  auto pht_alias_mult_h =
      Kokkos::create_mirror_view_and_copy(HostSpace(), photo_table.pht_alias_mult_1);

  const auto nw_ref       = read_int_vector(fixed["nw"])[0];
  const auto numj_ref     = read_int_vector(fixed["numj"])[0];
  const auto shape_ref    = read_int_vector(fixed["shape_of_rsf_tab"]);
  REQUIRE(shape_ref.size() == 5);
  const int nw_shape       = shape_ref[0];
  const int nump_shape     = shape_ref[1];
  const int numsza_shape   = shape_ref[2];
  const int numcolo3_shape = shape_ref[3];
  const int numalb_shape   = shape_ref[4];

  const auto sza_ref       = read_real_vector(fixed["sza"]);
  const auto del_sza_ref   = read_real_vector(fixed["del_sza"]);
  const auto alb_ref       = read_real_vector(fixed["alb"]);
  const auto del_alb_ref   = read_real_vector(fixed["del_alb"]);
  const auto colo3_ref     = read_real_vector(fixed["colo3"]);
  const auto o3rat_ref     = read_real_vector(fixed["o3rat"]);
  const auto del_o3rat_ref = read_real_vector(fixed["del_o3rat"]);
  const YAML::Node press_node = fixed["press"] ? fixed["press"] : fixed["pm"];
  REQUIRE(press_node);
  const auto press_ref = read_real_vector(press_node);
  const auto etfphot_ref   = read_real_vector(fixed["etfphot"]);
  const auto prs_ref       = read_real_vector(fixed["prs"]);
  const auto dprs_ref      = read_real_vector(fixed["dprs"]);
  const auto rsf_tab_2d    = read_real_vector(fixed["rsf_tab_2d"]);
  const auto xsqy_2d       = read_real_vector(fixed["xsqy_2d"]);

  SECTION("dimensions_match_expected_shapes") {
    REQUIRE(photo_table.nw == nw_ref);
    REQUIRE(photo_table.numj == numj_ref);
    REQUIRE(photo_table.nw == nw_shape);
    REQUIRE(photo_table.nump == nump_shape);
    REQUIRE(photo_table.numsza == numsza_shape);
    REQUIRE(photo_table.numcolo3 == numcolo3_shape);
    REQUIRE(photo_table.numalb == numalb_shape);

    REQUIRE(photo_table.sza.extent(0) == photo_table.numsza);
    REQUIRE(photo_table.alb.extent(0) == photo_table.numalb);
    REQUIRE(photo_table.colo3.extent(0) == photo_table.nump);
    REQUIRE(photo_table.o3rat.extent(0) == photo_table.numcolo3);
    REQUIRE(photo_table.prs.extent(0) == photo_table.np_xs);
    REQUIRE(photo_table.lng_indexer.extent(0) == mam4::mo_photo::phtcnt);
    REQUIRE(photo_table.pht_alias_mult_1.extent(0) == mam4::mo_photo::phtcnt);
  }

  SECTION("1d_tables_match_yaml_reference") {
    for (int i = 0; i < photo_table.numsza; ++i) {
      REQUIRE(nearly_equal(sza_h(i), sza_ref[i]));
    }
    for (int i = 0; i < photo_table.numalb; ++i) {
      REQUIRE(nearly_equal(alb_h(i), alb_ref[i]));
    }
    for (int i = 0; i < photo_table.nump; ++i) {
      REQUIRE(nearly_equal(colo3_h(i), colo3_ref[i]));
    }
    for (int i = 0; i < photo_table.numcolo3; ++i) {
      REQUIRE(nearly_equal(o3rat_h(i), o3rat_ref[i]));
    }
    for (int i = 0; i < photo_table.nump; ++i) {
      REQUIRE(nearly_equal(press_h(i), press_ref[i]));
    }
    for (int i = 0; i < photo_table.np_xs; ++i) {
      REQUIRE(nearly_equal(prs_h(i), prs_ref[i]));
    }
    for (int i = 0; i < photo_table.numsza - 1; ++i) {
      REQUIRE(nearly_equal(del_sza_h(i), del_sza_ref[i]));
    }
    for (int i = 0; i < photo_table.numalb - 1; ++i) {
      REQUIRE(nearly_equal(del_alb_h(i), del_alb_ref[i]));
    }
    for (int i = 0; i < photo_table.numcolo3 - 1; ++i) {
      REQUIRE(nearly_equal(del_o3rat_h(i), del_o3rat_ref[i]));
    }
    for (int i = 0; i < photo_table.np_xs - 1; ++i) {
      REQUIRE(nearly_equal(dprs_h(i), dprs_ref[i]));
    }
  }
  
  SECTION("rsf_corner_slice_matches_reference") {
    const int nw   = photo_table.nw;
    const int nump = photo_table.nump;
    int count = 0;
    for (int d2 = 0; d2 < nump; ++d2) {
      for (int d1 = 0; d1 < nw; ++d1) {
        const auto computed = rsf_tab_h(d1, d2, 0, 0, 0);;
        const auto expected = rsf_tab_2d[count];
        count++;
        INFO("rsf_tab mismatch at (i=" << d1 << ", j=" << d2
               << "), computed=" << computed << ", expected=" << expected);
        REQUIRE(nearly_equal(computed, expected,1e-6));
      }
    }
  }

  SECTION("xsqy_table_matches_yaml_reference") {
    REQUIRE(xsqy_2d.size() ==
            static_cast<std::size_t>(photo_table.numj) *
            static_cast<std::size_t>(photo_table.nw));
    int count = 0;
    for (int w = 0; w < photo_table.nw; ++w) {
      REQUIRE(nearly_equal(xsqy_h(0, w, 0, 0), xsqy_2d[count], 0, 0));
      ++count;
    }
  }

  SECTION("indexing_and_alias_arrays_are_initialized") {
    REQUIRE(lng_indexer_h(0) == 0);
    for (int i = 0; i < mam4::mo_photo::phtcnt; ++i) {
      REQUIRE(nearly_equal(pht_alias_mult_h(i), 1.0));
    }
  }

  SECTION("etfphot_is_finite") {
    REQUIRE(static_cast<int>(etfphot_ref.size()) == photo_table.nw);
    for (int i = 0; i < photo_table.nw; ++i) {
      REQUIRE(nearly_equal(etfphot_h(i), etfphot_ref[i]));
      REQUIRE(std::isfinite(etfphot_h(i)));
    }
  }

}

TEST_CASE("mam_photo_table_kernel_single_column_nlev72_regression",
          "[mam4][photo][kokkos]") {
  if constexpr (mam4::nlev != 72) return;
  constexpr int ncol = 1;
  constexpr int nlev = mam4::nlev;
  constexpr int nref = 1;
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);
  struct ScorpioGuard {
    explicit ScorpioGuard(const ekat::Comm& comm) : comm_(comm) {
      scorpio::init_subsystem(comm_);
    }
    ~ScorpioGuard() { scorpio::finalize_subsystem(); }
    const ekat::Comm& comm_;
  } scorpio_guard(comm);

  const std::string rsf_file =
      std::string(SCREAM_DATA_DIR) + "/mam4xx/photolysis/RSF_GT200nm_v3.0_c080811.nc";
  const std::string xs_long_file =
      std::string(SCREAM_DATA_DIR) + "/mam4xx/photolysis/temp_prs_GT200nm_JPL10_c130206.nc";
  const std::string input_yaml_file = "table_photo_input_ts_355.yaml";

  const YAML::Node root = YAML::LoadFile(input_yaml_file);
  REQUIRE(root["input"]);
  REQUIRE(root["input"]["fixed"]);
  const auto fixed = root["input"]["fixed"];

  const auto photo_table = scream::impl::read_photo_table(rsf_file, xs_long_file);
  const int work_len = mam4::mo_photo::get_photo_table_work_len(photo_table);
  const int npht     = mam4::mo_photo::phtcnt;

  REQUIRE(work_len > 0);
  REQUIRE(npht >= nref);

  // Allocate device views.
  view_2d work_photo_table("work_photo_table", ncol, work_len);
  view_2d pmid("pmid",     ncol, nlev);
  view_2d pdel("pdel",     ncol, nlev);
  view_2d temper("temper", ncol, nlev);
  view_2d o3col("o3col",   ncol, nlev);
  view_1d zen_angle("zen_angle", ncol);
  view_1d srf_alb("srf_alb",     ncol);
  view_2d qc("qc",   ncol, nlev);
  view_2d cld("cld", ncol, nlev);
  view_3d photo("photo", ncol, nlev, npht);

  Kokkos::deep_copy(work_photo_table, 0.0);
  Kokkos::deep_copy(photo, 0.0);

  // Read atmospheric-state reference data from YAML.
  const auto pmid_vals   = read_real_vector(fixed["pmid"]);
  const auto pdel_vals   = read_real_vector(fixed["pdel"]);
  const auto temper_vals = read_real_vector(fixed["temper"]);
  const auto o3col_vals  = read_real_vector(fixed["col_dens_1"]);
  const auto lwc_vals    = read_real_vector(fixed["lwc"]);
  const auto cloud_vals  = read_real_vector(fixed["clouds"]);
  const auto zen_vals    = read_real_vector(fixed["zen_angle"]);
  const auto alb_vals    = read_real_vector(fixed["srf_alb"]);
  const auto esfact_vals = read_real_vector(fixed["esfact"]);
  const auto photo_ref   = read_real_vector(fixed["photos"]);

  REQUIRE(pmid_vals.size()   >= static_cast<std::size_t>(nlev));
  REQUIRE(pdel_vals.size()   >= static_cast<std::size_t>(nlev));
  REQUIRE(temper_vals.size() >= static_cast<std::size_t>(nlev));
  REQUIRE(o3col_vals.size()  >= static_cast<std::size_t>(nlev));
  REQUIRE(lwc_vals.size()    >= static_cast<std::size_t>(nlev));
  REQUIRE(cloud_vals.size()  >= static_cast<std::size_t>(nlev));
  REQUIRE(zen_vals.size()    >= 1);
  REQUIRE(alb_vals.size()    >= 1);
  REQUIRE(esfact_vals.size() >= 1);

  const Real zen_val = zen_vals[0];
  const Real alb_val = alb_vals[0];
  const Real esfact  = esfact_vals[0];

  // Fill host mirrors and copy to device.
  auto pmid_h   = Kokkos::create_mirror_view(pmid);
  auto pdel_h   = Kokkos::create_mirror_view(pdel);
  auto temper_h = Kokkos::create_mirror_view(temper);
  auto o3col_h  = Kokkos::create_mirror_view(o3col);
  auto zen_h    = Kokkos::create_mirror_view(zen_angle);
  auto alb_h    = Kokkos::create_mirror_view(srf_alb);
  auto qc_h     = Kokkos::create_mirror_view(qc);
  auto cld_h    = Kokkos::create_mirror_view(cld);

  for (int k = 0; k < nlev; ++k) {
    pmid_h(0, k)   = pmid_vals[k];
    pdel_h(0, k)   = pdel_vals[k];
    temper_h(0, k) = temper_vals[k];
    o3col_h(0, k)  = o3col_vals[k];
    qc_h(0, k)     = lwc_vals[k];
    cld_h(0, k)    = cloud_vals[k];
  }
  zen_h(0) = zen_val;
  alb_h(0) = alb_val;

  Kokkos::deep_copy(pmid,      pmid_h);
  Kokkos::deep_copy(pdel,      pdel_h);
  Kokkos::deep_copy(temper,    temper_h);
  Kokkos::deep_copy(o3col,     o3col_h);
  Kokkos::deep_copy(zen_angle, zen_h);
  Kokkos::deep_copy(srf_alb,   alb_h);
  Kokkos::deep_copy(qc,        qc_h);
  Kokkos::deep_copy(cld,       cld_h);

  // Launch one-column photolysis kernel.
  TeamPolicy policy(ncol, Kokkos::AUTO());
  Kokkos::parallel_for(
    "unit_test_table_photo_nlev72", policy,
    KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();

      const auto work_icol = ekat::subview(work_photo_table, icol);
      mam4::mo_photo::PhotoTableWorkArrays photo_work_arrays;
      mam4::mo_photo::set_photo_table_work_arrays(photo_table, work_icol,
                                                  photo_work_arrays);
      team.team_barrier();

      mam4::mo_photo::table_photo(
        team,
        ekat::subview(photo,  icol),
        ekat::subview(pmid,   icol),
        ekat::subview(pdel,   icol),
        ekat::subview(temper, icol),
        ekat::subview(o3col,  icol),
        zen_angle(icol), srf_alb(icol),
        ekat::subview(qc,  icol),
        ekat::subview(cld, icol),
        esfact, photo_table, photo_work_arrays);
    });
  Kokkos::fence();

  auto photo_h = Kokkos::create_mirror_view_and_copy(HostSpace(), photo);

  SECTION("all_outputs_are_finite") {
    for (int k = 0; k < nlev; ++k) {
      for (int j = 0; j < nref; ++j) {
        INFO("Non-finite output at k=" << k << ", j=" << j
             << ", value=" << photo_h(0, k, j));
        REQUIRE(std::isfinite(photo_h(0, k, j)));
      }
    }
  }
  
  SECTION("compare_against_reference_when_available") {
  REQUIRE(photo_ref.size() == static_cast<std::size_t>(nlev * nref));

  int count = 0;
  for (int d2 = 0; d2 < nref; ++d2) {
    for (int d1 = 0; d1 < nlev; ++d1) {
      const auto computed = photo_h(0, d1, d2);
      const auto expected = photo_ref[count];
      count++;
      Real diff=computed - expected;
      Real rel = abs(diff)/expected;
      INFO("Reference mismatch at d1=" << d1 << ", d2=" << d2
              << ", computed=" << computed 
              << ", expected=" << expected);
      REQUIRE(nearly_equal(computed, expected, 1e-8, 1e-12));
    }
  }
  }
}
