#include <iomanip>

#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"Point Grid"});
  auto &pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

//-----------------------------------------------------------------------------------------------//
template <typename DeviceT>
void run(std::mt19937_64 &engine) {
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using view_1d    = typename KT::template view_1d<Real>;

  constexpr int num_levs = 33;
  constexpr Real g       = PC::gravit;
  constexpr Real macheps = PC::macheps;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 10;
  auto gm         = create_gm(comm, ncols, num_levs);

  // Input (randomized) views
  view_1d pseudo_density("pseudo_density", num_levs), qc("qc", num_levs),
      nc("nc", num_levs), qr("qr", num_levs), nr("nr", num_levs),
      qi("qi", num_levs), ni("ni", num_levs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qx(0.0, 1e-3), pdf_pd(1.0, 100.0);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();
  ekat::ParameterList params;

  REQUIRE_THROWS(
      diag_factory.create("NumberPath", comm, params));  // No 'Number Kind'
  params.set<std::string>("Number Kind", "Foo");
  REQUIRE_THROWS(diag_factory.create("NumberPath", comm,
                                     params));  // Invalid 'Number Kind'

  // Liquid
  params.set<std::string>("Number Kind", "Liq");
  auto diag_liq = diag_factory.create("NumberPath", comm, params);
  diag_liq->set_grids(gm);
  diags.emplace("lnp", diag_liq);
  // Ice
  params.set<std::string>("Number Kind", "Ice");
  auto diag_ice = diag_factory.create("NumberPath", comm, params);
  diag_ice->set_grids(gm);
  diags.emplace("inp", diag_ice);
  // Rain
  params.set<std::string>("Number Kind", "Rain");
  auto diag_rain = diag_factory.create("NumberPath", comm, params);
  diag_rain->set_grids(gm);
  diags.emplace("rnp", diag_rain);

  // Set the required fields for the diagnostic.
  std::map<std::string, Field> input_fields;
  for(const auto &dd : diags) {
    const auto &diag = dd.second;
    for(const auto &req : diag->get_required_field_requests()) {
      if(input_fields.find(req.fid.name()) == input_fields.end()) {
        Field f(req.fid);
        f.allocate_view();
        f.get_header().get_tracking().update_time_stamp(t0);
        input_fields.emplace(f.name(), f);
      }
      const auto &f = input_fields.at(req.fid.name());
      diag->set_required_field(f.get_const());
      REQUIRE_THROWS(diag->set_computed_field(f));
    }
    // Initialize the diagnostic
    diag->initialize(t0, RunType::Initial);
  }

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto &qc_f = input_fields.at("qc");
    const auto &qc_v = qc_f.get_view<Real **>();
    const auto &qc_h = qc_f.get_view<Real **, Host>();
    const auto &nc_f = input_fields.at("nc");
    const auto &nc_v = nc_f.get_view<Real **>();
    const auto &nc_h = nc_f.get_view<Real **, Host>();
    const auto &qi_f = input_fields.at("qi");
    const auto &qi_v = qi_f.get_view<Real **>();
    const auto &qi_h = qi_f.get_view<Real **, Host>();
    const auto &ni_f = input_fields.at("ni");
    const auto &ni_v = ni_f.get_view<Real **>();
    const auto &ni_h = ni_f.get_view<Real **, Host>();
    const auto &qr_f = input_fields.at("qr");
    const auto &qr_v = qr_f.get_view<Real **>();
    const auto &qr_h = qr_f.get_view<Real **, Host>();
    const auto &nr_f = input_fields.at("nr");
    const auto &nr_v = nr_f.get_view<Real **>();
    const auto &nr_h = nr_f.get_view<Real **, Host>();
    const auto &pd_f = input_fields.at("pseudo_density");
    const auto &pd_v = pd_f.get_view<Real **>();
    const auto &pd_h = pd_f.get_view<Real **, Host>();
    for(int icol = 0; icol < ncols; icol++) {
      const auto &qc_sub = ekat::subview(qc_v, icol);
      const auto &nc_sub = ekat::subview(nc_v, icol);
      const auto &qi_sub = ekat::subview(qi_v, icol);
      const auto &ni_sub = ekat::subview(ni_v, icol);
      const auto &qr_sub = ekat::subview(qr_v, icol);
      const auto &nr_sub = ekat::subview(nr_v, icol);
      const auto &dp_sub = ekat::subview(pd_v, icol);
      ekat::genRandArray(qc_sub, engine, pdf_qx);
      ekat::genRandArray(nc_sub, engine, pdf_qx);
      ekat::genRandArray(qi_sub, engine, pdf_qx);
      ekat::genRandArray(ni_sub, engine, pdf_qx);
      ekat::genRandArray(qr_sub, engine, pdf_qx);
      ekat::genRandArray(nr_sub, engine, pdf_qx);
      ekat::genRandArray(dp_sub, engine, pdf_pd);
    }
    // Grab views for each of the number path diagnostics
    const auto &lnp   = diags["lnp"]->get_diagnostic();
    const auto &inp   = diags["inp"]->get_diagnostic();
    const auto &rnp   = diags["rnp"]->get_diagnostic();
    const auto &lnp_h = lnp.get_view<const Real *, Host>();
    const auto &inp_h = inp.get_view<const Real *, Host>();
    const auto &rnp_h = rnp.get_view<const Real *, Host>();
    // Sync to host
    qc_f.sync_to_host();
    nc_f.sync_to_host();
    qi_f.sync_to_host();
    ni_f.sync_to_host();
    qr_f.sync_to_host();
    nr_f.sync_to_host();
    pd_f.sync_to_host();
    // Compute
    for(const auto &dd : diags) {
      dd.second->compute_diagnostic();
    }
    // Sync to host
    lnp.sync_to_host();
    inp.sync_to_host();
    rnp.sync_to_host();
    // Test manual calculation vs one provided by diags
    {
      for(int icol = 0; icol < ncols; icol++) {
        Real qndc_prod = 0.0;
        for(int ilev = 0; ilev < num_levs; ++ilev) {
          qndc_prod +=
              nc_h(icol, ilev) * qc_h(icol, ilev) * pd_h(icol, ilev) / g;
        }
        REQUIRE(std::abs(lnp_h(icol) - qndc_prod) < macheps);
        Real qndi_prod = 0.0;
        for(int ilev = 0; ilev < num_levs; ++ilev) {
          qndi_prod +=
              ni_h(icol, ilev) * qi_h(icol, ilev) * pd_h(icol, ilev) / g;
        }
        REQUIRE(std::abs(inp_h(icol) - qndi_prod) < macheps);
        Real qndr_prod = 0.0;
        for(int ilev = 0; ilev < num_levs; ++ilev) {
          qndr_prod +=
              nr_h(icol, ilev) * qr_h(icol, ilev) * pd_h(icol, ilev) / g;
        }
        REQUIRE(std::abs(rnp_h(icol) - qndr_prod) < macheps);
      }
    }
  }

  // Finalize the diagnostic
  for(const auto &dd : diags) {
    const auto &diag = dd.second;
    diag->finalize();
  }

}  // run()

TEST_CASE("number_path_test", "number_path_test]") {
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 10;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for(int irun = 0; irun < num_runs; ++irun) {
    run<Device>(engine);
  }
}

}  // namespace scream
