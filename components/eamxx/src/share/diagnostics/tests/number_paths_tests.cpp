#include <iomanip>

#include "catch2/catch.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

template <typename DeviceT>
void run()
{
  using PC = physics::Constants<Real>;

  constexpr Real g   = PC::gravit.value;
  constexpr Real tol = PC::macheps;

  ekat::Comm comm(MPI_COMM_WORLD);

  // For inputs randomization
  int seed = get_random_test_seed();

  // Create a grid
  const int ncols = 10;
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();
  ekat::ParameterList params;

  REQUIRE_THROWS(diag_factory.create("NumberPath", comm, params, grid));  // No 'number_kind'
  params.set<std::string>("number_kind", "Foo");
  REQUIRE_THROWS(diag_factory.create("NumberPath", comm, params, grid));  // Invalid 'number_kind'

  std::map<std::string,std::shared_ptr<AbstractDiagnostic>> diags;

  // Liquid
  params.set<std::string>("number_kind", "Liq");
  auto diag_liq = diag_factory.create("NumberPath", comm, params, grid);
  diags.emplace("lnp", diag_liq);
  // Ice
  params.set<std::string>("number_kind", "Ice");
  auto diag_ice = diag_factory.create("NumberPath", comm, params, grid);
  diags.emplace("inp", diag_ice);
  // Rain
  params.set<std::string>("number_kind", "Rain");
  auto diag_rain = diag_factory.create("NumberPath", comm, params, grid);
  diags.emplace("rnp", diag_rain);

  // Set the required fields for the diagnostic.
  std::map<std::string, Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    for(const auto &dd : diags) {
      const auto &diag = dd.second;
      for(const auto &fname : diag->get_input_fields_names()) {
        if(input_fields.find(fname) == input_fields.end()) {
          FieldIdentifier fid(fname, scalar3d, Pa, grid->name());
          Field f(fid);
          f.allocate_view();
          f.get_header().get_tracking().update_time_stamp(t0);
          input_fields.emplace(fname, f);
        }
        diag->set_input_field(input_fields.at(fname));
      }
      // Initialize the diagnostic
      diag->initialize();
    }
  }

  // Run tests
  {
    // Randomize input data
    auto& qc_f = input_fields.at("qc");
    auto& nc_f = input_fields.at("nc");
    auto& qi_f = input_fields.at("qi");
    auto& ni_f = input_fields.at("ni");
    auto& qr_f = input_fields.at("qr");
    auto& nr_f = input_fields.at("nr");
    auto& dp_f = input_fields.at("pseudo_density");
    randomize_uniform(qc_f,seed++,0,1e-3);
    randomize_uniform(nc_f,seed++,0,1e-3);
    randomize_uniform(qi_f,seed++,0,1e-3);
    randomize_uniform(ni_f,seed++,0,1e-3);
    randomize_uniform(qr_f,seed++,0,1e-3);
    randomize_uniform(nr_f,seed++,0,1e-3);
    randomize_uniform(dp_f,seed++,1,100);

    // Grab views for each of the number path diagnostics
    const auto &lnp   = diags["lnp"]->get_diagnostic();
    const auto &inp   = diags["inp"]->get_diagnostic();
    const auto &rnp   = diags["rnp"]->get_diagnostic();
    const auto &lnp_h = lnp.get_view<const Real *, Host>();
    const auto &inp_h = inp.get_view<const Real *, Host>();
    const auto &rnp_h = rnp.get_view<const Real *, Host>();

    // Compute & sync
    for(const auto &dd : diags) {
      dd.second->compute_diagnostic(t0);
      dd.second->get_diagnostic().sync_to_host();
    }

    // Test manual calculation vs one provided by diags
    auto qc_h = qc_f.get_view<const Real **, Host>();
    auto nc_h = nc_f.get_view<const Real **, Host>();
    auto qi_h = qi_f.get_view<const Real **, Host>();
    auto ni_h = ni_f.get_view<const Real **, Host>();
    auto qr_h = qr_f.get_view<const Real **, Host>();
    auto nr_h = nr_f.get_view<const Real **, Host>();
    auto dp_h = dp_f.get_view<const Real **, Host>();

    for(int icol = 0; icol < ncols; icol++) {
      auto scan = [&](auto lambda) {
        Real sum = 0;
        for(int ilev = 0; ilev < nlevs; ++ilev)
          sum += lambda(icol,ilev);
        return sum;
      };

      auto qndc = [&](int icol, int ilev) {
        return nc_h(icol, ilev) * qc_h(icol, ilev) * dp_h(icol, ilev) / g;
      };
      auto qndi = [&](int icol, int ilev) {
        return ni_h(icol, ilev) * qi_h(icol, ilev) * dp_h(icol, ilev) / g;
      };
      auto qndr = [&](int icol, int ilev) {
        return nr_h(icol, ilev) * qr_h(icol, ilev) * dp_h(icol, ilev) / g;
      };

      REQUIRE(std::abs(lnp_h(icol) - scan(qndc)) < tol);
      REQUIRE(std::abs(inp_h(icol) - scan(qndi)) < tol);
      REQUIRE(std::abs(rnp_h(icol) - scan(qndr)) < tol);
    }
  }
}

TEST_CASE("number_path_test", "number_path_test]") {
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 10;

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for(int irun = 0; irun < num_runs; ++irun) {
    run<Device>();
  }
}

}  // namespace scream
