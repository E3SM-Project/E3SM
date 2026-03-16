#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"point_grid"});
  auto& pl = gm_params.sublist("point_grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("cosp_diagnostic")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using namespace ekat::prefixes;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager
  constexpr int nlevs = 72;
  const int ncols = 2*comm.size();
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("physics");
  const auto& grid_name = grid->name();

  // Register diagnostics
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();

  // Create COSP diagnostic with parameters
  ekat::ParameterList params("cosp_test");
  params.set<int>("cosp_subcolumns", 10);
  auto diag = diag_factory.create("Cosp",comm,params);
  diag->set_grids(gm);

  // Verify it's multi-output
  REQUIRE(diag->is_multi_output());
  auto names = diag->get_diagnostic_names();
  REQUIRE(names.size() == 4);

  // Create input fields with physically reasonable constant values
  // and set them on the diagnostic
  Units nondim = Units::nondimensional();
  Units percent (nondim,"%");
  auto micron = micro*m;
  auto m2 = pow(m, 2);
  auto s2 = pow(s, 2);

  FieldLayout scalar2d     = grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = grid->get_3d_scalar_layout(false);

  auto make_field = [&](const std::string& name, const FieldLayout& layout,
                        const Units& units, Real fill_val) -> Field {
    FieldIdentifier fid(name, layout, units, grid_name);
    Field f(fid);
    f.allocate_view();
    f.deep_copy(fill_val);
    f.get_header().get_tracking().update_time_stamp(t0);
    return f;
  };

  // 2D fields
  auto surf_radiative_T = make_field("surf_radiative_T", scalar2d, K, 288.0);
  auto sunlit_mask      = make_field("sunlit_mask",      scalar2d, nondim, 1.0);
  auto phis             = make_field("phis",             scalar2d, m2/s2, 0.0);

  // 3D mid-level fields
  auto T_mid           = make_field("T_mid",           scalar3d_mid, K, 270.0);
  auto p_mid           = make_field("p_mid",           scalar3d_mid, Pa, 50000.0);
  auto pseudo_density   = make_field("pseudo_density",  scalar3d_mid, Pa, 100.0);
  auto cldfrac_rad      = make_field("cldfrac_rad",     scalar3d_mid, nondim, 0.5);
  auto dtau067          = make_field("dtau067",         scalar3d_mid, nondim, 1.0);
  auto dtau105          = make_field("dtau105",         scalar3d_mid, nondim, 1.0);
  auto eff_radius_qc    = make_field("eff_radius_qc",  scalar3d_mid, micron, 10.0);
  auto eff_radius_qi    = make_field("eff_radius_qi",  scalar3d_mid, micron, 25.0);

  // 3D interface field
  auto p_int = make_field("p_int", scalar3d_int, Pa, 50000.0);

  // Tracers (3D mid, with "tracers" group)
  auto qv = make_field("qv", scalar3d_mid, kg/kg, 0.01);
  auto qc = make_field("qc", scalar3d_mid, kg/kg, 1e-5);
  auto qi = make_field("qi", scalar3d_mid, kg/kg, 1e-6);

  // Set required fields on the diagnostic
  diag->set_required_field(surf_radiative_T.get_const());
  diag->set_required_field(sunlit_mask.get_const());
  diag->set_required_field(phis.get_const());
  diag->set_required_field(T_mid.get_const());
  diag->set_required_field(p_mid.get_const());
  diag->set_required_field(p_int.get_const());
  diag->set_required_field(pseudo_density.get_const());
  diag->set_required_field(cldfrac_rad.get_const());
  diag->set_required_field(dtau067.get_const());
  diag->set_required_field(dtau105.get_const());
  diag->set_required_field(eff_radius_qc.get_const());
  diag->set_required_field(eff_radius_qi.get_const());
  diag->set_required_field(qv.get_const());
  diag->set_required_field(qc.get_const());
  diag->set_required_field(qi.get_const());

  // Initialize (this calls CospFunc::initialize)
  diag->initialize(t0,RunType::Initial);

  // Compute the diagnostic
  diag->compute_diagnostic();

  // Check that all output fields are allocated and have valid timestamps
  for (const auto& oname : names) {
    auto out = diag->get_diagnostic(oname);
    REQUIRE(out.is_allocated());
    REQUIRE(out.get_header().get_tracking().get_time_stamp().is_valid());
  }

  // Check isccp_cldtot values: with sunlit_mask=1, all columns should
  // have valid (non-fill) values
  auto cldtot = diag->get_diagnostic("isccp_cldtot");
  cldtot.sync_to_host();
  auto cldtot_h = cldtot.get_view<const Real*, Host>();
  constexpr auto fill_value = constants::fill_value<Real>;
  for (int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
    // Sunlit columns should not have fill values
    REQUIRE(cldtot_h(icol) != fill_value);
    // Cloud total should be in [0, 100] percent range
    REQUIRE(cldtot_h(icol) >= 0.0);
    REQUIRE(cldtot_h(icol) <= 100.0);
  }

  // Test nighttime masking: set sunlit_mask to 0 and recompute
  sunlit_mask.deep_copy(0.0);
  // Advance timestamp so the diagnostic recomputes
  util::TimeStamp t1 ({2022,1,1},{0,0,1});
  for (auto* f : {&surf_radiative_T, &sunlit_mask, &phis, &T_mid, &p_mid, &p_int,
                   &pseudo_density, &cldfrac_rad, &dtau067, &dtau105,
                   &eff_radius_qc, &eff_radius_qi, &qv, &qc, &qi}) {
    f->get_header().get_tracking().update_time_stamp(t1);
  }

  diag->compute_diagnostic();

  cldtot.sync_to_host();
  cldtot_h = cldtot.get_view<const Real*, Host>();
  for (int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
    // Nighttime columns should have fill values
    REQUIRE(cldtot_h(icol) == fill_value);
  }

  // Finalize (this calls CospFunc::finalize)
  diag->finalize();
}

} // namespace scream
