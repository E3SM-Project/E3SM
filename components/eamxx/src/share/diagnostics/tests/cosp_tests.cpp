#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <cmath>

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

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  constexpr int nlevs = 72;
  const int ncols = 2*comm.size();
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("physics");
  const auto& grid_name = grid->name();
  const int nlcols = grid->get_num_local_dofs();

  // Register diagnostics
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();

  // Create COSP diagnostic
  ekat::ParameterList params("cosp_test");
  params.set<int>("cosp_subcolumns", 10);
  auto diag = diag_factory.create("Cosp",comm,params);
  diag->set_grids(gm);

  // Verify it's multi-output
  REQUIRE(diag->is_multi_output());
  auto names = diag->get_diagnostic_names();
  REQUIRE(names.size() == 4);

  // Units
  Units nondim = Units::nondimensional();
  Units percent (nondim,"%");
  auto micron = micro*m;
  auto m2 = pow(m, 2);
  auto s2 = pow(s, 2);

  // Layouts
  FieldLayout scalar2d     = grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = grid->get_3d_scalar_layout(false);

  // Helper to create and allocate a field
  auto make_field = [&](const std::string& name, const FieldLayout& layout,
                        const Units& units) -> Field {
    FieldIdentifier fid(name, layout, units, grid_name);
    Field f(fid);
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    return f;
  };

  // 2D fields
  auto surf_radiative_T = make_field("surf_radiative_T", scalar2d, K);
  auto sunlit_mask      = make_field("sunlit_mask",      scalar2d, nondim);
  auto phis             = make_field("phis",             scalar2d, m2/s2);

  surf_radiative_T.deep_copy(288.0);
  sunlit_mask.deep_copy(1.0);
  phis.deep_copy(0.0);

  // 3D fields — need realistic vertical profiles for the Fortran code.
  // Build a simple atmosphere: pressure decreasing with height from
  // ~1013 hPa at surface to ~1 hPa at top.
  auto T_mid          = make_field("T_mid",          scalar3d_mid, K);
  auto p_mid          = make_field("p_mid",          scalar3d_mid, Pa);
  auto pseudo_density = make_field("pseudo_density", scalar3d_mid, Pa);
  auto cldfrac_rad    = make_field("cldfrac_rad",    scalar3d_mid, nondim);
  auto dtau067        = make_field("dtau067",        scalar3d_mid, nondim);
  auto dtau105        = make_field("dtau105",        scalar3d_mid, nondim);
  auto eff_radius_qc  = make_field("eff_radius_qc", scalar3d_mid, micron);
  auto eff_radius_qi  = make_field("eff_radius_qi", scalar3d_mid, micron);
  auto p_int          = make_field("p_int",          scalar3d_int, Pa);
  auto qv             = make_field("qv",             scalar3d_mid, kg/kg);
  auto qc             = make_field("qc",             scalar3d_mid, kg/kg);
  auto qi             = make_field("qi",             scalar3d_mid, kg/kg);

  // Set up a realistic vertical profile on host
  {
    auto T_h    = T_mid.get_view<Real**, Host>();
    auto p_m_h  = p_mid.get_view<Real**, Host>();
    auto p_i_h  = p_int.get_view<Real**, Host>();
    auto pd_h   = pseudo_density.get_view<Real**, Host>();
    auto qv_h   = qv.get_view<Real**, Host>();

    const Real p_top  = 100.0;     // Top pressure [Pa]

    for (int icol = 0; icol < nlcols; ++icol) {
      // Vary surface pressure across columns; COSP Fortran internally
      // divides by (max_psurf - min_psurf) across columns, so identical
      // surface pressures cause division by zero.
      const Real p_surf = 101325.0 + icol * 50.0;

      // Interface pressures: evenly spaced in log(p)
      for (int k = 0; k <= nlevs; ++k) {
        Real frac = static_cast<Real>(k) / nlevs;
        p_i_h(icol, k) = p_top * std::exp(frac * std::log(p_surf / p_top));
      }
      // Mid-level pressure and temperature
      for (int k = 0; k < nlevs; ++k) {
        p_m_h(icol, k) = 0.5 * (p_i_h(icol, k) + p_i_h(icol, k + 1));
        // pseudo_density = dp / g  (approximate layer mass)
        pd_h(icol, k) = p_i_h(icol, k + 1) - p_i_h(icol, k);
        // Standard atmosphere temperature: ~220K at top, ~288K at surface
        Real frac = static_cast<Real>(k) / (nlevs - 1);
        T_h(icol, k) = 220.0 + 68.0 * frac;
        // Specific humidity decreasing with height
        qv_h(icol, k) = 0.01 * frac * frac;
      }
    }

    T_mid.sync_to_dev();
    p_mid.sync_to_dev();
    p_int.sync_to_dev();
    pseudo_density.sync_to_dev();
    qv.sync_to_dev();
  }

  // COSP-specific fields: constant values are fine for these
  cldfrac_rad.deep_copy(0.5);
  dtau067.deep_copy(1.0);
  dtau105.deep_copy(1.0);
  eff_radius_qc.deep_copy(10.0);
  eff_radius_qi.deep_copy(25.0);
  qc.deep_copy(1e-5);
  qi.deep_copy(1e-6);

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

  // Initialize and compute
  diag->initialize(t0, RunType::Initial);
  diag->compute_diagnostic();

  // Check that all output fields are allocated and have valid timestamps
  for (const auto& oname : names) {
    auto out = diag->get_diagnostic(oname);
    REQUIRE(out.is_allocated());
    REQUIRE(out.get_header().get_tracking().get_time_stamp().is_valid());
  }

  // Check isccp_cldtot values: with sunlit_mask=1, all columns should
  // have valid (non-fill) values in [0, 100]%
  constexpr auto fill_value = constants::fill_value<Real>;
  auto cldtot = diag->get_diagnostic("isccp_cldtot");
  cldtot.sync_to_host();
  auto cldtot_h = cldtot.get_view<const Real*, Host>();
  for (int icol = 0; icol < nlcols; ++icol) {
    REQUIRE(cldtot_h(icol) != fill_value);
    REQUIRE(cldtot_h(icol) >= 0.0);
    REQUIRE(cldtot_h(icol) <= 100.0);
  }

  // Test nighttime masking: set sunlit_mask to 0 and recompute
  sunlit_mask.deep_copy(0.0);
  util::TimeStamp t1 ({2022,1,1},{0,0,1});
  for (auto* f : {&surf_radiative_T, &sunlit_mask, &phis, &T_mid, &p_mid, &p_int,
                   &pseudo_density, &cldfrac_rad, &dtau067, &dtau105,
                   &eff_radius_qc, &eff_radius_qi, &qv, &qc, &qi}) {
    f->get_header().get_tracking().update_time_stamp(t1);
  }

  diag->compute_diagnostic();

  cldtot.sync_to_host();
  cldtot_h = cldtot.get_view<const Real*, Host>();
  for (int icol = 0; icol < nlcols; ++icol) {
    REQUIRE(cldtot_h(icol) == fill_value);
  }

  diag->finalize();
}

} // namespace scream
