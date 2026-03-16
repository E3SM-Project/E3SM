#include "eamxx_cosp.hpp"
#include "cosp_functions.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{
// =========================================================================================
Cosp::Cosp (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm, params)
{
  // How many subcolumns to use for COSP
  m_num_subcols = m_params.get<int>("cosp_subcolumns", 10);
}

// =========================================================================================
std::vector<std::string> Cosp::get_diagnostic_names () const
{
  return {"isccp_cldtot", "isccp_ctptau", "modis_ctptau", "misr_cthtau"};
}

// =========================================================================================
void Cosp::create_requests()
{
  using namespace ekat::units;
  using namespace ekat::prefixes;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  Units nondim = Units::nondimensional();
  Units percent (nondim,"%");
  auto micron = micro*m;
  auto m2 = pow(m, 2);
  auto s2 = pow(s, 2);

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this diagnostic

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);
  FieldLayout scalar4d_ctptau ( {COL,CMP,CMP},
                                {m_num_cols,m_num_tau,m_num_ctp},
                                {e2str(COL), "cosp_tau", "cosp_prs"});
  FieldLayout scalar4d_cthtau ( {COL,CMP,CMP},
                                {m_num_cols,m_num_tau,m_num_cth},
                                {e2str(COL), "cosp_tau", "cosp_cth"});

  // Set of fields used strictly as input
  add_field<Required>("surf_radiative_T", scalar2d    , K,      grid_name);
  add_field<Required>("sunlit_mask",       scalar2d    , nondim, grid_name);
  add_field<Required>("p_mid",             scalar3d_mid, Pa,     grid_name);
  add_field<Required>("p_int",             scalar3d_int, Pa,     grid_name);
  add_field<Required>("T_mid",            scalar3d_mid, K,      grid_name);
  add_field<Required>("phis",             scalar2d    , m2/s2,  grid_name);
  add_field<Required>("pseudo_density",   scalar3d_mid, Pa,     grid_name);
  add_field<Required>("cldfrac_rad",      scalar3d_mid, nondim, grid_name);
  add_tracer<Required>("qv", m_grid, kg/kg);
  add_tracer<Required>("qc", m_grid, kg/kg);
  add_tracer<Required>("qi", m_grid, kg/kg);
  // Optical properties, should be computed in radiation interface
  add_field<Required>("dtau067",     scalar3d_mid, nondim, grid_name); // 0.67 micron optical depth
  add_field<Required>("dtau105",     scalar3d_mid, nondim, grid_name); // 10.5 micron optical depth
  // Effective radii, should be computed in either microphysics or radiation interface
  add_field<Required>("eff_radius_qc",     scalar3d_mid, micron,      grid_name);
  add_field<Required>("eff_radius_qi",     scalar3d_mid, micron,      grid_name);

  // Allocate output fields in multi-output map (not via add_field<Computed>)
  auto mk_output = [&](const std::string& name, const FieldLayout& layout) {
    Field f(FieldIdentifier(name, layout, percent, grid_name));
    f.allocate_view();
    m_diagnostic_outputs[name] = f;
  };
  mk_output("isccp_cldtot", scalar2d);
  mk_output("isccp_ctptau", scalar4d_ctptau);
  mk_output("modis_ctptau", scalar4d_ctptau);
  mk_output("misr_cthtau",  scalar4d_cthtau);

  // Scratch fields for height computation
  m_z_mid = Field(FieldIdentifier("z_mid",scalar3d_mid,m,grid_name));
  m_z_int = Field(FieldIdentifier("z_int",scalar3d_int,m,grid_name));
  m_z_mid.allocate_view();
  m_z_int.allocate_view();
}

// =========================================================================================
void Cosp::initialize_impl (const RunType /* run_type */)
{
  // Initialize the COSP Fortran library
  CospFunc::initialize(m_num_cols, m_num_subcols, m_num_levs);

  // Set the mask field for each of the COSP output fields
  for (auto& [fname, field] : m_diagnostic_outputs) {
    field.get_header().set_extra_data("mask_field", get_field_in("sunlit_mask"));
    field.get_header().set_may_be_filled(true);
  }
}

// =========================================================================================
void Cosp::compute_diagnostic_impl ()
{
  // Get fields from field manager; note that we get host views because this
  // interface serves primarily as a wrapper to a C++ to F90 bridge for COSP.
  // All fields then need to be copied to layoutLeft views to permute the
  // indices for F90.

  // Ensure host data of input fields is up to date
  get_field_in("qv").sync_to_host();
  get_field_in("qc").sync_to_host();
  get_field_in("qi").sync_to_host();
  get_field_in("sunlit_mask").sync_to_host();
  get_field_in("surf_radiative_T").sync_to_host();
  get_field_in("T_mid").sync_to_host();
  get_field_in("p_mid").sync_to_host();
  get_field_in("p_int").sync_to_host();
  get_field_in("cldfrac_rad").sync_to_host();
  get_field_in("eff_radius_qc").sync_to_host();
  get_field_in("eff_radius_qi").sync_to_host();
  get_field_in("dtau067").sync_to_host();
  get_field_in("dtau105").sync_to_host();

  // Compute z_mid
  const auto T_mid_d = get_field_in("T_mid").get_view<const Real**>();
  const auto qv_d  = get_field_in("qv").get_view<const Real**>();
  const auto p_mid_d = get_field_in("p_mid").get_view<const Real**>();
  const auto phis_d  = get_field_in("phis").get_view<const Real*>();
  const auto pseudo_density_d = get_field_in("pseudo_density").get_view<const Real**>();
  const auto z_mid_d = m_z_mid.get_view<Real**>();
  const auto z_int_d = m_z_int.get_view<Real**>();
  const auto ncol = m_num_cols;
  const auto nlev = m_num_levs;

  using KT       = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;
  using PF       = scream::PhysicsFunctions<DefaultDevice>;

  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncol, nlev);
  const Real g = physics::Constants<Real>::gravit.value;
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
      const int i = team.league_rank();
      const auto p_mid_s = ekat::subview(p_mid_d, i);
      const auto T_mid_s = ekat::subview(T_mid_d, i);
      const auto qv_s = ekat::subview(qv_d, i);
      const auto z_int_s = ekat::subview(z_int_d, i);
      const auto z_mid_s = ekat::subview(z_mid_d, i);
      const Real z_surf  = phis_d(i) / g;
      const auto pseudo_density_s = ekat::subview(pseudo_density_d, i);

      // 1. Compute dz (recycle z_mid_s as a temporary)
      const auto dz_s = z_mid_s;
      PF::calculate_dz(team, pseudo_density_s, p_mid_s, T_mid_s, qv_s, dz_s);
      team.team_barrier();

      // 2. Compute z_int (vertical scan)
      PF::calculate_z_int(team,nlev,dz_s,z_surf,z_int_s);
      team.team_barrier();

      // 3. Compute z_mid (int->mid interpolation)
      PF::calculate_z_mid(team,nlev,z_int_s,z_mid_s);
      team.team_barrier();
  });
  Kokkos::fence();

  m_z_mid.sync_to_host();
  const auto z_mid_h = m_z_mid.get_view<const Real**,Host>();
  const auto T_mid_h   = get_field_in("T_mid").get_view<const Real**, Host>();
  const auto qv_h      = get_field_in("qv").get_view<const Real**, Host>();
  const auto p_mid_h   = get_field_in("p_mid").get_view<const Real**,Host>();
  const auto qc_h      = get_field_in("qc").get_view<const Real**, Host>();
  const auto qi_h      = get_field_in("qi").get_view<const Real**, Host>();
  const auto sunlit_h  = get_field_in("sunlit_mask").get_view<const Real*, Host>();
  const auto skt_h     = get_field_in("surf_radiative_T").get_view<const Real*, Host>();
  const auto p_int_h   = get_field_in("p_int").get_view<const Real**, Host>();
  const auto cldfrac_h = get_field_in("cldfrac_rad").get_view<const Real**, Host>();
  const auto reff_qc_h = get_field_in("eff_radius_qc").get_view<const Real**, Host>();
  const auto reff_qi_h = get_field_in("eff_radius_qi").get_view<const Real**, Host>();
  const auto dtau067_h = get_field_in("dtau067").get_view<const Real**, Host>();
  const auto dtau105_h = get_field_in("dtau105").get_view<const Real**, Host>();

  auto isccp_cldtot_h = m_diagnostic_outputs.at("isccp_cldtot").get_view<Real*, Host>();
  auto isccp_ctptau_h = m_diagnostic_outputs.at("isccp_ctptau").get_view<Real***, Host>();
  auto modis_ctptau_h = m_diagnostic_outputs.at("modis_ctptau").get_view<Real***, Host>();
  auto misr_cthtau_h  = m_diagnostic_outputs.at("misr_cthtau"). get_view<Real***, Host>();

  Real emsfc_lw = 0.99;
  CospFunc::main(
          m_num_cols, m_num_subcols, m_num_levs, m_num_tau, m_num_ctp, m_num_cth, emsfc_lw,
          sunlit_h, skt_h, T_mid_h, p_mid_h, p_int_h, z_mid_h, qv_h, qc_h, qi_h,
          cldfrac_h, reff_qc_h, reff_qi_h, dtau067_h, dtau105_h,
          isccp_cldtot_h, isccp_ctptau_h, modis_ctptau_h, misr_cthtau_h
  );
  // Mask night values
  constexpr auto fill_value = constants::fill_value<Real>;
  for (int i = 0; i < m_num_cols; i++) {
    if (sunlit_h(i) == 0) {
      // if night, set to fill val
      isccp_cldtot_h(i) = fill_value;
      for (int j = 0; j < m_num_tau; j++) {
        for (int k = 0; k < m_num_ctp; k++) {
          isccp_ctptau_h(i,j,k) = fill_value;
          modis_ctptau_h(i,j,k) = fill_value;
        }
        for (int k = 0; k < m_num_cth; k++) {
          misr_cthtau_h (i,j,k) = fill_value;
        }
      }
    }
  }

  // Make sure dev data is up to date
  m_diagnostic_outputs.at("isccp_cldtot").sync_to_dev();
  m_diagnostic_outputs.at("isccp_ctptau").sync_to_dev();
  m_diagnostic_outputs.at("modis_ctptau").sync_to_dev();
  m_diagnostic_outputs.at("misr_cthtau").sync_to_dev();
}

// =========================================================================================
void Cosp::finalize_impl()
{
  // Finalize COSP wrappers
  CospFunc::finalize();
}

} // namespace scream
