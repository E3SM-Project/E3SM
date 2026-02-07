#include "eamxx_cosp.hpp"
#include "cosp_functions.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{
// =========================================================================================
Cosp::Cosp (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
   // Determine how often to call COSP; units can be steps or hours
  m_cosp_frequency = m_params.get<Int>("cosp_frequency", 1);
  m_cosp_frequency_units = m_params.get<std::string>("cosp_frequency_units", "steps");
  EKAT_REQUIRE_MSG(
    (m_cosp_frequency_units == "steps") || (m_cosp_frequency_units == "hours"),
    "cosp_frequency_units " + m_cosp_frequency_units + " not supported"
  );

  // How many subcolumns to use for COSP
  m_num_subcols = m_params.get<Int>("cosp_subcolumns", 10);
}

// =========================================================================================
void Cosp::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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

  m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

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
  //                  Name in AD     Layout               Units   Grid       Group
  add_field<Required>("surf_radiative_T", scalar2d    , K,      grid_name);
  //add_field<Required>("surfelev",    scalar2d    , m,      grid_name);
  //add_field<Required>("landmask",    scalar2d    , nondim, grid_name);
  add_field<Required>("sunlit_mask",       scalar2d    , nondim, grid_name);
  add_field<Required>("p_mid",             scalar3d_mid, Pa,     grid_name);
  add_field<Required>("p_int",             scalar3d_int, Pa,     grid_name);
  //add_field<Required>("height_mid",  scalar3d_mid, m,      grid_name);
  //add_field<Required>("height_int",  scalar3d_int, m,      grid_name);
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
  // TODO: should these be meters or microns? Was meters before, but using "m" instead
  // of "micron" seemed to cause prim_model_finalize to throw error with the following:
  // ABORTING WITH ERROR: Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).
  // P3 defines this field with micron instead of meters units, so is this a unit conversion issue?
  add_field<Required>("eff_radius_qc",     scalar3d_mid, micron,      grid_name);
  add_field<Required>("eff_radius_qi",     scalar3d_mid, micron,      grid_name);
  // Set of fields used strictly as output
  // NOTE we set their corresponding masks in init impl
  add_field<Computed>("isccp_cldtot", scalar2d, percent, grid_name);
  add_field<Computed>("isccp_ctptau", scalar4d_ctptau, percent, grid_name, 1);
  add_field<Computed>("modis_ctptau", scalar4d_ctptau, percent, grid_name, 1);
  add_field<Computed>("misr_cthtau", scalar4d_cthtau, percent, grid_name, 1);

  // We can allocate these now
  m_z_mid = Field(FieldIdentifier("z_mid",scalar3d_mid,m,grid_name));
  m_z_int = Field(FieldIdentifier("z_int",scalar3d_int,m,grid_name));
  m_z_mid.allocate_view();
  m_z_int.allocate_view();

  // Add COSP dimension coordinate values and bounds as geometry data
  // These are needed for e3sm_diags to produce COSP diagnostics
  // Values match those hard-coded in components/eam/src/physics/cosp2/local/cosp_config.F90
  using namespace ShortFieldTagsNames;
  
  // COSP tau bins (optical depth) - standard ISCCP bins
  // From COSP: tau_binEdges = reshape((/0.0, 0.3, 0.3, 1.3, 1.3, 3.6, 3.6, 9.4, 9.4, 23.0, 23.0, 60.0, 60.0, 100000.0/))
  if (not m_grid->has_geometry_data("cosp_tau_bnds")) {
    FieldLayout tau_bnds_layout({CMP,CMP},{m_num_tau,2},{"cosp_tau","nbnd"});
    auto& tau_bnds = m_grid->create_geometry_data("cosp_tau_bnds", tau_bnds_layout, nondim);
    auto tau_bnds_h = tau_bnds.get_view<Real**,Host>();
    
    // Standard ISCCP tau bin edges from COSP Fortran source
    tau_bnds_h(0,0) = 0.0;   tau_bnds_h(0,1) = 0.3;
    tau_bnds_h(1,0) = 0.3;   tau_bnds_h(1,1) = 1.3;
    tau_bnds_h(2,0) = 1.3;   tau_bnds_h(2,1) = 3.6;
    tau_bnds_h(3,0) = 3.6;   tau_bnds_h(3,1) = 9.4;
    tau_bnds_h(4,0) = 9.4;   tau_bnds_h(4,1) = 23.0;
    tau_bnds_h(5,0) = 23.0;  tau_bnds_h(5,1) = 60.0;
    tau_bnds_h(6,0) = 60.0;  tau_bnds_h(6,1) = 100000.0;
    
    tau_bnds.sync_to_dev();
  }
  
  // COSP tau centers - use hard-coded values from COSP
  // From COSP: tau_binCenters = (/0.15, 0.80, 2.45, 6.5, 16.2, 41.5, 100.0/)
  if (not m_grid->has_geometry_data("cosp_tau")) {
    FieldLayout tau_layout({CMP},{m_num_tau},{"cosp_tau"});
    auto& tau = m_grid->create_geometry_data("cosp_tau", tau_layout, nondim);
    auto tau_h = tau.get_view<Real*,Host>();
    
    // Standard ISCCP tau bin centers from COSP Fortran source
    tau_h(0) = 0.15;
    tau_h(1) = 0.80;
    tau_h(2) = 2.45;
    tau_h(3) = 6.5;
    tau_h(4) = 16.2;
    tau_h(5) = 41.5;
    tau_h(6) = 100.0;
    
    tau.sync_to_dev();
  }
  
  // COSP pressure bins (cloud-top pressure in Pa) - standard ISCCP bins
  // From COSP: pres_binEdges = reshape((/100000, 80000, 80000, 68000, 68000, 56000, 56000, 44000, 44000, 31000, 31000, 18000, 18000, 0/))
  if (not m_grid->has_geometry_data("cosp_prs_bnds")) {
    FieldLayout prs_bnds_layout({CMP,CMP},{m_num_ctp,2},{"cosp_prs","nbnd"});
    auto& prs_bnds = m_grid->create_geometry_data("cosp_prs_bnds", prs_bnds_layout, Pa);
    auto prs_bnds_h = prs_bnds.get_view<Real**,Host>();
    
    // Standard ISCCP pressure bin edges from COSP Fortran source (in Pa)
    prs_bnds_h(0,0) = 100000.0; prs_bnds_h(0,1) = 80000.0;
    prs_bnds_h(1,0) = 80000.0;  prs_bnds_h(1,1) = 68000.0;
    prs_bnds_h(2,0) = 68000.0;  prs_bnds_h(2,1) = 56000.0;
    prs_bnds_h(3,0) = 56000.0;  prs_bnds_h(3,1) = 44000.0;
    prs_bnds_h(4,0) = 44000.0;  prs_bnds_h(4,1) = 31000.0;
    prs_bnds_h(5,0) = 31000.0;  prs_bnds_h(5,1) = 18000.0;
    prs_bnds_h(6,0) = 18000.0;  prs_bnds_h(6,1) = 0.0;
    
    prs_bnds.sync_to_dev();
  }
  
  // COSP pressure centers - use hard-coded values from COSP
  // From COSP: pres_binCenters = (/90000., 74000., 62000., 50000., 37500., 24500., 9000./)
  if (not m_grid->has_geometry_data("cosp_prs")) {
    FieldLayout prs_layout({CMP},{m_num_ctp},{"cosp_prs"});
    auto& prs = m_grid->create_geometry_data("cosp_prs", prs_layout, Pa);
    auto prs_h = prs.get_view<Real*,Host>();
    
    // Standard ISCCP pressure bin centers from COSP Fortran source (in Pa)
    prs_h(0) = 90000.0;
    prs_h(1) = 74000.0;
    prs_h(2) = 62000.0;
    prs_h(3) = 50000.0;
    prs_h(4) = 37500.0;
    prs_h(5) = 24500.0;
    prs_h(6) = 9000.0;
    
    prs.sync_to_dev();
  }
  
  // COSP height bins (cloud-top height in m) - standard MISR bins  
  // From COSP: hgt_binEdges = 1000*reshape((/-99.0, 0.0, 0.0, 0.5, 0.5, 1.0, ..., 17.0, 99.0/))
  // Note: first bin (-99 to 0) is "no retrieval" bin for cases with no cloud top height retrieval
  if (not m_grid->has_geometry_data("cosp_cth_bnds")) {
    FieldLayout cth_bnds_layout({CMP,CMP},{m_num_cth,2},{"cosp_cth","nbnd"});
    auto& cth_bnds = m_grid->create_geometry_data("cosp_cth_bnds", cth_bnds_layout, m);
    auto cth_bnds_h = cth_bnds.get_view<Real**,Host>();
    
    // Standard MISR height bin edges from COSP Fortran source (in meters)
    cth_bnds_h(0,0) = -99000.0;  cth_bnds_h(0,1) = 0.0;
    cth_bnds_h(1,0) = 0.0;       cth_bnds_h(1,1) = 500.0;
    cth_bnds_h(2,0) = 500.0;     cth_bnds_h(2,1) = 1000.0;
    cth_bnds_h(3,0) = 1000.0;    cth_bnds_h(3,1) = 1500.0;
    cth_bnds_h(4,0) = 1500.0;    cth_bnds_h(4,1) = 2000.0;
    cth_bnds_h(5,0) = 2000.0;    cth_bnds_h(5,1) = 2500.0;
    cth_bnds_h(6,0) = 2500.0;    cth_bnds_h(6,1) = 3000.0;
    cth_bnds_h(7,0) = 3000.0;    cth_bnds_h(7,1) = 4000.0;
    cth_bnds_h(8,0) = 4000.0;    cth_bnds_h(8,1) = 5000.0;
    cth_bnds_h(9,0) = 5000.0;    cth_bnds_h(9,1) = 7000.0;
    cth_bnds_h(10,0) = 7000.0;   cth_bnds_h(10,1) = 9000.0;
    cth_bnds_h(11,0) = 9000.0;   cth_bnds_h(11,1) = 11000.0;
    cth_bnds_h(12,0) = 11000.0;  cth_bnds_h(12,1) = 13000.0;
    cth_bnds_h(13,0) = 13000.0;  cth_bnds_h(13,1) = 15000.0;
    cth_bnds_h(14,0) = 15000.0;  cth_bnds_h(14,1) = 17000.0;
    cth_bnds_h(15,0) = 17000.0;  cth_bnds_h(15,1) = 99000.0;
    
    cth_bnds.sync_to_dev();
  }
  
  // COSP height centers - use hard-coded values from COSP
  // From COSP: hgt_binCenters = 1000*(/0., 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 6., 8., 10., 12., 14.5, 16., 18./)
  if (not m_grid->has_geometry_data("cosp_cth")) {
    FieldLayout cth_layout({CMP},{m_num_cth},{"cosp_cth"});
    auto& cth = m_grid->create_geometry_data("cosp_cth", cth_layout, m);
    auto cth_h = cth.get_view<Real*,Host>();
    
    // Standard MISR height bin centers from COSP Fortran source (in meters)
    cth_h(0) = 0.0;
    cth_h(1) = 250.0;
    cth_h(2) = 750.0;
    cth_h(3) = 1250.0;
    cth_h(4) = 1750.0;
    cth_h(5) = 2250.0;
    cth_h(6) = 2750.0;
    cth_h(7) = 3500.0;
    cth_h(8) = 4500.0;
    cth_h(9) = 6000.0;
    cth_h(10) = 8000.0;
    cth_h(11) = 10000.0;
    cth_h(12) = 12000.0;
    cth_h(13) = 14500.0;
    cth_h(14) = 16000.0;
    cth_h(15) = 18000.0;
    
    cth.sync_to_dev();
  }
}

// =========================================================================================
void Cosp::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  CospFunc::initialize(m_num_cols, m_num_subcols, m_num_levs);

  // Set the mask field for each of the cosp computed fields
  std::list<std::string> vnames = {"isccp_cldtot", "isccp_ctptau", "modis_ctptau", "misr_cthtau"};
  for (const auto& field_name : vnames) {
    // the mask here is just the sunlit mask, so set it
    get_field_out(field_name).get_header().set_extra_data("mask_field", get_field_in("sunlit_mask"));
    get_field_out(field_name).get_header().set_extra_data("mask_value", constants::fill_value<Real>);
    get_field_out(field_name).get_header().set_may_be_filled(true);
  }
}

// =========================================================================================
void Cosp::run_impl (const double dt)
{
  // Determine if we should update COSP this timestep
  // First get frequency in steps
  int cosp_freq_in_steps = 1;
  if (m_cosp_frequency_units == "steps") {
      cosp_freq_in_steps = m_cosp_frequency;
  } else if (m_cosp_frequency_units == "hours") {
      EKAT_REQUIRE_MSG((3600 % int(dt)) == 0, "cosp_frequency_units is hours but dt does not evenly divide 1 hour");
      cosp_freq_in_steps = 3600.0 * m_cosp_frequency / dt;
  } else {
      EKAT_ERROR_MSG("cosp_frequency_units " + m_cosp_frequency_units + " not supported");
  }

  // Make sure cosp frequency is multiple of rad frequency?

  // Compare frequency in steps with current timestep
  auto update_cosp = cosp_do(cosp_freq_in_steps, end_of_step_ts().get_num_steps());

  // Call COSP wrapper routines
  if (update_cosp) {
    // Get fields from field manager; note that we get host views because this
    // interface serves primarily as a wrapper to a c++ to f90 bridge for the COSP
    // all then need to be copied to layoutLeft views to permute the indices for
    // F90.

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
        const auto dz_s = z_mid_s; // 
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

    auto isccp_cldtot_h = get_field_out("isccp_cldtot").get_view<Real*, Host>();
    auto isccp_ctptau_h = get_field_out("isccp_ctptau").get_view<Real***, Host>();
    auto modis_ctptau_h = get_field_out("modis_ctptau").get_view<Real***, Host>();
    auto misr_cthtau_h  = get_field_out("misr_cthtau"). get_view<Real***, Host>();

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
    get_field_out("isccp_cldtot").sync_to_dev();
    get_field_out("isccp_ctptau").sync_to_dev();
    get_field_out("modis_ctptau").sync_to_dev();
    get_field_out("misr_cthtau").sync_to_dev();
  }
}

// =========================================================================================
void Cosp::finalize_impl()
{
  // Finalize COSP wrappers
  CospFunc::finalize();
}

} // namespace scream
