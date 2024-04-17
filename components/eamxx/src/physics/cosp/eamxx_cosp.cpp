#include "eamxx_cosp.hpp"
#include "cosp_functions.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include "share/field/field_utils.hpp"

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
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto nondim = Units::nondimensional();
  auto percent = Units::nondimensional();
  percent.set_string("%");
  auto micron = m / 1000000;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar2d_layout     { {COL},      {m_num_cols}              };
  FieldLayout scalar3d_layout_mid { {COL,LEV},  {m_num_cols,m_num_levs}   };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };
  FieldLayout scalar4d_layout_ctptau { {COL,ISCCPTAU,ISCCPPRS}, {m_num_cols,m_num_isccptau,m_num_isccpctp} };

  // Set of fields used strictly as input
  //                  Name in AD     Layout               Units   Grid       Group
  add_field<Required>("surf_radiative_T", scalar2d_layout    , K,      grid_name);
  //add_field<Required>("surfelev",    scalar2d_layout    , m,      grid_name);
  //add_field<Required>("landmask",    scalar2d_layout    , nondim, grid_name);
  add_field<Required>("sunlit",           scalar2d_layout    , nondim, grid_name);
  add_field<Required>("p_mid",             scalar3d_layout_mid, Pa,     grid_name);
  add_field<Required>("p_int",             scalar3d_layout_int, Pa,     grid_name);
  //add_field<Required>("height_mid",  scalar3d_layout_mid, m,      grid_name);
  //add_field<Required>("height_int",  scalar3d_layout_int, m,      grid_name);
  add_field<Required>("T_mid",            scalar3d_layout_mid, K,      grid_name);
  add_field<Required>("qv",               scalar3d_layout_mid, Q,      grid_name, "tracers");
  add_field<Required>("qc",               scalar3d_layout_mid, Q,      grid_name, "tracers");
  add_field<Required>("qi",               scalar3d_layout_mid, Q,      grid_name, "tracers");
  add_field<Required>("cldfrac_rad",      scalar3d_layout_mid, nondim, grid_name);
  // Optical properties, should be computed in radiation interface
  add_field<Required>("dtau067",     scalar3d_layout_mid, nondim, grid_name); // 0.67 micron optical depth
  add_field<Required>("dtau105",     scalar3d_layout_mid, nondim, grid_name); // 10.5 micron optical depth
  // Effective radii, should be computed in either microphysics or radiation interface
  // TODO: should these be meters or microns? Was meters before, but using "m" instead
  // of "micron" seemed to cause prim_model_finalize to throw error with the following:
  // ABORTING WITH ERROR: Error! prim_init_model_f90 was not called yet (or prim_finalize_f90 was already called).
  // P3 defines this field with micron instead of meters units, so is this a unit conversion issue?
  add_field<Required>("eff_radius_qc",     scalar3d_layout_mid, micron,      grid_name);
  add_field<Required>("eff_radius_qi",     scalar3d_layout_mid, micron,      grid_name);
  // Set of fields used strictly as output
  add_field<Computed>("isccp_cldtot", scalar2d_layout, percent, grid_name);
  add_field<Computed>("isccp_ctptau", scalar4d_layout_ctptau, percent, grid_name, 1);
  add_field<Computed>("isccp_mask"  , scalar2d_layout, nondim, grid_name);

}

// =========================================================================================
void Cosp::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  CospFunc::initialize(m_num_cols, m_num_subcols, m_num_levs);


  // Add note to output files about processing ISCCP fields that are only valid during
  // daytime. This can go away once I/O can handle masked time averages.
  using stratts_t = std::map<std::string,std::string>;
  std::list<std::string> vnames = {"isccp_cldtot", "isccp_ctptau"};
  for (const auto field_name : {"isccp_cldtot", "isccp_ctptau"}) {
      auto& f = get_field_out(field_name);
      auto& atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
      atts["note"] = "Night values are zero; divide by isccp_mask to get daytime mean";
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
  auto ts = timestamp();
  auto update_cosp = cosp_do(cosp_freq_in_steps, ts.get_num_steps());

  // Get fields from field manager; note that we get host views because this
  // interface serves primarily as a wrapper to a c++ to f90 bridge for the COSP
  // all then need to be copied to layoutLeft views to permute the indices for
  // F90.
  //
  // Need to make sure device data is synced to host?
  get_field_in("qv").sync_to_host();
  get_field_in("qc").sync_to_host();
  get_field_in("qi").sync_to_host();
  get_field_in("sunlit").sync_to_host();
  get_field_in("surf_radiative_T").sync_to_host();
  get_field_in("T_mid").sync_to_host();
  get_field_in("p_mid").sync_to_host();
  get_field_in("p_int").sync_to_host();
  get_field_in("cldfrac_rad").sync_to_host();
  get_field_in("eff_radius_qc").sync_to_host();
  get_field_in("eff_radius_qi").sync_to_host();
  get_field_in("dtau067").sync_to_host();
  get_field_in("dtau105").sync_to_host();

  auto qv      = get_field_in("qv").get_view<const Real**, Host>();
  auto qc      = get_field_in("qc").get_view<const Real**, Host>();
  auto qi      = get_field_in("qi").get_view<const Real**, Host>();
  auto sunlit  = get_field_in("sunlit").get_view<const Real*, Host>();
  auto skt     = get_field_in("surf_radiative_T").get_view<const Real*, Host>();
  auto T_mid   = get_field_in("T_mid").get_view<const Real**, Host>();
  auto p_mid   = get_field_in("p_mid").get_view<const Real**, Host>();
  auto p_int   = get_field_in("p_int").get_view<const Real**, Host>();
  auto cldfrac = get_field_in("cldfrac_rad").get_view<const Real**, Host>();
  auto reff_qc = get_field_in("eff_radius_qc").get_view<const Real**, Host>();
  auto reff_qi = get_field_in("eff_radius_qi").get_view<const Real**, Host>();
  auto dtau067 = get_field_in("dtau067").get_view<const Real**, Host>();
  auto dtau105 = get_field_in("dtau105").get_view<const Real**, Host>();
  auto isccp_cldtot = get_field_out("isccp_cldtot").get_view<Real*, Host>();
  auto isccp_ctptau = get_field_out("isccp_ctptau").get_view<Real***, Host>();
  auto isccp_mask   = get_field_out("isccp_mask"  ).get_view<Real*, Host>();  // Copy of sunlit flag with COSP frequency for proper averaging

  // Call COSP wrapper routines
  if (update_cosp) {
    Real emsfc_lw = 0.99;
    Kokkos::deep_copy(isccp_mask, sunlit);
    CospFunc::main(
            m_num_cols, m_num_subcols, m_num_levs, m_num_isccptau, m_num_isccpctp,
            emsfc_lw, sunlit, skt, T_mid, p_mid, p_int, qv,
            cldfrac, reff_qc, reff_qi, dtau067, dtau105,
            isccp_cldtot, isccp_ctptau
    );
    // Remask night values to ZERO since our I/O does not know how to handle masked/missing values
    // in temporal averages; this is all host data, so we can just use host loops like its the 1980s
    for (int i = 0; i < m_num_cols; i++) {
        if (sunlit(i) == 0) {
            isccp_cldtot(i) = 0;
            for (int j = 0; j < m_num_isccptau; j++) {
                for (int k = 0; k < m_num_isccpctp; k++) {
                    isccp_ctptau(i,j,k) = 0;
                }
            }
        }
    }
  } else {
    // If not updating COSP statistics, set these to ZERO; this essentially weights
    // the ISCCP cloud properties by the sunlit mask. What will be output for time-averages
    // then is the time-average mask-weighted statistics; to get true averages, we need to
    // divide by the time-average of the mask. I.e., if M is the sunlit mask, and X is the ISCCP
    // statistic, then
    //
    //     avg(X) = sum(M * X) / sum(M) = (sum(M * X)/N) / (sum(M)/N) = avg(M * X) / avg(M)
    //
    // TODO: mask this when/if the AD ever supports masked averages
    Kokkos::deep_copy(isccp_cldtot, 0.0);
    Kokkos::deep_copy(isccp_ctptau, 0.0);
    Kokkos::deep_copy(isccp_mask  , 0.0);
  }
  get_field_out("isccp_cldtot").sync_to_dev();
  get_field_out("isccp_ctptau").sync_to_dev();
  get_field_out("isccp_mask"  ).sync_to_dev();
}

// =========================================================================================
void Cosp::finalize_impl()
{
  // Finalize COSP wrappers
  CospFunc::finalize();
}
// =========================================================================================

} // namespace scream
