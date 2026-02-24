#include "eamxx_cld_fraction_process_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "physics/cld_fraction/cld_fraction_functions.hpp"

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#endif

#include <array>

namespace scream
{
  using namespace cld_fraction;
// =========================================================================================
CldFraction::CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void CldFraction::create_requests()
{
  using namespace ekat::units;
  using CldFractionFunc = cld_fraction::CldFractionFunctions<Real, DefaultDevice>;
  using Spack           = CldFractionFunc::Spack;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto nondim = Units::nondimensional();

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  auto scalar3d_layout_mid = m_grid->get_3d_scalar_layout(true);

  // Set of fields used strictly as input
  constexpr int ps = Spack::n;
  add_tracer<Required>("qi", m_grid, kg/kg, ps);
  add_field<Required>("cldfrac_liq", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used strictly as output
  add_field<Computed>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name,ps);
  add_field<Computed>("cldfrac_ice", scalar3d_layout_mid, nondim, grid_name,ps);
  // Note, we track two versions of the cloud fraction.  The versions below have "_for_analysis"
  // attached to the name because they're meant for use with fields that are exclusively
  // related to writing output.  This is an important distinction here because the internal ice
  // cloud fraction needs to be 100% whenever any ice at all is present in the cell (in order
  // for the model's ice processes to act on that cell). Folks evaluating cloud, on the other hand,
  // expect cloud fraction to represent cloud visible to the human eye (which corresponds to
  // ~1e-5 kg/kg).
  add_field<Computed>("cldfrac_tot_for_analysis", scalar3d_layout_mid, nondim, grid_name,ps);
  add_field<Computed>("cldfrac_ice_for_analysis", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.

  // Gather parameters for ice cloud thresholds from parameter list:
  m_icecloud_threshold = m_params.get<double>("ice_cloud_threshold",1e-12);  // Default = 1e-12
  m_icecloud_for_analysis_threshold = m_params.get<double>("ice_cloud_for_analysis_threshold",1e-5); // Default = 1e-5
}

// =========================================================================================
void CldFraction::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("cldfrac_ice"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_tot"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_ice_for_analysis"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_tot_for_analysis"),m_grid,0.0,1.0,false);
#ifdef EAMXX_HAS_PYTHON
  if (has_py_module()) {
    py_module_call("init");
  }
#endif
}

// =========================================================================================
void CldFraction::run_impl (const double /* dt */)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio.
  auto qi   = get_field_in("qi");
  auto liq_cld_frac = get_field_in("cldfrac_liq");
  auto ice_cld_frac = get_field_out("cldfrac_ice");
  auto tot_cld_frac = get_field_out("cldfrac_tot");
  auto ice_cld_frac_4out = get_field_out("cldfrac_ice_for_analysis");
  auto tot_cld_frac_4out = get_field_out("cldfrac_tot_for_analysis");

#ifdef EAMXX_HAS_PYTHON
  if (has_py_module()) {
    pybind11::array py_qi,
                    py_liq_cld_frac,
                    py_ice_cld_frac,
                    py_tot_cld_frac,
                    py_ice_cld_frac_4out,
                    py_tot_cld_frac_4out;

    if (m_params.get<std::string>("py_backend")=="device") {
      py_qi                = get_py_field_dev("qi");
      py_liq_cld_frac      = get_py_field_dev("cldfrac_liq");
      py_ice_cld_frac      = get_py_field_dev("cldfrac_ice");
      py_tot_cld_frac      = get_py_field_dev("cldfrac_tot");
      py_ice_cld_frac_4out = get_py_field_dev("cldfrac_ice_for_analysis");
      py_tot_cld_frac_4out = get_py_field_dev("cldfrac_tot_for_analysis");
    } else {
      qi.sync_to_host();
      liq_cld_frac.sync_to_host();
      py_qi                = get_py_field_host("qi");
      py_liq_cld_frac      = get_py_field_host("cldfrac_liq");
      py_ice_cld_frac      = get_py_field_host("cldfrac_ice");
      py_tot_cld_frac      = get_py_field_host("cldfrac_tot");
      py_ice_cld_frac_4out = get_py_field_host("cldfrac_ice_for_analysis");
      py_tot_cld_frac_4out = get_py_field_host("cldfrac_tot_for_analysis");
    }

    double ice_threshold      = m_params.get<double>("ice_cloud_threshold");
    double ice_4out_threshold = m_params.get<double>("ice_cloud_for_analysis_threshold");

    py_module_call("main",
                   ice_threshold,ice_4out_threshold,
                   py_qi,py_liq_cld_frac,
                   py_ice_cld_frac,py_tot_cld_frac,
                   py_ice_cld_frac_4out,py_tot_cld_frac_4out);

    if (m_params.get<std::string>("py_backend")=="host") {
      ice_cld_frac.sync_to_dev();
      tot_cld_frac.sync_to_dev();
      ice_cld_frac_4out.sync_to_dev();
      tot_cld_frac_4out.sync_to_dev();
    }
    return;
  }
#endif

  using CldFractionFunc = cld_fraction::CldFractionFunctions<Real, DefaultDevice>;
  using Spack = CldFractionFunc::Spack;

  auto qi_v                = qi.get_view<const Spack**>();
  auto liq_cld_frac_v      = liq_cld_frac.get_view<const Spack**>();
  auto ice_cld_frac_v      = ice_cld_frac.get_view<Spack**>();
  auto tot_cld_frac_v      = tot_cld_frac.get_view<Spack**>();
  auto ice_cld_frac_4out_v = ice_cld_frac_4out.get_view<Spack**>();
  auto tot_cld_frac_4out_v = tot_cld_frac_4out.get_view<Spack**>();

  CldFractionFunc::main(m_num_cols,m_num_levs,m_icecloud_threshold,m_icecloud_for_analysis_threshold,
    qi_v,liq_cld_frac_v,ice_cld_frac_v,tot_cld_frac_v,ice_cld_frac_4out_v,tot_cld_frac_4out_v);
}

// =========================================================================================
void CldFraction::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

} // namespace scream
