#include "eamxx_cld_frac_net_process_interface.hpp"

#include "share/atm_process/atmosphere_process_pyhelpers.hpp"

namespace scream
{

CldFracNet::CldFracNet (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  EKAT_REQUIRE_MSG( has_py_module(),
      "[CldFracNet] Error! Something went wrong while initializing the python module.\n");
  // Nothing to do here
}

void CldFracNet::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  const auto nondim = Units::nondimensional();
  const auto grid = grids_manager->get_grid("physics");
  const auto grid_name = grid->name();
  const auto layout = grid->get_3d_scalar_layout(true);

  // Input fields
  add_tracer<Required>("qi", grid, kg/kg);
  add_field<Required>("cldfrac_liq", layout, nondim, grid_name);

  // Output fields
  add_field<Computed>("cldfrac_tot", layout, nondim, grid_name);
  add_field<Computed>("cldfrac_ice", layout, nondim, grid_name);
}

void CldFracNet::initialize_impl (const RunType /* run_type */)
{
  py_module_call("init");
}

void CldFracNet::run_impl (const double /* dt */)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio.
  auto qi  = get_field_in("qi");
  auto liq = get_field_in("cldfrac_liq");
  auto ice = get_field_out("cldfrac_ice");
  auto tot = get_field_out("cldfrac_tot");

  } else {
  auto py_qi  = get_py_field_dev("qi");
  auto py_liq = get_py_field_dev("cldfrac_liq");
  auto py_ice = get_py_field_dev("cldfrac_ice");
  auto py_tot = get_py_field_dev("cldfrac_tot");

  double ice_threshold = m_params.get<double>("ice_cloud_threshold");

  py_module_call("forward",ice_threshold,py_qi,py_liq,py_ice,py_tot);
  }
}

void CldFracNet::finalize_impl()
{
  // Nothing to do here
}

} // namespace scream
