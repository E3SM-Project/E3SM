#include "eamxx_config.h"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_zm_process_interface.hpp"
#include "physics/share/physics_constants.hpp"

#include <ekat_assert.hpp>

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
/* Constructor for the ZMDeepConvection interface
* Inputs:
*     comm - an EKAT communication group
*     params - a parameter list of options for the process.
*/
ZMDeepConvection::ZMDeepConvection (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereProcess(comm,params)
{
  // params holds all runtime options - what do we need for ZM?
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // Specify which grid this process will act upon, typical options are "Dynamics" or "Physics".
  auto m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  const auto layout = m_grid->get_3d_scalar_layout(true);

  // retrieve local grid parameters
  m_ncols = m_grid->get_num_local_dofs();
  m_nlevs = m_grid->get_num_vertical_levels();

  constexpr int ps = Spack::n;

  const auto nondim = Units::nondimensional();
  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);

  FieldLayout scalar2d = m_grid->get_2d_scalar_layout();            // Layout for 2D variable
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);    // Layout for 3D variable at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);   // Layout for 3D variable at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);  // Layout for horiz_wind field

  // Input variables
  add_field<Required>("p_mid",          scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("p_int",          scalar3d_int, Pa,    grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("phis",           scalar2d    , m2/s2, grid_name, ps);
  add_field<Required>("omega",          scalar3d_mid, Pa/s,  grid_name, ps);

  // Input/Output variables
  add_field <Updated> ("T_mid",         scalar3d_mid, K,     grid_name, ps);
  add_field <Updated>("horiz_winds",    vector3d_mid, m/s,   grid_name, ps);
  add_tracer<Updated>("qv",             m_grid,       kg/kg,            ps);
  add_tracer<Updated>("qc",             m_grid,       kg/kg,            ps);

  // // Output variables
  // add_field<Computed>("???", scalar2d    , ???, grid_name);
  // add_field<Computed>("???", scalar3d_mid, ???, grid_name, ps);
  
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::initialize_impl (const RunType /* run_type */)
{
  // placeholder for initialization
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::run_impl (const double dt)
{

  // get fields
  const auto& T_mid    = get_field_out("T_mid").get_view<Spack**>();
  const auto& p_mid    = get_field_in("p_mid").get_view<const Spack**>();
  const auto& p_int    = get_field_in("p_int").get_view<const Spack**>();
  const auto& rho      = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega    = get_field_in("omega").get_view<const Spack**>();
  const auto& qc       = get_field_out("qc").get_view<Spack**>();
  const auto& qv       = get_field_out("qv").get_view<Spack**>();
  const auto& phis     = get_field_in("phis").get_view<const Real*>();
  
  // Run ZM
  
  // Update output fields
  
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::finalize_impl ()
{
  // placeholder for final cleanup
}

/*------------------------------------------------------------------------------------------------*/

} // namespace scream
