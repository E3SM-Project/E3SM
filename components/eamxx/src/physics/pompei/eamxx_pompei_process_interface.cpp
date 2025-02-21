/*
 * Add any #include statements pointing to headers needed by this process.
 * Typical headers would be
 *     eamxx_$PROCESS_process_interface.hpp  (the eamxx interface header file for this process)
 *     physics/share/physics_constants.hpp   (a header storing all physical constants in eamxx)
 *     $PROCESS.hpp                          (a header file for the process itself)
 */

#include "eamxx_pompei_process_interface.hpp"

#include "physics/share/physics_constants.hpp"
#include "pompei.hpp"

namespace scream
{

/*-----------------------------------------------------------------------------------------------
 * The Constructor for this interface
 *
 * Inputs:
 *     comm - an EKAT communication group
 *     params - a parameter list of options for the process.
 */

POMPEI::POMPEI (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereProcess(comm,params)
{
  // The 'params' variable will hold all runtime options.  
  // Note that `params.get<X>("Y") is the syntax to get the parameter labeled "Y" from the list
  // which is of type X. 
  m_eruption_start = util::str_to_time_stamp(params.get<std::string>("eruption_date"));
}

/*-----------------------------------------------------------------------------------------------
 * set_grids(grids_manager)
 *
 * set_grids establishes which grid this process will be run on.  
 * It is also where all fields (variables) that the process will need access to.
 *
 * Inputs:
 *     grids_manager - a grids manager object which stores all grids used in this simulation
 */
void POMPEI::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  // Some typical namespaces used in set_grids,
  //   using namespace ShortFieldTagsNames;
  //   using PC = scream::physics::Constants<Real>;

  // Specify which grid this process will act upon, typical options are "Dynamics" or "Physics".
  m_grid = grids_manager->get_grid("Physics");
  m_ncols = m_grid->get_num_local_dofs();
  m_nlevs = m_grid->get_num_vertical_levels();

  // When declaring a field for use we need:
  //
  // to set the units.  These can be found from ekat::units::Units
  // For example,
  constexpr auto Pa = ekat::units::Pa;
  constexpr auto kg = ekat::units::kg;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  // 
  // to set the layout, which is how the field is configured in space
  // For example,
  const auto layout_2d = m_grid->get_2d_scalar_layout();
  const auto layout_3d = m_grid->get_3d_scalar_layout(true);

  // We can now add all fields used by this process and designate how they will be used,
  // The format for declaring a field is:
  //   add_field<MODE>(NAME,LAYOUT,UNITS,GRID_NAME);
  // where,
  //   MODE is either "Required", "Computed" or "Updated".
  //     "Required" = input only  (process cannot change this field)
  //     "Computed" = output only (field doesn't have to exist, process will set values for this field)
  //     "Updated"  = input and output (process is able to update field which should already exist)
  //   NAME is the name of the field as a string
  //   LAYOUT is the physical layout of the field on the grid, e.g. 2D, 3D.  See layout example above.
  //   UNITS are the units of the field, see units example above
  //   GRID_NAME is the name of the grid this process is run on, typically "m_grid->name()"
  //   *Note, there is a special add_field call for tracers which has a different signature.  The format is,
  //     add_tracer<MODE>(NAME,GRID,UNITS)
  //   where,
  //     NAME and UNITS are as above
  //     GRID is the actual grid, typically "m_grid"
  
  // Required Fields
  add_field<Required>("pseudo_density",layout_3d,Pa,m_grid->name());
  // Computed Fields
  add_field<Computed>("ash_column_integrated",layout_2d,nondim,m_grid->name());
  // Updated Fields
  // Tracer Fields (advected by dynamics)
  add_tracer<Updated>("ash",m_grid,kg/kg);

}

/*-----------------------------------------------------------------------------------------------
 * intialize_impl(run_type)
 *
 * called once for each process at initialization of EAMxx.  This impl can be defined with any
 * actions or functions that are needed at initialization. 
 *
 * Inputs:
 *     run_type - an enum which describes the run type.  Initial or Restart
 *
 * can also be empty
 */
void POMPEI::initialize_impl (const RunType /* run_type */)
{
  // Define the a mask to control where ash is injected into the atmosphere.
  
  // Some universal constants that we will need
  using PC = scream::physics::Constants<Real>;
  constexpr Real deg2rad = PC::Pi / 180;
  constexpr Real r_earth = PC::r_earth;

  // We can create our helper "mask" field. It will be 1 where the volcano injection
  // in the atmosphere happens, and 0 elsewhere.
  const std::string mask_name = "emission_mask";
  constexpr auto nondim = ekat::units::Units::nondimensional();
  const auto layout_3d = m_grid->get_3d_scalar_layout(true);
  FieldIdentifier rate_fid(mask_name,layout_3d,nondim,m_grid->name());
  m_emission_mask = Field(rate_fid);
  m_emission_mask.allocate_view();
  m_emission_mask.deep_copy(0); // 0 means "not injecting here"

  // The grid already stores lat/lon, so go ahead and pre-compute the mask field
  auto lat = m_grid->get_geometry_data("lat"); // WATCH OUT: it's in degrees, NOT rad
  auto lon = m_grid->get_geometry_data("lon"); // WATCH OUT: it's in degrees, NOT rad

  // Target location of volcanic eruption
  auto volcano_lat = 40.8214 * deg2rad;
  auto volcano_lon = 14.4260 * deg2rad;
  auto emission_lev = 27;
  auto radius = m_params.get<double>("plume_radius_in_km");
  EKAT_REQUIRE_MSG (radius>0,
      "Error! Plume radius should be positive. Input value: " << radius << ".\n");

  // Extract Kokkos (device) views from the fields
  // NOTE: views are "just" multi-dimensional arrays, which are accessible on device
  auto emission_view = m_emission_mask.get_view<Real**>();
  auto lat_view = lat.get_view<const Real*>();
  auto lon_view = lon.get_view<const Real*>();

  // ------------------------------------------------------------------------------------------- //
  // We define a lambda function for computing the mask that we can then use with
  // Kokkos.
  // A lambda is just a function defined "on-the fly". KOKKOS_LAMBDA simply adds
  // some decoration for GPU execution, nothing to worry about for now
  auto compute_mask = KOKKOS_LAMBDA (const int icol) {
    auto lat_rad = lat_view(icol) * deg2rad;
    auto lon_rad = lon_view(icol) * deg2rad;
    auto delta_lat = lat_rad - volcano_lat;
    auto delta_lon = lon_rad - volcano_lon;

    auto dist = r_earth/1e3 * sqrt(delta_lat*delta_lat + delta_lon*delta_lon);

    // If the distance between this point and the center of the volcano is within the radius
    // set the MASK value to 1 (true).
    if (dist<radius) {
      emission_view(icol,emission_lev) = 1;
    }
  };
  // ------------------------------------------------------------------------------------------- //

  // A policy tells kokkos how to parallelize the loop. Here, we are doing
  // a single for loop over the range of indices [0,ncols)
  auto policy = Kokkos::RangePolicy<Field::device_t::execution_space>(0,m_ncols);

  // Execute the lambda in parallel according to the execution policy
  Kokkos::parallel_for("",policy,compute_mask);
}

/*-----------------------------------------------------------------------------------------------
 * run_impl(dt)
 *
 * The main run call for the process.  This is where most of the interface to the underlying process
 * takes place.  This impl is called every timestep. 
 *
 * Inputs:
 *     dt - the timestep for this run step.
 */
void POMPEI::run_impl (const double dt)
{
  using namespace pompei;

  auto t = timestamp() + dt;
  auto rate = ash_emission_rate (t.days_from(m_eruption_start));
  auto mass = dt*rate;

  auto qash = get_field_out("ash");
  auto rho = get_field_in("pseudo_density");

  // y.update(x,a,b) means y = b*y + a*x
  qash.scale(rho);
  qash.update(m_emission_mask,mass,1.0);
  qash.scale_inv(rho);

}

/*-----------------------------------------------------------------------------------------------
 * finalize_impl()
 *
 * Called at the end of the simulation, handles all finalization of the process.
 *
 * In most cases this is left blank, as EAMxx takes care of most finalization steps.
 * But just in case a process has specific needs the option is available.
 */
void POMPEI::finalize_impl ()
{
  printf("pompei has finished\n");
}
/*-----------------------------------------------------------------------------------------------*/

} // namespace scream
