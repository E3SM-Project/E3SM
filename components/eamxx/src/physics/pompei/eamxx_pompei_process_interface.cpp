#include "eamxx_pompei_process_interface.hpp"

#include "physics/share/physics_constants.hpp"
#include "pompei.hpp"

namespace scream {

/*-----------------------------------------------------------------------------------------------*/
PompeiEruption::PompeiEruption(const ekat::Comm &comm,
                               const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  // The 'params' variable will hold all runtime options.  In this case we
  // retrieve the parameter "eruption_date" from the parameter list and assign
  // it to a time stamp object. Note that `params.get<X>("Y") is the syntax to
  // get the parameter labeled "y" from the list which is of type X.
  m_eruption_start =
      util::str_to_time_stamp(params.get<std::string>("eruption_date"));
}

/*-----------------------------------------------------------------------------------------------*/
void PompeiEruption::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ShortFieldTagsNames;

  // Some constants
  using PC               = scream::physics::Constants<Real>;
  constexpr Real deg2rad = PC::Pi / 180;
  constexpr Real r_earth = PC::r_earth;

  constexpr auto Pa     = ekat::units::Pa;
  constexpr auto kg     = ekat::units::kg;
  constexpr auto s      = ekat::units::s;
  constexpr auto nondim = ekat::units::Units::nondimensional();

  // Specify which grid this process will act upon, typical options are
  // "Dynamics" or "Physics".
  auto m_grid       = grids_manager->get_grid("Physics");
  const auto layout = m_grid->get_3d_scalar_layout(true);

  // Declare the fields we need as input and/or output

  // Declare the field "ash" which will be updated each time step by this
  // process.
  add_tracer<Updated>("ash", m_grid, kg / kg);

  // Needed to compute tracer mix ratio: mass_of_ash/mass_of_air
  add_field<Required>("pseudo_density", layout, Pa, m_grid->name());

  // Number of columns/levels on this MPI rank
  int ncols = layout.dim(COL);
  int nlevs = layout.dim(LEV);

  // We can create our helper "mask" field. It will be 1 where the volcano
  // injection in the atmosphere happens, and 0 elsewhere.
  const std::string mask_name = "emission_mask";
  FieldIdentifier rate_fid(mask_name, layout, nondim, m_grid->name());
  m_emission_mask = Field(rate_fid);
  m_emission_mask.allocate_view();
  m_emission_mask.deep_copy(0);  // 0 means "not injecting here"

  // The grid already stores lat/lon, so go ahead and pre-compute the mask field
  auto lat =
      m_grid->get_geometry_data("lat");  // WATCH OUT: it's in degrees, NOT rad
  auto lon =
      m_grid->get_geometry_data("lon");  // WATCH OUT: it's in degrees, NOT rad

  // Target location of volcanic eruption
  auto volcano_lat  = 40.8214 * deg2rad;
  auto volcano_lon  = 14.4260 * deg2rad;
  auto emission_lev = 27;
  auto radius       = m_params.get<double>("plume_radius_in_km");
  EKAT_REQUIRE_MSG(radius > 0,
                   "Error! Plume radius should be positive. Input value: "
                       << radius << ".\n");

  // Extract Kokkos (device) views from the fields
  // NOTE: views are "just" multi-dimensional arrays, which are accessible on
  // device
  auto emission_view = m_emission_mask.get_view<Real **>();
  auto lat_view      = lat.get_view<const Real *>();
  auto lon_view      = lon.get_view<const Real *>();

  // A lambda is just a function defined "on-the fly". KOKKOS_LAMBDA simply adds
  // some decoration for GPU execution, nothing to worry for now
  auto compute_mask = KOKKOS_LAMBDA(const int icol) {
    auto lat_rad   = lat_view(icol) * deg2rad;
    auto lon_rad   = lon_view(icol) * deg2rad;
    auto delta_lat = lat_rad - volcano_lat;
    auto delta_lon = lon_rad - volcano_lon;

    auto dist =
        r_earth / 1e3 * sqrt(delta_lat * delta_lat + delta_lon * delta_lon);

    // If the distance between this point and the center of the volcano is
    // within the radius set the MASK value to 1 (true).
    if(dist < radius) {
      emission_view(icol, emission_lev) = 1;
    }
  };

  // A policy tells kokkos how to parallelize the loop. Here, we are doing
  // a single for loop over the range of indices [0,ncols)
  auto policy = Kokkos::RangePolicy<Field::device_t::execution_space>(0, ncols);

  // Execute the lambda in parallel according to the execution policy
  Kokkos::parallel_for("", policy, compute_mask);
}

/*-----------------------------------------------------------------------------------------------*/
void PompeiEruption::initialize_impl(const RunType /* run_type */) {
  // If pompei internally requires to initialize data, do that now.
  // NOTE: run_type tells us if this is an initial or restarted run,
  //       but this parametrization doesn't care
}

/*-----------------------------------------------------------------------------------------------*/
void PompeiEruption::run_impl(const double dt) {
  using namespace pompei;
  // Compute current emission rate and added mass
  // timestamp returns time at the *beginning* of the atm step.
  auto t    = timestamp() + dt;
  auto rate = ash_emission_rate(t.days_from(m_eruption_start));
  auto mass = dt * rate;

  // Update the output field: qash = (qash*rho + dt*injection_rate)/rho
  auto qash = get_field_out("ash");
  auto rho  = get_field_in("pseudo_density");

  // y.update(x,a,b) means y = b*y + a*x
  qash.scale(rho);
  qash.update(m_emission_mask, mass, 1.0);
  qash.scale_inv(rho);
}

/*-----------------------------------------------------------------------------------------------*/
void PompeiEruption::finalize_impl() {
  // If pompei internally requires to cleanup some data, do that now.
}
/*-----------------------------------------------------------------------------------------------*/

}  // namespace scream
