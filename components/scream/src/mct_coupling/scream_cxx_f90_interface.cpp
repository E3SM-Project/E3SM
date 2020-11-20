#include "scream_config.h"

#include "share/atm_process/atmosphere_process.hpp"
#include "control/atmosphere_driver.hpp"

#include "dynamics/register_dynamics.hpp"
#include "physics/register_physics.hpp"

#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/shoc/scream_shoc_interface.hpp"

#include "control/tests/dummy_grid.hpp"

#include "mct_coupling/ScreamContext.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/scream_session.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_assert.hpp"

using scream::Real;

// Anonymous namespace, for some utility functions
namespace {
// utility function, that wraps a lambda, and disable/restores
// fpes at entry/exit.
template<typename Lambda>
void fpe_guard_wrapper (const Lambda& f) {
  using namespace scream;

  // Disable all fpes we may have enabled, then enable scream default FPEs.
  // Store the mask, so we can restore before returning.
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes ();
  ekat::enable_default_fpes();

  // Execute wrapped function
  f();

  // Restore the FPE flag as it was when control was handed to us.
  ekat::disable_all_fpes();
  ekat::enable_fpes(fpe_mask);
}

scream::control::AtmosphereDriver& get_ad_nonconst () {
  auto& c = scream::ScreamContext::singleton();
  return c.getNonConst<scream::control::AtmosphereDriver>();
}

const scream::control::AtmosphereDriver& get_ad () {
  const auto& c = scream::ScreamContext::singleton();
  return c.get<scream::control::AtmosphereDriver>();
}

} // anonymous namespace

extern "C"
{

/*===============================================================================================*/
// WARNING: make sure input_yaml_file is a null-terminated string!
void scream_init (const MPI_Fint& f_comm,
                  const int& start_ymd,
                  const int& start_tod,
                  const char* input_yaml_file) {
                  // const int& compid) {
  using namespace scream;
  using namespace scream::control;

  fpe_guard_wrapper([&](){
    // First of all, initialize the scream session
    scream::initialize_scream_session();

    // Create the context
    auto& c = ScreamContext::singleton();

    // Create the C MPI_Comm from the Fortran one
    MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
    auto& atm_comm = c.create<ekat::Comm>(mpi_comm_c);


    // Create a parameter list for inputs
    printf("[scream] reading parameterr from yaml file: %s\n",input_yaml_file);
    ekat::ParameterList scream_params("Scream Parameters");
    parse_yaml_file (input_yaml_file, scream_params);
    scream_params.print();

    ekat::error::runtime_check(scream_params.isSublist("Atmosphere Driver"),
         "Error! Sublist 'Atmosphere Driver' not found inside '" +
         std::string(input_yaml_file) + "'.\n");

    auto& ad_params = scream_params.sublist("Atmosphere Driver");

    // Need to register products in the factories *before* we attempt to create any.
    // In particular, register all atm processes, and all grids managers.
    register_dynamics();
    register_physics();

    // Create the bare ad, then init it
    auto& ad = c.create<AtmosphereDriver>();

    // Recall that e3sm uses the int YYYYMMDD to store a date
    std::cout << "start_ymd: " << start_ymd << "\n";
    const int dd = start_ymd % 100;
    const int mm = (start_ymd / 100) % 100;
    const int yy = start_ymd / 10000;
    util::TimeStamp time (yy,mm,dd,start_tod);

    // Init and run (to finalize, wait till checks are completed,
    // or you'll clear the field repo!)
    ad.initialize(atm_comm,ad_params,time);
  });
}

void scream_setup_surface_coupling (
    const char*& x2a_names, const int& x2a_indices, double*& cpl_x2a_ptr, const int& num_imports,
    const char*& a2x_names, const int& a2x_indices, double*& cpl_a2x_ptr, const int& num_exports)
{
  fpe_guard_wrapper([&](){
    // Get the SurfaceCoupling from the AD, then register the pointer
    const auto& ad = get_ad();
    const auto& sc = ad.get_surface_coupling();
  });
}

/*===============================================================================================*/
void scream_run (const Real& dt) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  fpe_guard_wrapper([&](){
    // Get the AD, and run it
    auto& ad = get_ad_nonconst();
    ad.run(dt);
  });
}
/*===============================================================================================*/
void scream_finalize (/* args ? */) {
  fpe_guard_wrapper([&](){
    // Get the AD, and finalize it
    auto& ad = get_ad_nonconst();
    ad.finalize();

    // Clean up also P3 stuff
    scream::p3::P3GlobalForFortran::deinit();
  });
}

// Get the local (i.e., on current atm rank only) number of physics columns
int scream_get_num_local_cols () {
  int ncols;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& ref_grid = gm->get_reference_grid();

    ncols = ref_grid->get_num_local_dofs();
  });

  return ncols;
}

// Get the global (i.e., the whole earth) number of physics columns
int scream_get_num_global_cols () {
  int ncols;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& ref_grid = gm->get_reference_grid();

    ncols = ref_grid->get_num_global_dofs();
  });

  return ncols;
}

// Return the global ids of all physics column
void scream_get_local_cols_gids (void* const ptr) {
  fpe_guard_wrapper([&]() {
    auto gids_f = reinterpret_cast<int* const>(ptr);
    const auto& ad = get_ad();
    const auto& grid = ad.get_grids_manager()->get_reference_grid();

    auto gids_h = Kokkos::create_mirror_view(grid->get_dofs_gids());
    Kokkos::deep_copy(gids_h,grid->get_dofs_gids());

    for (int i=0; i<gids_h.extent_int(0); ++i) {
      gids_f[i] = gids_h(i);
    }
  });
}

// Retrieve the lat/lon coords of all physics FV dofs
void scream_get_cols_latlon (double* const& lat_ptr, double* const& lon_ptr) {
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& grid = ad.get_grids_manager()->get_reference_grid();
    const auto ncols = grid->get_num_local_dofs();

    auto lat_cxx = grid->get_geometry_data("lat");
    auto lon_cxx = grid->get_geometry_data("lon");

    using geo_view_f90 = ekat::Unmanaged<decltype(lat_cxx)::HostMirror>;
    geo_view_f90 lat_f90(lat_ptr, ncols);
    geo_view_f90 lon_f90(lon_ptr, ncols);
    Kokkos::deep_copy(lat_f90, lat_cxx);
    Kokkos::deep_copy(lon_f90, lon_cxx);
  });
}

// Retrieve the area of all physics FV dofs
void scream_get_cols_area (double* const& area_ptr) {
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& grid = ad.get_grids_manager()->get_reference_grid();
    const auto ncols = grid->get_num_local_dofs();

    auto area_cxx = grid->get_geometry_data("area");

    using geo_view_f90 = ekat::Unmanaged<decltype(area_cxx)::HostMirror>;
    geo_view_f90  area_f90 (area_ptr, ncols);
    Kokkos::deep_copy(area_f90,area_cxx);
  });
}

} // extern "C"
