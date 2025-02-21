#include "eamxx_config.h"

#include "share/atm_process/atmosphere_process.hpp"
#include "control/atmosphere_driver.hpp"
#include "control/surface_coupling_utils.hpp"

#include "dynamics/register_dynamics.hpp"
#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "control/register_surface_coupling.hpp"

#include "mct_coupling/ScreamContext.hpp"
#include "share/grid/point_grid.hpp"
#include "share/eamxx_session.hpp"
#include "share/eamxx_config.hpp"
#include "share/eamxx_types.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/logging/ekat_logger.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_assert.hpp"

#if defined(SCREAM_SYSTEM_WORKAROUND) && (SCREAM_SYSTEM_WORKAROUND == 1)
#include <hip/hip_runtime.h>
#endif

// Anonymous namespace, for some utility functions
namespace {

// Utility function, that wraps a lambda, and disable/restores
// fpes at entry/exit.
// When a function is callable from Fortran, its body should look like this:
//
// void foo ([args...]) {
//  fpe_guard_wrapper ([&]{
//    <actual body of foo here>
//  };
// }
//
// The fpe_guard_wrapper will store the current fpe mask, enable scream's
// default fpes, run the body of the function, the restore the original
// fpe mast before returning, ensuring both scream and the cpl use their
// own supported version of FPE mask.

template<typename Lambda>
void fpe_guard_wrapper (const Lambda& f) {
  using namespace scream;

  // Disable all fpes we may have enabled, then enable scream default FPEs.
  // Store the mask, so we can restore before returning.
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes ();
  ekat::enable_fpes(get_default_fpes());

  // Execute wrapped function
  try {
    f();
  } catch (std::exception &e) {
    // Print exception msg, then call MPI_Abort
    fprintf(stderr, "%s\n", e.what());

    // Get raw comm before cleaning up singleton
    auto& c = ScreamContext::singleton();
    auto raw_comm = c.get<ekat::Comm>().mpi_comm();
    c.clean_up();

    MPI_Abort (raw_comm,1);
  } catch (...) {

    // Get raw comm before cleaning up singleton
    auto& c = ScreamContext::singleton();
    auto raw_comm = c.get<ekat::Comm>().mpi_comm();
    c.clean_up();

    MPI_Abort (raw_comm,1);
  }

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
void scream_create_atm_instance (const MPI_Fint f_comm, const int atm_id,
                                 const char* input_yaml_file,
                                 const char* atm_log_file,
                                 const int run_type,
                                 const int run_start_ymd,
                                 const int run_start_tod,
                                 const int case_start_ymd,
                                 const int case_start_tod,
                                 const char* calendar_name,
                                 const char* caseid,
                                 const char* rest_caseid,
                                 const char* hostname,
                                 const char* username,
                                 const char* versionid)
{
  using namespace scream;
  using namespace scream::control;

  fpe_guard_wrapper([&](){

    // Create the context
    auto& c = ScreamContext::singleton();

    // Create the C MPI_Comm from the Fortran one
    MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
    auto& atm_comm = c.create<ekat::Comm>(mpi_comm_c);

    // Initialize the scream session.
    scream::initialize_eamxx_session(atm_comm.am_i_root());

    std::string cal = calendar_name;
    if (cal=="NO_LEAP") {
      scream::set_use_leap_year (false);
    } else if (cal=="GREGORIAN") {
      scream::set_use_leap_year (true);
    } else {
      EKAT_ERROR_MSG ("Error! Invalid/unsupported calendar name.\n"
          "   - input name : " + cal + "\n"
          "   - valid names: NO_LEAP, GREGORIAN\n");
    }

    // Create a parameter list for inputs
    ekat::ParameterList scream_params("Scream Parameters");
    parse_yaml_file (input_yaml_file, scream_params);

    scream_params.sublist("driver_options").set<std::string>("Atm Log File",atm_log_file);

    // Need to register products in the factories *before* we attempt to create any.
    // In particular, register all atm processes, grids managers, and diagnostics.
    register_dynamics();
    register_physics();
    register_diagnostics();
    register_surface_coupling();

    // Create the bare ad, then start the initialization sequence
    auto& ad = c.create<AtmosphereDriver>();

    // Recall that e3sm uses the int YYYYMMDD to store a date
    int yy,mm,dd,hr,min,sec;
    yy  = (run_start_ymd / 100) / 100;
    mm  = (run_start_ymd / 100) % 100;
    dd  =  run_start_ymd % 100;
    hr  = (run_start_tod / 60) / 60;
    min = (run_start_tod / 60) % 60;
    sec =  run_start_tod % 60;
    util::TimeStamp run_t0 (yy,mm,dd,hr,min,sec);

    yy  = (case_start_ymd / 100) / 100;
    mm  = (case_start_ymd / 100) % 100;
    dd  =  case_start_ymd % 100;
    hr  = (case_start_tod / 60) / 60;
    min = (case_start_tod / 60) % 60;
    sec =  case_start_tod % 60;
    util::TimeStamp case_t0 (yy,mm,dd,hr,min,sec);

    ad.set_comm(atm_comm);
    ad.set_params(scream_params);
    ad.init_scorpio(atm_id);
    ad.init_time_stamps(run_t0,case_t0,run_type);
    ad.set_provenance_data (caseid,rest_caseid,hostname,username,versionid);
    ad.create_output_managers ();
    ad.create_atm_processes ();
    ad.create_grids ();
    ad.create_fields ();
  });
}

void scream_setup_surface_coupling (const char*& import_field_names, int*& import_cpl_indices,
                                    double*& x2a_ptr,
#ifdef HAVE_MOAB
                                    double*& x2a_moab_ptr,
#endif
                                    int*& import_vector_components,
                                    double*& import_constant_multiple, bool*& do_import_during_init,
                                    const int& num_cpl_imports, const int& num_scream_imports, const int& import_field_size,
                                    char*& export_field_names, int*& export_cpl_indices,
                                    double*& a2x_ptr,
#ifdef HAVE_MOAB
                                    double*& a2x_moab_ptr,
#endif
                                    int*& export_vector_components,
                                    double*& export_constant_multiple, bool*& do_export_during_init,
                                    const int& num_cpl_exports, const int& num_scream_exports, const int& export_field_size)
{
  using namespace scream;

  fpe_guard_wrapper([&](){
    // Fortran gives a 1d array of 32char strings. So let's memcpy the input char
    // strings into 2d char arrays. Each string is null-terminated (atm_mct_mod
    // makes sure of that).
    using name_t = char[32];
    name_t* names_in  = new name_t[num_scream_imports];
    name_t* names_out = new name_t[num_scream_exports];
    std::memcpy(names_in,  import_field_names, num_scream_imports*32*sizeof(char));
    std::memcpy(names_out, export_field_names, num_scream_exports*32*sizeof(char));

    // Convert F90 -> C++ indexing
    using view_t = typename KokkosTypes<HostDevice>::template view_1d<int>;
    view_t import_cpl_indices_view(import_cpl_indices, num_scream_imports);
    view_t export_cpl_indices_view(export_cpl_indices, num_scream_exports);
    for (int f=0; f<num_scream_imports; ++f) import_cpl_indices_view(f) -= 1;
    for (int f=0; f<num_scream_exports; ++f) export_cpl_indices_view(f) -= 1;

    // Call the AD to setup surface coupling
    auto& ad = get_ad_nonconst();

    ad.setup_surface_coupling_data_manager(scream::SurfaceCouplingTransferType::Import,
                                           num_cpl_imports, num_scream_imports, import_field_size, x2a_ptr,
#ifdef HAVE_MOAB
                                           x2a_moab_ptr,
#endif
                                           names_in[0], import_cpl_indices, import_vector_components,
                                           import_constant_multiple, do_import_during_init);
    ad.setup_surface_coupling_data_manager(scream::SurfaceCouplingTransferType::Export,
                                           num_cpl_exports, num_scream_exports, export_field_size, a2x_ptr,
#ifdef HAVE_MOAB
                                           a2x_moab_ptr,
#endif
                                           names_out[0], export_cpl_indices, export_vector_components,
                                           export_constant_multiple, do_export_during_init);
  });
}

#if defined(SCREAM_SYSTEM_WORKAROUND) && (SCREAM_SYSTEM_WORKAROUND == 1)
void scream_init_hip_atm () {
    hipInit(0);
}
#endif

void scream_init_atm ()
{
  using namespace scream;
  using namespace scream::control;

  fpe_guard_wrapper([&](){
    // Get the ad, then complete initialization
    auto& ad = get_ad_nonconst();

    // Init all fields, atm processes, and output streams
    ad.initialize_fields ();
    ad.initialize_atm_procs ();
    // Do this before init-ing the output managers,
    // so the fields are valid if outputing at t=0
    ad.reset_accumulated_fields();
    ad.initialize_output_managers ();
  });
}

/*===============================================================================================*/
void scream_run (const int dt) {
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

    // Clean up the context
    // Note: doing the cleanup here via
    //   scream::ScreamContext::singleton().clean_up();
    // causes an ICE with C++17 on Summit/Ascent.
    // Wrapping it in a function seems to work though
    scream::cleanup_singleton();

    // Finalize scream session
    scream::finalize_eamxx_session();
  });
}

// Get the local (i.e., on current atm rank only) number of physics columns
int scream_get_num_local_cols () {
  int ncols = -1;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& phys_grid = gm->get_grid("Physics");

    ncols = phys_grid->get_num_local_dofs();
  });

  return ncols;
}

// Get the global (i.e., the whole earth) number of physics columns
int scream_get_num_global_cols () {
  int ncols = -1;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& phys_grid = gm->get_grid("Physics");

    ncols = phys_grid->get_num_global_dofs();
  });

  return ncols;
}

// Return the global ids of all physics column
void scream_get_local_cols_gids (void* const ptr) {
  using namespace scream;
  using gid_type = AbstractGrid::gid_type;
  fpe_guard_wrapper([&]() {
    auto gids_f = reinterpret_cast<int*>(ptr);
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("Physics");

    auto gids = phys_grid->get_dofs_gids();
    gids.sync_to_host();
    auto gids_h = gids.get_view<const gid_type*,Host>();

    for (int i=0; i<gids_h.extent_int(0); ++i) {
      gids_f[i] = gids_h(i);
    }
  });
}

// Retrieve the lat/lon coords of all physics FV dofs
void scream_get_cols_latlon (double* const& lat_ptr, double* const& lon_ptr) {
  using namespace scream;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("Physics");
    const auto ncols = phys_grid->get_num_local_dofs();

    auto lat_cxx = phys_grid->get_geometry_data("lat").get_view<const Real*, Host>();
    auto lon_cxx = phys_grid->get_geometry_data("lon").get_view<const Real*, Host>();

    using geo_view_f90 = ekat::Unmanaged<decltype(lat_cxx)::HostMirror>;
    geo_view_f90 lat_f90(lat_ptr, ncols);
    geo_view_f90 lon_f90(lon_ptr, ncols);
    Kokkos::deep_copy(lat_f90, lat_cxx);
    Kokkos::deep_copy(lon_f90, lon_cxx);
  });
}

// Retrieve the area of all physics FV dofs
void scream_get_cols_area (double* const& area_ptr) {
  using namespace scream;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("Physics");
    const auto ncols = phys_grid->get_num_local_dofs();

    auto area_cxx = phys_grid->get_geometry_data("area").get_view<const Real*, Host>();

    using geo_view_f90 = ekat::Unmanaged<decltype(area_cxx)::HostMirror>;
    geo_view_f90  area_f90 (area_ptr, ncols);
    Kokkos::deep_copy(area_f90,area_cxx);
  });
}

} // extern "C"
