/**
 * @file atm.cpp
 * @brief Atmosphere emulator component implementation.
 *
 * Stub implementation — fill in details for AI/ML inference,
 * coupling, and I/O.
 */

#include "eamxx_atm.hpp"
#include "emulator_c_api.hpp"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <ekat_fpe.hpp>

#include <mpi.h>

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

} // anonymous namespace

namespace eamxx {

Atm::Atm()
    : Emulator(EmulatorType::ATM_COMP, -1, "eamxx") {}

void Atm::create_instance(const MPI_Fint f_comm, const int atm_id,
                               const char* input_yaml_file,
                               const char* atm_log_file,
                               const int run_type,
                               const int run_start_ymd,
                               const int run_start_tod,
                               const int case_start_ymd,
                               const int case_start_tod,
                               const char* calendar_name) {
                             /*const char* caseid,
                               const char* rest_caseid,
                               const char* hostname,
                               const char* username,
                               const char* versionid) {*/
  using namespace scream;
  using namespace scream::control;

  fpe_guard_wrapper([&](){

    m_comm = config->comm;
    m_id = comp_id;  // set base class ID
    m_input_file = input_file;
    m_log_file = log_file;
    m_run_type = run_type;
    (void)start_ymd;
    (void)start_tod;

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
    // FIXME: can we infer provenance data from build/runtime information?
    // ad.set_provenance_data (caseid,rest_caseid,hostname,username,versionid);
    ad.init_scorpio(atm_id);
    ad.init_time_stamps(run_t0,case_t0,run_type);
    ad.create_output_managers();
    ad.create_atm_processes();
    ad.create_grids();
    ad.create_fields();
  });
}

void Atm::set_grid_data(const EmulatorGridDesc& grid) {
  m_nx = grid.nx;
  m_ny = grid.ny;
  m_num_local_cols = grid.num_local_cols;
  m_num_global_cols = grid.num_global_cols;

  m_col_gids.assign(grid.col_gids, grid.col_gids + grid.num_local_cols);
  m_lat.assign(grid.lat, grid.lat + grid.num_local_cols);
  m_lon.assign(grid.lon, grid.lon + grid.num_local_cols);
  m_area.assign(grid.area, grid.area + grid.num_local_cols);
}

void Atm::init_coupling_indices(
    const std::string &export_fields,
    const std::string &import_fields) {
  // TODO: Parse colon-separated MCT field lists and populate
  // m_coupling_idx with index positions.
  (void)export_fields;
  (void)import_fields;
}

void Atm::setup_coupling(const EmulatorCouplingDesc& cpl) {
  // FIXME: see scream_setup_surface_coupling
}

int Atm::get_num_local_cols() const {
  int ncols = -1;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& phys_grid = gm->get_grid("physics");

    ncols = phys_grid->get_num_local_dofs();
  });

  return ncols;
}

int Atm::get_num_global_cols() const {
  int ncols = -1;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& gm = ad.get_grids_manager();
    const auto& phys_grid = gm->get_grid("physics");

    ncols = phys_grid->get_num_global_dofs();
  });

  return ncols;
}

void Atm::get_local_col_gids(int *gids) const {
  using namespace scream;
  using gid_type = AbstractGrid::gid_type;
  fpe_guard_wrapper([&]() {
    auto gids_f = reinterpret_cast<int*>(ptr);
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("physics");

    auto gids = phys_grid->get_dofs_gids();
    gids.sync_to_host();
    auto gids_h = gids.get_view<const gid_type*,Host>();

    for (int i=0; i<gids_h.extent_int(0); ++i) {
      gids_f[i] = gids_h(i);
    }
  });
}

void Atm::get_cols_latlon(double *lat, double *lon) const {
  using namespace scream;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("physics");
    const auto ncols = phys_grid->get_num_local_dofs();

    auto lat_cxx = phys_grid->get_geometry_data("lat").get_view<const Real*, Host>();
    auto lon_cxx = phys_grid->get_geometry_data("lon").get_view<const Real*, Host>();

    using geo_view_f90 = ekat::Unmanaged<decltype(lat_cxx)::host_mirror_type>;
    geo_view_f90 lat_f90(lat_ptr, ncols);
    geo_view_f90 lon_f90(lon_ptr, ncols);
    Kokkos::deep_copy(lat_f90, lat_cxx);
    Kokkos::deep_copy(lon_f90, lon_cxx);
  });
}

void Atm::get_cols_area(double *area) const {
  using namespace scream;
  fpe_guard_wrapper([&]() {
    const auto& ad = get_ad();
    const auto& phys_grid = ad.get_grids_manager()->get_grid("physics");
    const auto ncols = phys_grid->get_num_local_dofs();

    auto area_cxx = phys_grid->get_geometry_data("area").get_view<const Real*, Host>();

    using geo_view_f90 = ekat::Unmanaged<decltype(area_cxx)::host_mirror_type>;
    geo_view_f90  area_f90 (area_ptr, ncols);
    Kokkos::deep_copy(area_f90,area_cxx);
  });
}

// =========================================================================
// Lifecycle implementations
// =========================================================================

void Atm::init_impl() {
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

void Atm::run_impl(int dt) {
  (void)dt;
  // FIXME: do steps 1, 2, 4 pertain? EAMxx uses only step 3.
  // 1. Import fields from coupler
  import_coupling_fields();

  // 2. Prepare model inputs
  prepare_inputs();

  // 3. Run the model
  fpe_guard_wrapper([&](){
    // Get the AD, and run it
    auto& ad = get_ad_nonconst();
    ad.run(dt);
  });

  // 4. Process AI outputs
  process_outputs();

  // 5. TODO: Diagnostic output

  // 6. Export fields to coupler
  export_coupling_fields();
}

void Atm::final_impl() {
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

// =========================================================================
// Coupling helpers
// =========================================================================

void Atm::import_coupling_fields() {
  // TODO: Transfer coupler import data → internal fields
}

void Atm::export_coupling_fields() {
  // TODO: Transfer internal fields → coupler export data
}

void Atm::prepare_inputs() {
  // TODO: Pack field data into m_fields.net_inputs tensor
  // for inference. Handle spatial_mode vs pointwise layout.
}

void Atm::process_outputs() {
  // TODO: Unpack m_fields.net_outputs tensor into field
  // vectors. Handle spatial_mode vs pointwise layout.
}

} // namespace eamxx
