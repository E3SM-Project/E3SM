#include "scream_config.h"

#include "share/atm_process/atmosphere_process.hpp"
#include "control/atmosphere_driver.hpp"
#include "control/surface_coupling.hpp"

#include "dynamics/register_dynamics.hpp"
#include "physics/register_physics.hpp"

#include "mct_coupling/ScreamContext.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/scream_session.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_assert.hpp"

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
void scream_create_atm_instance (const MPI_Fint& f_comm,
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

    if (atm_comm.am_i_root()) {
      scream_params.print();
    }

    ekat::error::runtime_check(scream_params.isSublist("Atmosphere Driver"),
         "Error! Sublist 'Atmosphere Driver' not found inside '" +
         std::string(input_yaml_file) + "'.\n");

    auto& ad_params = scream_params.sublist("Atmosphere Driver");

    // Need to register products in the factories *before* we attempt to create any.
    // In particular, register all atm processes, and all grids managers.
    register_dynamics();
    register_physics();

    // Create the bare ad, then start the initialization sequence
    auto& ad = c.create<AtmosphereDriver>();

    ad.set_comm(atm_comm);
    ad.set_params(ad_params);
    ad.create_atm_processes ();
    ad.create_grids ();
    ad.create_fields ();
    ad.set_surface_coupling (std::make_shared<SurfaceCoupling>(ad.get_ref_grid_field_mgr()));
  });
}

void scream_setup_surface_coupling (
    const char*& x2a_names, const int*& x2a_indices, double*& cpl_x2a_ptr, const int*& vec_comp_x2a,
    const int& num_cpl_imports, const int& num_scream_imports,
    const char*& a2x_names, const int*& a2x_indices, double*& cpl_a2x_ptr, const int*& vec_comp_a2x,
    const int& num_exports)
{
  fpe_guard_wrapper([&](){
    // Fortran gives a 1d array of 32char strings. So let's memcpy the input char
    // strings into 2d char arrays. Each string is null-terminated (atm_mct_mod
    // makes sure of that).
    using name_t = char[32];
    name_t* names_in = new name_t[num_cpl_imports];
    name_t* names_out = new name_t[num_exports];
    std::memcpy(names_in,x2a_names,num_cpl_imports*32*sizeof(char));
    std::memcpy(names_out,a2x_names,num_exports*32*sizeof(char));

    // Get the SurfaceCoupling from the AD, then register the fields
    const auto& ad = get_ad();
    const auto& sc = ad.get_surface_coupling();

    // Register import/export fields
    sc->set_num_fields(num_cpl_imports,num_scream_imports,num_exports);
    for (int i=0; i<num_cpl_imports; ++i) {
      sc->register_import(names_in[i],x2a_indices[i],vec_comp_x2a[i]);
    }
    for (int i=0; i<num_exports; ++i) {
      sc->register_export(names_out[i],a2x_indices[i],vec_comp_a2x[i]);
    }

    sc->registration_ends(cpl_x2a_ptr, cpl_a2x_ptr);
  });
}

void scream_init_atm (const int& start_ymd,
                      const int& start_tod)
{
  using namespace scream;
  using namespace scream::control;

  fpe_guard_wrapper([&](){
    // Get the ad, then complete initialization
    auto& ad = get_ad_nonconst();

    // Recall that e3sm uses the int YYYYMMDD to store a date
    std::cout << "start_ymd: " << start_ymd << "\n";
    const int dd = start_ymd % 100;
    const int mm = (start_ymd / 100) % 100;
    const int yy = start_ymd / 10000;
    util::TimeStamp t0 (yy,mm,dd,start_tod);

    // Init and run (to finalize, wait till checks are completed,
    // or you'll clear the field managers!)
    ad.initialize_fields (t0);
    ad.initialize_output_manager ();
    ad.initialize_atm_procs ();
  });
}

/*===============================================================================================*/
void scream_run (const scream::Real& dt) {
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
    scream::ScreamContext::singleton().clean_up();
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
