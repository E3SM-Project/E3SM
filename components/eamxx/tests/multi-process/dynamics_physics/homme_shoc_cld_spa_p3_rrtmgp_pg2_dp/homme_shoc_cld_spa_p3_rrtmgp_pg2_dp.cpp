#include "catch2/catch.hpp"

// The AD
#include "control/atmosphere_driver.hpp"

// Dynamics includes
#include "dynamics/register_dynamics.hpp"

// Physics includes
#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"

// Surface coupling includes
#include "control/register_surface_coupling.hpp"
#include "control/atmosphere_surface_coupling_importer.hpp"

// EKAT headers
#include "ekat/kokkos/ekat_kokkos_types.hpp"

TEST_CASE("scream_homme_physics", "scream_homme_physics") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  ad_params.print();

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  // Register all atm procs and the grids manager in the respective factories
  register_dynamics();
  register_physics();
  register_diagnostics();
  register_surface_coupling();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  // NOTE: Kokkos is finalize in ekat_catch_main.cpp, and YAKL is finalized
  //       during RRTMGPRatiation::finalize_impl, after RRTMGP has deallocated
  //       all its arrays.
  ad.set_comm(atm_comm);
  ad.set_params(ad_params);
  ad.init_scorpio ();
  ad.init_time_stamps (t0, t0);
  ad.create_output_managers ();
  ad.create_atm_processes ();
  ad.create_grids ();
  ad.create_fields ();

  // Setup surface coupler import to be NaNs for fields IOP should overwrite
  const int ncols = ad.get_grids_manager()->get_grid("Physics")->get_num_local_dofs();
  static constexpr int num_imports = 4;
  char import_names[num_imports][32];
  std::strcpy(import_names[0], "surf_radiative_T");
  std::strcpy(import_names[1], "surf_lw_flux_up");
  std::strcpy(import_names[2], "surf_sens_flux");
  std::strcpy(import_names[3], "surf_evap");
  KokkosTypes<HostDevice>::view_2d<Real> import_data("import_data",ncols, num_imports);
  Kokkos::deep_copy(import_data, std::nan(""));
  KokkosTypes<HostDevice>::view_1d<int> import_cpl_indices("import_cpl_indices", num_imports);
  std::iota(import_cpl_indices.data(), import_cpl_indices.data()+num_imports, 0);
  KokkosTypes<HostDevice>::view_1d<int> import_vec_comps("import_vec_comps", num_imports);
  Kokkos::deep_copy(import_vec_comps, -1);
  KokkosTypes<HostDevice>::view_1d<Real> import_constant_multiple("import_constant_multiple", num_imports);
  Kokkos::deep_copy(import_constant_multiple, 1);
  KokkosTypes<HostDevice>::view_1d<bool> do_import_during_init("do_import_during_init", num_imports);
  Kokkos::deep_copy(do_import_during_init, false);
  do_import_during_init(2) = true; do_import_during_init(3) = true;

  ad.setup_surface_coupling_data_manager(SurfaceCouplingTransferType::Import,
                                         4, 4, ncols, import_data.data(), import_names[0], import_cpl_indices.data(),
                                         import_vec_comps.data(), import_constant_multiple.data(), do_import_during_init.data());
  ad.initialize_fields ();
  ad.initialize_output_managers ();
  ad.initialize_atm_procs ();

  if (atm_comm.am_i_root()) {
    printf("Start time stepping loop...       [  0%%]\n");
  }
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    if (atm_comm.am_i_root()) {
      std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
      std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
    }
  }
  ad.finalize();


  // If we got here, we were able to run without errors.
  REQUIRE (true);
}
