#include "catch2/catch.hpp"

#include "nudging_tests_helpers.hpp"

using namespace scream;

TEST_CASE("create_nudging_data") {

  // Generate nudging data on a grid with 12 global cols and 20 levs.
  // We create two copies of the data set. One with some filled values,
  // and one without any filled values

  ekat::Comm comm(MPI_COMM_WORLD);

  const int nlevs  = nlevs_data;
  const int ngcols = ngcols_data;

  const int dt     = dt_data;
  const int nsteps = nsteps_data;
  const auto t0    = get_t0();

  // Initialize the pio_subsystem for this test:
  scorpio::eam_init_pio_subsystem(comm);

  // Create a grids manager
  const auto gm = create_gm(comm,ngcols,nlevs);
  const auto grid = gm->get_grid("Point Grid");

  // Create a field manager, and init fields (since OM's need t0 values)
  const auto fm1 = create_fm(grid);
  const auto fm2 = create_fm(grid);
  const auto fm3 = create_fm(grid);
  update_fields(fm1,t0,0);
  update_fields(fm2,t0,nlevs_filled);
  update_fields(fm3,t0,0);

  // Create output manager
  const auto om1 = create_om("nudging_data",fm1,gm,t0,comm);
  const auto om2 = create_om("nudging_data_filled",fm2,gm,t0,comm);
  const auto om3 = create_om("nudging_data_nonconst_p",fm3,gm,t0,comm);

  auto time = t0;
  for (int istep=1; istep<=nsteps; ++istep) {
    time += dt;

    // Compute fields, but keep p_mid constnat in fm1 and fm2, to avoid vinterp
    update_fields(fm1,time,0,false);
    update_fields(fm2,time,nlevs_filled,false);
    update_fields(fm3,time,0);

    om1->run(time);
    om2->run(time);
    om3->run(time);
  }
  om1->finalize();
  om2->finalize();
  om3->finalize();

  scorpio::eam_pio_finalize();
}
