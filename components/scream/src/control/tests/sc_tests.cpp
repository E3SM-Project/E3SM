#include "share/grid/point_grid.hpp"
#include "control/surface_coupling.hpp"

#include <ekat/util/ekat_test_utils.hpp>

#include <catch2/catch.hpp>
#include <memory>

TEST_CASE ("surface_coupling")
{
  /*
   * This test performs an export-import sequence, and it verifies that
   * the resulting fields (after the import) match the fields that were
   * originally exported. If successful, the test "ensures" that the
   * import/export kernels in the SurfaceCoupling class are correct.
   */

  // Some namespaces/aliases
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FL = FieldLayout;
  using FID = FieldIdentifier;
  using FR = FieldRequest;
  using GR = GroupRequest;
  using RPDF = std::uniform_real_distribution<Real>;
  using rngAlg = std::mt19937_64;

  // Some constants
  constexpr int ncols = 4;
  constexpr int nlevs = 8;
  constexpr int nruns = 10;

  // Create a comm
  ekat::Comm comm (MPI_COMM_WORLD);

  // The random numbers generator
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  if (comm.am_i_root()) {
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  }
  rngAlg engine(seed);
  RPDF pdf(0.0,1.0);

  // Create a grid
  auto grid = create_point_grid("my grid",ncols*comm.size(), nlevs, comm);

  // Create some field ids, and register them in a field manager
  // Note: we create two fms, so we can compare outputs with inputs
  FID s2d_id("s2d",FL{{COL},{ncols}},Pa,grid->name());
  FID s3d_id("s3d",FL{{COL,LEV},{ncols,nlevs}},Pa,grid->name());
  FID v2d_id("v2d",FL{{COL,CMP},{ncols,2}},m/s,grid->name());
  FID v3d_id("v3d",FL{{COL,CMP,LEV},{ncols,2,nlevs}},m/s,grid->name());

  FID sub_s2d1_id("sub_s2d1",FL{{COL},{ncols}},Pa,grid->name());
  FID sub_s2d2_id("sub_s2d2",FL{{COL},{ncols}},Pa,grid->name());
  FID sub_s3d1_id("sub_s3d1",FL{{COL,LEV},{ncols,nlevs}},Pa,grid->name());
  FID sub_s3d2_id("sub_s3d2",FL{{COL,LEV},{ncols,nlevs}},Pa,grid->name());


  // NOTE: if you add fields above, you will have to modify these counters too.
  const int num_s2d = 1;
  const int num_s3d = 1;
  const int num_v2d = 1;
  const int num_v3d = 1;
  const int num_sub_s3d = 2;
  const int num_fields = num_s2d+num_s3d+2*(num_v2d+num_v3d)+num_sub_s3d;

  // Keep two separate fms, so we can compare original and final fields.
  auto fm_in = std::make_shared<FieldManager<Real>> (grid);
  auto fm_out = std::make_shared<FieldManager<Real>> (grid);
  fm_in->registration_begins();
  fm_in->register_field(s2d_id);
  fm_in->register_field(s3d_id);
  fm_in->register_field(v2d_id);
  fm_in->register_field(v3d_id);
  fm_in->register_field(FR{sub_s3d1_id,"group3d"});
  fm_in->register_field(FR{sub_s3d2_id,"group3d"});
  fm_in->register_group(GR("group3d",grid->name(),Bundling::Required));
  fm_in->registration_ends();

  fm_out->registration_begins();
  fm_out->register_field(s2d_id);
  fm_out->register_field(s3d_id);
  fm_out->register_field(v2d_id);
  fm_out->register_field(v3d_id);
  fm_out->register_field(FR{sub_s3d1_id,"group3d"});
  fm_out->register_field(FR{sub_s3d2_id,"group3d"});
  fm_out->register_group(GR("group3d",grid->name(),Bundling::Required));
  fm_out->registration_ends();

  // Create two SC objects, to import and export
  control::SurfaceCoupling importer(fm_in);
  control::SurfaceCoupling exporter(fm_out);

  importer.set_num_fields(num_fields,0); // Recall that SC counts *scalar* fields, so vector3d counts as 2 fields
  exporter.set_num_fields(0,num_fields);

  // Register fields in the importer/exporter
  // Note: the 1st integer is the field "idx" (the idx used by the component cpl to retrieve it),
  //       which here is not really used (though it still needs to be passed). The 2nd integer
  //       is needed for vector fields, to tell the importer/exporter which component of the
  //       vector field is imported/exported.
  importer.register_import("s2d",0);
  importer.register_import("s3d",1);
  importer.register_import("v2d",2,0);
  importer.register_import("v2d",3,1);
  importer.register_import("v3d",4,0);
  importer.register_import("v3d",5,1);
  importer.register_import("sub_s3d1",6);
  importer.register_import("sub_s3d2",7);

  exporter.register_export("s2d",0);
  exporter.register_export("s3d",1);
  exporter.register_export("v2d",2,0);
  exporter.register_export("v2d",3,1);
  exporter.register_export("v3d",4,0);
  exporter.register_export("v3d",5,1);
  exporter.register_export("sub_s3d1",6);
  exporter.register_export("sub_s3d2",7);

  // Create a raw array big enough to contain all the 2d data for import/export
  double* raw_data = new double[ncols*num_fields];

  // Complete setup of importer/exporter
  importer.registration_ends(raw_data,nullptr);
  exporter.registration_ends(nullptr,raw_data);

  // Repeat experiment N times: fill export fields, export, import, check import fields
  auto s2d_exp = fm_out->get_field(s2d_id);
  auto s3d_exp = fm_out->get_field(s3d_id);
  auto v2d_exp = fm_out->get_field(v2d_id);
  auto v3d_exp = fm_out->get_field(v3d_id);
  auto sub_s3d1_exp = fm_out->get_field(sub_s3d1_id);
  auto sub_s3d2_exp = fm_out->get_field(sub_s3d2_id);
  auto G3d_exp = fm_out->get_field(fm_out->get_field_group("group3d").m_bundle->get_header().get_identifier().name());

  auto s2d_imp = fm_in->get_field(s2d_id);
  auto s3d_imp = fm_in->get_field(s3d_id);
  auto v2d_imp = fm_in->get_field(v2d_id);
  auto v3d_imp = fm_in->get_field(v3d_id);
  auto sub_s3d1_imp = fm_in->get_field(sub_s3d1_id);
  auto sub_s3d2_imp = fm_in->get_field(sub_s3d2_id);
  auto G3d_imp = fm_in->get_field(fm_in->get_field_group("group3d").m_bundle->get_header().get_identifier().name());

  auto s2d_exp_d = s2d_exp.get_reshaped_view<Real*>();
  auto s3d_exp_d = s3d_exp.get_reshaped_view<Real**>();
  auto v2d_exp_d = v2d_exp.get_reshaped_view<Real**>();
  auto v3d_exp_d = v3d_exp.get_reshaped_view<Real***>();
  auto sub_s3d1_exp_d = sub_s3d1_exp.get_reshaped_view<Real**>();
  auto sub_s3d2_exp_d = sub_s3d2_exp.get_reshaped_view<Real**>();

  auto s2d_imp_d = s2d_imp.get_reshaped_view<Real*>();
  auto s3d_imp_d = s3d_imp.get_reshaped_view<Real**>();
  auto v2d_imp_d = v2d_imp.get_reshaped_view<Real**>();
  auto v3d_imp_d = v3d_imp.get_reshaped_view<Real***>();
  auto sub_s3d1_imp_d = sub_s3d1_imp.get_reshaped_view<Real**>();
  auto sub_s3d2_imp_d = sub_s3d2_imp.get_reshaped_view<Real**>();

  auto s2d_exp_h = s2d_exp.get_reshaped_view<Real*,Host>();
  auto s3d_exp_h = s3d_exp.get_reshaped_view<Real**,Host>();
  auto v2d_exp_h = v2d_exp.get_reshaped_view<Real**,Host>();
  auto v3d_exp_h = v3d_exp.get_reshaped_view<Real***,Host>();
  auto sub_s3d1_exp_h = sub_s3d1_exp.get_reshaped_view<Real**,Host>();
  auto sub_s3d2_exp_h = sub_s3d2_exp.get_reshaped_view<Real**,Host>();

  auto s2d_imp_h = s2d_imp.get_reshaped_view<Real*,Host>();
  auto s3d_imp_h = s3d_imp.get_reshaped_view<Real**,Host>();
  auto v2d_imp_h = v2d_imp.get_reshaped_view<Real**,Host>();
  auto v3d_imp_h = v3d_imp.get_reshaped_view<Real***,Host>();
  auto sub_s3d1_imp_h = sub_s3d1_imp.get_reshaped_view<Real**,Host>();
  auto sub_s3d2_imp_h = sub_s3d2_imp.get_reshaped_view<Real**,Host>();

  for (int i=0; i<nruns; ++i) {
    // Fill export fields
    ekat::genRandArray(s2d_exp_d,engine,pdf);
    ekat::genRandArray(s3d_exp_d,engine,pdf);
    ekat::genRandArray(v2d_exp_d,engine,pdf);
    ekat::genRandArray(v3d_exp_d,engine,pdf);
    ekat::genRandArray(G3d_exp.get_view(),engine,pdf);

    // Set all raw_data to -1 (might be helpful for debugging)
    std::fill_n(raw_data,4*ncols,-1);

    // Perform export
    exporter.do_export();

    // Perform import
    importer.do_import();

    // Check f_imported==f_exported (on surface only)
    s2d_exp.sync_to_host();
    s3d_exp.sync_to_host();
    v2d_exp.sync_to_host();
    v3d_exp.sync_to_host();
    sub_s3d1_exp.sync_to_host();
    sub_s3d2_exp.sync_to_host();
    s2d_imp.sync_to_host();
    s3d_imp.sync_to_host();
    v2d_imp.sync_to_host();
    v3d_imp.sync_to_host();
    sub_s3d1_imp.sync_to_host();
    sub_s3d2_imp.sync_to_host();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE (s2d_exp_h(icol)==s2d_imp_h(icol));
      REQUIRE (s3d_exp_h(icol,nlevs-1)==s3d_imp_h(icol,nlevs-1));
      REQUIRE (v2d_exp_h(icol,0)==v2d_imp_h(icol,0));
      REQUIRE (v2d_exp_h(icol,1)==v2d_imp_h(icol,1));
      REQUIRE (v3d_exp_h(icol,0,nlevs-1)==v3d_imp_h(icol,0,nlevs-1));
      REQUIRE (v3d_exp_h(icol,1,nlevs-1)==v3d_imp_h(icol,1,nlevs-1));
      REQUIRE (sub_s3d1_exp_h(icol,nlevs-1)==sub_s3d1_imp_h(icol,nlevs-1));
      REQUIRE (sub_s3d2_exp_h(icol,nlevs-1)==sub_s3d2_imp_h(icol,nlevs-1));
    }
  }

  // Clean up
  delete[] raw_data;
}


TEST_CASE ("recreate_mct_coupling")
{
  /*
   * This test aims to recreate the import/export of fields inside the mct-coupling
   * for CIME runs. Much of the code inside surface coupling is hard-coded for these runs.
   */

  // Some namespaces/aliases
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FL     = FieldLayout;
  using FID    = FieldIdentifier;
  using FR     = FieldRequest;
  using GR     = GroupRequest;
  using RPDF   = std::uniform_real_distribution<Real>;
  using rngAlg = std::mt19937_64;

  // Some constants
  constexpr int ncols = 4;
  constexpr int nlevs = 8;
  constexpr int nruns = 10;

  // Create a comm
  ekat::Comm comm (MPI_COMM_WORLD);

  // The random numbers generator
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  if (comm.am_i_root()) {
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  }
  rngAlg engine(seed);
  RPDF pdf(0.0,1.0);

  // Create a grid
  auto grid = create_point_grid("my_grid",ncols*comm.size(), nlevs, comm);
  const auto grid_name = grid->name();

  // Layouts matching those in AD
  FL scalar2d_layout{ {COL          }, {ncols          } };
  FL scalar3d_layout{ {COL, LEV     }, {ncols,    nlevs} };
  FL vector3d_layout{ {COL, CMP, LEV}, {ncols, 2, nlevs} };

  // Create import fields
  FID surf_latent_flux_id ("surf_latent_flux", scalar2d_layout, W/(m*m), grid_name);
  FID surf_sens_flux_id   ("surf_sens_flux",   scalar2d_layout, W/(m*m), grid_name);
  FID surf_u_mom_flux_id  ("surf_u_mom_flux",  scalar2d_layout, W/(m*m), grid_name);
  FID surf_v_mom_flux_id  ("surf_v_mom_flux",  scalar2d_layout, W/(m*m), grid_name);

  // Create necessary fields for export. Tracers qc and qr are unnecessary, but
  // are included to verify that subviewed fields (qv) are correctly handled
  const auto nondim = Units::nondimensional();
  FID T_mid_id           ("T_mid",           scalar3d_layout, K,      grid_name);
  FID p_mid_id           ("p_mid",           scalar3d_layout, Pa,     grid_name);
  FID z_mid_id           ("z_mid",           scalar3d_layout, m,      grid_name);
  FID horiz_winds_id     ("horiz_winds",     vector3d_layout, m/s,    grid_name);
  FID pseudo_density_id  ("pseudo_density",  scalar3d_layout, Pa,     grid_name);
  FID qv_id              ("qv",              scalar3d_layout, nondim, grid_name);
  FID precip_liq_surf_id ("precip_liq_surf", scalar2d_layout, m/s,    grid_name);

  // NOTE: if you add fields above, you will have to modify these counters too.
  const int num_total_imports  = 21;
  const int num_scream_imports = 4;
  const int num_exports        = 13;

  // Register fields and tracer group in a FieldManager
  auto fm = std::make_shared<FieldManager<Real>> (grid);
  fm->registration_begins();
  fm->register_field(FR{surf_latent_flux_id});
  fm->register_field(FR{surf_sens_flux_id});
  fm->register_field(FR{surf_u_mom_flux_id});
  fm->register_field(FR{surf_v_mom_flux_id});
  fm->register_field(FR{T_mid_id});
  fm->register_field(FR{p_mid_id});
  fm->register_field(FR{z_mid_id});
  fm->register_field(FR{horiz_winds_id});
  fm->register_field(FR{pseudo_density_id});
  fm->register_field(FR{qv_id,"tracers"});
  fm->register_field(FR{precip_liq_surf_id});

  fm->register_group(GR("tracers", grid_name ,Bundling::Required));
  fm->registration_ends();

  // Create alias to field views
  auto surf_latent_flux_f = fm->get_field(surf_latent_flux_id);
  auto surf_sens_flux_f   = fm->get_field(surf_sens_flux_id);
  auto surf_u_mom_flux_f  = fm->get_field(surf_u_mom_flux_id);
  auto surf_v_mom_flux_f  = fm->get_field(surf_v_mom_flux_id);
  auto T_mid_f            = fm->get_field(T_mid_id);
  auto p_mid_f            = fm->get_field(p_mid_id);
  auto z_mid_f            = fm->get_field(z_mid_id);
  auto horiz_winds_f      = fm->get_field(horiz_winds_id);
  auto pseudo_density_f   = fm->get_field(pseudo_density_id);
  auto qv_f               = fm->get_field(qv_id);
  auto precip_liq_surf_f  = fm->get_field(precip_liq_surf_id);

  auto group = fm->get_field_group("tracers");
  const auto& Q_name = group.m_bundle->get_header().get_identifier().name();
  auto Q = fm->get_field(Q_name);

  auto surf_latent_flux_d = surf_latent_flux_f.get_reshaped_view<Real*>();
  auto surf_sens_flux_d   = surf_sens_flux_f.get_reshaped_view<Real*>();
  auto surf_u_mom_flux_d  = surf_u_mom_flux_f.get_reshaped_view<Real*>();
  auto surf_v_mom_flux_d  = surf_v_mom_flux_f.get_reshaped_view<Real*>();
  auto T_mid_d            = T_mid_f.get_reshaped_view<Real**>();
  auto p_mid_d            = p_mid_f.get_reshaped_view<Real**>();
  auto z_mid_d            = z_mid_f.get_reshaped_view<Real**>();
  auto horiz_winds_d      = horiz_winds_f.get_reshaped_view<Real***>();
  auto pseudo_density_d   = pseudo_density_f.get_reshaped_view<Real**>();
  auto qv_d               = qv_f.get_reshaped_view<Real**>();
  auto precip_liq_surf_d  = precip_liq_surf_f.get_reshaped_view<Real*>();

  auto surf_latent_flux_h = surf_latent_flux_f.get_reshaped_view<Real*,Host>();
  auto surf_sens_flux_h   = surf_sens_flux_f.get_reshaped_view<Real*,Host>();
  auto surf_u_mom_flux_h  = surf_u_mom_flux_f.get_reshaped_view<Real*,Host>();
  auto surf_v_mom_flux_h  = surf_v_mom_flux_f.get_reshaped_view<Real*,Host>();
  auto T_mid_h            = T_mid_f.get_reshaped_view<Real**,Host>();
  auto p_mid_h            = p_mid_f.get_reshaped_view<Real**,Host>();
  auto z_mid_h            = z_mid_f.get_reshaped_view<Real**,Host>();
  auto pseudo_density_h   = pseudo_density_f.get_reshaped_view<Real**,Host>();
  auto horiz_winds_h      = horiz_winds_f.get_reshaped_view<Real***,Host>();
  auto qv_h               = qv_f.get_reshaped_view<Real**,Host>();
  auto precip_liq_surf_h  = precip_liq_surf_f.get_reshaped_view<Real*,Host>();

  // Create SC object and set number of import/export fields
  control::SurfaceCoupling coupler(fm);
  coupler.set_num_fields(num_total_imports, num_scream_imports, num_exports);

  // Register fields in the coupler. These match the scr_names_x2a/a2x from
  // scream_cpl_indices.F90. When radiation is fully implemented, RRTMGP
  // fields will be replaced with correct fields.
  coupler.register_import("surf_latent_flux", 0);
  coupler.register_import("surf_sens_flux",   1);
  coupler.register_import("unused",           2);
  coupler.register_import("surf_u_mom_flux",  3);
  coupler.register_import("surf_v_mom_flux",  4);
  coupler.register_import("RRTMGP",           5);
  coupler.register_import("RRTMGP",           6);
  coupler.register_import("RRTMGP",           7);
  coupler.register_import("RRTMGP",           8);
  coupler.register_import("RRTMGP",           9);
  coupler.register_import("unused",           10);
  coupler.register_import("unused",           11);
  coupler.register_import("unused",           12);
  coupler.register_import("unused",           13);
  coupler.register_import("unused",           14);
  coupler.register_import("unused",           15);
  coupler.register_import("RRTMGP",           16);
  coupler.register_import("unused",           17);
  coupler.register_import("RRTMGP",           18);
  coupler.register_import("unused",           19);
  coupler.register_import("unused",           20);

  coupler.register_export("T_mid",            0);
  coupler.register_export("Sa_ptem",          1);
  coupler.register_export("z_mid",            2);
  coupler.register_export("Sa_u",             3);
  coupler.register_export("Sa_v",             4);
  coupler.register_export("p_mid",            5);
  coupler.register_export("Sa_dens",          6);
  coupler.register_export("qv",               7);
  coupler.register_export("set_zero",         8);
  coupler.register_export("precip_liq_surf",  9);
  coupler.register_export("set_zero",         10);
  coupler.register_export("set_zero",         11);
  coupler.register_export("set_zero",         12);

  // Complete setup of importer/exporter, providing raw_data
  double* import_raw_data = new double[ncols*num_total_imports];
  double* export_raw_data = new double[ncols*num_exports];
  coupler.registration_ends(import_raw_data, export_raw_data);

  for (int i=0; i<nruns; ++i) {

    // Set import field views to 0
    surf_latent_flux_f.deep_copy(0.0);
    surf_sens_flux_f.deep_copy(0.0);
    surf_u_mom_flux_f.deep_copy(0.0);
    surf_v_mom_flux_f.deep_copy(0.0);

    // Fill views needed in the export with random values
    ekat::genRandArray(T_mid_d,engine,pdf);
    ekat::genRandArray(p_mid_d,engine,pdf);
    ekat::genRandArray(z_mid_d,engine,pdf);
    ekat::genRandArray(horiz_winds_d,engine,pdf);
    ekat::genRandArray(pseudo_density_d,engine,pdf);
    ekat::genRandArray(precip_liq_surf_d,engine,pdf);
    ekat::genRandArray(Q.get_view(),engine,pdf);

    // Fill import_raw_data with random values
    for (int i=0; i<ncols*num_total_imports; ++i) {
      import_raw_data[i] = pdf(engine);
    }

    // Set all export_raw_data to -1 (might be helpful for debugging)
    std::fill_n(export_raw_data, num_exports*ncols, -1);

    // Perform import
    coupler.do_import();

    // Perform export
    coupler.do_export();

    // Sync host to device
    surf_latent_flux_f.sync_to_host();
    surf_sens_flux_f.sync_to_host();
    surf_u_mom_flux_f.sync_to_host();
    surf_v_mom_flux_f.sync_to_host();
    T_mid_f.sync_to_host();
    p_mid_f.sync_to_host();
    z_mid_f.sync_to_host();
    horiz_winds_f.sync_to_host();
    pseudo_density_f.sync_to_host();
    precip_liq_surf_f.sync_to_host();
    Q.sync_to_host();

    // Check values
    for (int icol=0; icol<ncols; ++icol) {

      // Imports

      REQUIRE (surf_latent_flux_h(icol) == import_raw_data[0 + icol*num_total_imports]); // 1st scream import
      REQUIRE (surf_sens_flux_h  (icol) == import_raw_data[1 + icol*num_total_imports]); // 2nd scream import
      REQUIRE (surf_u_mom_flux_h (icol) == import_raw_data[3 + icol*num_total_imports]); // 3rd scream import (4th total import)
      REQUIRE (surf_v_mom_flux_h (icol) == import_raw_data[4 + icol*num_total_imports]); // 4th scream import (5th total import)

      // Exports

      // These exports are direct values from a scream field
      REQUIRE (export_raw_data[0 + icol*num_exports] == T_mid_h          (icol,    nlevs-1)); // 1st export
      REQUIRE (export_raw_data[2 + icol*num_exports] == z_mid_h          (icol,    nlevs-1)); // 3rd export
      REQUIRE (export_raw_data[3 + icol*num_exports] == horiz_winds_h    (icol, 0, nlevs-1)); // 4th export
      REQUIRE (export_raw_data[4 + icol*num_exports] == horiz_winds_h    (icol, 1, nlevs-1)); // 5th export
      REQUIRE (export_raw_data[5 + icol*num_exports] == p_mid_h          (icol,    nlevs-1)); // 6th export
      REQUIRE (export_raw_data[7 + icol*num_exports] == qv_h             (icol,    nlevs-1)); // 8th export
      REQUIRE (export_raw_data[9 + icol*num_exports] == precip_liq_surf_h(icol            )); // 10th export

      // These exports should be set to 0
      const auto eps = std::numeric_limits<Real>::epsilon();
      REQUIRE (abs(export_raw_data[8  + icol*num_exports]) < eps); // 9th export
      REQUIRE (abs(export_raw_data[10 + icol*num_exports]) < eps); // 11th export
      REQUIRE (abs(export_raw_data[11 + icol*num_exports]) < eps); // 12th export
      REQUIRE (abs(export_raw_data[12 + icol*num_exports]) < eps); // 13th export
    }
  }

  // Clean up
  delete[] export_raw_data;
  delete[] import_raw_data;
}
