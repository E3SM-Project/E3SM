#include "share/grid/point_grid.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "control/surface_coupling.hpp"
#include "physics/share/physics_constants.hpp"

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

  // Some constants
  constexpr int ncols = 4;
  constexpr int nlevs = 8;
  constexpr int nruns = 10;

  // Create a comm
  ekat::Comm comm (MPI_COMM_WORLD);

  // The random numbers generator
  auto engine = setup_random_test(&comm);
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
  auto fm_in = std::make_shared<FieldManager> (grid);
  auto fm_out = std::make_shared<FieldManager> (grid);
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

  auto s2d_exp_d = s2d_exp.get_view<Real*>();
  auto s3d_exp_d = s3d_exp.get_view<Real**>();
  auto v2d_exp_d = v2d_exp.get_view<Real**>();
  auto v3d_exp_d = v3d_exp.get_view<Real***>();
  auto sub_s3d1_exp_d = sub_s3d1_exp.get_view<Real**>();
  auto sub_s3d2_exp_d = sub_s3d2_exp.get_view<Real**>();

  auto s2d_imp_d = s2d_imp.get_view<Real*>();
  auto s3d_imp_d = s3d_imp.get_view<Real**>();
  auto v2d_imp_d = v2d_imp.get_view<Real**>();
  auto v3d_imp_d = v3d_imp.get_view<Real***>();
  auto sub_s3d1_imp_d = sub_s3d1_imp.get_view<Real**>();
  auto sub_s3d2_imp_d = sub_s3d2_imp.get_view<Real**>();

  auto s2d_exp_h = s2d_exp.get_view<Real*,Host>();
  auto s3d_exp_h = s3d_exp.get_view<Real**,Host>();
  auto v2d_exp_h = v2d_exp.get_view<Real**,Host>();
  auto v3d_exp_h = v3d_exp.get_view<Real***,Host>();
  auto sub_s3d1_exp_h = sub_s3d1_exp.get_view<Real**,Host>();
  auto sub_s3d2_exp_h = sub_s3d2_exp.get_view<Real**,Host>();

  auto s2d_imp_h = s2d_imp.get_view<Real*,Host>();
  auto s3d_imp_h = s3d_imp.get_view<Real**,Host>();
  auto v2d_imp_h = v2d_imp.get_view<Real**,Host>();
  auto v3d_imp_h = v3d_imp.get_view<Real***,Host>();
  auto sub_s3d1_imp_h = sub_s3d1_imp.get_view<Real**,Host>();
  auto sub_s3d2_imp_h = sub_s3d2_imp.get_view<Real**,Host>();

  for (int i=0; i<nruns; ++i) {
    // Fill export fields
    ekat::genRandArray(s2d_exp_d,engine,pdf);
    ekat::genRandArray(s3d_exp_d,engine,pdf);
    ekat::genRandArray(v2d_exp_d,engine,pdf);
    ekat::genRandArray(v3d_exp_d,engine,pdf);
    auto G3d_size = G3d_exp.get_header().get_alloc_properties().get_num_scalars();
    ekat::genRandArray(G3d_exp.get_internal_view_data<Real,Host>(),G3d_size,engine,pdf);

    // Set all raw_data to -1 (might be helpful for debugging)
    std::fill_n(raw_data,ncols*num_fields,-1);

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
  using C = scream::physics::Constants<Real>;

  // Some constants
  constexpr int ncols = 4;
  constexpr int nlevs = 72;
  constexpr int nruns = 10;

  // Create a comm
  ekat::Comm comm (MPI_COMM_WORLD);

  // The random numbers generator
  auto engine = setup_random_test(&comm);
  RPDF pdf(0.0,1.0);

  // Create a grid
  auto grid = create_point_grid("my_grid",ncols*comm.size(), nlevs, comm);
  const auto grid_name = grid->name();

  // Layouts matching those in AD
  FL scalar2d_layout{ {COL          }, {ncols          } };
  FL vector2d_layout{ {COL, CMP     }, {ncols, 2       } };
  FL scalar3d_layout{ {COL, LEV     }, {ncols,    nlevs} };
  FL vector3d_layout{ {COL, CMP, LEV}, {ncols, 2, nlevs} };

  // Create import fields
  const auto nondim = Units::nondimensional();
  FID surf_latent_flux_id ("surf_latent_flux",   scalar2d_layout, W/(m*m), grid_name);
  FID surf_sens_flux_id   ("surf_sens_flux",     scalar2d_layout, W/(m*m), grid_name);
  FID surf_mom_flux_id    ("surf_mom_flux",      vector2d_layout, W/(m*m), grid_name);
  FID sfc_alb_dir_vis_id  ("sfc_alb_dir_vis",    scalar2d_layout, nondim,  grid_name);
  FID sfc_alb_dir_nir_id  ("sfc_alb_dir_nir",    scalar2d_layout, nondim,  grid_name);
  FID sfc_alb_dif_vis_id  ("sfc_alb_dif_vis",    scalar2d_layout, nondim,  grid_name);
  FID sfc_alb_dif_nir_id  ("sfc_alb_dif_nir",    scalar2d_layout, nondim,  grid_name);
  FID surf_lw_flux_up_id  ("surf_lw_flux_up",    scalar2d_layout, nondim,  grid_name);

  // Create necessary fields for export. Tracers qc and qr are unnecessary, but
  // are included to verify that subviewed fields (qv) are correctly handled
  FID T_mid_id           ("T_mid",           scalar3d_layout, K,      grid_name);
  FID p_mid_id           ("p_mid",           scalar3d_layout, Pa,     grid_name);
  FID horiz_winds_id     ("horiz_winds",     vector3d_layout, m/s,    grid_name);
  FID pseudo_density_id  ("pseudo_density",  scalar3d_layout, Pa,     grid_name);
  FID qv_id              ("qv",              scalar3d_layout, nondim, grid_name);
  FID precip_liq_surf_id ("precip_liq_surf", scalar2d_layout, m/s,    grid_name);
  FID precip_ice_surf_id ("precip_ice_surf", scalar2d_layout, m/s,    grid_name);
  FID sfc_flux_dir_nir_id("sfc_flux_dir_nir", scalar2d_layout, W/(m*m), grid_name);
  FID sfc_flux_dir_vis_id("sfc_flux_dir_vis", scalar2d_layout, W/(m*m), grid_name);
  FID sfc_flux_dif_nir_id("sfc_flux_dif_nir", scalar2d_layout, W/(m*m), grid_name);
  FID sfc_flux_dif_vis_id("sfc_flux_dif_vis", scalar2d_layout, W/(m*m), grid_name);
  FID sfc_flux_sw_net_id ("sfc_flux_sw_net",  scalar2d_layout, W/(m*m), grid_name);
  FID sfc_flux_lw_dn_id  ("sfc_flux_lw_dn",   scalar2d_layout, W/(m*m), grid_name);

  // NOTE: if you add fields above, you will have to modify these counters too.
  const int num_cpl_imports    = 30;
  const int num_scream_imports = 9;
  const int num_cpl_exports    = 36;

  // Register fields and tracer group in a FieldManager
  auto fm = std::make_shared<FieldManager> (grid);
  fm->registration_begins();
  fm->register_field(FR{surf_latent_flux_id});
  fm->register_field(FR{surf_sens_flux_id});
  fm->register_field(FR{surf_mom_flux_id});
  fm->register_field(FR{sfc_alb_dir_vis_id});
  fm->register_field(FR{sfc_alb_dir_nir_id});
  fm->register_field(FR{sfc_alb_dif_vis_id});
  fm->register_field(FR{sfc_alb_dif_nir_id});
  fm->register_field(FR{surf_lw_flux_up_id});
  fm->register_field(FR{T_mid_id});
  fm->register_field(FR{p_mid_id});
  fm->register_field(FR{horiz_winds_id});
  fm->register_field(FR{pseudo_density_id});
  fm->register_field(FR{qv_id,"tracers"});
  fm->register_field(FR{precip_liq_surf_id});
  fm->register_field(FR{precip_ice_surf_id});
  fm->register_field(FR{sfc_flux_dir_nir_id});
  fm->register_field(FR{sfc_flux_dir_vis_id});
  fm->register_field(FR{sfc_flux_dif_nir_id});
  fm->register_field(FR{sfc_flux_dif_vis_id});
  fm->register_field(FR{sfc_flux_sw_net_id});
  fm->register_field(FR{sfc_flux_lw_dn_id});

  fm->register_group(GR("tracers", grid_name ,Bundling::Required));
  fm->registration_ends();

  // Create alias to field views
  auto surf_latent_flux_f = fm->get_field(surf_latent_flux_id);
  auto surf_sens_flux_f   = fm->get_field(surf_sens_flux_id);
  auto surf_mom_flux_f    = fm->get_field(surf_mom_flux_id);
  auto sfc_alb_dir_vis_f  = fm->get_field(sfc_alb_dir_vis_id);
  auto sfc_alb_dir_nir_f  = fm->get_field(sfc_alb_dir_nir_id);
  auto sfc_alb_dif_vis_f  = fm->get_field(sfc_alb_dif_vis_id);
  auto sfc_alb_dif_nir_f  = fm->get_field(sfc_alb_dif_nir_id);
  auto surf_lw_flux_up_f  = fm->get_field(surf_lw_flux_up_id);
  auto T_mid_f            = fm->get_field(T_mid_id);
  auto p_mid_f            = fm->get_field(p_mid_id);
  auto horiz_winds_f      = fm->get_field(horiz_winds_id);
  auto pseudo_density_f   = fm->get_field(pseudo_density_id);
  auto qv_f               = fm->get_field(qv_id);
  auto precip_liq_surf_f  = fm->get_field(precip_liq_surf_id);
  auto precip_ice_surf_f  = fm->get_field(precip_ice_surf_id);
  auto sfc_flux_dir_nir_f = fm->get_field(sfc_flux_dir_nir_id);
  auto sfc_flux_dir_vis_f = fm->get_field(sfc_flux_dir_vis_id);
  auto sfc_flux_dif_nir_f = fm->get_field(sfc_flux_dif_nir_id);
  auto sfc_flux_dif_vis_f = fm->get_field(sfc_flux_dif_vis_id);
  auto sfc_flux_sw_net_f  = fm->get_field(sfc_flux_sw_net_id);
  auto sfc_flux_lw_dn_f   = fm->get_field(sfc_flux_lw_dn_id);

  auto group = fm->get_field_group("tracers");
  const auto& Q_name = group.m_bundle->get_header().get_identifier().name();
  auto Q = fm->get_field(Q_name);

  auto surf_latent_flux_d = surf_latent_flux_f.get_view<Real*>();
  auto surf_sens_flux_d   = surf_sens_flux_f.get_view<Real*>();
  auto surf_mom_flux_d    = surf_mom_flux_f.get_view<Real**>();
  auto sfc_alb_dir_vis_d  = sfc_alb_dir_vis_f.get_view<Real*>();
  auto sfc_alb_dir_nir_d  = sfc_alb_dir_nir_f.get_view<Real*>();
  auto sfc_alb_dif_vis_d  = sfc_alb_dif_vis_f.get_view<Real*>();
  auto sfc_alb_dif_nir_d  = sfc_alb_dif_nir_f.get_view<Real*>();
  auto surf_lw_flux_up_d  = surf_lw_flux_up_f.get_view<Real*>();
  auto T_mid_d            = T_mid_f.get_view<Real**>();
  auto p_mid_d            = p_mid_f.get_view<Real**>();
  auto horiz_winds_d      = horiz_winds_f.get_view<Real***>();
  auto pseudo_density_d   = pseudo_density_f.get_view<Real**>();
  auto qv_d               = qv_f.get_view<Real**>();
  auto precip_liq_surf_d  = precip_liq_surf_f.get_view<Real*>();
  auto precip_ice_surf_d  = precip_ice_surf_f.get_view<Real*>();
  auto sfc_flux_dir_nir_d = sfc_flux_dir_nir_f.get_view<Real*>();
  auto sfc_flux_dir_vis_d = sfc_flux_dir_vis_f.get_view<Real*>();
  auto sfc_flux_dif_nir_d = sfc_flux_dif_nir_f.get_view<Real*>();
  auto sfc_flux_dif_vis_d = sfc_flux_dif_vis_f.get_view<Real*>();
  auto sfc_flux_sw_net_d  = sfc_flux_sw_net_f.get_view<Real*>();
  auto sfc_flux_lw_dn_d   = sfc_flux_lw_dn_f.get_view<Real*>();

  auto surf_latent_flux_h = surf_latent_flux_f.get_view<Real*,Host>();
  auto surf_sens_flux_h   = surf_sens_flux_f.get_view<Real*,Host>();
  auto surf_mom_flux_h    = surf_mom_flux_f.get_view<Real**,Host>();
  auto sfc_alb_dir_vis_h  = sfc_alb_dir_vis_f.get_view<Real*,Host>();
  auto sfc_alb_dir_nir_h  = sfc_alb_dir_nir_f.get_view<Real*,Host>();
  auto sfc_alb_dif_vis_h  = sfc_alb_dif_vis_f.get_view<Real*,Host>();
  auto sfc_alb_dif_nir_h  = sfc_alb_dif_nir_f.get_view<Real*,Host>();
  auto surf_lw_flux_up_h  = surf_lw_flux_up_f.get_view<Real*,Host>();
  auto T_mid_h            = T_mid_f.get_view<Real**,Host>();
  auto p_mid_h            = p_mid_f.get_view<Real**,Host>();
  auto pseudo_density_h   = pseudo_density_f.get_view<Real**,Host>();
  auto horiz_winds_h      = horiz_winds_f.get_view<Real***,Host>();
  auto qv_h               = qv_f.get_view<Real**,Host>();
  auto precip_liq_surf_h  = precip_liq_surf_f.get_view<Real*,Host>();
  auto precip_ice_surf_h  = precip_ice_surf_f.get_view<Real*,Host>();
  auto sfc_flux_dir_nir_h = sfc_flux_dir_nir_f.get_view<Real*,Host>();
  auto sfc_flux_dir_vis_h = sfc_flux_dir_vis_f.get_view<Real*,Host>();
  auto sfc_flux_dif_nir_h = sfc_flux_dif_nir_f.get_view<Real*,Host>();
  auto sfc_flux_dif_vis_h = sfc_flux_dif_vis_f.get_view<Real*,Host>();
  auto sfc_flux_sw_net_h  = sfc_flux_sw_net_f.get_view<Real*,Host>();
  auto sfc_flux_lw_dn_h   = sfc_flux_lw_dn_f.get_view<Real*,Host>();

  // Create SC object and set number of import/export fields
  control::SurfaceCoupling coupler(fm);
  coupler.set_num_fields(num_cpl_imports, num_scream_imports, num_cpl_exports);

  // Register fields in the coupler. These match the scr_names_x2a/a2x from
  // scream_cpl_indices.F90. When radiation is fully implemented, RRTMGP
  // fields will be replaced with correct fields. For fluxes, we flip the sign.
  coupler.register_import("unused",           0);
  coupler.register_import("unused",           1);
  coupler.register_import("unused",           2);
  coupler.register_import("sfc_alb_dir_vis",  3);
  coupler.register_import("sfc_alb_dir_nir",  4);
  coupler.register_import("sfc_alb_dif_vis",  5);
  coupler.register_import("sfc_alb_dif_nir",  6);
  coupler.register_import("unused",           7);
  coupler.register_import("unused",           8);
  coupler.register_import("unused",           9);
  coupler.register_import("unused",           10);
  coupler.register_import("unused",           11);
  coupler.register_import("unused",           12);
  coupler.register_import("unused",           13);
  coupler.register_import("unused",           14);
  coupler.register_import("unused",           15);
  coupler.register_import("unused",           16);
  coupler.register_import("unused",           17);
  coupler.register_import("unused",           18);
  coupler.register_import("surf_mom_flux",    19, 0);
  coupler.register_import("surf_mom_flux",    20, 1);
  coupler.register_import("unused",           21);
  coupler.register_import("surf_sens_flux",   22);
  coupler.register_import("surf_lw_flux_up",  23);
  coupler.register_import("surf_latent_flux", 24);
  coupler.register_import("unused",           25);
  coupler.register_import("unused",           26);
  coupler.register_import("unused",           27);
  coupler.register_import("unused",           28);
  coupler.register_import("unused",           29);

  coupler.register_export("Sa_z",            0);
  coupler.register_export("set_zero",        1);
  coupler.register_export("horiz_winds",     2, 0);
  coupler.register_export("horiz_winds",     3, 1);
  coupler.register_export("T_mid",           4);
  coupler.register_export("Sa_ptem",         5);
  coupler.register_export("qv",              6);
  coupler.register_export("p_mid",           7);
  coupler.register_export("Sa_dens",         8);
  coupler.register_export("set_zero",        9);
  coupler.register_export("set_zero",        10);
  coupler.register_export("set_zero",        11);
  coupler.register_export("set_zero",        12);
  coupler.register_export("set_zero",        13);
  coupler.register_export("set_zero",        14);
  coupler.register_export("Faxa_rainl",      15);  
  coupler.register_export("Faxa_snowl",      16);
  coupler.register_export("sfc_flux_lw_dn",  17);
  coupler.register_export("sfc_flux_dir_nir",18);
  coupler.register_export("sfc_flux_dir_vis",19);
  coupler.register_export("sfc_flux_dif_nir",20);
  coupler.register_export("sfc_flux_dif_vis",21);
  coupler.register_export("sfc_flux_sw_net", 22);
  coupler.register_export("set_zero",        23);
  coupler.register_export("set_zero",        24);
  coupler.register_export("set_zero",        25);
  coupler.register_export("set_zero",        26);
  coupler.register_export("set_zero",        27);
  coupler.register_export("set_zero",        28);
  coupler.register_export("set_zero",        29);
  coupler.register_export("set_zero",        30);
  coupler.register_export("set_zero",        31);
  coupler.register_export("set_zero",        32);
  coupler.register_export("set_zero",        33);
  coupler.register_export("set_zero",        34);
  coupler.register_export("set_zero",        35);


  // Complete setup of importer/exporter, providing raw_data
  double* import_raw_data = new double[ncols*num_cpl_imports];
  double* export_raw_data = new double[ncols*num_cpl_exports];
  coupler.registration_ends(import_raw_data, export_raw_data);

  for (int i=0; i<nruns; ++i) {

    // Set import field views to 0
    surf_latent_flux_f.deep_copy(0.0);
    surf_sens_flux_f.deep_copy(0.0);
    surf_mom_flux_f.deep_copy(0.0);
    sfc_alb_dir_vis_f.deep_copy(0.0);
    sfc_alb_dir_nir_f.deep_copy(0.0);
    sfc_alb_dif_vis_f.deep_copy(0.0);
    sfc_alb_dif_nir_f.deep_copy(0.0);
    surf_lw_flux_up_f.deep_copy(0.0);

    // Fill views needed in the export with random values
    ekat::genRandArray(T_mid_d,engine,pdf);
    ekat::genRandArray(p_mid_d,engine,pdf);
    ekat::genRandArray(horiz_winds_d,engine,pdf);
    ekat::genRandArray(pseudo_density_d,engine,pdf);
    ekat::genRandArray(precip_liq_surf_d,engine,pdf);
    ekat::genRandArray(precip_ice_surf_d,engine,pdf);
    ekat::genRandArray(sfc_flux_lw_dn_d,engine,pdf);
    ekat::genRandArray(sfc_flux_dir_nir_d,engine,pdf);
    ekat::genRandArray(sfc_flux_dir_vis_d,engine,pdf);
    ekat::genRandArray(sfc_flux_dif_nir_d,engine,pdf);
    ekat::genRandArray(sfc_flux_dif_vis_d,engine,pdf);
    ekat::genRandArray(sfc_flux_sw_net_d,engine,pdf);
    auto Q_size = Q.get_header().get_alloc_properties().get_num_scalars();
    ekat::genRandArray(Q.get_internal_view_data<Real,Host>(),Q_size,engine,pdf);

    // Fill import_raw_data with random values
    for (int icol=0; icol<ncols*num_cpl_imports; ++icol) {
      import_raw_data[icol] = pdf(engine);
    }

    // Set all export_raw_data to -1 (might be helpful for debugging)
    std::fill_n(export_raw_data, num_cpl_exports*ncols, -1);

    // Perform import
    coupler.do_import();

    // Perform export
    coupler.do_export();

    // Sync host to device
    surf_latent_flux_f.sync_to_host();
    surf_sens_flux_f.sync_to_host();
    surf_mom_flux_f.sync_to_host();
    sfc_alb_dir_vis_f.sync_to_host();
    sfc_alb_dir_nir_f.sync_to_host();
    sfc_alb_dif_vis_f.sync_to_host();
    sfc_alb_dif_nir_f.sync_to_host();
    surf_lw_flux_up_f.sync_to_host();
    T_mid_f.sync_to_host();
    p_mid_f.sync_to_host();
    horiz_winds_f.sync_to_host();
    pseudo_density_f.sync_to_host();
    precip_liq_surf_f.sync_to_host();
    precip_ice_surf_f.sync_to_host();
    sfc_flux_lw_dn_f.sync_to_host();
    sfc_flux_dir_nir_f.sync_to_host();
    sfc_flux_dir_vis_f.sync_to_host();
    sfc_flux_dif_nir_f.sync_to_host();
    sfc_flux_dif_vis_f.sync_to_host();
    sfc_flux_sw_net_f.sync_to_host();
    Q.sync_to_host();

    // Check values
    for (int icol=0; icol<ncols; ++icol) {

      // Imports

      REQUIRE (sfc_alb_dir_vis_h (icol)    ==  import_raw_data[3 + icol*num_cpl_imports]);  // 1st scream import (4th cpl import)
      REQUIRE (sfc_alb_dir_nir_h (icol)    ==  import_raw_data[4 + icol*num_cpl_imports]);  // 2nd scream import (5th cpl import)
      REQUIRE (sfc_alb_dif_vis_h (icol)    ==  import_raw_data[5 + icol*num_cpl_imports]);  // 3rd scream import (6th cpl import)
      REQUIRE (sfc_alb_dif_nir_h (icol)    ==  import_raw_data[6 + icol*num_cpl_imports]);  // 4th scream import (7th cpl import)
      // These imports should have opposite signs
      REQUIRE (surf_mom_flux_h   (icol, 0) == -import_raw_data[19 + icol*num_cpl_imports]); // 5th scream import (20th cpl import)
      REQUIRE (surf_mom_flux_h   (icol, 1) == -import_raw_data[20 + icol*num_cpl_imports]); // 6th scream import (21st cpl import)
      REQUIRE (surf_sens_flux_h  (icol)    == -import_raw_data[22 + icol*num_cpl_imports]); // 7th scream import (23rd cpl import)
      REQUIRE (surf_lw_flux_up_h (icol)    == -import_raw_data[23 + icol*num_cpl_imports]); // 8th scream import (24th cpl import)
      REQUIRE (surf_latent_flux_h(icol)    == -import_raw_data[24 + icol*num_cpl_imports]); // 9th scream import (24th cpl import)

      // Exports

      // These exports are direct values from a scream field
      REQUIRE (export_raw_data[2 + icol*num_cpl_exports]  == horiz_winds_h    (icol, 0, nlevs-1)); // 3rd export
      REQUIRE (export_raw_data[3 + icol*num_cpl_exports]  == horiz_winds_h    (icol, 1, nlevs-1)); // 4th export
      REQUIRE (export_raw_data[4 + icol*num_cpl_exports]  == T_mid_h          (icol,    nlevs-1)); // 5th export
      REQUIRE (export_raw_data[6 + icol*num_cpl_exports]  == qv_h             (icol,    nlevs-1)); // 7th export
      REQUIRE (export_raw_data[7 + icol*num_cpl_exports]  == p_mid_h          (icol,    nlevs-1)); // 8th export
      REQUIRE (export_raw_data[15 + icol*num_cpl_exports] == C::RHO_H2O*precip_liq_surf_h(icol));  // 16th export
      REQUIRE (export_raw_data[16 + icol*num_cpl_exports] == C::RHO_H2O*precip_ice_surf_h(icol));  // 17th export
      REQUIRE (export_raw_data[17 + icol*num_cpl_exports] == sfc_flux_lw_dn_h(icol));              // 18th export
      REQUIRE (export_raw_data[18 + icol*num_cpl_exports] == sfc_flux_dir_nir_h(icol));            // 19th export
      REQUIRE (export_raw_data[19 + icol*num_cpl_exports] == sfc_flux_dir_vis_h(icol));            // 20th export
      REQUIRE (export_raw_data[20 + icol*num_cpl_exports] == sfc_flux_dif_nir_h(icol));            // 21st export
      REQUIRE (export_raw_data[21 + icol*num_cpl_exports] == sfc_flux_dif_vis_h(icol));            // 22nd export
      REQUIRE (export_raw_data[22 + icol*num_cpl_exports] == sfc_flux_sw_net_h(icol));             // 23rd export

      // These exports should be set to 0
      REQUIRE (export_raw_data[1  + icol*num_cpl_exports] == 0); // 2nd export
      REQUIRE (export_raw_data[9  + icol*num_cpl_exports] == 0); // 10th export
      REQUIRE (export_raw_data[10 + icol*num_cpl_exports] == 0); // 11th export
      REQUIRE (export_raw_data[11 + icol*num_cpl_exports] == 0); // 12th export
      REQUIRE (export_raw_data[12 + icol*num_cpl_exports] == 0); // 13th export
      REQUIRE (export_raw_data[13 + icol*num_cpl_exports] == 0); // 14th export
      REQUIRE (export_raw_data[14 + icol*num_cpl_exports] == 0); // 15th export
      REQUIRE (export_raw_data[23 + icol*num_cpl_exports] == 0); // 24th export
      REQUIRE (export_raw_data[24 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[25 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[26 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[27 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[28 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[29 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[30 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[31 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[32 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[33 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[34 + icol*num_cpl_exports] == 0);
      REQUIRE (export_raw_data[35 + icol*num_cpl_exports] == 0); // 35th export
    }
  }

  // Clean up
  delete[] export_raw_data;
  delete[] import_raw_data;
}


TEST_CASE ("do_initial_export")
{
  /*
   * This test performs 2 exports, one with init_phase=true, to
   * test that feature.
   */

  // Some namespaces/aliases
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FL = FieldLayout;
  using FID = FieldIdentifier;
  using RPDF = std::uniform_real_distribution<Real>;

  // Some constants
  constexpr int ncols = 4;
  constexpr int nlevs = 8;
  constexpr int nruns = 10;

  // Create a comm
  ekat::Comm comm (MPI_COMM_WORLD);

  // The random numbers generator
  auto engine = setup_random_test(&comm);
  RPDF pdf(0.0,1.0);

  // Create a grid
  auto grid = create_point_grid("my grid",ncols*comm.size(), nlevs, comm);

  // Create some field ids, and register them in a field manager
  FID f1_id("f1",FL{{COL},{ncols}},Pa,grid->name());
  FID f2_id("f2",FL{{COL},{ncols}},Pa,grid->name());

  // NOTE: if you add fields above, you will have to modify these counters too.
  const int num_fields = 2;

  // Keep two separate fms, so we can compare original and final fields.
  auto fm = std::make_shared<FieldManager> (grid);
  fm->registration_begins();
  fm->register_field(f1_id);
  fm->register_field(f2_id);
  fm->registration_ends();

  // Create a raw array big enough to contain all the 2d data for import/export
  double* raw_data = new double[ncols*num_fields];

  // Repeat experiment N times: fill export fields, export init fields, check values, export all, check values
  auto f1_exp = fm->get_field(f1_id);
  auto f2_exp = fm->get_field(f2_id);
  auto f1_exp_d = f1_exp.get_view<Real*>();
  auto f2_exp_d = f2_exp.get_view<Real*>();
  auto f1_exp_h = f1_exp.get_view<Real*,Host>();
  auto f2_exp_h = f2_exp.get_view<Real*,Host>();

  for (int i=0; i<nruns; ++i) {
    // Create two SC objects, to import and export
    control::SurfaceCoupling exporter(fm);
    exporter.set_num_fields(0,num_fields);

    // Register fields in the exporter. Set f2 to not export during init phase
    exporter.register_export("f1",0);
    exporter.register_export("f2",1,-1,false);

    // Set all raw_data to -1
    std::fill_n(raw_data,ncols*num_fields,-1);

    // Complete setup of exporter. This needs to be done in the run loop for this test since
    // this is the last place which
    exporter.registration_ends(nullptr,raw_data);

    // Fill export fields
    ekat::genRandArray(f1_exp_d,engine,pdf);
    ekat::genRandArray(f2_exp_d,engine,pdf);

    // Perform export with init_phase=true
    exporter.do_export(true);

    // Check that only f1 was exported
    f1_exp.sync_to_host();
    f2_exp.sync_to_host();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE (raw_data[0 + icol*num_fields] == f1_exp_h(icol));
      REQUIRE (raw_data[1 + icol*num_fields] == -1);
    }

    // Set all raw_data back to -1
    std::fill_n(raw_data,ncols*num_fields,-1);

    // Perform export with init_phase=false (default)
    exporter.do_export();

    // Check that both f1 and f2 were exported
    f1_exp.sync_to_host();
    f2_exp.sync_to_host();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE (raw_data[0 + icol*num_fields] == f1_exp_h(icol));
      REQUIRE (raw_data[1 + icol*num_fields] == f2_exp_h(icol));
    }
  }

  // Clean up
  delete[] raw_data;
}
