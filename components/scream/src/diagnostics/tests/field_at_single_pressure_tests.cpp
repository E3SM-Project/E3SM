#include "catch2/catch.hpp"

#include "ekat/ekat_pack_utils.hpp"

#include "diagnostics/field_at_single_pressure.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("field_at_single_pressure")
{

  // Test that output at a single pressure level works as expected.
  // For this test we set a field "M" to be defined as 100*i + k,
  // where i=column and k=level
  //
  // We then set the pressure levels to be 100*k where again k=level.
  //
  // Lastly we define a pressure_at_single_level to be for p=150, thus
  // the output should be (M_1+M_2)/2, or halfway between level 1 and level
  // 2 of the data.  Given the formula above the output should be exactly
  //   100*i + (1+2)/2 = 100*i + 1.5

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = scream::setup_random_test(&comm);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;  // Note, we need at least 3 levels for the test to work
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("Point Grid");

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create input fields
  const auto units = ekat::units::Units::invalid();

  FieldIdentifier fid_mid ("V_mid",FL({COL,LEV},{ncols,nlevs}),units,grid->name());
  Field f_mid (fid_mid);
  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_mid.allocate_view();
  f_mid.get_header().get_tracking().update_time_stamp(t0);

  FieldIdentifier fid_int ("V_int",FL({COL,ILEV},{ncols,nlevs+1}),units,grid->name());
  Field f_int (fid_int);
  f_int.get_header().get_alloc_properties().request_allocation(packsize);
  f_int.allocate_view();
  f_int.get_header().get_tracking().update_time_stamp(t0);

  ekat::ParameterList params_mid, params_int;
  params_mid.set("Field Name",f_mid.name());
  params_mid.set("Field Units",fid_mid.get_units());
  params_mid.set("Field Layout",fid_mid.get_layout());
  params_mid.set("Grid Name",fid_mid.get_grid_name());
  params_mid.set<int>("Field Target Pressure",150);
  params_int.set("Field Name",f_int.name());
  params_int.set("Field Units",fid_int.get_units());
  params_int.set("Field Layout",fid_int.get_layout());
  params_int.set("Grid Name",fid_int.get_grid_name());
  params_int.set<int>("Field Target Pressure",150);

  auto diag_mid = std::make_shared<FieldAtSinglePressure>(comm,params_mid);
  diag_mid->set_grids(gm);
  diag_mid->set_required_field(f_mid);
  auto diag_int = std::make_shared<FieldAtSinglePressure>(comm,params_int);
  diag_int->set_grids(gm);
  diag_int->set_required_field(f_int);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag_mid->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(packsize);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag_mid->set_required_field(f);
    input_fields.emplace(name,f);
  }
  for (const auto& req : diag_int->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(packsize);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag_int->set_required_field(f);
    input_fields.emplace(name,f);
  }

  Field p_mid_f = input_fields["p_mid"];
  Field p_int_f = input_fields["p_int"];
  //Fill data to interpolate
  auto f_mid_v   = f_mid.get_view<Real**>();
  auto p_mid_v   = p_mid_f.get_view<Real**>();
  auto f_mid_v_h = Kokkos::create_mirror_view(f_mid_v);
  auto p_mid_v_h = Kokkos::create_mirror_view(p_mid_v);
  auto f_int_v   = f_int.get_view<Real**>();
  auto p_int_v   = p_int_f.get_view<Real**>();
  auto f_int_v_h = Kokkos::create_mirror_view(f_int_v);
  auto p_int_v_h = Kokkos::create_mirror_view(p_int_v);
  for (int ilev=0; ilev<nlevs; ilev++){
    for (int icol=0; icol<ncols; icol++){
      f_mid_v_h(icol,ilev) = icol*100 + ilev;
      p_mid_v_h(icol,ilev) = 100*ilev;
      f_int_v_h(icol,ilev) = icol*100 + ilev;
      p_int_v_h(icol,ilev) = 100*ilev;
      f_int_v_h(icol,ilev+1) = icol*100 + ilev+1;
      p_int_v_h(icol,ilev+1) = 100*(ilev+1);
    }
  }
  f_mid.sync_to_dev();
  p_mid_f.sync_to_dev();
  f_int.sync_to_dev();
  p_int_f.sync_to_dev();
  

  diag_mid->initialize(t0,RunType::Initial);
  diag_int->initialize(t0,RunType::Initial);
  
  // Run diagnostics
  diag_mid->compute_diagnostic();
  diag_int->compute_diagnostic();

  auto d_mid = diag_mid->get_diagnostic();
  d_mid.sync_to_host();
  auto d_mid_v = d_mid.get_view<const Real*,Host>();
  auto d_int = diag_int->get_diagnostic();
  d_int.sync_to_host();
  auto d_int_v = d_int.get_view<const Real*,Host>();
  
  for (int icol=0; icol<ncols; ++icol) {
    REQUIRE (d_mid_v(icol)==icol*100 + 1.5);
    REQUIRE (d_int_v(icol)==icol*100 + 1.5);
  }
  
}

} // namespace scream
