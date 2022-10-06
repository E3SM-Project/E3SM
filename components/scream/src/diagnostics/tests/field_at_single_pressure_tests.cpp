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
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = scream::setup_random_test(&comm);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("Point Grid");

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create input fields
  const auto units = ekat::units::Units::invalid();

  FieldIdentifier fid_mid ("M",FL({COL,LEV},{ncols,nlevs}),units,grid->name());
  Field f_mid (fid_mid);
  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_mid.allocate_view();
  f_mid.get_header().get_tracking().update_time_stamp(t0);

  //Fill data to interpolate
  auto f_mid_v = f_mid.get_view<Real**,Host>();
  for (int icol=0; icol<ncols; icol++){
    for (int ilev=0; ilev<nlevs; ilev++){
      f_mid_v(icol,ilev) = icol*100 + ilev;
      //f_p_v(icol,ilev) = 50*ilev;
    }
  }
  
  ekat::ParameterList params_mid;
  params_mid.set("Field Name",f_mid.name());
  std::cout<<"f_mid.name(): "<<f_mid.name()<<std::endl;
  params_mid.set("Field Units",fid_mid.get_units());
  params_mid.set("Field Layout",fid_mid.get_layout());
  params_mid.set("Grid Name",fid_mid.get_grid_name());
  params_mid.set<std::string>("Field Target Pressure","500");
  //std::string location_pressure_file = "/usr/workspace/rebassoo/Climate/press_tgt_levels.txt";
  //params_mid.set<std::string>("Field Pressure file",location_pressure_file);

  auto diag_mid = std::make_shared<FieldAtSinglePressure>(comm,params_mid);
  diag_mid->set_grids(gm);

  diag_mid->set_required_field(f_mid);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag_mid->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(packsize);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    //diag_mid->set_required_field(f.get_const());
    diag_mid->set_required_field(f);
    REQUIRE_THROWS(diag_mid->set_computed_field(f));
    input_fields.emplace(name,f);
  }


  //Fill pressure fields
  Field p_mid_f = input_fields["p_mid"];
  auto p_mid_v  = p_mid_f.get_view<Real**>();

  for (int icol=0; icol<ncols; icol++){
    for (int ilev=0; ilev<nlevs; ilev++){
      //f_mid_v(icol,ilev) = icol*100 + ilev;
      p_mid_v(icol,ilev) = 25+50*ilev;
    }
  }
  
  diag_mid->initialize(t0,RunType::Initial);
  
  // Run diagnostics
  diag_mid->compute_diagnostic();

  auto d_mid = diag_mid->get_diagnostic();
  d_mid.sync_to_host();

  auto d_mid_v = d_mid.get_view<const Real**,Host>();
  
  //auto d_int_v = d_int.get_view<const Real*,Host>();
  for (int icol=0; icol<ncols; ++icol) {
    std::cout<<"d_mid_v("<<icol<<"): "<<d_mid_v(icol,0)<<std::endl;
    std::cout<<"d_mid_v("<<icol<<"): "<<d_mid_v(icol,0)<<std::endl;
    //REQUIRE (d_mid_v(icol)==f_mid_v(icol,lev));
    // REQUIRE (d_int_v(icol)==f_int_v(icol,lev));
    }
  
}

} // namespace scream
