#include "catch2/catch.hpp"
#include "physics/nudging/atmosphere_nudging.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

using namespace scream;

//namespace {

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


TEST_CASE("nudging") {
  std::cout<<"Get into nudging test"<<std::endl;

  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT = KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Pack>;
  using view_1d_const = typename KT::template view_1d<const Pack>;

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;
  std::cout<<"packsize: "<<packsize<<std::endl;

  ekat::Comm comm(MPI_COMM_WORLD);
  
  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;  // Note, we need at least 3 levels for the test to work
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("Physics");

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});
  
  const auto units = ekat::units::Units::invalid();

  std::cout<<"Columns outside: "<<ncols<<std::endl;
  std::cout<<"Levs outside: "<<nlevs<<std::endl;

  FieldIdentifier fid_mid ("T_mid",FL({COL,LEV},{ncols,nlevs}),units,grid->name());
  Field f_mid (fid_mid);
  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_mid.allocate_view();
  f_mid.get_header().get_tracking().update_time_stamp(t0);

  ekat::ParameterList params_mid;
  auto nudging_mid = std::make_shared<NUDGING>(comm,params_mid);
  nudging_mid->set_grids(gm);
  nudging_mid->set_required_field(f_mid);

  //fill data
  auto f_mid_v_h   = f_mid.get_view<Real**, Host>();
  for (int ilev=0; ilev<nlevs; ilev++){ 
    for (int icol=0; icol<ncols; icol++){
      f_mid_v_h(icol,ilev) = (-icol)*100 + ilev;
    }
  }
  f_mid.sync_to_dev();

  //Fill data to interpolate   

  //initialize
  nudging_mid->initialize(t0,RunType::Initial);

  //run
  nudging_mid->run(100);
  f_mid.sync_to_host();

  //Now check that I was able to nudge it by a value of 1
  for (int ilev=0; ilev<nlevs; ilev++){ 
    for (int icol=0; icol<ncols; icol++){
      REQUIRE(f_mid_v_h(icol,ilev) == (-icol)*100 + ilev + 1);
    }
  }

  
}


//} // empty namespace
