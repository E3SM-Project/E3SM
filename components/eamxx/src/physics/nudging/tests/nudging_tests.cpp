#include "catch2/catch.hpp"
#include "physics/nudging/atmosphere_nudging.hpp"

namespace {

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
  
  FieldIdentifier fid_mid ("T_mid",FL({COL,LEV},{ncols,nlevs}),units,grid->name());
  Field f_mid (fid_mid);
  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_mid.allocate_view();
  f_mid.get_header().get_tracking().update_time_stamp(t0);

  ekat::ParameterList params_mid;
  auto nudging_mid = std::make_shared<NUDGING>(comm,params_mid);
  nudging_mid->set_grids(gm);
  //nudging_mid->set_required_field(f_mid);

  //fill data

  //initialize
  nudging_mid->initialize();

  //run
  nudging_mid->run();

}


} // empty namespace
