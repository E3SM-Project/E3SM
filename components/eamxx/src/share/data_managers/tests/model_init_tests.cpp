#include <catch2/catch.hpp>

#include "share/data_managers/model_init.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_parameter_list.hpp>

namespace scream {

TEST_CASE ("model_init_constants_and_copy", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);

  const std::string gn = "point_grid";
  const int ncols = 4*comm.size();
  const int nlevs = 3;
  const int ncmp  = 2;

  auto grid = create_point_grid(gn,ncols,nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>(grid);

  auto fm = std::make_shared<FieldManager>(gm);

  FieldIdentifier a_id ("A", {{COL,LEV},    {ncols,nlevs}},       none, gn);
  FieldIdentifier z_id ("Z", {{COL,LEV},    {ncols,nlevs}},       none, gn);
  // Rank-2 vector field, with CMP as the fastest-striding dim (idim=1):
  // exercises the "fast path" branch of set_constant_field.
  FieldIdentifier v_id ("V", {{COL,CMP},    {ncols,ncmp}},        none, gn);
  // Rank-3 vector field: exercises the generic per-component branch.
  FieldIdentifier w_id ("W", {{COL,CMP,LEV},{ncols,ncmp,nlevs}},  none, gn);

  fm->register_field(FR{a_id});
  fm->register_field(FR{z_id});
  fm->register_field(FR{v_id});
  fm->register_field(FR{w_id});
  fm->register_group(GroupRequest("STARTUP",gn));
  fm->registration_ends();

  for (const auto& name : {"A","Z","V","W"}) {
    fm->add_to_group(name,gn,"STARTUP");
  }

  ekat::ParameterList params("initial_conditions");
  params.set<double>("A",1.5);
  params.set<std::string>("Z","A");
  params.set<std::vector<double>>("V",{2.0,3.0});
  params.set<std::vector<double>>("W",{4.0,5.0});

  ModelInit model_init(params,comm);

  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Initial);

  auto a = fm->get_field("A",gn);
  auto z = fm->get_field("Z",gn);
  auto v = fm->get_field("V",gn);
  auto w = fm->get_field("W",gn);

  // A is a constant, Z is a copy of A, so they should be equal
  REQUIRE (a.get_header().get_tracking().get_time_stamp()==t0);
  REQUIRE (z.get_header().get_tracking().get_time_stamp()==t0);
  REQUIRE (views_are_equal(a,z));

  Field a_check(a_id); a_check.allocate_view(); a_check.deep_copy(1.5);
  REQUIRE (views_are_equal(a,a_check));

  Field v_check(v_id); v_check.allocate_view();
  v_check.get_component(0).deep_copy(2.0);
  v_check.get_component(1).deep_copy(3.0);
  REQUIRE (views_are_equal(v,v_check));
  REQUIRE (v.get_header().get_tracking().get_time_stamp()==t0);

  Field w_check(w_id); w_check.allocate_view();
  w_check.get_component(0).deep_copy(4.0);
  w_check.get_component(1).deep_copy(5.0);
  REQUIRE (views_are_equal(w,w_check));
  REQUIRE (w.get_header().get_tracking().get_time_stamp()==t0);
}

} // namespace scream
