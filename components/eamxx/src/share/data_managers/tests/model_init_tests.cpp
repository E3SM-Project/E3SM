#include <catch2/catch.hpp>

#include "share/data_managers/model_init.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_parameter_list.hpp>

#include <any>

namespace scream {

namespace {

struct Fixture {
  std::shared_ptr<PointGrid>    grid;
  std::shared_ptr<FieldManager> fm;
  std::string                   gn;
};

// Create a field manager on a point grid, with ncols columns *per rank*, and the given fields. If add_to_groups=true,
// all fields are added to all the given groups (as the atm driver does with the
// STARTUP/RESTART groups). Otherwise, the groups are created empty.
Fixture make_fixture (const ekat::Comm& comm,
                      const int ncols, const int nlevs,
                      const std::vector<FieldIdentifier>& fids,
                      const std::vector<std::string>& groups,
                      const bool add_to_groups = true)
{
  auto grid = create_point_grid("point_grid",ncols*comm.size(),nlevs,comm);
  const auto gn = grid->name();
  auto gm = std::make_shared<LibraryGridsManager>(grid);
  auto fm = std::make_shared<FieldManager>(gm);

  for (const auto& fid : fids) {
    fm->register_field(FieldRequest(fid));
  }
  for (const auto& g : groups) {
    fm->register_group(GroupRequest(g,gn));
  }
  fm->registration_ends();
  if (add_to_groups) {
    for (const auto& fid : fids) {
      for (const auto& g : groups) {
        fm->add_to_group(fid.name(),gn,g);
      }
    }
  }
  return { grid, fm, gn };
}

// Write a fresh netcdf file, containing (COL,LEV) fields with 'ncol'/'lev' on-file
// dimension names, with every entry of each field set to the corresponding value.
// The ncol dimension is decomposed linearly among the ranks, with ncols columns per rank.
// Optionally, also write some global attributes.
void write_file (const std::string& filename,
                 const std::map<std::string,Real>& vars,
                 const ekat::Comm& comm,
                 const int ncols, const int nlevs,
                 const std::map<std::string,int>& int_atts = {},
                 const std::map<std::string,std::string>& str_atts = {})
{
  using namespace ShortFieldTagsNames;

  scorpio::register_file(filename,scorpio::Write);
  scorpio::define_dim(filename,"ncol",ncols*comm.size());
  scorpio::define_dim(filename,"lev",nlevs);
  for (const auto& [name,value] : vars) {
    scorpio::define_var(filename,name,{"ncol","lev"},"real",false);
  }
  for (const auto& [name,value] : int_atts) {
    scorpio::set_attribute<int>(filename,"GLOBAL",name,value);
  }
  for (const auto& [name,value] : str_atts) {
    scorpio::set_attribute<std::string>(filename,"GLOBAL",name,value);
  }
  scorpio::set_dim_decomp(filename,"ncol",comm.rank()*ncols,ncols);
  scorpio::enddef(filename);

  for (const auto& [name,value] : vars) {
    FieldIdentifier fid(name,{{COL,LEV},{ncols,nlevs}},ekat::units::none,"");
    Field f(fid);
    f.allocate_view();
    f.deep_copy(value);
    f.sync_to_host();
    scorpio::write_var(filename,name,f.get_internal_view_data<const Real,Host>());
  }
  scorpio::release_file(filename);
}

// Check that all entries of f are equal to value
void check_constant (const Field& f, const Real value)
{
  Field f_check(f.get_header().get_identifier());
  f_check.allocate_view();
  f_check.deep_copy(value);
  REQUIRE (views_are_equal(f,f_check));
}

} // anonymous namespace

TEST_CASE ("model_init_constants", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 3;
  const int ncmp  = 2;

  FieldIdentifier a_id ("A", {{COL,LEV},    {ncols,nlevs}},      none, "point_grid");
  FieldIdentifier b_id ("B", {{COL,LEV},    {ncols,nlevs}},      none, "point_grid");
  FieldIdentifier c_id ("C", {{COL,LEV},    {ncols,nlevs}},      none, "point_grid");
  // Rank-2 vector field, with CMP as the fastest-striding dim
  FieldIdentifier v_id ("V", {{COL,CMP},    {ncols,ncmp}},       none, "point_grid");
  // Rank-3 vector field
  FieldIdentifier w_id ("W", {{COL,CMP,LEV},{ncols,ncmp,nlevs}}, none, "point_grid");

  auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{a_id,b_id,c_id,v_id,w_id},{"STARTUP"});

  ekat::ParameterList params("initial_conditions");
  auto& overrides = params.sublist("constant_fields").sublist("overrides");
  auto& defaults  = params.sublist("constant_fields").sublist("defaults");
  overrides.set<double>("A",1.5);
  overrides.set<std::vector<double>>("B",{2.5});
  overrides.set<std::vector<double>>("V",{2.0,3.0});
  // 'overrides' wins over 'defaults'
  overrides.set<double>("C",4.0);
  defaults.set<double>("C",5.0);
  // No IC file: the default is used
  defaults.set<std::vector<double>>("W",{6.0,7.0});

  util::TimeStamp t0(2000,1,1,0,0,0);
  ModelInit model_init(params);
  model_init.run(fm,t0,RunType::Initial);

  for (const auto& name : {"A","B","C","V","W"}) {
    REQUIRE (fm->get_field(name,gn).get_header().get_tracking().get_time_stamp()==t0);
  }

  check_constant(fm->get_field("A",gn),1.5);
  check_constant(fm->get_field("B",gn),2.5);
  check_constant(fm->get_field("C",gn),4.0);

  auto v = fm->get_field("V",gn);
  Field v_check(v_id); v_check.allocate_view();
  v_check.get_component(0).deep_copy(2.0);
  v_check.get_component(1).deep_copy(3.0);
  REQUIRE (views_are_equal(v,v_check));

  auto w = fm->get_field("W",gn);
  Field w_check(w_id); w_check.allocate_view();
  w_check.get_component(0).deep_copy(6.0);
  w_check.get_component(1).deep_copy(7.0);
  REQUIRE (views_are_equal(w,w_check));
}

TEST_CASE ("model_init_constants_errors", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 3;

  std::map<std::string,FieldIdentifier> ids = {
    {"A", FieldIdentifier("A", {{COL,LEV}, {ncols,nlevs}}, none, "point_grid")},
    {"I", FieldIdentifier("I", {{COL},     {ncols}},       none, "point_grid", DataType::IntType)},
    {"V", FieldIdentifier("V", {{COL,CMP}, {ncols,2}},     none, "point_grid")}
  };

  util::TimeStamp t0(2000,1,1,0,0,0);

  // Run ModelInit on a field manager with the given field, which is set to a constant
  auto run_with_override = [&](const std::string& name, const std::vector<double>& values) {
    auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{ids.at(name)},{"STARTUP"});
    ekat::ParameterList params("initial_conditions");
    params.sublist("constant_fields").sublist("overrides").set<std::vector<double>>(name,values);
    ModelInit model_init(params);
    model_init.run(fm,t0,RunType::Initial);
    return fm;
  };

  // Integer fields can be set to integral values only
  auto fm = run_with_override("I",{3.0});
  auto i = fm->get_field("I","point_grid");
  i.sync_to_host();
  auto i_h = i.get_view<const int*,Host>();
  for (int icol=0; icol<ncols; ++icol) {
    REQUIRE (i_h(icol)==3);
  }
  REQUIRE_THROWS (run_with_override("I",{2.5}));

  // Wrong number of components
  REQUIRE_THROWS (run_with_override("V",{1.0,2.0,3.0}));
  REQUIRE_THROWS (run_with_override("A",{1.0,2.0}));
}

TEST_CASE ("model_init_startup_from_file", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4;
  const int nlevs = 3;

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");
  FieldIdentifier d1_id("D1",    {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");
  FieldIdentifier d2_id("D2",    {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");
  FieldIdentifier o_id ("O",     {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");

  auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{t_id,d1_id,d2_id,o_id},{"STARTUP"});

  // D1 and O are in the file, D2 is not.
  const std::string filename = "model_init_startup_np" + std::to_string(comm.size()) + ".nc";
  write_file(filename,{{"T_mid",300.0},{"D1",5.0},{"O",1.0}},comm,ncols,nlevs);

  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("filename",filename);
  auto& overrides = params.sublist("constant_fields").sublist("overrides");
  auto& defaults  = params.sublist("constant_fields").sublist("defaults");
  defaults.set<double>("D1",7.0);   // in the file: default is NOT used
  defaults.set<double>("D2",9.0);   // not in the file: default is used
  overrides.set<double>("O",2.0);   // in the file: file is ignored

  util::TimeStamp t0(2000,1,1,0,0,0);
  ModelInit model_init(params);
  model_init.run(fm,t0,RunType::Initial);

  for (const auto& name : {"T_mid","D1","D2","O"}) {
    REQUIRE (fm->get_field(name,gn).get_header().get_tracking().get_time_stamp()==t0);
  }
  check_constant(fm->get_field("T_mid",gn),300.0);
  check_constant(fm->get_field("D1",gn),5.0);
  check_constant(fm->get_field("D2",gn),9.0);
  check_constant(fm->get_field("O",gn),2.0);

  scorpio::finalize_subsystem();
}

TEST_CASE ("model_init_startup_missing_file", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 3;

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");
  auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{t_id},{"STARTUP"});

  // T_mid needs to be read from file, but there is no file
  ekat::ParameterList params("initial_conditions");
  ModelInit model_init(params);
  REQUIRE_THROWS (model_init.run(fm,util::TimeStamp(2000,1,1,0,0,0),RunType::Initial));
}

TEST_CASE ("model_init_subfields", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 3;

  FieldIdentifier hw_id("horiz_winds", {{COL,CMP,LEV},{ncols,2,nlevs}}, none, "point_grid");
  auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{hw_id},{"STARTUP"},false);

  // As the atm driver does, add U/V subfields, and put them (rather than their parent)
  // in the STARTUP group.
  auto hw = fm->get_field("horiz_winds",gn);
  for (int icmp=0; icmp<2; ++icmp) {
    const std::string name = icmp==0 ? "U" : "V";
    fm->add_field(hw.get_component(icmp).alias(name));
    fm->add_to_group(name,gn,"STARTUP");
  }

  util::TimeStamp t0(2000,1,1,0,0,0);

  // Cannot set a constant for a field with subfields
  {
    ekat::ParameterList params("initial_conditions");
    params.sublist("constant_fields").sublist("overrides").set<double>("horiz_winds",1.0);
    ModelInit model_init(params);
    REQUIRE_THROWS (model_init.run(fm,t0,RunType::Initial));
  }

  // But we can for each of the subfields. Once they are both inited, so is the parent.
  ekat::ParameterList params("initial_conditions");
  auto& overrides = params.sublist("constant_fields").sublist("overrides");
  overrides.set<double>("U",1.0);
  REQUIRE (not hw.get_header().get_tracking().get_time_stamp().is_valid());
  overrides.set<double>("V",2.0);

  ModelInit model_init(params);
  model_init.run(fm,t0,RunType::Initial);

  REQUIRE (fm->get_field("U",gn).get_header().get_tracking().get_time_stamp()==t0);
  REQUIRE (fm->get_field("V",gn).get_header().get_tracking().get_time_stamp()==t0);
  REQUIRE (hw.get_header().get_tracking().get_time_stamp()==t0);

  check_constant(fm->get_field("U",gn),1.0);
  check_constant(fm->get_field("V",gn),2.0);
}

TEST_CASE ("model_init_restart", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4;
  const int nlevs = 3;

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, "point_grid");

  auto [grid,fm,gn] = make_fixture(comm,ncols,nlevs,{t_id},{"STARTUP","RESTART"});

  const std::string filename = "model_init_restart_np" + std::to_string(comm.size()) + ".nc";
  write_file(filename,{{"T_mid",300.0}},comm,ncols,nlevs,{{"nsteps",42}},{{"label","restarted"}});

  // Constants must be ignored on a restart run
  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("filename",filename);
  params.sublist("constant_fields").sublist("overrides").set<double>("T_mid",1.0);

  auto nsteps = std::make_shared<std::any>(0);
  auto label  = std::make_shared<std::any>(std::string());

  ModelInit model_init(params);
  model_init.set_restart_extra_data({{"nsteps",nsteps},{"label",label}});

  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Restart);

  auto f = fm->get_field("T_mid",gn);
  REQUIRE (f.get_header().get_tracking().get_time_stamp()==t0);
  check_constant(f,300.0);

  REQUIRE (std::any_cast<int>(*nsteps)==42);
  REQUIRE (std::any_cast<std::string>(*label)=="restarted");

  scorpio::finalize_subsystem();
}

} // namespace scream
