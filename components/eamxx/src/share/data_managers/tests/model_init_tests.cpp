#include <catch2/catch.hpp>

#include "share/data_managers/model_init.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/physics_constants.hpp"
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
                      const bool add_to_groups = true,
                      const std::string& grid_name = "point_grid")
{
  auto grid = create_point_grid(grid_name,ncols*comm.size(),nlevs,comm);
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

TEST_CASE ("model_init_topography", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4;
  const int nlevs = 3;
  const std::string gn = "physics_gll";

  FieldIdentifier phis_id ("phis", {{COL},     {ncols}},       none, gn);
  FieldIdentifier sgh_id  ("sgh",  {{COL},     {ncols}},       none, gn);
  FieldIdentifier t_id    ("T_mid",{{COL,LEV}, {ncols,nlevs}}, none, gn);

  // The topography file uses 'ncol_d' (not 'ncol') for the GLL grid's column dimension,
  // and PHIS_d (not phis) as the name of the surface geopotential.
  const std::string topo_file = "model_init_topo_np" + std::to_string(comm.size()) + ".nc";
  {
    Field f(phis_id); f.allocate_view(); f.deep_copy(42.0); f.sync_to_host();

    scorpio::register_file(topo_file,scorpio::Write);
    scorpio::define_dim(topo_file,"ncol_d",ncols*comm.size());
    scorpio::define_var(topo_file,"PHIS_d",{"ncol_d"},"real",false);
    scorpio::set_dim_decomp(topo_file,"ncol_d",comm.rank()*ncols,ncols);
    scorpio::enddef(topo_file);
    scorpio::write_var(topo_file,"PHIS_d",f.get_internal_view_data<const Real,Host>());
    scorpio::release_file(topo_file);
  }
  const std::string ic_file = "model_init_topo_ic_np" + std::to_string(comm.size()) + ".nc";
  write_file(ic_file,{{"T_mid",300.0},{"phis",-1.0}},comm,ncols,nlevs);

  util::TimeStamp t0(2000,1,1,0,0,0);

  // phis is read from the topography file, NOT from the IC file (which has a decoy),
  // and its default is not used. sgh is overridden, so there is no need to read it.
  {
    auto [grid,fm,_] = make_fixture(comm,ncols,nlevs,{phis_id,sgh_id,t_id},{"STARTUP"},true,gn);

    ekat::ParameterList params("initial_conditions");
    params.set<std::string>("filename",ic_file);
    params.set<std::string>("topography_filename",topo_file);
    params.sublist("constant_fields").sublist("defaults").set<double>("phis",10.0);
    params.sublist("constant_fields").sublist("overrides").set<double>("sgh",2.0);

    ModelInit model_init(params);
    model_init.run(fm,t0,RunType::Initial);

    for (const auto& name : {"phis","sgh","T_mid"}) {
      REQUIRE (fm->get_field(name,gn).get_header().get_tracking().get_time_stamp()==t0);
    }
    check_constant(fm->get_field("phis",gn),42.0);
    check_constant(fm->get_field("sgh",gn),2.0);
    check_constant(fm->get_field("T_mid",gn),300.0);
  }

  // Topography is requested, but no topography file is given
  {
    auto [grid,fm,_] = make_fixture(comm,ncols,nlevs,{phis_id},{"STARTUP"},true,gn);

    ekat::ParameterList params("initial_conditions");
    ModelInit model_init(params);
    REQUIRE_THROWS (model_init.run(fm,t0,RunType::Initial));
  }

  scorpio::finalize_subsystem();
}

namespace {

// Create a fixture on a physics_gll grid with hyam/hybm set so that the reference pressure
// increases with level index: with perturbation_minimum_pressure=1000mb, only the bottom
// two (of four) levels end up perturbed.
Fixture make_perturbation_fixture (const ekat::Comm& comm, const int ncols, const int nlevs,
                                   const std::vector<FieldIdentifier>& fids,
                                   const std::vector<std::string>& groups)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto fix = make_fixture(comm,ncols,nlevs,fids,groups,true,"physics_gll");

  constexpr auto ps0 = physics::Constants<Real>::P0.value;
  Field hyam(FieldIdentifier("hyam",{{LEV},{nlevs}},none,fix.gn)); hyam.allocate_view();
  Field hybm(FieldIdentifier("hybm",{{LEV},{nlevs}},none,fix.gn)); hybm.allocate_view();
  auto hyam_h = hyam.get_view<Real*,Host>();
  auto hybm_h = hybm.get_view<Real*,Host>();
  for (int k=0; k<nlevs; ++k) {
    const Real pref_mb = 500.0 + k*300.0; // 500, 800, 1100, 1400
    hyam_h(k) = hybm_h(k) = pref_mb*100/ps0/2; // hyam*ps0 + hybm*ps0 = pref_mb*100 (Pa)
  }
  hyam.sync_to_dev();
  hybm.sync_to_dev();
  fix.grid->set_geometry_data(hyam);
  fix.grid->set_geometry_data(hybm);
  return fix;
}

} // anonymous namespace

TEST_CASE ("model_init_perturbation", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 4;
  const std::string gn = "physics_gll";

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  FieldIdentifier q_id ("qv",    {{COL,LEV}, {ncols,nlevs}}, none, gn);
  auto [grid,fm,_] = make_perturbation_fixture(comm,ncols,nlevs,{t_id,q_id},{"STARTUP"});

  ekat::ParameterList params("initial_conditions");
  auto& overrides = params.sublist("constant_fields").sublist("overrides");
  overrides.set<double>("T_mid",300.0);
  overrides.set<double>("qv",0.01);
  params.set<std::vector<std::string>>("perturbed_fields",{"T_mid"});
  params.set<double>("perturbation_limit",0.05);
  params.set<int>("perturbation_random_seed",123);
  params.set<double>("perturbation_minimum_pressure",1000.0);

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Initial);

  // Only fields in perturbed_fields are perturbed
  check_constant(fm->get_field("qv",gn),0.01);

  auto f = fm->get_field("T_mid",gn);
  f.sync_to_host();
  auto v = f.get_view<const Real**,Host>();

  // Levels 0,1 (pref 500,800mb) are below the pressure threshold: untouched.
  for (int i=0; i<ncols; ++i) {
    for (int k=0; k<2; ++k) {
      REQUIRE (v(i,k)==300.0);
    }
  }

  // Levels 2,3 (pref 1100,1400mb) are perturbed: within [1-lim,1+lim]*300,
  // and (with overwhelming probability) not exactly equal to the unperturbed
  // value for at least one column/level.
  bool any_perturbed = false;
  for (int i=0; i<ncols; ++i) {
    for (int k=2; k<nlevs; ++k) {
      REQUIRE (v(i,k)>=300.0*0.95);
      REQUIRE (v(i,k)<=300.0*1.05);
      any_perturbed |= (v(i,k)!=300.0);
    }
  }
  REQUIRE (any_perturbed);
}

TEST_CASE ("model_init_perturbation_errors", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4;
  const int nlevs = 4;
  const std::string gn = "physics_gll";

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  FieldIdentifier q_id ("qv",    {{COL,LEV}, {ncols,nlevs}}, none, gn);
  util::TimeStamp t0(2000,1,1,0,0,0);

  // Perturbing a field that is not inited by ModelInit
  {
    auto [grid,fm,_] = make_perturbation_fixture(comm,ncols,nlevs,{t_id,q_id},{"STARTUP"});
    ekat::ParameterList params("initial_conditions");
    params.sublist("constant_fields").sublist("overrides").set<double>("T_mid",300.0);
    params.set<std::vector<std::string>>("perturbed_fields",{"qv"});
    ModelInit model_init(params);
    REQUIRE_THROWS (model_init.run(fm,t0,RunType::Initial));
  }

  // generate_perturbation_random_seed and perturbation_random_seed are mutually exclusive
  {
    auto [grid,fm,_] = make_perturbation_fixture(comm,ncols,nlevs,{t_id},{"STARTUP"});
    ekat::ParameterList params("initial_conditions");
    params.sublist("constant_fields").sublist("overrides").set<double>("T_mid",300.0);
    params.set<std::vector<std::string>>("perturbed_fields",{"T_mid"});
    params.set<bool>("generate_perturbation_random_seed",true);
    params.set<int>("perturbation_random_seed",1);
    ModelInit model_init(params);
    REQUIRE_THROWS (model_init.run(fm,t0,RunType::Initial));
  }
}

TEST_CASE ("model_init_no_perturbation_on_restart", "")
{
  // A restarted field is already whatever it was at the end of the previous
  // run, so ModelInit must not perturb it again.
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4;
  const int nlevs = 4;
  const std::string gn = "physics_gll";

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  auto [grid,fm,_] = make_perturbation_fixture(comm,ncols,nlevs,{t_id},{"RESTART"});

  const std::string filename = "model_init_no_perturb_restart_np" + std::to_string(comm.size()) + ".nc";
  write_file(filename,{{"T_mid",300.0}},comm,ncols,nlevs);

  // Even though a 'perturbed_fields' entry is present, it must be ignored on a restart run.
  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("filename",filename);
  params.set<std::vector<std::string>>("perturbed_fields",{"T_mid"});
  params.set<double>("perturbation_limit",0.05);

  ModelInit model_init(params);
  model_init.run(fm,util::TimeStamp(2000,1,1,0,0,0),RunType::Restart);

  check_constant(fm->get_field("T_mid",gn),300.0);

  scorpio::finalize_subsystem();
}

} // namespace scream
