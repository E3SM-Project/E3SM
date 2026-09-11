#include <catch2/catch.hpp>

#include "share/data_managers/model_init.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_parameter_list.hpp>

namespace scream {

namespace {

// A grid+field_manager pair, set up on a single grid named grid_name, with
// ncols columns and nlevs levels (partitioned evenly across comm's ranks).
struct Fixture {
  std::shared_ptr<PointGrid>   grid;
  std::shared_ptr<FieldManager> fm;
};

Fixture make_fixture (const std::string& grid_name, int ncols, int nlevs, const ekat::Comm& comm)
{
  auto grid = create_point_grid(grid_name,ncols,nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>(grid);
  return { grid, std::make_shared<FieldManager>(gm) };
}

// Write a single (real-valued) variable to an already-open (and enddef'd)
// scorpio file.
void write_field_to_file (const std::string& filename, Field& f)
{
  f.sync_to_host();
  scorpio::write_var(filename, f.name(),
                     f.get_internal_view_data<const Real,Host>());
}

// Write a fresh, single-variable netcdf file, containing a (COL,LEV)
// field named varname, with 'ncol'/'lev' on-file dimension names, and every
// entry set to value.
void write_constant_2d_field (const std::string& filename,
                              const std::string& varname,
                              int ncols, int nlevs, Real value)
{
  using namespace ShortFieldTagsNames;
  FieldIdentifier fid(varname,{{COL,LEV},{ncols,nlevs}},ekat::units::none,"");
  Field f(fid); f.allocate_view(); f.deep_copy(value);

  scorpio::register_file(filename,scorpio::Write);
  scorpio::define_dim(filename,"ncol",ncols);
  scorpio::define_dim(filename,"lev",nlevs);
  scorpio::define_var(filename,varname,{"ncol","lev"},"real",false);
  scorpio::enddef(filename);

  write_field_to_file(filename,f);
  scorpio::release_file(filename);
}

} // anonymous namespace

TEST_CASE ("model_init_constants_and_copy", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4*comm.size();
  const int nlevs = 3;
  const int ncmp  = 2;

  auto [grid,fm] = make_fixture("point_grid",ncols,nlevs,comm);
  const auto& gn = grid->name();

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

  ModelInit model_init(params);

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

TEST_CASE ("model_init_startup_from_file", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4*comm.size();
  const int nlevs = 3;

  auto [grid,fm] = make_fixture("point_grid",ncols,nlevs,comm);
  const auto& gn = grid->name();

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  fm->register_field(FR{t_id});
  fm->register_group(GroupRequest("STARTUP",gn));
  fm->registration_ends();
  fm->add_to_group("T_mid",gn,"STARTUP");

  const std::string filename = "model_init_startup_np" + std::to_string(comm.size()) + ".nc";
  write_constant_2d_field(filename,"T_mid",ncols,nlevs,300.0);

  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("filename",filename);

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Initial);

  auto f = fm->get_field("T_mid",gn);
  REQUIRE (f.get_header().get_tracking().get_time_stamp()==t0);

  Field f_check(t_id); f_check.allocate_view(); f_check.deep_copy(300.0);
  REQUIRE (views_are_equal(f,f_check));

  scorpio::finalize_subsystem();
}

TEST_CASE ("model_init_group_selection_rules", "")
{
  // This test does not perform any file I/O: it checks the selection rules
  // that get_leaf_fields/get_fields use to decide which fields of a
  // monolithically allocated group need to be read/inited, and that
  // fixup_parents_time_stamp correctly propagates the time stamp to the
  // parent once all its children (which do NOT auto-update their parent's
  // time stamp) are inited.
  //
  // Since these are protected, this test uses a trivial derived class to
  // expose them.
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;
  using SL = std::list<std::string>;

  struct TestModelInit : public ModelInit {
    using ModelInit::ModelInit;
    using ModelInit::get_leaf_fields;
    using ModelInit::get_fields;
    using ModelInit::fixup_parents_time_stamp;
  };

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4*comm.size();
  const int nlevs = 3;

  auto [grid,fm] = make_fixture("point_grid",ncols,nlevs,comm);
  const auto& gn = grid->name();

  FieldIdentifier qv_id ("qv", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  FieldIdentifier qc_id ("qc", {{COL,LEV}, {ncols,nlevs}}, none, gn);

  fm->register_field(FR{qv_id,SL{"tracers"}});
  fm->register_field(FR{qc_id,SL{"tracers"}});
  fm->register_group(GroupRequest("tracers",gn,MonolithicAlloc::Required));
  fm->register_group(GroupRequest("STARTUP",gn));
  fm->register_group(GroupRequest("RESTART",gn));
  fm->registration_ends();

  // Mimic what the driver does: put the group's monolithic field (only) in
  // STARTUP/RESTART, since that's what a restart file stores; the IC file
  // stores the leaves instead, which is why get_leaf_fields must expand
  // parents.
  fm->add_to_group("tracers",gn,"STARTUP");
  fm->add_to_group("tracers",gn,"RESTART");

  ekat::ParameterList params("initial_conditions");
  TestModelInit model_init(params);

  // STARTUP: parent has children, so it must be expanded into its leaves.
  auto startup_fields = model_init.get_leaf_fields(fm,"STARTUP",gn);
  REQUIRE (startup_fields.size()==2);
  std::map<std::string,bool> found = {{"qv",false},{"qc",false}};
  for (const auto& f : startup_fields) {
    REQUIRE (found.count(f.name())==1);
    found[f.name()] = true;
  }
  REQUIRE (found.at("qv"));
  REQUIRE (found.at("qc"));

  // RESTART: the parent's children are also in the group (transitively,
  // since only the parent was added), but the parent itself must be
  // returned as-is (not expanded), since the restart file stores it whole.
  auto restart_fields = model_init.get_fields(fm,"RESTART",gn);
  REQUIRE (restart_fields.size()==1);
  REQUIRE (restart_fields[0].name()=="tracers");

  // Now manually stamp the two leaves (as init_startup_fields would, after
  // reading them from file), and check that fixup_parents_time_stamp
  // propagates the stamp to the (still un-stamped) parent.
  util::TimeStamp t0(2000,1,1,0,0,0);
  auto qv = fm->get_field("qv",gn);
  auto qc = fm->get_field("qc",gn);
  qv.get_header().get_tracking().update_time_stamp(t0);
  REQUIRE (not fm->get_field("tracers",gn).get_header().get_tracking().get_time_stamp().is_valid());

  qc.get_header().get_tracking().update_time_stamp(t0);
  REQUIRE (not fm->get_field("tracers",gn).get_header().get_tracking().get_time_stamp().is_valid());

  model_init.fixup_parents_time_stamp(fm,"STARTUP",gn,t0);
  REQUIRE (fm->get_field("tracers",gn).get_header().get_tracking().get_time_stamp()==t0);
}

TEST_CASE ("model_init_topography", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4*comm.size();

  auto [grid,fm] = make_fixture("physics_gll",ncols,3,comm);
  const auto& gn = grid->name();

  FieldIdentifier phis_id ("phis", {{COL}, {ncols}}, none, gn);
  fm->register_field(FR{phis_id});
  fm->register_group(GroupRequest("TOPOGRAPHY",gn));
  fm->registration_ends();
  fm->add_to_group("phis",gn,"TOPOGRAPHY");

  const std::string filename = "model_init_topo_np" + std::to_string(comm.size()) + ".nc";
  {
    Field f(phis_id); f.allocate_view();
    f.deep_copy(42.0);

    // The topography file uses 'ncol_d' (not 'ncol') for the GLL grid's
    // column dimension: see ModelInit::get_tag_rename. The field is also
    // named 'PHIS_d' on file: see ModelInit::get_topography_file_names.
    scorpio::register_file(filename,scorpio::Write);
    scorpio::define_dim(filename,"ncol_d",ncols);
    scorpio::define_var(filename,"PHIS_d",{"ncol_d"},"real",false);
    scorpio::enddef(filename);

    auto f_alias = f.alias("PHIS_d");
    write_field_to_file(filename,f_alias);
    scorpio::release_file(filename);
  }

  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("topography_filename",filename);

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Initial);

  auto phis = fm->get_field("phis",gn);
  REQUIRE (phis.get_header().get_tracking().get_time_stamp()==t0);

  Field phis_check(phis_id); phis_check.allocate_view(); phis_check.deep_copy(42.0);
  REQUIRE (views_are_equal(phis,phis_check));

  scorpio::finalize_subsystem();
}

TEST_CASE ("model_init_topography_not_loaded_on_restart", "")
{
  // Restart files store the full model state, so ModelInit must NOT attempt
  // to load topography on a restart run (which would also throw, since a
  // 'topography_filename' need not even be given on restart).
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4*comm.size();

  auto [grid,fm] = make_fixture("physics_gll",ncols,3,comm);
  const auto& gn = grid->name();

  FieldIdentifier phis_id ("phis", {{COL}, {ncols}}, none, gn);
  fm->register_field(FR{phis_id});
  fm->register_group(GroupRequest("TOPOGRAPHY",gn));
  fm->registration_ends();
  fm->add_to_group("phis",gn,"TOPOGRAPHY");

  // Note: no 'topography_filename' entry, so init_topography_fields would
  // throw if it were (erroneously) invoked.
  ekat::ParameterList params("initial_conditions");

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  REQUIRE_NOTHROW (model_init.run(fm,t0,RunType::Restart));

  REQUIRE (not fm->get_field("phis",gn).get_header().get_tracking().get_time_stamp().is_valid());
}

TEST_CASE ("model_init_perturbation", "")
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 4*comm.size();
  const int nlevs = 4;

  // Perturbation only ever applies to the physics_gll grid.
  auto [grid,fm] = make_fixture("physics_gll",ncols,nlevs,comm);
  const auto& gn = grid->name();

  // Set up hyam/hybm so that the reference pressure increases with level
  // index: with perturbation_minimum_pressure=1000mb below, only the
  // bottom two (of four) levels end up perturbed.
  constexpr auto ps0 = physics::Constants<Real>::P0.value;
  Field hyam(FieldIdentifier("hyam",{{LEV},{nlevs}},none,gn)); hyam.allocate_view();
  Field hybm(FieldIdentifier("hybm",{{LEV},{nlevs}},none,gn)); hybm.allocate_view();
  auto hyam_h = hyam.get_view<Real*,Host>();
  auto hybm_h = hybm.get_view<Real*,Host>();
  for (int k=0; k<nlevs; ++k) {
    const Real pref_mb = 500.0 + k*300.0; // 500, 800, 1100, 1400
    hyam_h(k) = hybm_h(k) = pref_mb*100/ps0/2; // hyam*ps0 + hybm*ps0 = pref_mb*100 (Pa)
  }
  hyam.sync_to_dev();
  hybm.sync_to_dev();
  grid->set_geometry_data(hyam);
  grid->set_geometry_data(hybm);

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  fm->register_field(FR{t_id});
  fm->register_group(GroupRequest("STARTUP",gn));
  fm->registration_ends();
  fm->add_to_group("T_mid",gn,"STARTUP");

  ekat::ParameterList params("initial_conditions");
  params.set<double>("T_mid",300.0);
  params.set<std::vector<std::string>>("perturbed_fields",{"T_mid"});
  params.set<double>("perturbation_limit",0.05);
  params.set<int>("perturbation_random_seed",123);
  params.set<double>("perturbation_minimum_pressure",1000.0);

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Initial);

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
  // and (with overwhelming probability, given a real RNG draw) not exactly
  // equal to the unperturbed value for at least one column/level.
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

TEST_CASE ("model_init_no_perturbation_on_restart", "")
{
  // A restarted field is already whatever it was at the end of the previous
  // run, so ModelInit must not perturb it again.
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR = FieldRequest;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int ncols = 4*comm.size();
  const int nlevs = 4;

  auto [grid,fm] = make_fixture("physics_gll",ncols,nlevs,comm);
  const auto& gn = grid->name();

  FieldIdentifier t_id ("T_mid", {{COL,LEV}, {ncols,nlevs}}, none, gn);
  fm->register_field(FR{t_id});
  fm->register_group(GroupRequest("RESTART",gn));
  fm->registration_ends();
  fm->add_to_group("T_mid",gn,"RESTART");

  const std::string filename = "model_init_no_perturb_restart_np" + std::to_string(comm.size()) + ".nc";
  write_constant_2d_field(filename,"T_mid",ncols,nlevs,300.0);

  // Even though a 'perturbed_fields' entry is present, it must be ignored
  // on a restart run.
  ekat::ParameterList params("initial_conditions");
  params.set<std::string>("filename",filename);
  params.set<std::vector<std::string>>("perturbed_fields",{"T_mid"});
  params.set<double>("perturbation_limit",0.05);

  ModelInit model_init(params);
  util::TimeStamp t0(2000,1,1,0,0,0);
  model_init.run(fm,t0,RunType::Restart);

  auto f = fm->get_field("T_mid",gn);
  Field f_check(t_id); f_check.allocate_view(); f_check.deep_copy(300.0);
  REQUIRE (views_are_equal(f,f_check));

  scorpio::finalize_subsystem();
}

} // namespace scream
