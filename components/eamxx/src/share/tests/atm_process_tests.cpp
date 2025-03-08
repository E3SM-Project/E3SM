#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/remap/inverse_remapper.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_scalar_traits.hpp"

namespace scream {

ekat::ParameterList create_test_params ()
{
  using strvec_t = std::vector<std::string>;

  // Create a parameter list for inputs
  ekat::ParameterList params ("Atmosphere Processes");

  params.set<std::string>("schedule_type","Sequential");
  params.set<strvec_t>("atm_procs_list",{"Foo","BarBaz"});

  auto& p0 = params.sublist("Foo");
  p0.set<std::string>("Type", "Foo");
  p0.set<std::string>("Grid Name", "Point Grid");

  auto& p1 = params.sublist("BarBaz");
  p1.set<strvec_t>("atm_procs_list",{"Bar","Baz"});
  p1.set<std::string>("Type", "Group");
  p1.set<std::string>("schedule_type","Sequential");

  auto& p1_0 = p1.sublist("Bar");
  p1_0.set<std::string>("Type", "Bar");
  p1_0.set<std::string>("Grid Name", "Point Grid");

  auto& p1_1 = p1.sublist("Baz");
  p1_1.set<std::string>("Type", "Baz");
  p1_1.set<std::string>("Grid Name", "Point Grid");

  return params;
}

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm) {

  const int num_local_elems = 4;
  const int np = 4;
  const int nlevs = 10;
  const int num_local_cols = 13;
  const int num_global_cols = num_local_cols*comm.size();

  auto gm = create_mesh_free_grids_manager(comm,num_local_elems,np,nlevs,num_global_cols);
  gm->build_grids();

  return gm;
}

// =============================== Diagnostics ========================== //
// A dummy diagnostic
class DummyDiag : public AtmosphereDiagnostic {

public:
  DummyDiag (const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereDiagnostic(comm, params)
  {
    m_name = params.name();
    m_grid_name = params.get<std::string> ("Grid Name");
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

protected:

  void compute_diagnostic_impl () {}

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize_impl (const RunType /* run_type */ ) {}

  // Clean up
  void finalize_impl ( /* inputs */ ) {}

  std::string m_name;
  std::string m_grid_name;

};

class DiagFail : public DummyDiag
{
public:
  DiagFail (const ekat::Comm& comm, const ekat::ParameterList& params)
    : DummyDiag(comm,params)
  {
    // Nothing to do here
  }

  std::string name() const { return "Failure Dianostic"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_2d_scalar_layout ();

    add_field<Required>("Field A",lt,K,m_grid_name);
    add_field<Computed>("Field B",lt,K,m_grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, K, m_grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
  }
protected:
  void compute_diagnostic_impl () {
    // Do nothing, this diagnostic should fail.
  }
};

class DiagIdentity : public DummyDiag
{
public:
  DiagIdentity (const ekat::Comm& comm, const ekat::ParameterList& params)
    : DummyDiag(comm,params)
  {
    // Nothing to do here
  }

  std::string name() const { return "Identity Dianostic"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_2d_scalar_layout ();

    add_field<Required>("Field A",lt,K,m_grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, K, m_grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
  }
protected:
  void compute_diagnostic_impl () {
    auto f = get_field_in("Field A", m_grid_name);
    auto v_A = f.get_view<const Real*,Host>();
    auto v_me = m_diagnostic_output.get_view<Real*,Host>();
    for (size_t i=0; i<v_me.size(); ++i) {
      v_me[i] = v_A[i];
    }
  }
};

class DiagSum : public DummyDiag
{
public:
  DiagSum (const ekat::Comm& comm, const ekat::ParameterList&params)
    : DummyDiag(comm,params)
  {
    // Nothing to do here
  }

  std::string name() const { return "Summation Dianostic"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_2d_scalar_layout ();

    add_field<Required>("Field A",lt,K,m_grid_name);
    add_field<Required>("Field B",lt,K,m_grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, K, m_grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
  }

protected:
  void compute_diagnostic_impl () {
    auto f_A = get_field_in("Field A", m_grid_name);
    auto f_B = get_field_in("Field B", m_grid_name);
    auto v_A = f_A.get_view<const Real*,Host>();
    auto v_B = f_B.get_view<const Real*,Host>();
    auto v_me = m_diagnostic_output.get_view<Real*,Host>();
    for (size_t i=0; i<v_me.size(); ++i) {
      v_me[i] = v_A[i]+v_B[i];
    }
  }
};
// =============================== Processes ========================== //
// A dummy atm proc
class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const ekat::Comm& comm,const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params)
  {
    m_name = params.name();
    m_grid_name = params.get<std::string> ("Grid Name");
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

protected:

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize_impl (const RunType /* run_type */ ) {}

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run_impl (const double /* dt */) {}

  // Clean up
  void finalize_impl ( /* inputs */ ) {}

  std::string m_name;
  std::string m_grid_name;
};

class Foo : public DummyProcess
{
public:

  Foo (const ekat::Comm& comm,const ekat::ParameterList& params)
   : DummyProcess(comm,params)
  {
    // Nothing to do here
  }

  // The type of the atm proc
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_3d_scalar_layout(true);

    add_field<Required>("Temperature tendency",lt,K/s,m_grid_name);
    add_field<Computed>("Temperature",lt,K,m_grid_name);
  }
};

class Bar : public DummyProcess
{
public:
  Bar (const ekat::Comm& comm,const ekat::ParameterList& params)
   : DummyProcess(comm,params)
  {
    // Nothing to do here
  }

  // The type of the atm proc
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_3d_scalar_layout (true);

    add_field<Required>("Temperature",lt,K,m_grid_name);
    add_field<Computed>("Concentration A",lt,kg/pow(m,3),m_grid_name);
  }
};

class Baz : public DummyProcess
{
public:
  Baz (const ekat::Comm& comm,const ekat::ParameterList& params)
   : DummyProcess(comm,params)
  {
    // Nothing to do here
  }

  // The type of the atm proc
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto phys_lt = grid->get_3d_scalar_layout (true);

    add_field<Required>("Temperature",phys_lt,K,m_grid_name);
    add_field<Required>("Concentration A",phys_lt,kg/pow(m,3),m_grid_name);

    add_field<Computed>("Temperature tendency",phys_lt,K/s,m_grid_name);
  }
};

class AddOne : public DummyProcess
{
public:
  AddOne (const ekat::Comm& comm,const ekat::ParameterList& params)
   : DummyProcess(comm,params)
  {
    // Nothing to do here
  }

  // The type of the atm proc
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_2d_scalar_layout ();

    add_field<Updated>("Field A",lt,K,m_grid_name);
  }
protected:
    void run_impl (const double /* dt */) {
    auto v = get_field_out("Field A", m_grid_name).get_view<Real*,Host>();

    for (int i=0; i<v.extent_int(0); ++i) {
      v[i] += Real(1.0);
    }
  }
};

// ================================ TESTS ============================== //

TEST_CASE("process_factory", "") {
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create then factory, and register constructors
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("Foo",&create_atmosphere_process<Foo>);
  factory.register_product("Bar",&create_atmosphere_process<Bar>);
  factory.register_product("Baz",&create_atmosphere_process<Baz>);
  factory.register_product("grouP",&create_atmosphere_process<AtmosphereProcessGroup>);
  factory.register_product("DiagIdentity",&create_atmosphere_process<DiagIdentity>);

  // Load ad parameter list
  std::string fname = "atm_process_tests_named_procs.yaml";
  ekat::ParameterList params ("Atmosphere Processes");
  parse_yaml_file(fname,params);

  // Create the processes
  std::shared_ptr<AtmosphereProcess> atm_process (factory.create("group",comm,params));

  // CHECKS
  auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_process);

  // 1) Must be a group
  REQUIRE (static_cast<bool>(group));

  // 2) Must store 2 processes: a Physics and a Group
  REQUIRE (group->get_num_processes()==2);
  REQUIRE (group->get_process(0)->type()==AtmosphereProcessType::Dynamics);
  REQUIRE (group->get_process(1)->type()==AtmosphereProcessType::Group);

  // 3) The group must store two physics
  auto group_2 = std::dynamic_pointer_cast<const AtmosphereProcessGroup>(group->get_process(1));
  REQUIRE (static_cast<bool>(group_2));
  REQUIRE (group_2->get_num_processes()==2);
  REQUIRE (group_2->get_process(0)->type()==AtmosphereProcessType::Physics);
  REQUIRE (group_2->get_process(1)->type()==AtmosphereProcessType::Physics);
}

TEST_CASE("atm_proc_dag", "") {
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create then factory, and register constructors
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("Foo",&create_atmosphere_process<Foo>);
  factory.register_product("Bar",&create_atmosphere_process<Bar>);
  factory.register_product("Baz",&create_atmosphere_process<Baz>);
  factory.register_product("grouP",&create_atmosphere_process<AtmosphereProcessGroup>);

  // Create a grids manager
  auto gm = create_gm(comm);

  auto create_fields = [](const AtmosphereProcess& ap)
    -> std::map<std::string,Field>
  {
    std::map<std::string,Field> fields;
    for (auto r : ap.get_required_field_requests()) {
      fields[r.fid.name()] = Field(r.fid);
      fields[r.fid.name()].allocate_view();
    }
    for (auto r : ap.get_computed_field_requests()) {
      fields[r.fid.name()] = Field(r.fid);
      fields[r.fid.name()].allocate_view();
    }
    return fields;
  };

  auto create_and_set_fields = [&] (AtmosphereProcess& ap)
  {
    auto fields = create_fields(ap);
    for (auto r : ap.get_required_field_requests()) {
      ap.set_required_field(fields.at(r.fid.name()));
    }
    for (auto r : ap.get_computed_field_requests()) {
      ap.set_computed_field(fields.at(r.fid.name()));
    }
  };

  // Test a case where the dag has no unmet deps
  SECTION ("working") {
    auto params = create_test_params ();

    // Create the processes
    std::shared_ptr<AtmosphereProcess> atm_process (factory.create("group",comm,params));

    // Set the grids, so the remappers in the group are not empty
    atm_process->set_grids(gm);

    create_and_set_fields (*atm_process);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_process));
    dag.write_dag("working_atm_proc_dag.dot",4);

    REQUIRE (not dag.has_unmet_dependencies());
  }

  SECTION ("broken") {

    using strvec_t = std::vector<std::string>;
    auto params = create_test_params();
    auto& p1 = params.sublist("BarBaz");

    // Make sure there's a missing piece (whatever Baz computes);
    p1.set<strvec_t>("atm_procs_list",{"Bar"});
    std::shared_ptr<AtmosphereProcess> broken_atm_group (factory.create("group",comm,params));
    broken_atm_group->set_grids(gm);

    create_and_set_fields (*broken_atm_group);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(broken_atm_group));
    dag.write_dag("broken_atm_proc_dag.dot",4);

    REQUIRE (dag.has_unmet_dependencies());
  }
}

TEST_CASE("field_checks", "") {
  using namespace scream;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  auto gm = create_gm(comm);
  auto grid = gm->get_grid("Point Grid");

  // Create a parameter list
  ekat::ParameterList params ("Atmosphere Processes");
  params.set<std::string>("Grid Name", "Point Grid");

  const auto lt = grid->get_3d_scalar_layout(true);
  FieldIdentifier fid_T_tend("Temperature tendency",lt,K/s,"Point Grid");
  FieldIdentifier fid_T("Temperature",lt,K,"Point Grid");
  Field T(fid_T), T_tend(fid_T_tend);
  T.allocate_view();
  T_tend.allocate_view();
  T_tend.deep_copy(-1.0);
  T.deep_copy(-1.0);
  util::TimeStamp t0(1,1,1,1,1,1);

  constexpr auto Warning = CheckFailHandling::Warning;
  auto pos_check_pre = std::make_shared<FieldLowerBoundCheck>(T_tend,grid,0,false);
  auto pos_check_post = std::make_shared<FieldLowerBoundCheck>(T,grid,0,false);
  for (bool allow_failure : {true,false}) {
    for (bool check_pre : {true, false}) {
      for (bool check_post : {true, false}) {

        params.set("enable_precondition_checks",check_pre);
        params.set("enable_postcondition_checks",check_post);

        // Create the process
        auto foo = create_atmosphere_process<Foo>(comm,params);
        foo->set_grids(gm);

        foo->set_required_field(T_tend);
        foo->set_computed_field(T);
        foo->initialize(t0,RunType::Initial);

        if (allow_failure) {
          foo->add_precondition_check(pos_check_pre,Warning);
          foo->add_postcondition_check(pos_check_post,Warning);
        } else {
          // Default CheckFailHandling is Fatal
          foo->add_precondition_check(pos_check_pre);
          foo->add_postcondition_check(pos_check_post);
        }

        if (not allow_failure && (check_pre || check_post)) {
          REQUIRE_THROWS (foo->run(1));
        } else {
          foo->run(1);
        }
      }
    }
  }
}

TEST_CASE ("subcycling") {
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager
  auto gm = create_gm(comm);

  ekat::ParameterList params, params_sub;
  params.set<std::string>("Grid Name", "Point Grid");
  params_sub.set<std::string>("Grid Name", "Point Grid");
  params_sub.set<int>("number_of_subcycles", 5);

  // Create and init two atm procs, one subcycled and one not subcycled
  auto ap     = std::make_shared<AddOne>(comm,params);
  auto ap_sub = std::make_shared<AddOne>(comm,params_sub);

  ap->set_grids(gm);
  ap_sub->set_grids(gm);

  // Create fields (should be just one) and set it in the atm procs
  for(const auto& req : ap->get_required_field_requests()) {
    Field f(req.fid);
    f.allocate_view();
    f.deep_copy(0);
    f.get_header().get_tracking().update_time_stamp(t0);
    ap->set_required_field(f.get_const());
    ap->set_computed_field(f);

    Field f_sub(req.fid);
    f_sub.allocate_view();
    f_sub.deep_copy(0);
    f_sub.get_header().get_tracking().update_time_stamp(t0);
    ap_sub->set_required_field(f_sub.get_const());
    ap_sub->set_computed_field(f_sub);
  }

  ap->initialize(t0,RunType::Initial);
  ap_sub->initialize(t0,RunType::Initial);

  // Now run both procs for dt=5.
  const int dt = 5;
  ap->run(dt);
  ap_sub->run(dt);

  // Now, ap_sub should have added one 5 times, while ap only once
  auto v = ap->get_fields_in().front().get_view<const Real*,Host>();
  auto v_sub = ap_sub->get_fields_in().front().get_view<const Real*,Host>();

  // Safety check
  REQUIRE (v.size()==v_sub.size());
  for (size_t i=0; i<v.size(); ++i) {
    REQUIRE (v_sub[i]==5*v[i]);
  }
}

TEST_CASE ("diagnostics") {

  //TODO: This test needs a field manager so that changes in Field A are seen everywhere.
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager
  auto gm = create_gm(comm);

  // Create the identity diagnostic
  ekat::ParameterList params_identity("DiagIdentity");
  params_identity.set<std::string>("Grid Name", "Point Grid");
  auto diag_identity = std::make_shared<DiagIdentity>(comm,params_identity);
  diag_identity->set_grids(gm);

  // Create the sum diagnostic
  ekat::ParameterList params_sum("DiagSum");
  params_sum.set<std::string>("Grid Name", "Point Grid");
  auto diag_sum = std::make_shared<DiagSum>(comm,params_sum);
  diag_sum->set_grids(gm);

  // Create the fail diagnostic
  ekat::ParameterList params_fail("DiagFail");
  params_fail.set<std::string>("Grid Name", "Point Grid");
  auto diag_fail = std::make_shared<DiagFail>(comm,params_fail);
  diag_fail->set_grids(gm);

  std::map<std::string,Field> input_fields;
  for (const auto& req : diag_sum->get_required_field_requests()) {
    Field f(req.fid);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag_sum->set_required_field(f.get_const());
    REQUIRE_THROWS(diag_fail->set_computed_field(f));
    if (name == "Field A") {
      diag_identity->set_required_field(f.get_const());
      f.deep_copy<double,Host>(1.0);
    } else {
      f.deep_copy<double,Host>(2.0);
    } 
    input_fields.emplace(name,f);
  }
  auto f_A        = input_fields["Field A"];
  auto f_B        = input_fields["Field B"];
  auto v_A        = f_A.get_view<Real*,Host>();
  auto v_B        = f_B.get_view<Real*,Host>();

  diag_identity->initialize(t0,RunType::Initial);
  diag_sum->initialize(t0,RunType::Initial);

  // Run the diagnostics
  diag_identity->compute_diagnostic();
  diag_sum->compute_diagnostic();

  // Get diagnostics outputs
  const auto& f_identity = diag_identity->get_diagnostic();
  const auto& f_sum      = diag_sum->get_diagnostic();

  // For the identity diagnostic check that the fields match
  auto v_identity = f_identity.get_view<const Real*,Host>();
  auto v_sum      = f_sum.get_view<const Real*,Host>();
  REQUIRE (v_A.size()==v_identity.size());
  REQUIRE (v_A.size()==v_sum.size());
  REQUIRE (v_B.size()==v_sum.size());
  for (size_t i=0; i<v_A.size(); ++i) {
    REQUIRE (v_identity[i]==v_A[i]);
    REQUIRE (v_sum[i]==v_A[i]+v_B[i]);
    REQUIRE (v_sum[i]==3);
  }
}

} // empty namespace
