#include <catch2/catch.hpp>

#include "ekat/ekat_parameter_list.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

ekat::ParameterList create_test_params ()
{
  // Create a parameter list for inputs
  ekat::ParameterList params ("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Foo");
  p0.set<std::string>("Grid Name", "Point Grid");

  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Group");
  p1.set("Number of Entries",2);
  p1.set<std::string>("Schedule Type","Sequential");

  auto& p1_0 = p1.sublist("Process 0");
  p1_0.set<std::string>("Process Name", "Bar");
  p1_0.set<std::string>("Grid Name", "Point Grid");

  auto& p1_1 = p1.sublist("Process 1");
  p1_1.set<std::string>("Process Name", "Baz");
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

  ekat::ParameterList gm_params;
  gm_params.set<std::string>("Reference Grid", "Point Grid");
  gm_params.sublist("Mesh Free").set<int>("Number of Global Columns", num_global_cols);
  gm_params.sublist("Mesh Free").set<int>("Number of Local Elements", num_local_elems);
  gm_params.sublist("Mesh Free").set<int>("Number of Vertical Levels", nlevs);
  gm_params.sublist("Mesh Free").set<int>("Number of Gauss Points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_all_grids();

  return gm;
}

// A dummy atm proc
class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const ekat::Comm& comm,const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params)
  {
    m_name = params.get<std::string> ("Process Name");
    m_grid_name = params.get<std::string> ("Grid Name");
  }

  // The type of grids on which the process is defined
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_grid_name);
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

protected:

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize_impl (const RunType /* run_type */ ) {}

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run_impl (const int /* dt */) {}

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

// ================================ TESTS ============================== //

TEST_CASE("process_factory", "") {
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "atm_process_tests.yaml";
  ekat::ParameterList params ("Atmosphere Processes");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,params) );

  // Create then factory, and register constructors
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("Foo",&create_atmosphere_process<Foo>);
  factory.register_product("Bar",&create_atmosphere_process<Bar>);
  factory.register_product("Baz",&create_atmosphere_process<Baz>);
  factory.register_product("grouP",&create_atmosphere_process<AtmosphereProcessGroup>);

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

  // Test a case where the dag has no unmet deps
  SECTION ("working") {
    auto params = create_test_params ();

    // Create the processes
    std::shared_ptr<AtmosphereProcess> atm_process (factory.create("group",comm,params));

    // Set the grids, so the remappers in the group are not empty
    atm_process->set_grids(gm);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_process),{});
    dag.write_dag("working_atm_proc_dag.dot",4);

    REQUIRE (not dag.has_unmet_dependencies());
  }

  SECTION ("broken") {

    auto params = create_test_params();
    auto p1 = params.sublist("Process 1");

    // Make it look like we forgot to request MyPhysicsB
    p1.set("Number of Entries", 1);
    // p1_0.set<std::string>("Process Name", "MyPhysicsB");
    std::shared_ptr<AtmosphereProcess> broken_atm_group (factory.create("group",comm,params));
    broken_atm_group->set_grids(gm);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(broken_atm_group),{});
    dag.write_dag("broken_atm_proc_dag.dot",4);

    REQUIRE (dag.has_unmet_dependencies());
  }
}

} // empty namespace
