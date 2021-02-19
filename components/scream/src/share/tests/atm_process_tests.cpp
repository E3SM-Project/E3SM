#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/remap/inverse_remapper.hpp"
#include "share/tests/dummy_se_point_remapper.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

template<AtmosphereProcessType PType>
class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const ekat::Comm& comm,const ekat::ParameterList& params)
   : m_comm(comm)
  {
    m_name = params.get<std::string> ("Process Name");
    m_grid_name = params.get<std::string> ("Grid Name");
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return PType; }

  // The type of grids on which the process is defined
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_grid_name);
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& /* field_repo */) const {}

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_fids_in; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_fids_out; }

protected:

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize_impl (const util::TimeStamp& /* t0 */ ) {}

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run_impl (const Real /* dt */) {}

  // Clean up
  void finalize_impl ( /* inputs */ ) {}

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real>& /* f */) {}
  void set_computed_field_impl (const Field<      Real>& /* f */) {}

  std::set<FieldIdentifier> m_fids_in;
  std::set<FieldIdentifier> m_fids_out;

  std::string m_name;
  std::string m_grid_name;

  ekat::Comm    m_comm;
};

class MyDynamics : public DummyProcess<AtmosphereProcessType::Dynamics>
{
public:
  using base = DummyProcess<AtmosphereProcessType::Dynamics>;

  MyDynamics (const ekat::Comm& comm,const ekat::ParameterList& params)
   : base(comm,params)
  {
    // Nothing to do here
  }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto dyn_lt = grid->get_3d_scalar_layout(true);

    m_fids_in.emplace("Temperature tendency",dyn_lt,K/s,m_grid_name);

    m_fids_out.emplace("Temperature",dyn_lt,K,m_grid_name);
  }
};

class MyPhysicsA : public DummyProcess<AtmosphereProcessType::Physics>
{
public:
  using base = DummyProcess<AtmosphereProcessType::Physics>;

  MyPhysicsA (const ekat::Comm& comm,const ekat::ParameterList& params)
   : base(comm,params)
  {
    // Nothing to do here
  }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto phys_lt = grid->get_3d_scalar_layout (true);

    m_fids_in.emplace("Temperature",phys_lt,K,m_grid_name);

    m_fids_out.emplace("Concentration A",phys_lt,kg/pow(m,3),m_grid_name);
  }
};

class MyPhysicsB : public DummyProcess<AtmosphereProcessType::Physics>
{
public:
  using base = DummyProcess<AtmosphereProcessType::Physics>;

  MyPhysicsB (const ekat::Comm& comm,const ekat::ParameterList& params)
   : base(comm,params)
  {
    // Nothing to do here
  }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto phys_lt = grid->get_3d_scalar_layout (true);

    m_fids_in.emplace("Temperature",phys_lt,K,m_grid_name);
    m_fids_in.emplace("Concentration A",phys_lt,kg/pow(m,3),m_grid_name);

    m_fids_out.emplace("Temperature tendency",phys_lt,K/s,m_grid_name);
  }
};

std::shared_ptr<UserProvidedGridsManager>
setup_upgm (const int ne) {

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_local_elems = 4;
  const int num_global_elems = 4*comm.size();
  const int np = 4;
  const int nvl = 128;
  const int ncols = 6*ne*ne*9 + 2;

  // Note: our test does not use actual dof info, but we need to set these
  //       views in the SEGrid's, so that the num local dofs is set
  AbstractGrid::dofs_list_type dyn_dofs("",num_local_elems*np*np);
  AbstractGrid::dofs_list_type phys_dofs("",ncols);

  AbstractGrid::lid_to_idx_map_type dyn_dofs_map ("",num_local_elems*np*np,3);
  AbstractGrid::lid_to_idx_map_type phys_dofs_map ("",ncols,1);

  // Create a grids manager
  auto upgm = std::make_shared<UserProvidedGridsManager>();
  auto dummy_dyn_grid  = std::make_shared<SEGrid>("Dynamics",num_global_elems,num_local_elems,np,nvl);
  auto dummy_phys_grid = create_point_grid("Physics",ncols,nvl,comm);
  dummy_dyn_grid->set_dofs(dyn_dofs,dyn_dofs_map);
  upgm->set_grid(dummy_dyn_grid);
  upgm->set_grid(dummy_phys_grid);
  upgm->set_reference_grid(dummy_phys_grid->name());

  using dummy_remapper = DummySEPointRemapper<Real>;
  using inverse_remapper = InverseRemapper<Real>;
  auto dummy_phys_dyn_remapper = std::make_shared<dummy_remapper>(dummy_phys_grid,dummy_dyn_grid);
  auto dummy_phys_dyn_remapper2 = std::make_shared<dummy_remapper>(dummy_phys_grid,dummy_dyn_grid);
  auto dummy_dyn_phys_remapper = std::make_shared<inverse_remapper>(dummy_phys_dyn_remapper2);
  upgm->set_remapper(dummy_phys_dyn_remapper);
  upgm->set_remapper(dummy_dyn_phys_remapper);

  return upgm;
}

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
  factory.register_product("MyPhysicsA",&create_atmosphere_process<MyPhysicsA>);
  factory.register_product("mYphysicsb",&create_atmosphere_process<MyPhysicsB>);
  factory.register_product("mYdynAmics",&create_atmosphere_process<MyDynamics>);
  factory.register_product("grouP",&create_atmosphere_process<AtmosphereProcessGroup>);

  // Create the processes
  std::shared_ptr<AtmosphereProcess> atm_process (factory.create("group",comm,params));

  // CHECKS
  auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_process);

  // 1) Must be a group
  REQUIRE (static_cast<bool>(group));

  // 2) Must store 2 processes: a dynamics and a group
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

  constexpr int ne = 4;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a parameter list for inputs
  ekat::ParameterList params ("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "MyDynamics");
  p0.set<std::string>("Grid Name", "Dynamics");

  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Group");
  p1.set("Number of Entries",2);
  p1.set<std::string>("Schedule Type","Sequential");

  auto& p1_0 = p1.sublist("Process 0");
  p1_0.set<std::string>("Process Name", "MyPhysicsA");
  p1_0.set<std::string>("Grid Name", "Physics");

  auto& p1_1 = p1.sublist("Process 1");
  p1_1.set<std::string>("Process Name", "MyPhysicsB");
  p1_1.set<std::string>("Grid Name", "Physics");

  // Create then factory, and register constructors
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("MyPhysicsA",&create_atmosphere_process<MyPhysicsA>);
  factory.register_product("mYphysicsb",&create_atmosphere_process<MyPhysicsB>);
  factory.register_product("mYdynAmics",&create_atmosphere_process<MyDynamics>);
  factory.register_product("grouP",&create_atmosphere_process<AtmosphereProcessGroup>);

  // Test a case where the dag has no unmet deps
  SECTION ("working") {
    // Create the processes
    std::shared_ptr<AtmosphereProcess> atm_process (factory.create("group",comm,params));

    // Create a grids manager
    auto upgm = setup_upgm (ne);

    // Set the grids, so the remappers in the group are not empty
    atm_process->set_grids(upgm);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_process),nullptr);
    dag.write_dag("working_atm_proc_dag.dot",4);

    // Clean up
    upgm->clean_up();
  }

  SECTION ("broken") {
    auto upgm = setup_upgm(ne);

    // Make it look like we forgot to request MyPhysicsA
    p1.set("Number of Entries", 1);
    p1_0.set<std::string>("Process Name", "MyPhysicsB");
    std::shared_ptr<AtmosphereProcess> broken_atm_group (factory.create("group",comm,params));
    broken_atm_group->set_grids(upgm);

    // Create the dag
    AtmProcDAG dag;
    dag.create_dag(*std::dynamic_pointer_cast<AtmosphereProcessGroup>(broken_atm_group),nullptr);
    dag.write_dag("broken_atm_proc_dag.dot",4);

    upgm->clean_up();
  }
}

} // empty namespace
