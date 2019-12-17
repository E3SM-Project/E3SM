#include <catch2/catch.hpp>
#include "share/atmosphere_process.hpp"
#include "share/atmosphere_process_group.hpp"

namespace scream {

template<AtmosphereProcessType PType>
class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const Comm& comm,const ParameterList& /* params */)
   : m_comm(comm)
  {
    // Nothing to do
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return PType; }

  // The type of grids on which the process is defined
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(e2str(GridType::Undefined));
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return e2str(PType); }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_comm; }

  void set_grids (const std::shared_ptr<const GridsManager> /* grids_manager */) {
    // Do nothing
  }

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize (const util::TimeStamp& /* t0 */ ) {}

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run (const Real /* dt */) {}

  // Clean up
  void finalize ( /* inputs */ ) {}

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& /* field_repo */) const {}

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_dummy_fids; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_dummy_fids; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& /* f */) {}
  void set_computed_field_impl (const Field<      Real, device_type>& /* f */) {}

  std::set<FieldIdentifier> m_dummy_fids;

  Comm    m_comm;
};

TEST_CASE("process_factory", "") {
  using namespace scream;

  // A world comm
  Comm comm(MPI_COMM_WORLD);

  // Create a parameter list for inputs
  ParameterList params ("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Dummy Dynamics");

  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Group");
  p1.set("Number of Entries",2);
  p1.set<std::string>("Schedule Type","Sequential");

  auto& p1_0 = p1.sublist("Process 0");
  p1_0.set<std::string>("Process Name", "Dummy Physics");

  auto& p1_1 = p1.sublist("Process 1");
  p1_1.set<std::string>("Process Name", "Dummy Physics");

  // Create then factory, and register constructors
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("duMmy pHySics",&create_atmosphere_process<DummyProcess<AtmosphereProcessType::Physics>>);
  factory.register_product("dummY dynAmics",&create_atmosphere_process<DummyProcess<AtmosphereProcessType::Dynamics>>);
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

} // empty namespace
