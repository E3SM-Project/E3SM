#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"

#include "share/grid/point_grid.hpp"
#include "share/data_managers/library_grids_manager.hpp"

#include <ekat_parameter_list.hpp>

namespace scream {

ekat::ParameterList create_test_params ()
{
  using strvec_t = std::vector<std::string>;

  // Create a parameter list for inputs
  ekat::ParameterList params ("Atmosphere Processes");

  params.set<std::string>("schedule_type","sequential");
  params.set<strvec_t>("atm_procs_list",{"Foo","BarBaz"});

  auto& p0 = params.sublist("Foo");
  p0.set<std::string>("type", "Foo");
  p0.set<std::string>("grid_name", "point_grid");

  auto& p1 = params.sublist("BarBaz");
  p1.set<strvec_t>("atm_procs_list",{"Bar","Baz"});
  p1.set<std::string>("type", "group");
  p1.set<std::string>("schedule_type","sequential");

  auto& p1_0 = p1.sublist("Bar");
  p1_0.set<std::string>("type", "Bar");
  p1_0.set<std::string>("grid_name", "point_grid");

  auto& p1_1 = p1.sublist("Baz");
  p1_1.set<std::string>("type", "Baz");
  p1_1.set<std::string>("grid_name", "point_grid");

  return params;
}

// =============================== Processes ========================== //
// A dummy atm proc
class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const ekat::Comm& comm,const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params)
  {
    m_name = params.name();
    m_grid_name = params.get<std::string> ("grid_name");
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
  const int nlcols = 3;
  const int nlevs = 10;
  auto grid = create_point_grid ("point_grid",nlcols*comm.size(),nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>(grid);

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

} // empty namespace
