#include <catch2/catch.hpp>
#include "share/atmosphere_process.hpp"
#include "control/atmosphere_driver.hpp"
#include "control/surface_coupling.hpp"

namespace scream {

template<AtmosphereProcessType PType>
class DummyProcess : public scream::AtmosphereProcess {
public:

  explicit DummyProcess (const ParameterList& /* params */) {}

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return PType; }

  // Return some sort of name, linked to PType
  std::string name () const { return e2str(PType); }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_comm; }

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize (const Comm& comm) { m_comm = comm; }

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run ( /* inputs? */ ) {}

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

template<AtmosphereProcessType PType>
AtmosphereProcess* create_dummy_process (const ParameterList& p) {
  return new DummyProcess<PType>(p);
}

TEST_CASE("process_factory", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a parameter list for inputs
  ParameterList ad_params("Atmosphere Driver");
  auto& params = ad_params.sublist("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Dummy Physics");

  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Dummy Dynamics");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& factory = AtmosphereProcessFactory::instance();
  factory.register_product("duMmy pHySics",&create_dummy_process<AtmosphereProcessType::Physics>);
  factory.register_product("dummY dynAmics",&create_dummy_process<AtmosphereProcessType::Dynamics>);
  factory.register_product("Surface Coupling",&create_surface_coupling);

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  AtmosphereDriver ad;
  ad.initialize(atm_comm,ad_params);

  // If we get here, the ad initialized correctly
  REQUIRE (true);
}

} // empty namespace
